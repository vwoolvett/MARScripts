# Module to find and reduce calibrations upon new scan arrival for watchdog
# at APEX

def __watchdog__findCalScans(fe='LFA', verbose=False):
    assert fe=='LFA' or fe=='HFA', 'fe must be LFA or HFA'

    # Create Blacklist if non-existent
    if not os.path.exists('.watchdog_blacklist.txt'):
            with open('.watchdog_blacklist.txt', 'w') as f:
                    f.write('# List of calibration scans that have timed out\n')
                    f.write('# in their reduction to avoid reducing them again\n')
                    f.write('# Example:\n')
                    f.write('# 12345,LFA\n')
                    f.write('# 12346,LFA\n')
                    f.write('# 12347,HFA\n')

    # Read blacklist
    with open('.watchdog_blacklist.txt', 'r') as f:
        lines = f.readlines()
        blacklist = [int(line.split(',')[0]) for line in lines if line.split(',')[0].isdigit() and line.split(',')[1]==fe]

    # no HFA-only mode, so:
    FeBedict = {'LFA': 'AMKID870-AMKID870BE', 'HFA':'AMKID350-AMKID350BE'}
    febe = FeBedict[fe]

    calscanlist = []
    fsweep_list = []
    wscan_list = []
    messages = []
    files = os.listdir('../obslogs')
    c=0
    for file in files:
        fullfilename = '../obslogs/' + file
        f = open(fullfilename,'r')
        lines = f.readlines()
        index = 0
        start = False
        keys = []
        for index in range(len(lines)):
            line = lines[index]
            if line[0:4]=='<th>':
                keys.append(line[4:-6])
                index+=1
            elif line[0:4]=='<tr>':
                start=True
                index+=1    
            elif line[0:5]=='</tr>':
                index+=1
            else:
                index+=1
            if start:
                message=''
                scan=0 
                for key in keys:
                    line=lines[index]
                    index+=1 
                    if key == 'Scan':                               # 0
                        scan_int = int(line[4:-6])
                        scan = (line[4:-6].ljust(6) + ' | ')

                    if key == 'Source':                             # 1 
                        src = (line[4:-6].ljust(12) + ' | ') 

                    if key == 'Scan type':                          # 2
                        scantype = (line[4:-6].ljust(12) + ' | ')

                    if key == 'Observ. mode':                       # 3
                        mode = (line[4:-6].ljust(12) + ' | ')

                    if key == 'Frontend-backend':                   # 4
                        thisFeBe = line[4:-6]
                        if FeBedict['HFA'] in thisFeBe:
                            thisfe = ('LFA + HFA'.ljust(9) + ' | ')
                        elif FeBedict['LFA'] in thisFeBe:
                            thisfe = ('LFA'.ljust(9) + ' | ')
                        else:
                            thisfe = ('NOT AMKID!'.ljust(9) + ' | ')

                    if key == 'Command':
                        command = (line[4:-6].ljust(12) + ' | ')
                        # only add first 12 characters
                        command = command[:12] + ' | '              # 5

                    if key == 'Scan duration':
                        duration = (line[4:-6].ljust(8) + ' | ')    # 6

                    if key == 'Scan status':
                        status = (line[4:-6].ljust(12) + ' | ')     # 7

                    if key == 'Comment':                            # last
                        comment = (line[4:-6])[0:85]       
                        message += scan + src + scantype + mode + thisfe +\
                                   command + duration + status

                start = False

                # If it is a calibration
                if (('calibrate(' in command)\
                    or ('beamscan(' in command))\
                    and 'go(' not in command\
                    and 'load(' not in command\
                    and fe in thisfe\
                    and 'OK' in status:

                    # and is finished and not in blacklist, consider it
                    if '-999' not in duration and scan_int not in blacklist:
                        message += 'SCAN CONSIDERED'.ljust(15) + ' | '
                        message += comment
                        calscanlist.append(scan_int)
                        if 'calibrate(' in command:
                            fsweep_list.append(scan_int)
                        if 'beamscan(' in command:
                            wscan_list.append(scan_int)

                    # if is finished but blacklisted for timeout, mention
                    elif '-999' not in duration and scan_int not in blacklist:
                        message += 'REDUCTION TIMEOUT'.ljust(15) + ' | '
                        message += comment

                    # if ongoing, don't add
                    else:
                        message += 'SCAN ONGOING'.ljust(15) + ' | '
                        message += comment

                # not a calibration
                else:
                    message += 'SCAN DISCARDED'.ljust(15) + ' | '
                    message += comment

                if len(keys) != 0:
                    keys_toshow = keys
                    c+=1
                messages.append(message)

    calscanarr = np.array(calscanlist)
    order = calscanarr.argsort()
    messages_sorted = [messages[i] for i in order]
    if verbose:
        for message in messages_sorted:
            print(message)
        print('============')
        print('OBSLOG KEYS:')
        print('============')
        print(keys_toshow)

    calscanlist.sort()
    fsweep_list.sort()
    wscan_list.sort()
    return calscanlist, fsweep_list, wscan_list



def __watchdog__timeout_handler(signum, frame):
    raise TimeoutError("Reduction timed out!")



def __watchdog__write_to_blacklist(scan, fe, caltype, exception, timeout=0.0):
    # Create Blacklist if non-existent
    if not os.path.exists('.watchdog_blacklist.txt'):
        with open('.watchdog_blacklist.txt', 'w') as f:
            f.write('# List of calibration scans that have timed out\n')
            f.write('# in their reduction to avoid reducing them again\n')
            f.write('# Example:\n')
            f.write('# 12345,LFA,Fsweep,TimeoutError,180\n')
            f.write('# 12346,LFA,WireScanner,TimeoutError,180\n')
            f.write('# 12347,HFA,Fsweep,OtherError,0\n')

    with open('.watchdog_blacklist.txt', 'r') as f:
        lines = f.readlines()
        blacklist = [int(line.split(',')[0]) for line in lines if line.split(',')[0].isdigit() and line.split(',')[1]==fe]
    
    with open('.watchdog_blacklist.txt', 'a') as f:
        if scan not in blacklist:
            f.write("%i,%s,%s,%s,%.3f\n"%(scan, fe, caltype, type(exception).__name__, timeout))



def __reduce_Fsweep_list(fsweep_list, fe='LFA', overwrite=False, timeout=180.):
    import signal as temp
    # track reductions
    total = len(fsweep_list)
    already_reduced = 0

    # Iterate over fsweep list
    for fsweep in fsweep_list:
        # check if there is a reduction
        globlist = glob('CalFiles/fsweep_%s_%s.csv'%(fsweep, fe))
        if len(globlist)==0:
            if scanIsFsweep(fsweep, fe=fe):
                temp.signal(temp.SIGALRM, __watchdog__timeout_handler)
                try:
                    temp.alarm(timeout)
                    reduceFsweep(fsweep, fe=fe)
                except TimeoutError as e:
                    __watchdog__write_to_blacklist(fsweep,fe,'Fsweep',e,timeout)
                except Exception as e:
                    __watchdog__write_to_blacklist(fsweep,fe,'Fsweep',e,0.)
                finally:
                    temp.alarm(0)
                
        if len(globlist)==1 and overwrite==True:
            if scanIsFsweep(fsweep, fe=fe):
                temp.signal(temp.SIGALRM, __watchdog__timeout_handler)
                try:
                    temp.alarm(timeout)
                    reduceFsweep(fsweep, fe=fe, overwrite=True)
                except TimeoutError as e:
                    __watchdog__write_to_blacklist(fsweep,fe,'WireScanner',e,timeout)
                except Exception as e:
                    __watchdog__write_to_blacklist(fsweep,fe,'WireScanner',e,0.)
                finally:
                    temp.alarm(0)
        else:
            already_reduced += 1

    if already_reduced==total:
        info('No new WireScans to reduce...')



def __reduce_WireScan_list(wscan_list, fe='LFA', doPlot=True, overwrite=False, timeout=180.):
    import signal as temp
    if doPlot==True:
        figure_initiated = False

    # track reductions
    total = len(wscan_list)
    already_reduced = 0

    # Iterate over wscan list
    for wscan in wscan_list:
        # check if there is a reduction
        globlist = glob('CalFiles/beam_scanner_%s_%s.csv'%(wscan, fe))
        if len(globlist)==0:
            if scanIsBeamscan(wscan, fe=fe):
                if doPlot==True and figure_initiated==False:
                    figure = plt.figure(figsize=(10,10))
                    figure_initiated=True
                temp.signal(temp.SIGALRM, __watchdog__timeout_handler)
                try:
                    temp.alarm(timeout)
                    reduceWireScan(wscan, fe=fe, doPlot=True, fig=figure)
                except TimeoutError as e:
                    __watchdog__write_to_blacklist(wscan,fe,e,timeout)
                except Exception as e:
                    __watchdog__write_to_blacklist(wscan,fe,e,0.)
                finally:
                    temp.alarm(0)

        if len(globlist)==1 and overwrite==True:
            if scanIsBeamscan(wscan, fe=fe):
                if doPlot==True and figure_initiated==False:
                    figure = plt.figure(figsize=(10,10))
                    figure_initiated=True
                temp.signal(temp.SIGALRM, __watchdog__timeout_handler)
                try:
                    temp.alarm(timeout)
                    reduceWireScan(wscan, fe=fe, overwrite=True, doPlot=True, fig=figure)
                except TimeoutError as e:
                    __watchdog__write_to_blacklist(wscan,fe,e,timeout)
                except Exception as e:
                    __watchdog__write_to_blacklist(wscan,fe,e,0.)
                finally:
                    temp.alarm(0)
        else:
            already_reduced += 1

    if already_reduced==total:
        info('No new WireScans to reduce...')


def __watchdog_run(fe='LFA', mode='loop', timeout=180.):

    if mode=='single':
        info('Watchdog in single-run mode.')
        calscan_list, fsweep_list, wscan_list = __watchdog__findCalScans(fe=fe)
        __reduce_Fsweep_list(fsweep_list, fe=fe, timeout=timeout)
        __reduce_WireScan_list(wscan_list, fe=fe, timeout=timeout)

    if mode=='loop':
        info('Watchdog in loop mode, CTRL+F to exit!')
        time_0 = time.time()
        calscan_list, fsweep_list, wscan_list = __watchdog__findCalScans(fe=fe)
        __reduce_Fsweep_list(fsweep_list, fe=fe, timeout=timeout)
        __reduce_WireScan_list(wscan_list, fe=fe, timeout=timeout)
        while True:
            time_current = time.time()
            if time_current - time_0 > 10:  # at least 10 seconds
                __watchdog_run(fe=fe, mode='loop')