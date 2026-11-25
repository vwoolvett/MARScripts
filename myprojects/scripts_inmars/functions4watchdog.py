# Module to find and reduce calibrations upon new scan arrival for watchdog
# at APEX

def __watchdog__timeout_handler(signum, frame):
    raise TimeoutError("Reduction timed out!")



def __watchdog__checkblacklist():
    if not os.path.exists('.watchdog_blacklist.txt'):
        with open('.watchdog_blacklist.txt', 'w') as f:
            f.write('# List of calibration scans that have failed\n')
            f.write('# in their reduction to avoid reducing them again\n')
            f.write('# Scan, Array, CalType, Timeout, Error, Errordescription\n')



def __watchdog__write_to_blacklist(scan, fe, caltype, error, timeout=0):
    fe = str.upper(fe)
    timeout = int(timeout)
    # Create Blacklist if non-existent
    __watchdog__checkblacklist()

    with open('.watchdog_blacklist.txt', 'r') as f:
        lines = f.readlines()
        blacklist = __watchdog_readblacklist(fe=fe)
    
    with open('.watchdog_blacklist.txt', 'a') as f:
        if scan not in blacklist:
            if type(error)==str:
                errorstr = error
                errordesc = 'void'
            else:
                errorstr = type(error).__name__
                errordesc = str(error)
            f.write("%i,%s,%s,%s,%s,%i\n"%(scan, fe, caltype, timeout, errorstr, errordesc))



def __watchdog_readblacklist(fe='LFA', returnDict=False):
    fe = str.upper(fe)
    # Create Blacklist if non-existent
    __watchdog__checkblacklist()

    # Read blacklist
    with open('.watchdog_blacklist.txt', 'r') as f:
        blacklistdict = {}
        lines = f.readlines()
        for line in lines:
            if line.split(',')[0].isdigit() and str.upper(line.split(',')[1])==fe:
                calscan = int(line.split(',')[0])
                calscanfe = line.split(',')[1]
                calscantype = line.split(',')[2]
                timeout = line.split(',')[3]
                calscanerror = line.split(',')[4]
                calscanerrordesc = line.split(',')[5]

                blacklistdict[calscan] = {'scan':calscan,
                                          'fe':calscanfe,
                                          'type':calscantype,
                                          'timeout':timeout,
                                          'error':calscanerror,
                                          'errordesc':calscanerrordesc
                                          }

        if returnDict:
            return blacklistdict

        else:
            blacklist = [blacklisted for blacklisted in blacklistdict.keys()]
            blacklist.sort()
            return blacklist



def __watchdog_flash_figure(fig, textstr='', fontsize=15, flashcolor='red', flashtime=0.3):
    assert type(textstr)==str
    fig.clf()
    fig.text(0.5, 0.5, textstr, fontdict={'fontsize':fontsize}, horizontalalignment='center', color='white')
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.pause(flashtime)
    frames = np.arange(10)
    for frame in frames:
        if frame % 2 == 0:
            color = flashcolor
        else:
            color = 'black'
        fig.clf()
        fig.set_facecolor(color)
        fig.text(0.5, 0.5, textstr, fontdict={'fontsize':fontsize}, horizontalalignment='center', color='white')
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(flashtime)



def __watchdog__findCalScans(fe='LFA', verbose=False):
    fe = str.upper(fe)
    assert fe=='LFA' or fe=='HFA', 'fe must be LFA or HFA'

    # Create Blacklist if non-existent
    __watchdog__checkblacklist()

    # Read blacklist
    blacklist = __watchdog_readblacklist(fe=fe)

    # no HFA-only mode, so:
    FeBedict = {'LFA': 'AMKID870-AMKID870BE', 'HFA':'AMKID350-AMKID350BE'}
    febe = FeBedict[fe]

    scanlist = []
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
                    elif '-999' not in duration and scan_int in blacklist:
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
                scanlist.append(scan_int)
                messages.append(message)

    scanarr = np.array(scanlist)
    order = scanarr.argsort()
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



def __reduce_Fsweep_list(fsweep_list, fe='LFA', overwrite=False, timeout=60):
    import signal as temp
    fe = str.upper(fe)

    # info print
    info('Processing Fsweeps (%s)'%fe)

    # track reductions
    total = len(fsweep_list)
    already_reduced = 0

    # Iterate over fsweep list
    for fsweep in fsweep_list:
        temp.signal(temp.SIGALRM, __watchdog__timeout_handler)
        try:
            temp.alarm(timeout)
            # check if there is a reduction
            globlist = glob('CalFiles/fsweep_%s_%s.csv'%(fsweep, fe))

            if len(globlist)==0:
                fD = reduceFsweep(fsweep, fe=fe, doReturn=True)
                if fD == None:  # when not Fsweep
                    raise AssertionError
                
            if len(globlist)==1 and overwrite==True:
                fD = reduceFsweep(fsweep, fe=fe, doReturn=True, overwrite=True)
                if fD == None:  # when not Fsweep
                    raise AssertionError

            else:
                already_reduced += 1

        except AssertionError:
            __watchdog__write_to_blacklist(fsweep,fe,'Fsweep','NotFsweep',0)
        except TimeoutError as e:
            __watchdog__write_to_blacklist(fsweep,fe,'Fsweep',e,timeout)
        except Exception as e:
            __watchdog__write_to_blacklist(fsweep,fe,'Fsweep',e,0)
        finally:
            temp.alarm(0)

    if already_reduced==total:
        info('No new Fsweeps to reduce.')



def __reduce_WireScan_list(wscan_list, fe='LFA', doPlot=True, fig=None, overwrite=False, timeout=60):
    import signal as temp
    fe = str.upper(fe)

    # info print
    info('Processing WireScans (%s)'%fe)

    # track reductions
    total = len(wscan_list)
    already_reduced = 0

    # Iterate over wscan list
    for wscan in wscan_list:
        temp.signal(temp.SIGALRM, __watchdog__timeout_handler)
        try:
            temp.alarm(timeout)
            # check if there is a reduction
            globlist = glob('CalFiles/beam_scanner_%s_%s.csv'%(wscan, fe))
            
            if len(globlist)==0:
                if doPlot and fig is not None:
                    mytext = 'Now processing WireScan %s (%s)!'%(wscan, fe)
                    __watchdog_flash_figure(fig=fig, textstr=mytext)
                wD = reduceWireScan(wscan, fe=fe, doPlot=doPlot, fig=fig, doReturn=True)
                if wD == {}:  # when not WireScan
                    raise AssertionError
                if doPlot and fig is not None:
                    fig.canvas.draw()
                    fig.canvas.flush_events()
                    plt.pause(0.1)

            if len(globlist)==1 and overwrite==True:
                if doPlot and fig is not None:
                    mytext = 'Now processing WireScan %s (%s)!'%(wscan, fe)
                    __watchdog_flash_figure(fig=fig, textstr=mytext)
                wD = reduceWireScan(wscan, fe=fe, overwrite=True, doPlot=doPlot, fig=fig, doReturn=True)
                if wD == {}:  # when not WireScan
                    raise AssertionError
                if doPlot and fig is not None:
                    fig.canvas.draw()
                    fig.canvas.flush_events()
                    plt.pause(0.1)
                
            else:
                already_reduced += 1
            
        except AssertionError:
            __watchdog__write_to_blacklist(wscan,fe,'WireScanner','NotWireScan',0)
        except TimeoutError as e:
            __watchdog__write_to_blacklist(wscan,fe,'WireScanner',e,timeout)
        except Exception as e:
            __watchdog__write_to_blacklist(wscan,fe,'WireScanner',e,0)
        finally:
            temp.alarm(0)

    if already_reduced==total:
        info('No new WireScans to reduce.')


__watchdog_figure = None
def __watchdog_run(fe='LFA', mode='loop', doPlot=True, timeout='default'):
    assert str.upper(fe) in ['LFA', 'HFA'], 'fe must be LFA or HFA'
    fe = str.upper(fe)
    global __watchdog_figure
    timeout = 60 if fe=='LFA' else 180  # 3x KIDs in HFA

    if doPlot:
        if __watchdog_figure==None:
            __watchdog_figure = plt.figure(figsize=(10,10))
            __watchdog_figure.show()

    timeout = int(timeout)

    if mode=='single':
        info('Watchdog in single-run mode for %s.'%fe)
        calscan_list, fsweep_list, wscan_list = __watchdog__findCalScans(fe=fe)
        __reduce_Fsweep_list(fsweep_list, fe=fe, timeout=timeout)
        __reduce_WireScan_list(wscan_list, fe=fe, timeout=timeout, fig=__watchdog_figure)

    if mode=='loop':
        info('Watchdog in loop mode for %s, CTRL+C to exit!'%fe)
        calscan_list, fsweep_list, wscan_list = __watchdog__findCalScans(fe=fe)
        info('Processing Fsweeps...')
        __reduce_Fsweep_list(fsweep_list, fe=fe, timeout=timeout)
        info('Processing WireScans...')
        __reduce_WireScan_list(wscan_list, fe=fe, timeout=timeout, fig=__watchdog_figure)
        info('Checking back in 15 seconds...')
        time_0 = time.time()
        print('')
        while True:
            time_current = time.time()
            if time_current - time_0 > 15:  # at least 15 seconds
                __watchdog_run(fe=fe, mode='loop', doPlot=doPlot, timeout=timeout)