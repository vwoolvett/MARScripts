# Module to find and reduce calibrations upon new scan arrival for watchdog
# at APEX

def __watchdog__findCalScans(fe='LFA', verbose=False):
    assert fe=='LFA' or fe=='HFA', 'fe must be LFA or HFA'
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

                    # and is finished, consider it
                    if '-999' not in duration:
                        message += 'SCAN CONSIDERED'.ljust(15) + ' | '
                        message += comment
                        calscanlist.append(scan_int)
                        if 'calibrate(' in command:
                            fsweep_list.append(scan_int)
                        if 'beamscan(' in command:
                            wscan_list.append(scan_int)

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



def __reduce_Fsweep_list(fsweep_list, fe='LFA', overwrite=False):

    # track reductions
    total = len(wscan_list)
    already_reduced = 0

    # Iterate over fsweep list
    for fsweep in fsweep_list:
        # check if there is a reduction
        globlist = glob('CalFiles/fsweep_%s_%s.csv'%(fsweep, fe))
        if len(globlist)==0:
            if scanIsFsweep(fsweep, fe=fe):
                reduceFsweep(fsweep, fe=fe)
        if len(globlist)==1 and overwrite==True:
            if scanIsFsweep(fsweep, fe=fe):
                reduceFsweep(fsweep, fe=fe, overwrite=True)
        else:
            already_reduced += 1

    if len(already_reduced)==total:
            info('No new WireScans to reduce...')



def __reduce_WireScan_list(wscan_list, fe='LFA', doPlot=True, overwrite=False):
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
                reduceWireScan(wscan, fe=fe, doPlot=True, fig=figure)

        if len(globlist)==1 and overwrite==True:
            if scanIsBeamscan(wscan, fe=fe):
                if doPlot==True and figure_initiated==False:
                    figure = plt.figure(figsize=(10,10))
                    figure_initiated=True
                reduceWireScan(wscan, fe=fe, overwrite=True, doPlot=True, fig=figure)
        else:
            already_reduced += 1

    if len(already_reduced)==total:
        info('No new WireScans to reduce...')


def __watchdog_run(fe='LFA', mode='loop'):
    if mode=='single':
        info('Watchdog in single-run mode.')
        calscan_list, fsweep_list, wscan_list = __watchdog__findCalScans(fe=fe)
        __reduce_Fsweep_list(fsweep_list, fe=fe)
        __reduce_WireScan_list(wscan_list, fe=fe)

    if mode=='loop':
        info('Watchdog in loop mode, CTRL+F to exit!')
        time_0 = time.time()
        calscan_list, fsweep_list, wscan_list = __watchdog__findCalScans(fe=fe)
        __reduce_Fsweep_list(fsweep_list, fe=fe)
        __reduce_WireScan_list(wscan_list, fe=fe)
        while True:
            time_current = time.time()
            if time_current - time_0 > 10:  # at least 10 seconds
                __watchdog_run(fe=fe, mode='loop')