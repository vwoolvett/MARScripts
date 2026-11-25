# Module to find and reduce calibrations upon new scan arrival for watchdog
# at APEX

# Initialize figure variable
__watchdog_figure = None


def __watchdog_timeout_handler(signum, frame):
    raise TimeoutError('Reduction timed out')



def __watchdog_checkblacklist():
    # Create watchdog folder for blacklist and others
    if os.path.exists('./.watchdog')==False:
        os.makedirs('./.watchdog')

    # Create csv header
    if not os.path.exists('./.watchdog/watchdog_blacklist.csv'):
        with open('./.watchdog/watchdog_blacklist.csv', 'w') as f:
            f.write('# List of calibration scans that have failed\n')
            f.write('# in their reduction to avoid reducing them again\n')
            f.write('# scan , fe , type           , timeout, error                    , errordescription\n')
            f.close()



def __watchdog_write_to_blacklist(scan, fe, caltype, error, timeout=0):
    fe = str.upper(fe)
    # Create Blacklist if non-existent
    __watchdog_checkblacklist()

    # Array-specific blacklist
    blacklist = __watchdog_readblacklist(fe=fe)

    if scan not in blacklist:
        with open('./.watchdog/watchdog_blacklist.csv', 'a') as f:
            if type(error)==str:
                errorstr = error
                errordesc = 'No description, Exception is not raised but was catched'
            else:
                errorstr = type(error).__name__
                errordesc = error.message
                if errordesc == '':
                    errordesc = str(error)
                    if len(errordesc)!=0:
                        if errordesc[0]=="'":
                            errordesc = errordesc[1:]
                        if errordesc[-1]=="'":
                            errordesc = errordesc[:-1]
            
            caltype = caltype.ljust(15)
            timeout = str(int(timeout)).ljust(7)
            errorstr = errorstr.ljust(25)

            f.write("  %s, %s, %s, %s, %s, %s\n"%(scan, fe, caltype, timeout, errorstr, errordesc))
            f.close()



def __watchdog_readblacklist(fe='LFA', returnDict=False):
    fe = str.upper(fe)

    # Create Blacklist if non-existent
    __watchdog_checkblacklist()

    # Read blacklist
    with open('./.watchdog/watchdog_blacklist.csv', 'r') as f:
        blacklistdict = {}
        lines = f.readlines()
        for line in lines:
            if line[0] != '#' and str.upper(line.split(',')[1].split()[0])==fe:
                calscan = int(line.split(',')[0])
                calscanfe = line.split(',')[1]
                calscantype = line.split(',')[2]
                timeout = line.split(',')[3]
                error = line.split(',')[4]
                errordesc = line.split(',')[5]  # might have commas but is last CSV value

                blacklistdict[calscan] = {'scan':calscan,
                                          'fe':calscanfe,
                                          'type':calscantype,
                                          'timeout':timeout,
                                          'error':error,
                                          'errordescription':errordesc
                                          }

        if returnDict:
            return blacklistdict

        else:
            blacklist = [blacklisted for blacklisted in blacklistdict.keys()]
            blacklist.sort()
            return blacklist



def __watchdog_figure_update(fig=__watchdog_figure, pausetime=0.1):
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.pause(pausetime)



def __watchdog_figure_flash(fig, textstr='', fontsize=20, flashcolor='yellow', flashtime=0.5, backtoblack=True):
    assert type(textstr)==str
    bbox=dict(
        facecolor='black',
        edgecolor='black',
        linewidth=1,
        boxstyle='square,pad=0.3'
    )
    frames = np.arange(5)
    for frame in frames:
        if frame % 2 == 0:
            color = flashcolor
        else:
            color = 'black'
        fig.clf()
        fig.set_facecolor(color)
        fig.text(0.5, 0.5, textstr, bbox=bbox, fontdict={'fontsize':fontsize}, horizontalalignment='center', color='white')
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(flashtime)

    if backtoblack:
        fig.set_facecolor('black')



def __watchdog_figure_text(fig, textstr='', fontsize=20, facecolor='green', backtoblack=True, duration=10):
    assert type(textstr)==str
    bbox=dict(
        facecolor='black',
        edgecolor='black',
        linewidth=1,
        boxstyle='square,pad=0.3'
    )
    fig.clf()
    fig.set_facecolor(facecolor)
    fig.text(0.5, 0.5, textstr, bbox=bbox, fontdict={'fontsize':fontsize}, horizontalalignment='center', color='white')
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.pause(duration)
    if backtoblack:
        fig.set_facecolor('black')



def __watchdog_findCalScans(fe='LFA', verbose=False):
    fe = str.upper(fe)
    assert fe=='LFA' or fe=='HFA', 'fe must be LFA or HFA'

    # Create Blacklist if non-existent
    __watchdog_checkblacklist()

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



def __reduce_Fsweep_list(fsweep_list, fe='LFA', timeout=60, overwrite=False, doPlot=True, fig=None):
    import signal as temp
    fe = str.upper(fe)

    # info print
    info('Processing Fsweeps (%s)'%fe)

    # track reductions
    total = len(fsweep_list)
    already_reduced = 0

    # Iterate over fsweep list
    for fsweep in fsweep_list:
        temp.signal(temp.SIGALRM, __watchdog_timeout_handler)
        try:
            temp.alarm(timeout)
            # check if there is a reduction
            #globlist = glob('CalFiles/fsweep_%s_%s.csv'%(fsweep, fe))
            #if len(globlist)==0:
            reductionexists = os.path.exists('CalFiles/fsweep_%s_%s.csv'%(fsweep, fe))
            if reductionexists==False:
                if doPlot and fig is not None:
                    mytext = 'Now processing Fsweep %s (%s)!'%(fsweep, fe)
                    __watchdog_figure_flash(fig=fig, textstr=mytext)
                if scanIsFsweep(fsweep, fe=fe)==False:
                    raise AssertionError
                reduceFsweep(fsweep, fe=fe, doReturn=True)
                if doPlot and fig is not None:
                    mytext  = 'Reduction of Fsweep %s (%s) successful!\n'%(fsweep, fe)
                    mytext += '\n'
                    mytext += 'Statistics (to be implemented):'
                    __watchdog_figure_flash(fig=fig, textstr=mytext, flashcolor='green')
                    __watchdog_figure_text(fig=fig, textstr=mytext, facecolor='green')
            #if len(globlist)==1 and overwrite==True:
            elif reductionexists==True and overwrite==True:
                if doPlot and fig is not None:
                    mytext = 'Now processing Fsweep %s (%s)!'%(fsweep, fe)
                    __watchdog_figure_flash(fig=fig, textstr=mytext)
                if scanIsFsweep(fsweep, fe=fe)==False:
                    raise AssertionError
                reduceFsweep(fsweep, fe=fe, doReturn=True, overwrite=True)
                if doPlot and fig is not None:
                    mytext  = 'Reduction of Fsweep %s (%s) successful!\n'%(fsweep, fe)
                    mytext += '\n'
                    mytext += 'Statistics (to be implemented):'
                    __watchdog_figure_flash(fig=fig, textstr=mytext, flashcolor='green')
                    __watchdog_figure_text(fig=fig, textstr=mytext, facecolor='green')


            else:
                already_reduced += 1

        except AssertionError:
            temp.alarm(0)
            __watchdog_write_to_blacklist(fsweep,fe,'Fsweep','NotFsweep',0)
            if doPlot and fig is not None:
                textstr = 'Reduction of Fsweep %i (%s) failed:\nNotFsweep'%(fsweep, fe)
                __watchdog_figure_flash(fig=fig, textstr=textstr, flashcolor='red')
        except TimeoutError as e:
            temp.alarm(0)
            __watchdog_write_to_blacklist(fsweep,fe,'Fsweep',e,timeout)
            if doPlot and fig is not None:
                textstr = 'Reduction of Fsweep %i (%s) failed:\n%s'%(fsweep, fe, type(e).__name__)
                __watchdog_figure_flash(fig=fig, textstr=textstr, flashcolor='red')
        except Exception as e:
            temp.alarm(0)
            __watchdog_write_to_blacklist(fsweep,fe,'Fsweep',e,0)
            if doPlot and fig is not None:
                textstr = 'Reduction of Fsweep %i (%s) failed:\n%s'%(fsweep, fe, type(e).__name__)
                __watchdog_figure_flash(fig=fig, textstr=textstr, flashcolor='red')
        finally:
            temp.alarm(0)

    if already_reduced==total:
        info('No new Fsweeps to reduce.')



def __reduce_WireScan_list(wscan_list, fe='LFA', timeout=60, overwrite=False, doPlot=True, fig=None):
    import signal as temp
    fe = str.upper(fe)

    # info print
    info('Processing WireScans (%s)'%fe)

    # track reductions
    total = len(wscan_list)
    already_reduced = 0

    # Iterate over wscan list
    for wscan in wscan_list:
        temp.signal(temp.SIGALRM, __watchdog_timeout_handler)
        try:
            temp.alarm(timeout)
            # check if there is a reduction
            #globlist = glob('CalFiles/beam_scanner_%s_%s.csv'%(wscan, fe))
            #if len(globlist)==0:
            reductionexists = os.path.exists('CalFiles/beam_scanner_%s_%s.csv'%(wscan, fe))
            if reductionexists==False:
                if doPlot and fig is not None:
                    mytext = 'Now processing WireScan %s (%s)!'%(wscan, fe)
                    __watchdog_figure_flash(fig=fig, textstr=mytext)
                if scanIsBeamscan(wscan, fe=fe)==False:
                    raise AssertionError
                reduceWireScan(wscan, fe=fe, doPlot=False, fig=fig)
                if doPlot and fig is not None:
                    mytext = 'Reduction of WireScan %s (%s) successful!'%(wscan, fe)
                    __watchdog_figure_flash(fig=fig, textstr=mytext, flashcolor='green')
                fig.set_facecolor('black')
                reduceWireScan(wscan, fe=fe, doPlot=doPlot, fig=fig)
                if doPlot and fig is not None:
                    __watchdog_figure_update(fig=fig)
            #if len(globlist)==1 and overwrite==True:
            if reductionexists==True and overwrite==True:
                if doPlot and fig is not None:
                    mytext = 'Now processing WireScan %s (%s)!'%(wscan, fe)
                    __watchdog_figure_flash(fig=fig, textstr=mytext)
                if scanIsBeamscan(wscan, fe=fe)==False:
                    raise AssertionError
                reduceWireScan(wscan, fe=fe, overwrite=True, doPlot=False, fig=fig)
                if doPlot and fig is not None:
                    mytext = 'Reduction of WireScan %s (%s) successful!'%(wscan, fe)
                    __watchdog_figure_flash(fig=fig, textstr=mytext, flashcolor='green')
                fig.set_facecolor('black')
                reduceWireScan(wscan, fe=fe, doPlot=doPlot, fig=fig)
                if doPlot and fig is not None:
                    __watchdog_figure_update(fig=fig)
                
            else:
                already_reduced += 1
            
        except AssertionError:
            temp.alarm(0)
            __watchdog_write_to_blacklist(wscan,fe,'WireScanner','NotWireScan',0)
            if doPlot and fig is not None:
                textstr = 'Reduction of WireScanner %i (%s) failed:\nNotWireScan'%(wscan, fe)
                __watchdog_figure_flash(fig=fig, textstr=textstr, flashcolor='red')
        except TimeoutError as e:
            temp.alarm(0)
            __watchdog_write_to_blacklist(wscan,fe,'WireScanner',e,timeout)
            if doPlot and fig is not None:
                textstr = 'Reduction of WireScanner %i (%s) failed:\n%s'%(wscan, fe, type(e).__name__)
                __watchdog_figure_flash(fig=fig, textstr=textstr, flashcolor='red')
        except Exception as e:
            temp.alarm(0)
            __watchdog_write_to_blacklist(wscan,fe,'WireScanner',e,0)
            if doPlot and fig is not None:
                textstr = 'Reduction of WireScanner %i (%s) failed:\n%s'%(wscan, fe, type(e).__name__)
                __watchdog_figure_flash(fig=fig, textstr=textstr, flashcolor='red')
        finally:
            temp.alarm(0)

    if already_reduced==total:
        info('No new WireScans to reduce.')



def __watchdog_run(fe='LFA', mode='loop', doPlot=True, timeout='default', firstloop=True):
    assert str.upper(fe) in ['LFA', 'HFA'], 'fe must be LFA or HFA'

    # Check current directory is a "Reduction" directory, otherwise ask if continue
    if firstloop==True:
        cwd = os.getcwd()
        if str.upper(cwd.split('/')[-1])!='REDUCTION':
            message  = '-----------------------------------------------------\n'
            message += '-----------            Warning!           -----------\n'
            message += '-----------------------------------------------------\n'
            message += 'Current directory:\n'
            message += '%s\n'%cwd
            message += 'is not a Reduction folder. You may miss the required\n'
            message += 'blacklist file for Watchdog. It is recommended to run\n'
            message += 'Watchdog in the Reduction folder of the project.\n'
            message += '\nContinue anyways? (y/N)'

            user_input = raw_input(message)
            user_input = str(user_input)

            if str.upper(user_input) in ['', 'N', 'NO']:
                info('Quitting Watchdog')
                return
            else:
                firstloop=False
                pass


    # check blacklist
    __watchdog_checkblacklist()

    global __watchdog_figure
    fe = str.upper(fe)

    if timeout=='default':
        timeout = 60 if fe=='LFA' else 180  # 3x KIDs in HFA
    else:
        timeout = int(timeout)
        assert timeout < 600, 'Timeout is larger than 10 minutes!'

    # Initialize figure and assign to global variable
    if doPlot:
        # if unassigned
        if __watchdog_figure==None:
            __watchdog_figure = plt.figure(figsize=(10,10))
            __watchdog_figure.show()

        else:
            # if Window closed, reopen
            if __watchdog_figure.canvas.manager.window==None:
                __watchdog_figure = None
                __watchdog_figure = plt.figure(figsize=(10,10))
                __watchdog_figure.show()
            else:
                # keep using same figure with opened window
                pass


    if mode=='single':
        info('Watchdog in single-run mode for %s.'%fe)
        calscan_list, fsweep_list, wscan_list = __watchdog_findCalScans(fe=fe)
        __reduce_Fsweep_list(fsweep_list, fe=fe, timeout=timeout, doPlot=doPlot, fig=__watchdog_figure)
        __reduce_WireScan_list(wscan_list, fe=fe, timeout=timeout, doPlot=doPlot, fig=__watchdog_figure)

    if mode=='loop':
        info('Watchdog in loop mode for %s, CTRL+C to exit!'%fe)
        calscan_list, fsweep_list, wscan_list = __watchdog_findCalScans(fe=fe)
        __reduce_Fsweep_list(fsweep_list, fe=fe, timeout=timeout, doPlot=doPlot, fig=__watchdog_figure)
        __reduce_WireScan_list(wscan_list, fe=fe, timeout=timeout, doPlot=doPlot, fig=__watchdog_figure)
        info('Checking back in 5 seconds...')
        time_0 = time.time()
        print('')
        while True:
            time_current = time.time()
            if time_current - time_0 > 5:  # at least 5 seconds
                __watchdog_run(fe=fe, mode='loop', doPlot=doPlot, timeout=timeout, firstloop=firstloop)