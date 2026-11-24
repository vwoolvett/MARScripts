def findspikes_Jy_speed(windowtime=5., sig=5.5, expspikefree=75, crosstones=0., ignoreblinds=True, full_output=False, doplots=False, debug=False):
    '''
    ** VERSION 1 - 10.08.2026 **

    Same as findspikes_IQBT but intended to be used at Jy scale

    @param windowtime:      size of windows in seconds
    @type windowtime:       float
    @param sig:             how many STDs away from mean tone speed to flag as spikes
    @type sig:              float
    @param expspikefree:    percentage of the timelines that is expected to be spike-free to determine "usual tone behavior"
    @type expspikefree:     float
    @param crosstones:      for any time window, what percentage of tones must be spiked to consider the window as really spiked. Default 25% is one chip.
    @type crosstones:       float
    @param ignoreblinds:    whether to ignore (True) or consider (False) blindtones. Useful if blindtones are spiked too.
    @type ignoreblinds:     bool
    @param full_output:     if False, only mask for data.Data is returned. Otherwise, the mask (0) and the 1-D array (1) of
                            fraction of timeline lost to spikes for all used tones (dependent of ignoreblinds!) are returned
    @type full_output:      bool
    @param doplots:         whether to do an example plot of the spike detection
    @type doplots:          bool
    @param debug:           Whether to plot debug info for spike cross-tone check
    @type debug:            bool
    '''
    
    info('Beginning spike/instability detection (Jy/b)...')
    # retrieve flagged data mask
    isflagged = data.FlagHandler._aFlags > 0
    
    # data array
    F = copy.deepcopy(data.Data)  # is now flux in Jy/b
    
    # set flagged data to NaN in this temporary variable to avoid flagging source model as spike
    F[isflagged] = np.NaN
    
    # time array
    scannum = data.ScanParam.ScanNum
    info('Retrieving time axis data...')
    time = (data.ScanParam.MJD - data.ScanParam.MJD[0]) * 24 * 3600
    totaltime = time[-1] - time[0]
    dt = np.nanmedian(np.diff(time))
    
    # set ntones and other stuff
    ntones = F.shape[1]
    febe = data.BolometerArray.FeBe
    if febe == 'AMKID870-AMKID870BE':
        fe = 'LFA'
        nkids = 880 * 4
    elif febe == 'AMKID350-AMKID350BE':
        fe = 'HFA'
        nkids = 800 * 5 * 4
    else:
        warn('This is not AMKID data!')
        return
    
    if ignoreblinds:
        F = F[:, :nkids]
        nused = nkids
    else:
        nused = ntones

    # Ensure at least 100 windows
    nwindows = totaltime / windowtime
    if nwindows<50.0:
        nwindows = 50
        warn('Adjusting windowtime to %1.2f seconds to ensure 50 windows...'%(totaltime/nwindows))
        windowtime = totaltime / nwindows

    # Derivatives
    dFdt = np.diff(F, axis=0) / dt  # Jy/b / s
    speeds = np.abs(dFdt) # Jy/b / s
    auxtime = time[:-1] + dt/2

    # define windows
    windows_tstart = np.arange(0, totaltime, windowtime)
    windows_time = windows_tstart + windowtime/2
    info('Windowing: %.2f s -> %i windows of %.1f s'%(totaltime, len(windows_tstart), windowtime))

    # RMS of windows
    # info('Computing spike detection metrics...')
    windows_speed_mean = []
    windows_speed_std = []
    windows_speed_max = []

    for t_start in windows_tstart:
        t_end = t_start + windowtime
        valid = (auxtime >= t_start) & (auxtime < t_end)
        speeds_window = speeds[valid, :nused]

        if speeds_window.shape[0] == 0:
            windows_speed_mean.append(np.full(nused, np.nan))
            windows_speed_std.append(np.full(nused, np.nan))
            windows_speed_max.append(np.full(nused, np.nan))
        else:
            windows_speed_mean.append(np.nanmean(speeds_window, axis=0))
            windows_speed_std.append(np.nanstd(speeds_window, axis=0))
            windows_speed_max.append(np.nanmax(speeds_window, axis=0))

    # convert to arrays
    windows_speed_mean = np.array(windows_speed_mean)
    windows_speed_std = np.array(windows_speed_std)
    windows_speed_max = np.array(windows_speed_max)

    # initialize window flagging array
    windowflag = np.zeros_like(windows_speed_std, dtype=bool)
    if doplots:
        # initialize used values list
        # has as many values as tones
        tone_floor_speedMEANs = []
        tone_floor_speedSTDs = []
        tone_thresholds_speed = []

    info('Detecting spikes and instabilities...')
    # fill windows flagging array - has True for each window timestamp and tone if window has rms above a carefully chosen threshold
    for toneidx in range(nused):
        # only for non-nan tones!
        if not(np.all(np.isnan(data.Data[:, toneidx]))):
            # Tag windows to determine Flux-speed mean and STD floors excluding the 100-expspikefree% most noisy windows
            # (hard-coded max spike proportion of timelines)
            thresh_tonebased = np.nanpercentile(windows_speed_std[:, toneidx], expspikefree)
            valid = windows_speed_std[:, toneidx] <= thresh_tonebased

            # Compute the spike-free mean of speed means for tone (MEAN floor)
            floor_speedMEAN = np.nanmean(windows_speed_mean[valid, toneidx])

            # compute the spike-free mean of speed STDs for tone (STD floor)
            floor_speedSTD = np.nanmean(windows_speed_std[valid, toneidx])

            # compute the mean + sig STD threshold
            threshold_speed = floor_speedMEAN + sig * floor_speedSTD
            
            if doplots:
                tone_floor_speedMEANs.append(floor_speedMEAN)
                tone_floor_speedSTDs.append(floor_speedSTD)
                tone_thresholds_speed.append(threshold_speed)
                
            # Tone likely has spike if ANY speed sample in window surpasses threshold
            # or equivalently if the maximum speed does (less comparisons)
            windowflag[:, toneidx] = windows_speed_max[:, toneidx] >= threshold_speed
            
        else:
            if doplots:
                # Tone is only NaNs, then fill masks to leave un-edited (no spikes)
                tone_floor_speedMEANs.append(np.nan)
                tone_floor_speedSTDs.append(np.nan)
                tone_thresholds_speed.append(np.nan)
            windowflag[:, toneidx] = False
    
    if crosstones > 0.0:        
        # copy original windowflag
        windowflag2 = copy.deepcopy(windowflag)
        
        if ignoreblinds==True:
            N_TOT = len(data.BolometerArray.checkChanList([]))
        else:
            N_TOT = len(data.BolometerArray.checkChanList([])) +\
                    len([tone for tone in data.tone_dict.keys() if data.tone_dict[tone]['type']=='REF' and data.tone_dict[tone]['is_good']==True])
        N_THRESH = int(crosstones/100. * float(N_TOT))
        info('Cross-checking tones...')
        for windowidx in range(len(windows_tstart)):
            flaggedtones_thiswindow = np.sum(windowflag[windowidx, :])
            if debug:
                print("DEBUG: WIN=%i | time=%.1f s | N_FLAGGED=%i | N_THRESH=%i | N_TOT=%i"%(windowidx+1, windows_time[windowidx], flaggedtones_thiswindow, N_THRESH, N_TOT))

            # if this window fulfills the criterium
            should_this_window_be_flagged = (flaggedtones_thiswindow >= N_THRESH)

            # otherwise, it's not a spike and we should not flag
            if not should_this_window_be_flagged:
                # Not real spike
                windowflag2[windowidx, :] = False

        # re-create windowflag and save memory
        windowflag = windowflag2
        windowflag2 = None
    
    # inizialize data flagging array
    # if blindtones are ignored in process, flag is false for all blindtones so they will not be affected
    flagmask = np.zeros_like(data.Data, dtype=bool)
    
    # fill data mask array
    for toneidx in range(nused):
        # if tone is spiked more than 30% of the time, simply flag it
        #if float(np.sum(windowflag[:, toneidx]))/float(len(windowflag[:, toneidx])) > 0.30:
        #    windowflag[:, toneidx] = True
        flagged_windows = np.where(windowflag[:, toneidx])[0]
        for winidx in flagged_windows:
            t_start = windows_tstart[winidx]
            t_end = t_start + windowtime
            flagmask[:, toneidx] |= (time >= t_start) & (time < t_end)

    if doplots:
        # convert to arrays
        tone_floor_speedMEANs = np.array(tone_floor_speedMEANs)
        tone_floor_speedSTDs = np.array(tone_floor_speedSTDs)
        tone_thresholds_speed = np.array(tone_thresholds_speed)

    # final info
    usedkids_idx = data.BolometerArray.checkChanList([]) - 1
    spikedfraction_usedkids = np.sum(windowflag[:, usedkids_idx], axis=0) / float(np.shape(windowflag)[0])
    spikedfraction_avg = np.nanmean(spikedfraction_usedkids)
    if np.all(spikedfraction_usedkids==0):
        info('Scan is free of spikes!')
    else:
        info('%.3f'%(spikedfraction_avg*100) + '% (avg.) of the KID timelines is lost')

    if doplots:
        info('Initializing de-spiking process interactive plot...')

        # =====================
        # Interactive plot
        # =====================

        # Message to interact with user
        msg = '\n------------------------------------------------\n'
        msg += 'Plot next Tone:                           <Enter>\n'
        msg += 'Plot previous Tone:                     - <Enter>\n'
        msg += 'Go to Tone Number:               <number> <Enter>\n'
        msg += 'Exit script execution:                  q <Enter>\n'
        msg += '-------------------------------------------------\n'
        msg += 'Enter choice: '

        # Create figure
        fig, ax = plt.subplots(2, 2, figsize=(9, 6))

        usedtones = data.BolometerArray.checkChanList([])
        testtone = usedtones[0]

        while True:
            # define data array tone index
            toneidx = testtone - 1

            # extract Flux speed of tone
            speed_testtone = speeds[:, toneidx]

            # create mask applied to speed array
            auxflagmask = np.zeros_like(speed_testtone, dtype=bool)
            for winidx in range(len(windows_time)):
                if windowflag[winidx, toneidx]:
                    # flag all data points that fall in this window
                    t_start = windows_time[winidx] - windowtime/2
                    t_end = windows_time[winidx] + windowtime/2
                    auxflagmask[(auxtime >= t_start) & (auxtime < t_end)] = True

            # now we make nan all values in the original data that fall in a spiked window
            despiked_F_testtone = F[:, toneidx].copy()
            despiked_F_testtone[flagmask[:, toneidx]] = np.nan

            # and for the speed
            despiked_speed_testtone = speed_testtone.copy()
            despiked_speed_testtone[auxflagmask] = np.nan

            # and the window speed mean array
            despiked_windows_speed_mean_testtone = windows_speed_mean[:, toneidx].copy()
            despiked_windows_speed_mean_testtone[windowflag[:, toneidx]] = np.nan

            # and the window speed std array
            despiked_windows_speed_std_testtone = windows_speed_std[:, toneidx].copy()
            despiked_windows_speed_std_testtone[windowflag[:, toneidx]] = np.nan

            # Flux
            ax[0, 0].plot(time, F[:, toneidx], label='Spiked flux')
            ax[0, 0].plot(time, despiked_F_testtone, label='Spike-masked flux', zorder=1e9)
            ax[0, 0].set_xlabel('Time (s)')
            ax[0, 0].set_ylabel('Flux (Jy/b)')
            ax[0, 0].legend(framealpha=1)

            # Flux speed (normalized)
            ax[0, 1].plot(auxtime, speed_testtone)
            ax[0, 1].plot(auxtime, despiked_speed_testtone)
            ax[0, 1].axhline(tone_floor_speedMEANs[toneidx], c='green', lw=2, label='Average spike-free tone speed')
            ax[0, 1].axhline(tone_thresholds_speed[toneidx], c='red', lw=2, label='Mean + %s Sigma'%sig)
            try:
                ax[0, 1].set_ylim(0, tone_thresholds_speed[toneidx]*3.)
            except:
                pass
            ax[0, 1].set_xlabel('Time (s)')
            ax[0, 1].set_ylabel('Flux speed (Jy/b/s)')
            ax[0, 1].legend(framealpha=1)

            # flux speed window mean
            ax[1, 0].plot(windows_time, windows_speed_mean[:, testtone-1])
            ax[1, 0].plot(windows_time, despiked_windows_speed_mean_testtone)
            ax[1, 0].axhline(0)
            ax[1, 0].axhline(tone_floor_speedMEANs[toneidx], c='green', lw=2, label='Average spike-free tone speed')
            ax[1, 0].legend(framealpha=1)
            ax[1, 0].set_xlabel('Time (s)')
            ax[1, 0].set_ylabel('Windowed flux-speed mean (Jy/b/s)')

            # flux speed window std
            ax[1, 1].plot(windows_time, windows_speed_std[:, testtone-1])
            ax[1, 1].plot(windows_time, despiked_windows_speed_std_testtone)
            ax[1, 1].axhline(0)
            ax[1, 1].axhline(tone_floor_speedSTDs[toneidx], c='red', lw=2, label='Average spike-free tone speed STD')
            ax[1, 1].legend(framealpha=1)
            ax[1, 1].set_xlabel('Time (s)')
            ax[1, 1].set_ylabel('Windowed flux-speed STD (Jy/b/s)')

            fig.suptitle('Spike/instability detection for tone %i in scan %i'%(testtone, scannum))

            plt.show()

            userInput = raw_input(msg)
            userInput=str(userInput)
            if userInput == '':
                current = np.where(np.array(usedtones)==testtone)[0][0]
                testtone = usedtones[current+1]
            elif userInput == '-':
                current = np.where(np.array(usedtones)==testtone)[0][0]
                testtone = usedtones[current-1]
            elif userInput == 'q':
                plt.close(fig)
                break
            else:    
                try:
                    oldtesttone = testtone 
                    testtone = int(userInput)
                except:
                    print('Tone number must be a number.')
                    continue

            if testtone not in usedtones:
                print('Tone %s not in used tones.'%testtone)
                testtone = oldtesttone

            ax[0, 0].cla()
            ax[0, 1].cla()
            ax[1, 0].cla()
            ax[1, 1].cla()

    if full_output:
        return flagmask, spikedfraction_alltones
    else:
        return flagmask



def findspikes_Jy_flux(windowtime=2., sig=5., expspikefree=75, crosstones=40., ignoreblinds=True, full_output=False, doplots=False, debug=False):
    '''
    ** VERSION 1 - 10.08.2026 **

    Same as findspikes_IQBT but intended to be used at Jy scale

    @param windowtime:      size of windows in seconds
    @type windowtime:       float
    @param sig:             how many STDs away from mean tone flux to flag as spikes
    @type sig:              float
    @param expspikefree:    percentage of the timelines that is expected to be spike-free to determine "usual tone behavior"
    @type expspikefree:     float
    @param crosstones:      for any time window, what percentage of tones must be spiked to consider the window as really spiked. Default 25% is one chip.
    @type crosstones:       float
    @param ignoreblinds:    whether to ignore (True) or consider (False) blindtones. Useful if blindtones are spiked too.
    @type ignoreblinds:     bool
    @param full_output:     if False, only mask for data.Data is returned. Otherwise, the mask (0) and the 1-D array (1) of
                            fraction of timeline lost to spikes for all used tones (dependent of ignoreblinds!) are returned
    @type full_output:      bool
    @param doplots:         whether to do an example plot of the spike detection
    @type doplots:          bool
    @param debug:           Whether to plot debug info for spike cross-tone check
    @type debug:            bool
    '''
    
    info('Beginning spike/instability detection (Jy/b)...')
    # retrieve flagged data mask
    isflagged = data.FlagHandler._aFlags > 0
    
    # data array
    F = copy.deepcopy(data.Data)  # is now flux in Jy/b
    
    # set flagged data to NaN in this temporary variable to avoid flagging source model as spike
    F[isflagged] = np.NaN
    
    # time array
    scannum = data.ScanParam.ScanNum
    info('Retrieving time axis data...')
    time = (data.ScanParam.MJD - data.ScanParam.MJD[0]) * 24 * 3600
    totaltime = time[-1] - time[0]
    dt = np.nanmedian(np.diff(time))
    
    # set ntones and other stuff
    ntones = F.shape[1]
    febe = data.BolometerArray.FeBe
    if febe == 'AMKID870-AMKID870BE':
        fe = 'LFA'
        nkids = 880 * 4
    elif febe == 'AMKID350-AMKID350BE':
        fe = 'HFA'
        nkids = 800 * 5 * 4
    else:
        warn('This is not AMKID data!')
        return
    
    if ignoreblinds:
        F = F[:, :nkids]
        nused = nkids
    else:
        nused = ntones

    # Ensure at least 100 windows
    nwindows = totaltime / windowtime
    if nwindows<50.0:
        nwindows = 50
        warn('Adjusting windowtime to %1.2f seconds to ensure 50 windows...'%(totaltime/nwindows))
        windowtime = totaltime / nwindows

    # Others
    absF = np.abs(F) # Jy/b / s

    # define windows
    windows_tstart = np.arange(0, totaltime, windowtime)
    windows_time = windows_tstart + windowtime/2
    info('Windowing: %.2f s -> %i windows of %.1f s'%(totaltime, len(windows_tstart), windowtime))

    # RMS of windows
    # info('Computing spike detection metrics...')
    windows_flux_mean = []
    windows_flux_std = []
    windows_absflux_max = []

    for t_start in windows_tstart:
        t_end = t_start + windowtime
        valid = (time >= t_start) & (time < t_end)
        flux_window = F[valid, :nused]
        absflux_window = absF[valid, :nused]

        if flux_window.shape[0] == 0:
            windows_flux_mean.append(np.full(nused, np.nan))
            windows_flux_std.append(np.full(nused, np.nan))
            windows_absflux_max.append(np.full(nused, np.nan))
        else:
            windows_flux_mean.append(np.nanmean(flux_window, axis=0))
            windows_flux_std.append(np.nanstd(flux_window, axis=0))
            windows_absflux_max.append(np.nanmax(absflux_window, axis=0))

    # convert to arrays
    windows_flux_mean = np.array(windows_flux_mean)
    windows_flux_std = np.array(windows_flux_std)
    windows_absflux_max = np.array(windows_absflux_max)

    # initialize window flagging array
    windowflag = np.zeros_like(windows_flux_std, dtype=bool)
    if doplots:
        # initialize used values list
        # has as many values as tones
        tone_floor_fluxMEANs = []
        tone_floor_fluxSTDs = []
        tone_thresholds_absflux = []
        tone_thresholds_absmean = []

    info('Detecting spikes and instabilities...')
    # fill windows flagging array - has True for each window timestamp and tone if window has rms above a carefully chosen threshold
    for toneidx in range(nused):
        # only for non-nan tones!
        if not(np.all(np.isnan(data.Data[:, toneidx]))):
            # Tag windows to determine flux mean and STD floors excluding the 100-expspikefree% most noisy windows
            # (hard-coded max spike proportion of timelines)
            thresh_tonebased = np.nanpercentile(windows_flux_std[:, toneidx], expspikefree)
            valid = windows_flux_std[:, toneidx] <= thresh_tonebased

            # Compute the spike-free mean flux for tone (MEAN floor)
            #floor_fluxMEAN = np.nanmean(windows_flux_mean[valid, toneidx]) # should be zero
            floor_fluxMEAN = 0.

            # compute the spike-free std flux for tone (STD floor)
            floor_fluxSTD = np.nanmean(windows_flux_std[valid, toneidx])
            
            # compute max theoretical mean
            ptsperwin = windowtime / dt
            sigma_mean = floor_fluxSTD / (ptsperwin)**0.5
            threshold_absmean = sig * sigma_mean

            # compute the mean + sig STD threshold
            threshold_absflux = floor_fluxMEAN + sig * floor_fluxSTD
            
            if doplots:
                tone_floor_fluxMEANs.append(floor_fluxMEAN)
                tone_floor_fluxSTDs.append(floor_fluxSTD)
                tone_thresholds_absflux.append(threshold_absflux)
                tone_thresholds_absmean.append(threshold_absmean)
                
            # Tone likely has spike if ANY absflux sample in window surpasses threshold
            # or equivalently if the maximum absflux does (less comparisons)
            windowflag[:, toneidx] = windows_absflux_max[:, toneidx] >= threshold_absflux
            
            # or, as is the case for LFA 4, if the abs mean flux in the window exceeds the expected
            # mean error (from central limit theorem SIG / sqrt(N)) away from nominal mean=0 
            windowflag[:, toneidx] |= np.abs(windows_flux_mean[:, toneidx]) >= threshold_absmean
            
        else:
            if doplots:
                # Tone is only NaNs, then fill masks to leave un-edited (no spikes)
                tone_floor_fluxMEANs.append(np.nan)
                tone_floor_fluxSTDs.append(np.nan)
                tone_thresholds_absflux.append(np.nan)
                tone_thresholds_absmean.append(np.nan)
            windowflag[:, toneidx] = False
    
    if crosstones > 0.0:        
        # copy original windowflag
        windowflag2 = copy.deepcopy(windowflag)
        
        if ignoreblinds==True:
            N_TOT = len(data.BolometerArray.checkChanList([]))
        else:
            N_TOT = len(data.BolometerArray.checkChanList([])) +\
                    len([tone for tone in data.tone_dict.keys() if data.tone_dict[tone]['type']=='REF' and data.tone_dict[tone]['is_good']==True])
        N_THRESH = int(crosstones/100. * float(N_TOT))
        info('Cross-checking tones...')
        for windowidx in range(len(windows_tstart)):
            flaggedtones_thiswindow = np.sum(windowflag[windowidx, :])
            if debug:
                print("DEBUG: WIN=%i | time=%.1f s | N_FLAGGED=%i | N_THRESH=%i | N_TOT=%i"%(windowidx+1, windows_time[windowidx], flaggedtones_thiswindow, N_THRESH, N_TOT))

            # if this window fulfills the criterium
            should_this_window_be_flagged = (flaggedtones_thiswindow >= N_THRESH)

            # otherwise, it's not a spike and we should not flag
            if not should_this_window_be_flagged:
                # Not real spike
                windowflag2[windowidx, :] = False

        # re-create windowflag and save memory
        windowflag = windowflag2
        windowflag2 = None
    
    # inizialize data flagging array
    # if blindtones are ignored in process, flag is false for all blindtones so they will not be affected
    flagmask = np.zeros_like(data.Data, dtype=bool)
    
    # fill data mask array
    for toneidx in range(nused):
        # if tone is spiked more than 30% of the time, simply flag it
        #if float(np.sum(windowflag[:, toneidx]))/float(len(windowflag[:, toneidx])) > 0.30:
        #    windowflag[:, toneidx] = True
        flagged_windows = np.where(windowflag[:, toneidx])[0]
        for winidx in flagged_windows:
            t_start = windows_tstart[winidx]
            t_end = t_start + windowtime
            flagmask[:, toneidx] |= (time >= t_start) & (time < t_end)

    if doplots:
        # convert to arrays
        tone_floor_fluxMEANs = np.array(tone_floor_fluxMEANs)
        tone_floor_fluxSTDs = np.array(tone_floor_fluxSTDs)
        tone_thresholds_absflux = np.array(tone_thresholds_absflux)
        tone_thresholds_absmean = np.array(tone_thresholds_absmean)

    # final info
    usedkids_idx = data.BolometerArray.checkChanList([]) - 1
    spikedfraction_usedkids = np.sum(windowflag[:, usedkids_idx], axis=0) / float(np.shape(windowflag)[0])
    spikedfraction_avg = np.nanmean(spikedfraction_usedkids)
    if np.all(spikedfraction_usedkids==0):
        info('Scan is free of spikes!')
    else:
        info('%.3f'%(spikedfraction_avg*100) + '% (avg.) of the KID timelines is lost')

    if doplots:
        info('Initializing de-spiking process interactive plot...')

        # =====================
        # Interactive plot
        # =====================

        # Message to interact with user
        msg = '\n------------------------------------------------\n'
        msg += 'Plot next Tone:                           <Enter>\n'
        msg += 'Plot previous Tone:                     - <Enter>\n'
        msg += 'Go to Tone Number:               <number> <Enter>\n'
        msg += 'Exit script execution:                  q <Enter>\n'
        msg += '-------------------------------------------------\n'
        msg += 'Enter choice: '

        # Create figure
        fig, ax = plt.subplots(2, 2, figsize=(9, 6))

        usedtones = data.BolometerArray.checkChanList([])
        testtone = usedtones[0]

        while True:
            # define data array tone index
            toneidx = testtone - 1

            # extract absflux of testtone
            absflux_testtone = absF[:, toneidx]

            # now we make nan all values in the original data that fall in a spiked window
            despiked_F_testtone = F[:, toneidx].copy()
            despiked_F_testtone[flagmask[:, toneidx]] = np.nan

            # and for the absflux
            despiked_absflux_testtone = absflux_testtone.copy()
            despiked_absflux_testtone[flagmask[:, toneidx]] = np.nan

            # and the window flux mean array
            despiked_windows_flux_mean_testtone = windows_flux_mean[:, toneidx].copy()
            despiked_windows_flux_mean_testtone[windowflag[:, toneidx]] = np.nan

            # and the window flux std array
            despiked_windows_flux_std_testtone = windows_flux_std[:, toneidx].copy()
            despiked_windows_flux_std_testtone[windowflag[:, toneidx]] = np.nan

            # Flux
            ax[0, 0].plot(time, F[:, toneidx], label='Spiked flux')
            ax[0, 0].plot(time, despiked_F_testtone, label='Spike-masked flux', zorder=1e9)
            ax[0, 0].set_xlabel('Time (s)')
            ax[0, 0].set_ylabel('Flux (Jy/b)')
            ax[0, 0].legend(framealpha=1)

            # absFlux
            ax[0, 1].plot(time, absflux_testtone)
            ax[0, 1].plot(time, despiked_absflux_testtone)
            ax[0, 1].axhline(tone_thresholds_absflux[toneidx], c='red', lw=2, label='Mean + %s Sigma'%sig)
            try:
                ax[0, 1].set_ylim(0, tone_thresholds_absflux[toneidx]*3.)
            except:
                pass
            ax[0, 1].set_xlabel('Time (s)')
            ax[0, 1].set_ylabel('abs(Flux) (Jy/b)')
            ax[0, 1].legend(framealpha=1)

            # flux window mean
            ax[1, 0].plot(windows_time, windows_flux_mean[:, testtone-1])
            ax[1, 0].plot(windows_time, despiked_windows_flux_mean_testtone)
            ax[1, 0].axhline(0)
            ax[1, 0].axhline(tone_floor_fluxMEANs[toneidx], c='green', lw=2, label='Average spike-free flux')
            ax[1, 0].axhline(tone_thresholds_absmean[toneidx], c='magenta', lw=2, label='Expected mean envelope')
            ax[1, 0].axhline(-tone_thresholds_absmean[toneidx], c='magenta', lw=2)
            ax[1, 0].set_ylim(tone_floor_fluxMEANs[toneidx]-3*tone_thresholds_absmean[toneidx], tone_floor_fluxMEANs[toneidx]+3*tone_thresholds_absmean[toneidx])
            ax[1, 0].legend(framealpha=1)
            ax[1, 0].set_xlabel('Time (s)')
            ax[1, 0].set_ylabel('Windowed flux mean (Jy/b)')

            # flux window std
            ax[1, 1].plot(windows_time, windows_flux_std[:, testtone-1])
            ax[1, 1].plot(windows_time, despiked_windows_flux_std_testtone)
            ax[1, 1].axhline(0)
            ax[1, 1].axhline(tone_floor_fluxSTDs[toneidx], c='red', lw=2, label='Average spike-free flux STD')
            ax[1, 1].legend(framealpha=1)
            ax[1, 1].set_xlabel('Time (s)')
            ax[1, 1].set_ylabel('Windowed flux STD (Jy/b)')

            fig.suptitle('Spike/instability detection for tone %i in scan %i'%(testtone, scannum))

            plt.show()

            userInput = raw_input(msg)
            userInput=str(userInput)
            if userInput == '':
                current = np.where(np.array(usedtones)==testtone)[0][0]
                testtone = usedtones[current+1]
            elif userInput == '-':
                current = np.where(np.array(usedtones)==testtone)[0][0]
                testtone = usedtones[current-1]
            elif userInput == 'q':
                plt.close(fig)
                break
            else:    
                try:
                    oldtesttone = testtone 
                    testtone = int(userInput)
                except:
                    print('Tone number must be a number.')
                    continue

            if testtone not in usedtones:
                print('Tone %s not in used tones.'%testtone)
                testtone = oldtesttone

            ax[0, 0].cla()
            ax[0, 1].cla()
            ax[1, 0].cla()
            ax[1, 1].cla()

    if full_output:
        return flagmask, spikedfraction_alltones
    else:
        return flagmask


def findspikes_Jy(windowtime=5., sig=5.5, expspikefree=75, crosstones=0., ignoreblinds=True, full_output=False, doplots=False, debug=False):
    '''Alias for findspikes_Jy_speed'''
    return findspikes_Jy_speed(windowtime=windowtime, sig=sig, expspikefree=expspikefree, crosstones=crosstones,\
                               ignoreblinds=ignoreblinds, full_output=full_output, doplots=doplots, debug=debug)

def findspikes_IQBT(windowtime=10., sig=5.5, expspikefree=75, crosstones=0., ignoreblinds=True, full_output=False, doplots=False, debug=False):
    '''
    ** VERSION 2 - 06.08.2026 **

    Finds spiked windowtime-long windows in BT-corrected I and Q data for each tone based
    on the statistics of the IQ-speed of the tone. Then cross-checks whether each spiked window
    in a tone is also spiked in other (at least crosstones%) tones and decides whether it is a spike or not.
    Outputs a mask for data.Data of whether a timestamp should be used for each tone.
    
    Timestamps within a spiked window are all flagged for a given tone, even if not all
    timestamps are spiked. This is to ensure at least <windowtime> seconds of stable data.
    Whether all non-spiked windows have the same KID state is unknown: the algorithm can
    only detect jumps, not if the tone went back to the original KID and KID state for which
    it was calibrated.

    @param windowtime:      size of windows in seconds
    @type windowtime:       float
    @param sig:             how many STDs away from mean tone speed to flag as spikes
    @type sig:              float
    @param expspikefree:    percentage of the timelines that is expected to be spike-free to determine "usual tone behavior"
    @type expspikefree:     float
    @param crosstones:      for any time window, what percentage of tones must be spiked to consider the window as really spiked. Default 25% is one chip.
    @type crosstones:       float
    @param ignoreblinds:    whether to ignore (True) or consider (False) blindtones. Useful if blindtones are spiked too.
    @type ignoreblinds:     bool
    @param full_output:     if False, only mask for data.Data is returned. Otherwise, the mask (0) and the 1-D array (1) of
                            fraction of timeline lost to spikes for all used tones (dependent of ignoreblinds!) are returned
    @type full_output:      bool
    @param doplots:         whether to do an example plot of the spike detection
    @type doplots:          bool
    @param debug:           Whether to plot debug info for spike cross-tone check
    @type debug:            bool
    '''
    
    info('Beginning spike/instability detection (IQBT)...')
    # data array
    scannum = data.ScanParam.ScanNum
    Z = data.Data
    ntones = Z.shape[1]
    febe = data.BolometerArray.FeBe
    if febe == 'AMKID870-AMKID870BE':
        fe = 'LFA'
        nkids = 880 * 4
    elif febe == 'AMKID350-AMKID350BE':
        fe = 'HFA'
        nkids = 800 * 5 * 4
    else:
        warn('This is not AMKID data!')
        return
    
    if ignoreblinds:
        Z = Z[:, :nkids]
        nused = nkids
    else:
        nused = ntones

    # time array
    info('Retrieving time axis data...')
    time = (data.ScanParam.MJD - data.ScanParam.MJD[0]) * 24 * 3600
    totaltime = time[-1] - time[0]
    dt = np.nanmedian(np.diff(time))

    # Ensure at least 100 windows
    nwindows = totaltime / windowtime
    if nwindows<50.0:
        nwindows = 50
        warn('Adjusting windowtime to %1.2f seconds to ensure 50 windows...'%(totaltime/nwindows))
        windowtime = totaltime / nwindows

    # Derivatives
    dZdt = np.diff(Z, axis=0) / dt * 1000  # mV / s
    speeds = np.abs(dZdt)
    auxtime = time[:-1] + dt/2

    # define windows
    windows_tstart = np.arange(0, totaltime, windowtime)
    windows_time = windows_tstart + windowtime/2
    info('Windowing: %.2f s -> %i windows of %.1f s'%(totaltime, len(windows_tstart), windowtime))

    # RMS of windows
    # info('Computing spike detection metrics...')
    windows_speed_mean = []
    windows_speed_std = []
    windows_speed_max = []

    for t_start in windows_tstart:
        t_end = t_start + windowtime
        valid = (auxtime >= t_start) & (auxtime < t_end)
        speeds_window = speeds[valid, :nused]

        if speeds_window.shape[0] == 0:
            windows_speed_mean.append(np.full(nused, np.nan))
            windows_speed_std.append(np.full(nused, np.nan))
            windows_speed_max.append(np.full(nused, np.nan))
        else:
            windows_speed_mean.append(np.nanmean(speeds_window, axis=0))
            windows_speed_std.append(np.nanstd(speeds_window, axis=0))
            windows_speed_max.append(np.nanmax(speeds_window, axis=0))

    # convert to arrays
    windows_speed_mean = np.array(windows_speed_mean)
    windows_speed_std = np.array(windows_speed_std)
    windows_speed_max = np.array(windows_speed_max)

    # initialize window flagging array
    windowflag = np.zeros_like(windows_speed_std, dtype=bool)
    if doplots:
        # initialize used values list
        # has as many values as tones
        tone_floor_speedMEANs = []
        tone_floor_speedSTDs = []
        tone_thresholds_speed = []

    info('Detecting spikes and instabilities...')
    # fill windows flagging array - has True for each window timestamp and tone if window has rms above a carefully chosen threshold
    for toneidx in range(nused):
        # only for non-nan tones!
        if not(np.all(np.isnan(data.Data[:, toneidx]))):
            # Tag windows to determine IQBT-speed mean and STD floors excluding the 100-expspikefree% most noisy windows
            # (hard-coded max spike proportion of timelines)
            thresh_tonebased = np.nanpercentile(windows_speed_std[:, toneidx], expspikefree)
            valid = windows_speed_std[:, toneidx] <= thresh_tonebased

            # Compute the spike-free mean of speed means for tone (MEAN floor)
            floor_speedMEAN = np.nanmean(windows_speed_mean[valid, toneidx])

            # compute the spike-free mean of speed STDs for tone (STD floor)
            floor_speedSTD = np.nanmean(windows_speed_std[valid, toneidx])

            # compute the mean + sig STD threshold
            threshold_speed = floor_speedMEAN + sig * floor_speedSTD
            
            if doplots:
                tone_floor_speedMEANs.append(floor_speedMEAN)
                tone_floor_speedSTDs.append(floor_speedSTD)
                tone_thresholds_speed.append(threshold_speed)
                
            # Tone likely has spike if ANY speed sample in window surpasses threshold
            # or equivalently if the maximum speed does (less comparisons)
            windowflag[:, toneidx] = windows_speed_max[:, toneidx] >= threshold_speed
            
        else:
            if doplots:
                # Tone is only NaNs, then fill masks to leave un-edited (no spikes)
                tone_floor_speedMEANs.append(np.nan)
                tone_floor_speedSTDs.append(np.nan)
                tone_thresholds_speed.append(np.nan)
            windowflag[:, toneidx] = False
    
    if crosstones > 0.0:        
        # copy original windowflag
        windowflag2 = copy.deepcopy(windowflag)
        
        if ignoreblinds==True:
            N_TOT = len(data.BolometerArray.checkChanList([]))
        else:
            N_TOT = len(data.BolometerArray.checkChanList([])) +\
                    len([tone for tone in data.tone_dict.keys() if data.tone_dict[tone]['type']=='REF' and data.tone_dict[tone]['is_good']==True])
        N_THRESH = int(crosstones/100. * float(N_TOT))
        info('Cross-checking tones...')
        for windowidx in range(len(windows_tstart)):
            flaggedtones_thiswindow = np.sum(windowflag[windowidx, :])
            if debug:
                print("DEBUG: WIN=%i | N_FLAGGED=%i | N_THRESH=%i | N_TOT=%i"%(windowidx+1, flaggedtones_thiswindow, N_THRESH, N_TOT))

            # if this window fulfills the criterium
            should_this_window_be_flagged = (flaggedtones_thiswindow >= N_THRESH)

            # or the previous (if exists)
            #if windowidx != 0:
            #    flaggedtones_prevwindow = np.sum(windowflag[windowidx - 1, :])
            #    should_this_window_be_flagged |= (flaggedtones_prevwindow >= N_THRESH)

            # or the next (if exists)
            #if windowidx != len(windows_tstart) - 1:
            #    flaggedtones_nextwindow = np.sum(windowflag[windowidx + 1, :])
            #    should_this_window_be_flagged |= (flaggedtones_nextwindow >= N_THRESH)

            # fulfill the criterium, then this window should be flagged
            # to allow propagation of spike across array

            # otherwise, it's not a spike and we should not flag
            if not should_this_window_be_flagged:
                # Not real spike
                windowflag2[windowidx, :] = False

        # re-create windowflag and save memory
        windowflag = windowflag2
        windowflag2 = None
    
    # inizialize data flagging array
    # if blindtones are ignored in process, flag is false for all blindtones so they will not be affected
    flagmask = np.zeros_like(data.Data, dtype=bool)
    
    # fill data mask array
    for toneidx in range(nused):
        flagged_windows = np.where(windowflag[:, toneidx])[0]
        for winidx in flagged_windows:
            t_start = windows_tstart[winidx]
            t_end = t_start + windowtime
            flagmask[:, toneidx] |= (time >= t_start) & (time < t_end)

    if doplots:
        # convert to arrays
        tone_floor_speedMEANs = np.array(tone_floor_speedMEANs)
        tone_floor_speedSTDs = np.array(tone_floor_speedSTDs)
        tone_thresholds_speed = np.array(tone_thresholds_speed)

    # final info
    usedkids_idx = data.BolometerArray.checkChanList([]) - 1
    spikedfraction_usedkids = np.sum(windowflag[:, usedkids_idx], axis=0) / float(np.shape(windowflag)[0])
    spikedfraction_avg = np.nanmean(spikedfraction_usedkids)
    if np.all(spikedfraction_usedkids==0):
        info('Scan is free of spikes!')
    else:
        info('%.3f'%(spikedfraction_avg*100) + '% (avg.) of the KID timelines is lost')

    if doplots:
        info('Initializing de-spiking process interactive plot...')

        # =====================
        # Interactive plot
        # =====================

        # Message to interact with user
        msg = '\n------------------------------------------------\n'
        msg += 'Plot next Tone:                           <Enter>\n'
        msg += 'Plot previous Tone:                     - <Enter>\n'
        msg += 'Go to Tone Number:               <number> <Enter>\n'
        msg += 'Exit script execution:                  q <Enter>\n'
        msg += '-------------------------------------------------\n'
        msg += 'Enter choice: '

        # Create figure
        fig, ax = plt.subplots(2, 2, figsize=(12, 6))

        usedtones = list(data.tone_dict.keys())
        testtone = usedtones[0]

        while True:
            # define data array tone index
            toneidx = testtone - 1

            # extract IQBT-plane speed of tone
            speed_testtone = speeds[:, toneidx]

            # create mask applied to speed array
            auxflagmask = np.zeros_like(speed_testtone, dtype=bool)
            for winidx in range(len(windows_time)):
                if windowflag[winidx, toneidx]:
                    # flag all data points that fall in this window
                    t_start = windows_time[winidx] - windowtime/2
                    t_end = windows_time[winidx] + windowtime/2
                    auxflagmask[(auxtime >= t_start) & (auxtime < t_end)] = True

            # now we make nan all values in the original data that fall in a spiked window
            despiked_Z_testtone = data.Data[:, toneidx].copy()
            despiked_Z_testtone[flagmask[:, toneidx]] = np.nan

            # and for the speed
            despiked_speed_testtone = speed_testtone.copy()
            despiked_speed_testtone[auxflagmask] = np.nan

            # and the window speed mean array
            despiked_windows_speed_mean_testtone = windows_speed_mean[:, toneidx].copy()
            despiked_windows_speed_mean_testtone[windowflag[:, toneidx]] = np.nan

            # and the window speed std array
            despiked_windows_speed_std_testtone = windows_speed_std[:, toneidx].copy()
            despiked_windows_speed_std_testtone[windowflag[:, toneidx]] = np.nan

            # IQBT
            ax[0, 0].plot(data.Data.real[:, toneidx], data.Data.imag[:, toneidx], label='Spiked IQBT')
            ax[0, 0].plot(despiked_Z_testtone.real, despiked_Z_testtone.imag, label='Spike-masked IQBT', zorder=1e9)
            ax[0, 0].set_xlabel('I (V)')
            ax[0, 0].set_ylabel('Q (V)')
            ax[0, 0].axis('equal')
            ax[0, 0].legend(framealpha=1)

            # IQBT speed
            ax[0, 1].plot(auxtime, speed_testtone)
            ax[0, 1].plot(auxtime, despiked_speed_testtone)
            ax[0, 1].axhline(tone_floor_speedMEANs[toneidx], c='green', lw=2, label='Average spike-free tone speed')
            ax[0, 1].axhline(tone_thresholds_speed[toneidx], c='red', lw=2, label='Mean + %s Sigma'%sig)
            try:
                ax[0, 1].set_ylim(0, tone_thresholds_speed[toneidx]*3.)
            except:
                pass
            ax[0, 1].set_xlabel('Time (s)')
            ax[0, 1].set_ylabel('IQ-plane speed (mV/s)')
            ax[0, 1].legend(framealpha=1)

            # IQBT speed window mean
            ax[1, 0].plot(windows_time, windows_speed_mean[:, testtone-1])
            ax[1, 0].plot(windows_time, despiked_windows_speed_mean_testtone)
            ax[1, 0].axhline(0)
            ax[1, 0].axhline(tone_floor_speedMEANs[toneidx], c='green', lw=2, label='Average spike-free tone speed')
            ax[1, 0].legend(framealpha=1)
            ax[1, 0].set_xlabel('Time (s)')
            ax[1, 0].set_ylabel('Windowed IQBT-speed mean (mV/s)')

            # IQBT speed window mean
            ax[1, 1].plot(windows_time, windows_speed_std[:, testtone-1])
            ax[1, 1].plot(windows_time, despiked_windows_speed_std_testtone)
            ax[1, 1].axhline(0)
            ax[1, 1].axhline(tone_floor_speedSTDs[toneidx], c='red', lw=2, label='Average spike-free tone speed STD')
            ax[1, 1].legend(framealpha=1)
            ax[1, 1].set_xlabel('Time (s)')
            ax[1, 1].set_ylabel('Windowed IQBT-speed STD (mV/s)')

            fig.suptitle('Spike/instability detection for tone %i in scan %i'%(testtone, scannum))

            plt.show()

            userInput = raw_input(msg)
            userInput=str(userInput)
            if userInput == '':
                current = np.where(np.array(usedtones)==testtone)[0][0]
                testtone = usedtones[current+1]
            elif userInput == '-':
                current = np.where(np.array(usedtones)==testtone)[0][0]
                testtone = usedtones[current-1]
            elif userInput == 'q':
                plt.close(fig)
                break
            else:    
                try:
                    oldtesttone = testtone 
                    testtone = int(userInput)
                except:
                    print('Tone number must be a number.')
                    continue

            if testtone not in usedtones:
                print('Tone %s not in used tones.'%testtone)
                testtone = oldtesttone

            ax[0, 0].cla()
            ax[0, 1].cla()
            ax[1, 0].cla()
            ax[1, 1].cla()

    if full_output:
        return flagmask, spikedfraction_alltones
    else:
        return flagmask
