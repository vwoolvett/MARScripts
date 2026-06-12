def findspikes_IQBT(windowtime=10., sig=4.5, expspikefree=75., crosstones=10., ignoreblinds=True, full_output=False, doplots=False, debug=False):
    '''
    ** VERSION 4.1 - 11.06.2026 **

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
    @param crosstones:      for any time window, what percentage of tones must be spiked to consider them as really spiked
    @type crosstones:       float
    @param ignoreblinds:    whether to ignore (True) or consider (False) blindtones. Useful if blindtones are spiked too.
    @type ignoreblinds:     bool
    @param full_output:     if False, only mask for data.Data is returned. Otherwise, the mask (0) and the 1-D array (1) of
                            fraction of timeline lost to spikes for all used tones (dependent of ignoreblinds!) are returned
    @type full_output:      bool
    @param doplots:         whether to do an example plot of the spike detection
    @type doplots:          bool
    @param testtone:        tone to use for example plot
    @type testtone:         int
    '''
    
    info('** Beginning spike detection...')
    # data array
    scannum = data.ScanParam.ScanNum
    Z = data.Data
    ntones = Z.shape[1]
    febe = data.BolometerArray.FeBe
    if febe == 'AMKID870-AMKID870BE':
        fe = 'LFA'
    elif febe == 'AMKID350-AMKID350BE':
        fe = 'HFA'
    else:
        warn('This is not AMKID data!')
        return
    
    if scanIsBeamscan(scannum, fe):
        warn('Running de-spiking on WireScanner. Changing crosstones to 50%...')
        crosstones = 50
    
    _ , chains, kidsPerChain = getFebe(fe)
    nkids = len(chains) * kidsPerChain
    
    if ignoreblinds:
        Z = Z[:, :nkids]
        nused = nkids
    else:
        nused = ntones

    # time array
    info('Retrieving time axis data and computing IQ speeds...')
    time = (data.ScanParam.MJD - data.ScanParam.MJD[0]) * 24 * 3600
    totaltime = time[-1] - time[0]
    dt = np.nanmedian(np.diff(time))

    # Derivatives
    dZdt = np.diff(Z, axis=0) / dt * 1000  # mV / s
    speeds = np.abs(dZdt)
    auxtime = time[:-1] + dt/2

    # Ensure at least 50 windows
    nwindows = totaltime / windowtime
    if nwindows < 50.0:
        nwindows = 50
        warn('Default windowtime of %1.2f seconds is too large for only %.2f seconds of data.'%(windowtime, totaltime))
        warn('Changing windowtime to %1.2f seconds to ensure 50 windows...'%(totaltime/nwindows))
        windowtime = totaltime / nwindows

    # define windows
    windows_tstart = np.arange(0, totaltime, windowtime)
    windows_time = windows_tstart + windowtime/2

    # print final window information
    info('%.2f seconds of data: %i windows of %.1f seconds'%(totaltime, len(windows_tstart), windowtime))

    # RMS of windows
    info('Computing spike detection metrics...')
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

    info('Detecting spikes...')
    # fill windows flagging array - has True for each window timestamp and tone if window has rms above a carefully chosen threshold
    for toneidx in range(nused):
        # only for non-nan tones!
        if not(np.all(np.isnan(data.Data[:, toneidx]))):
            # Tag windows to determine IQBT-speed mean and STD floors excluding the 100-expspikefree% most noisy windows
            # (hard-coded max spike proportion of timelines)
            thresh_tonebased = np.nanpercentile(windows_speed_std[:, toneidx], expspikefree)
            valid = windows_speed_std[:, toneidx] <= thresh_tonebased

            # Compute the spike-free mean of window speed means for tone (MEAN floor)
            floor_speedMEAN = np.nanmean(windows_speed_mean[valid, toneidx])
            
            # compute the spike-free mean of window speed STDs for tone (STD floor)
            floor_speedSTD = np.nanmean(windows_speed_std[valid, toneidx])
            
            # compute the mean + sig STD threshold
            threshold_speed = floor_speedMEAN + sig * floor_speedSTD

            if doplots:
                # save values
                tone_floor_speedMEANs.append(floor_speedMEAN)
                tone_floor_speedSTDs.append(floor_speedSTD)
                tone_thresholds_speed.append(threshold_speed)

            # Tone likely has spike if ANY speed sample in window surpasses threshold
            # or equivalently if the maximum speed does (less comparisons)
            windowflag[:, toneidx] = windows_speed_max[:, toneidx] >= threshold_speed

        else:
            # Tone is only NaNs, then fill masks to leave un-edited (no spikes)
            tone_floor_speedMEANs.append(np.nan)
            tone_floor_speedSTDs.append(np.nan)
            tone_thresholds_speed.append(np.nan)
            windowflag[:, toneidx] = False

    
    # NOTE: OLD METHOD COMPARE ALL KIDS (ACTUALLY WORKS BETTER THAN CHAIN-BASED)
    # Spikes are only real if they appear in the same time window as at least
    # crosstones% of all the tones.
    # info('Cross-checking tones...')
    # for windowidx in range(len(windows_tstart)):
    #     flaggedtones_thiswindow = np.sum(windowflag[windowidx, :])
    #     if debug:
    #         print("DEBUG: WIN=%i | NFLAGGED=%i | NTHRESH=%.3f"%(windowidx, flaggedtones_thiswindow, int(crosstones/100. * float(nused))))
    #     if flaggedtones_thiswindow <= int(crosstones/100. * float(nused)):
    #         # Not real spike
    #         windowflag[windowidx, :] = False


    # Spikes are only real if they appear in the same time window as at least
    # crosstones% of the tones, chain-wise.
    # NOTE: NEW METHOD COMPARE PER CHAIN (ACTUALLY DOES NOT WORK AS WELL AS OLD METHOD)
    # info('Cross-checking tones per chain...')
    # _, chains, kidsPerChain = getFebe(fe)
    # nkids = len(chains) * kidsPerChain
    # kididx_in_chain = np.array([np.arange(kidsPerChain*(chain-1), kidsPerChain*(chain)) for chain in chains])
    # for chain in chains:
    #     chainidx = chain - 1
    #     kididx_here = kididx_in_chain[chainidx]
    #     for windowidx in range(len(windows_tstart)):
    #         flaggedtones_thischain_thiswindow = np.sum(windowflag[windowidx, kididx_here])
    #         if debug:
    #             print("DEBUG: CHAIN=%i, WIN=%i | NFLAGGED=%i | NTHRESH=%.3f"%(chain, windowidx, flaggedtones_thischain_thiswindow, int(crosstones/100. * float(kidsPerChain))))
    #         if flaggedtones_thischain_thiswindow <= int(crosstones/100. * float(kidsPerChain)):
    #             # Not real spike
    #             windowflag[windowidx, kididx_here] = False


    windowflag2 = copy.deepcopy(windowflag)

    # NOTE: TAKE 3: COMPARE ADJACENT WINDOWS TOO (WORKS BEST)
    info('Cross-checking tones...')
    for windowidx in range(len(windows_tstart)):
        flaggedtones_thiswindow = np.sum(windowflag[windowidx, :])
        if debug:
            print("DEBUG: WIN=%i | NFLAGGED=%i | NTHRESH=%.3f"%(windowidx, flaggedtones_thiswindow, int(crosstones/100. * float(nused))))

        # if this window
        should_this_window_be_flagged = (flaggedtones_thiswindow >= int(crosstones/100. * float(nused)))

        # or the previous (if exists)
        if windowidx != 0:
            flaggedtones_prevwindow = np.sum(windowflag[windowidx - 1, :])
            should_this_window_be_flagged |= (flaggedtones_prevwindow >= int(crosstones/100. * float(nused)))

        # or the next (if exists)
        if windowidx != len(windows_tstart) - 1:
            flaggedtones_nextwindow = np.sum(windowflag[windowidx + 1, :])
            should_this_window_be_flagged |= (flaggedtones_nextwindow >= int(crosstones/100. * float(nused)))

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
    spikedfraction_alltones = np.sum(windowflag, axis=0) / float(np.shape(windowflag)[0])
    spikedfraction_alltones_nonzero = spikedfraction_alltones[spikedfraction_alltones != 0]
    if len(spikedfraction_alltones_nonzero) != 0:
        spikedfraction_avg = np.nanmean(spikedfraction_alltones_nonzero)
        spikedtime_avg = spikedfraction_avg * totaltime
        info('On average for spiked tones, %.2f percent (%1.1f / %1.1f seconds) of the timelines is lost.'%(spikedfraction_avg*100, spikedtime_avg, totaltime))
        scanisspikefree = False
    else:
        info('Scan is spike-free!')
        scanisspikefree = True


    if doplots:
        info('Initializing de-spiking process histograms and interactive plot...')
        # ==========
        # Histograms
        # ==========
        fighist, axhist = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

        # speed MEAN floor histogram
        avg = np.nanmean(tone_floor_speedMEANs)
        nonan = ~np.isnan(tone_floor_speedMEANs)
        axhist[0].hist(np.log10(tone_floor_speedMEANs[nonan]), bins=100)
        axhist[0].axvline(np.log10(avg), c='green', lw=2, label='Mean: %.2f mV / s'%avg)
        axhist[0].set_xlabel('log10(Average tone speed in mV/s)')
        axhist[0].set_ylabel('Number of tones')
        axhist[0].legend(framealpha=1)
        
        # speed STD floor histogram
        avg = np.nanmean(tone_floor_speedSTDs)
        nonan = ~np.isnan(tone_floor_speedSTDs)
        axhist[1].hist(np.log10(tone_floor_speedSTDs[nonan]), bins=100)
        axhist[1].axvline(np.log10(avg), c='red', lw=2, label='Mean: %.2f mV / s'%avg)
        axhist[1].set_xlabel('log10(Average tone speed STD in mV/s)')
        axhist[1].legend(framealpha=1)

        # spike percentage histogram
        if not scanisspikefree:
            meanpercentlost = np.nanmean(spikedfraction_alltones_nonzero*100)
            axhist[2].hist(spikedfraction_alltones_nonzero*100, bins=100, range=(0, 100))
            axhist[2].axvline(meanpercentlost, c='magenta', lw=3, label='Mean: %.2f percent'%meanpercentlost)
            axhist[2].set_xlim(0, np.nanmax(spikedfraction_alltones_nonzero*100)+2)
            axhist[2].set_xlabel('Percentage of timeline lost to spikes (spiked tones only!)')
            axhist[2].legend(framealpha=1)
        
        else:
            meanpercentlost = np.nanmean(spikedfraction_alltones*100)
            axhist[2].hist(spikedfraction_alltones*100, bins=100, range=(0, 100))
            axhist[2].axvline(meanpercentlost, c='magenta', lw=3, label='Mean: %.2f percent'%meanpercentlost)
            axhist[2].set_xlim(0, np.nanmax(spikedfraction_alltones*100)+2)
            axhist[2].set_xlabel('Percentage of timeline lost to spikes (scan is spike-free!)')
            axhist[2].legend(framealpha=1)

        fighist.suptitle('Spike detection histograms for scan %i'%scannum)
        
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
        fig, ax = plt.subplots(2, 2, figsize=(20, 10))

        usedtones = np.arange(1, nused+1)
        testtone = 1

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

            # testplot
            

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

            fig.suptitle('Spike detection for test tone %i in scan %i'%(testtone, scannum))

            plt.show()

            userInput = raw_input(msg)
            userInput=str(userInput)
            if userInput == '':
                testtone += 1
            elif userInput == '-':
                testtone -=1
            elif userInput == 'q':
                plt.close(fig)
                plt.close(fighist)
                break
            else:    
                try: 
                    testtone = int(userInput)
                except:
                    print('Tone number must be a number.')
                    continue

            if testtone not in usedtones:
                print('Resulting tone out of bounds. Going back to tone 1.')
                testtone = 1

            ax[0, 0].cla()
            ax[0, 1].cla()
            ax[1, 0].cla()
            ax[1, 1].cla()

    if full_output:
        return flagmask, spikedfraction_alltones
    else:
        return flagmask