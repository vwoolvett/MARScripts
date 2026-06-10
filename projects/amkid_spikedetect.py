def findspikes_IQBT(windowtime=10, sig=5, expspikefree=75, ignoreblinds=True, full_output=False, doplots=False, testtone=2642):
    '''

    VERSION 4.0 - 10.06.2026

    Finds spiked windowtime-long windows in BT-corrected I and Q data for each tone based
    on the statistics of the IQ-speed of the tone. Then cross-checks whether each spiked window
    in a tone is also spiked in other tones and decides whether it is a spike or not.

    Outputs a mask for data.Data of whether a timestamp should be used for each tone.
    
    Timestamps within a spiked window are all flagged for a given tone, even if not all
    timestamps are spiked. This is to ensure at least <windowtime> seconds of stable data.
    Whether all non-spiked windows have the same KID state is unknown: the algorithm can
    only detect jumps, not if the tone went back to the original KID and KID state for which
    it was calibrated.
    '''
    
    info('Beginning spike detection...')
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

    # Ensure at least 50 windows
    nwindows = totaltime / windowtime
    if nwindows < 50.0:
        nwindows = 50
        warn('Default windowtime of %1.2f seconds yields too few windows. Changing windowtime to %1.2f seconds to ensure 50 windows...'%(windowtime, totaltime/nwindows))
        windowtime = totaltime / nwindows

    # Derivatives
    dZdt = np.diff(Z, axis=0) / dt * 1000  # mV / s
    speeds = np.abs(dZdt)
    auxtime = time[:-1] + dt/2

    # define windows
    windows_tstart = np.arange(0, totaltime, windowtime)
    windows_time = windows_tstart + windowtime/2
    info('The %.2f seconds of data will be divided in %i windows of %.1f seconds each...'%(totaltime, len(windows_tstart), windowtime))

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
            flagmask[:, toneidx] = False
            windowflag[:, toneidx] = False

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
    spikedfraction_avg = np.nanmean(spikedfraction_alltones)
    spikedtime_avg = spikedfraction_avg * totaltime
    info('On average, %.2f percent (%1.1f / %1.1f seconds) of the timelines is lost due to spikes'%(spikedfraction_avg*100, spikedtime_avg, totaltime))

    if doplots:
        if not(testtone in range(1, nused+1)):
            info('Tone %s not in analized array, displaying tone 2642 (default for LFA4)'%testtone)
            testtone = 2642
        else:
            info('Plotting spike detection process for tone %s'%testtone)
         
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
        fig, ax = plt.subplots(2, 2, figsize=(20, 10))

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
        ax[0, 1].axhline(tone_floor_speedMEANs[testtone - 1], c='green', lw=2, label='Average spike-free tone speed')
        ax[0, 1].axhline(tone_thresholds_speed[testtone - 1], c='red', lw=2, label='Mean + %s Sigma'%sig)
        ax[0, 1].set_xlabel('Time (s)')
        ax[0, 1].set_ylabel('IQ-plane speed (mV/s)')
        ax[0, 1].legend(framealpha=1)

        # IQBT speed window mean
        ax[1, 0].plot(windows_time, windows_speed_mean[:, testtone-1])
        ax[1, 0].plot(windows_time, despiked_windows_speed_mean_testtone)
        ax[1, 0].axhline(0)
        ax[1, 0].axhline(tone_floor_speedMEANs[testtone - 1], c='green', lw=2, label='Average spike-free tone speed')
        ax[1, 0].legend(framealpha=1)
        ax[1, 0].set_xlabel('Time (s)')
        ax[1, 0].set_ylabel('Windowed IQBT-speed mean (mV/s)')

        # IQBT speed window mean
        ax[1, 1].plot(windows_time, windows_speed_std[:, testtone-1])
        ax[1, 1].plot(windows_time, despiked_windows_speed_std_testtone)
        ax[1, 1].axhline(0)
        ax[1, 1].axhline(tone_floor_speedSTDs[testtone - 1], c='red', lw=2, label='Average spike-free tone speed STD')
        ax[1, 1].legend(framealpha=1)
        ax[1, 1].set_xlabel('Time (s)')
        ax[1, 1].set_ylabel('Windowed IQBT-speed STD (mV/s)')

        fig.suptitle('Spike detection for test tone %i in scan %i'%(testtone, scannum))

        # Histograms
        fig, ax = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

        # speed MEAN floor histogram
        avg = np.nanmean(tone_floor_speedMEANs)
        nonan = ~np.isnan(tone_floor_speedMEANs)
        ax[0].hist(np.log10(tone_floor_speedMEANs[nonan]), bins=100)
        ax[0].axvline(np.log10(avg), c='green', lw=2, label='Mean: %.2f mV / s'%avg)
        ax[0].set_xlabel('log10(Average tone speed in mV/s)')
        ax[0].set_ylabel('Number of tones')
        ax[0].legend(framealpha=1)
        
        # speed STD floor histogram
        avg = np.nanmean(tone_floor_speedSTDs)
        nonan = ~np.isnan(tone_floor_speedSTDs)
        ax[1].hist(np.log10(tone_floor_speedSTDs[nonan]), bins=100)
        ax[1].axvline(np.log10(avg), c='red', lw=2, label='Mean: %.2f mV / s'%avg)
        ax[1].set_xlabel('log10(Average tone speed STD in mV/s)')
        ax[1].legend(framealpha=1)

        # spike percentage histogram
        meanpercentlost = np.nanmean(spikedfraction_alltones*100)
        ax[2].hist(spikedfraction_alltones*100, bins=100, range=(0, 100))
        ax[2].axvline(meanpercentlost, c='magenta', lw=3, label='Mean: %.2f percent'%meanpercentlost)
        ax[2].set_xlim(0, np.nanmax(spikedfraction_alltones*100)+2)
        ax[2].set_xlabel('Percentage of timeline lost to spikes')
        ax[2].legend(framealpha=1)

        fig.suptitle('Spike detection histograms for scan %i'%scannum)

    if full_output:
        return flagmask, spikedfraction_alltones
    else:
        return flagmask