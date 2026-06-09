def test_newreduceFsweep(fsweep ,fe='LFA', chain=None, wirescan=None):
    """
    Show diagnostics for individual KIDS for the given fsweep.
    
    @param fsweep: Scan number of the fsweep.
    @type fsweep: int
    @param fe: array under investigation ('LFA' or 'HFA')
    @type fsweep: string
    @param wirescan: Scan number of the wirescan (used to retrieve NET).
    @type fsweep: int    
    @rtype: void
    
    """
    tone_dict = reduceFsweep(fsweep,fe=fe, overwrite=False,doReturn=True)
    if not tone_dict: 
        print("Error: tone dictionary was not created")
        return

    #We need to get the sweep trace
    getIQBT(fsweep,fe=fe,sub=1)
    sweepData=copy.copy(data.Data)

    # get the frequency steps from HK data
    febe, chains, kidsPerChain = getFebe(fe)      
    data.read(fsweep,febe=febe+'_HK',subscans=(1,))
    f = np.nanmean(data.Data[11:,:], axis=1)

    #check if there is a noise trace:
    scan_dir = getScanDirFromScanNum(fsweep)
    num_subscans = len(tuple(getSubDirectories(scan_dir, '?')))
    if num_subscans==2:
        #get the noise trace
        getIQBT(fsweep,fe=fe,sub=2)
        #tone_dict = calibIQ(data.Data, data.tone_dict, is_fsweep=True, keep_iq=False)
        noiseData = copy.copy(data.Data)
        dt = np.nanmedian(data.ScanParam.get('deltat'))
        f_samp = 1.0/dt
    else:
        noiseData = np.array([],'f')
        print("INFO:Skipping noise trace analysis")

    aux = 30 - 10 * np.log10(50)
    # If a wirescan is indicated, get data from dictionary

    wire_dict={}
    if wirescan:
        wire_dict = reduceWireScan(wirescan,fe=fe)

    #Message to interact with user
    msg = '\n------------------------------------------------\n'
    msg += 'Plot next KID:                            <Enter>\n'
    msg += 'Plot previous KID:                      - <Enter>\n'
    msg += 'Go to KID Number:                <number> <Enter>\n'
    msg += 'Exit script execution:                  q <Enter>\n'
    msg += '-------------------------------------------------\n'
    msg += 'Enter choice: '
    
    fig,ax=plt.subplots(2,2) 
    fig.set_size_inches(12,8)
    ax11b=ax[1,1].twinx()

    kid_list=sorted(tone_dict.keys())
    if chain==None:
        kid_index=0
    else:
        for k in kid_list:         
            if k > ((chain-1)*kidsPerChain):
                kid=k
                kid_index=kid_list.index(kid)
                break
            if k>chain*kidsPerChain:
                print("Error: No Kids in chain %i"%chain)
                return
    kid=kid_list[kid_index]
     
    while True:
        #get NET from wire scanner (if wirescan!=None)
        if kid in wire_dict:
            NET = '%.1f'%(wire_dict[kid]['NET'])
            ID = 'fsweep %i, wire %i, NET=%s[mk*sqrt(s)]'%(fsweep,wirescan,NET)
        else:
            NET = '--'
            ID = 'fsweep %i'%(fsweep)
            
        #get Parameters from fsweep dictionary
        Z_0 =  tone_dict[kid]['circleI'] + 1j*tone_dict[kid]['circleQ']
        r   =  tone_dict[kid]['circleR']
        Z_tune = tone_dict[kid]['z_tuningpoint']
        df = tone_dict[kid]['df']
        sweep_width = tone_dict[kid]['sweep_width']
        freq = tone_dict[kid]['freq']
        resonance_freq = tone_dict[kid]['actual_freq']
        power = tone_dict[kid]['dBm']
        Q = tone_dict[kid]['Qfactor']
        responsivity = tone_dict[kid]['responsivity']
        PSD_PHI = tone_dict[kid]['KID_phase_noise']
        PSD_AMP = tone_dict[kid]['KID_amp_noise']
        dfs = freq + f 
        status=int(tone_dict[kid]['overdriven'])
        Status=['under','under','under','normal',
                'over','over','over']
        chain=int(kid/kidsPerChain)+1
        ## FFT frequencies across fsweep
        chanRes = 2200./32768.
        if tone_dict[kid]['band'] == 'LO':
            band = 0
        elif tone_dict[kid]['band'] == 'HI':
            band = 1
        f_lo3 = tone_dict[kid]['f_lo3']
        chan = int(np.round((f_lo3 - freq) / (1 * (1 - 2*band)) / chanRes))
        f0ffts = f_lo3 - 1 * (1 - 2*band) * chan * 2200. / 32768.    
        fftFreq=f0ffts+np.array(range(-10,11))*chanRes
        fftmask= np.where( (np.array(fftFreq) >= np.nanmin(dfs)) *  (np.array(fftFreq) <= np.nanmax(dfs)) )
        fftFreq=fftFreq[fftmask]
        #Figure title        
        fig.suptitle('KID %i, %s-%i, '%(kid,fe,chain)+ID)
        
        #subplot with circle information
        N=len(sweepData[:,kid-1])
        circle=plt.Circle((1E3*np.real(Z_0),1E3*np.imag(Z_0)),1E3*r,fill=False)
        ax[0,0].add_patch(circle) 
        ax[0,0].plot(1E3*np.real(sweepData[0:int(N/3+1),kid-1]),
                     1E3*np.imag(sweepData[0:int(N/3+1),kid-1]),
                     'gray', linewidth=1,label=None)
        ax[0,0].plot(1E3*np.real(sweepData[int(N/3):int(2*N/3),kid-1]),
                     1E3*np.imag(sweepData[int(N/3):int(2*N/3),kid-1]),
                     color='red', linewidth=1,
                     label='sweep trace')
        ax[0,0].plot(1E3*np.real(sweepData[int(2*N/3-1):N,kid-1]),
                     1E3*np.imag(sweepData[int(2*N/3-1):N,kid-1]),
                     'gray', linewidth=1,label=None)
        
        if noiseData.any():                 
            ax[0,0].plot(1E3*np.real(noiseData[:,kid-1]),1E3*np.imag(noiseData[:,kid-1]),
                         color='green', linewidth=4, marker='*',label='noise trace')   
        ax[0,0].plot([1E3*np.real(Z_0), 1E3*np.real(Z_tune)],[1E3*np.imag(Z_0),1E3*np.imag(Z_tune)], marker='*',color='white',linewidth=1)
        ax[0,0].set(xlabel='I [mV], BT corrected', ylabel='Q [mV], BT corrected')
        ax[0,0].axis('equal')
        ax[0,0].legend(loc='lower left')
        
        #plot resonance in dB  
        trace=20*np.log10(np.abs(sweepData[:,kid-1])) + aux

        ax[0,1].plot(dfs,trace,label='Q=%i, %s'%(Q,Status[status+3])) 
        if resonance_freq == -1:
            ax[0,1].vlines([freq],min(trace),max(trace),colors=['red','salmon'], 
                       label='Frequency shift=-[KHz]')
        else: 
            ax[0,1].vlines([freq],min(trace),max(trace),colors =['red'], 
                           label='current freq=%8.3f[MHz]'%((freq)))
            ax[0,1].vlines([resonance_freq],min(trace),max(trace),colors=['lightgreen'],linestyles=['dashed'], 
                           label='optimal freq=%8.3f[MHz]'%((resonance_freq)))

            
        ax[0,1].text(freq, min(trace), 'freq offset %i [kHz]'%((resonance_freq-freq)*1e3), ha='center',va='top')
        ax[0,1].set(xlabel='Frequency [MHz]', ylabel='Power [dBm]')
        ax[0,1].set(ylim=[min(trace)-3,max(trace+3)])
        ax[0,1].legend(loc='upper right')
        
        #subplot with phase and responsivity of sweep trace:
        Z=sweepData[:,kid-1]-Z_0
        PHI = -np.unwrap(np.angle(Z)-np.angle(Z_tune-Z_0))
        ax[1,1].plot(dfs,PHI,label='Phi')
        ax[1,1].vlines([freq],min(PHI),max(PHI),colors =['red'])
        ax[1,1].vlines([resonance_freq],min(PHI),max(PHI),colors=['lightgreen'],linestyles=['dashed'])

        ax[1,1].vlines(fftFreq,min(PHI),max(PHI),colors=['lightgray'],linestyles=['dotted'])
        ax[1,1].set(xlabel='Frequency [MHz]', ylabel='Phi [rad]')   
        ax[1,1].set(xlim=[min(dfs),max(dfs)])
        
        d_PHI=(PHI[1:]-PHI[:-1])/df
        aux_dfs=(dfs[1:]+dfs[:-1])/2
        
        ax11b.plot(aux_dfs,d_PHI,label='dPhi/dF= %.1f'%responsivity, color='yellow',marker='*')
        ax11b.set(ylim=[min(d_PHI),max(d_PHI)])
        ax11b.set_ylabel('Responsivity [rad/MHz]',color='yellow')
        ax11b.legend(loc='upper left')
        
        #subplot with noise statistics
        if noiseData.any():
            #getPhi of noise trace:
            Z=noiseData[:,kid-1]-Z_0
            PHI = -np.unwrap(np.angle(Z)-np.angle(Z_tune-Z_0))
            df = np.abs(df) # reverse sweep have a negative df
            if np.sum(np.isnan(PHI)) > 0:
                print("Warning: Noise trace has missed datapoints")
                #Noise trace reparation 
                index_nan=np.argwhere(np.isnan(PHI)>0)[0]
                print(index_nan)
                PHI=PHI[:index_nan]
            N=len(PHI)
            # Calculate power spectral distribution
            P=np.abs(np.fft.ifft(r*PHI))**2 #Signal power per bin
            P=P/df #Signal power per Hz
            power_corrected=20*np.log10(np.abs(Z_0)+r) 
            P=10*np.log10(P)  - power_corrected
            f_samples=np.linspace(0,f_samp,N)
            ax[1,0].plot(f_samples[0:N/2],P[0:N/2],color='yellow',
                         label='mean PSD_phase= %.1f'%(PSD_PHI-power_corrected))   
   
            #getAMP of noise trace:
            Z=noiseData[:,kid-1]-Z_0
            AMP = np.abs(Z)
        
            # Calculate power spectral distribution
            P=np.abs(np.fft.ifft(AMP))**2 #Signal power per bin
            P=P/df #Signal power per Hz
            P=10*np.log10(P) - power_corrected   
            ax[1,0].plot(f_samples[0:N/2],P[0:N/2],color='green',
                         label='mean PSD_AMP= %.1f'%(PSD_AMP-power_corrected))    
            ax[1,0].set(xlabel='Frequency [Hz]', ylabel='PSD [dBc/Hz]')
            ax[1,0].legend()           

        plt.show()  
        
        userInput = raw_input(msg)
        userInput=str(userInput)

        if userInput == '':
            kid_index+=1
            kid_input=kid_list[kid_index]

        elif userInput == '-':
            kid_index-=1
            kid_input=kid_list[kid_index]

        elif userInput == 'q':
            plt.close(fig)
    
        else:
            try: 
                kid_input=int(userInput)
            except:
                print('KID number must be a number')
                continue
        
        if kid_input not in kid_list:
            print('INFO: KID not in sweep dictionary')
        else:
            kid=kid_input
            kid_index=kid_list.index(kid)
        
        ax[0,0].cla()
        ax[1,0].cla()
        ax[0,1].cla()
        ax[1,1].cla()
        ax11b.cla()