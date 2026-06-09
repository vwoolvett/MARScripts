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
    ax01b=ax[0,1].twinx()

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
            
        # get Parameters from fsweep dictionary
        Z_0 =  tone_dict[kid]['circleI'] + 1j*tone_dict[kid]['circleQ']
        r   =  tone_dict[kid]['circleR']
        Z_tune = tone_dict[kid]['z_tuningpoint']
        df = tone_dict[kid]['df']
        freq = tone_dict[kid]['freq']
        resonance_freq = tone_dict[kid]['actual_freq']
        Q = tone_dict[kid]['Qfactor']
        responsivity = tone_dict[kid]['responsivity']
        PSD_PHI = tone_dict[kid]['KID_phase_noise']
        PSD_AMP = tone_dict[kid]['KID_amp_noise']
        dfs = f # offsets from sweepcenter
        freqs = freq + dfs  # absolute frequencies
        dfData = tone_dict[kid]['dfData']
        phData = tone_dict[kid]['phData']
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
        fftmask= np.where( (np.array(fftFreq) >= np.nanmin(freqs)) *  (np.array(fftFreq) <= np.nanmax(freqs)) )
        fftFreq=fftFreq[fftmask]

        # NOW COMPUTE NEW TONE PLACING
        Z = sweepData[:, kid-1]
        dZdf = np.diff(Z) / df
        dIdf = np.real(dZdf)
        dQdf = np.imag(dZdf)
        absspeed = np.abs(dZdf)

        #Figure title        
        fig.suptitle('KID %i, %s-%i, '%(kid,fe,chain)+ID)
        
        #subplot with circle information
        N = len(sweepData[:,kid-1])
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
     
        ax[0,0].plot([1E3*np.real(Z_0), 1E3*np.real(Z_tune)],[1E3*np.imag(Z_0),1E3*np.imag(Z_tune)], marker='*',color='white',linewidth=1)
        ax[0,0].set(xlabel='I [mV], BT corrected', ylabel='Q [mV], BT corrected')
        ax[0,0].axis('equal')
        ax[0,0].legend(loc='lower left')
        
        #subplot with phase and responsivity of sweep trace:
        # VWO: almost same as reduceFsweep
        Z_norm = (Z - Z_0)/(Z_tune-Z_0)
        PHI = - _correctPhases(np.angle(Z_norm))
        # Extract absolute radian shift from reduction:
        # dfs and PHI are larger arrays than dfData and phData, but
        # at df=0 the phase arrays are shifted because findPropertyForKid
        # removes the necessary 2pi multiples to bring the point of minimum
        # amplitude as close as possible to phi=0. Then we ensure that both
        # arrays coincide at df=0 by applying this constant shift:
        df0_reduction_index = np.argmin(np.abs(dfData))
        phase_df0_reduction = phData[df0_reduction_index]
        df0_here_index = np.argmin(np.abs(dfs))
        phase_df0_here = PHI[df0_here_index]
        excess = phase_df0_here - phase_df0_reduction
        PHI -= excess

        ax[0,1].plot(dfs, PHI,label='Phi')
        ax[0,1].vlines([0], min(PHI),max(PHI),colors =['red'])
        ax[0,1].vlines([resonance_freq - freq], min(PHI),max(PHI),colors=['lightgreen'],linestyles=['dashed'])
        ax[0,1].vlines(fftFreq - freq, min(PHI),max(PHI),colors=['lightgray'],linestyles=['dotted'])
        ax[0,1].set(xlabel='Frequency [MHz]', ylabel='Phi [rad]')   
        ax[0,1].set(xlim=[min(dfs),max(dfs)])
        
        d_PHI=(PHI[1:]-PHI[:-1])/df
        aux_dfs=(dfs[1:]+dfs[:-1])/2
        ax01b.plot(aux_dfs, d_PHI, label='dPhi/dF= %.1f'%responsivity, color='yellow',marker='*')
        ax01b.set(ylim=[min(d_PHI),max(d_PHI)])
        ax01b.set_ylabel('Responsivity [rad/MHz]',color='yellow')
        ax01b.legend(loc='upper left')

        # subplot with new considered IQBT trace

        # subplot with IQBT-plane speeds
        ax[1, 1].plot(dfs[:-1], absspeed)



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
            break
    
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
        ax01b.cla()