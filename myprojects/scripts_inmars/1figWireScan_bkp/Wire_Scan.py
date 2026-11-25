# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 16:23:11 2021
Modified to operate in Mars Sept 15
@author: Reyes
"""
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as stat
from matplotlib import patches
from scipy.optimize import minimize
from operator import itemgetter
from itertools import groupby  

#Code to reduce the data from wire scan and from noise scan.
#12.10.2023 add command to decorrelate data (correlChains())
#18.10.2023 rename data to dataset in several functions, noisespectra function to work with data.data
#14.2.2024  get number of kids and chains from getFebe(), 
#           give option of using beamscanner to estimate noise.
#           avoid the use of copy.deepcopy, script must be compatible with MARS enviroment
#           add KID resonance frequency as parameter for the analysis
#           comments and cleaning on makeRCP()
#5.7.2024   add option to save and reload results from file, add elevation to analysis
#11.7.2024  results are now stored in a dictionary, deprecating of unused functions
#           Spectra() and plot_Chain_NoiseSpectra() merged in single function.
#14.7.2025  changing calculation of NET, new NET will be higuher by a factor sqrt(2)!
#30.9.2025  plot_NoiseSpectra new version. Avoid normalization and better control of plot.
#15.12.2025 changed NET computation (VWO)
#16.12.2025 NR: change the code to select trace in X+,X-,Y+ and Y-
#17.2.2026  NR: add CompareRCPtoWire()



#This set of global variables will be initialized in the WireScan_full_Analisys function
freq=0  
T_sky=0 
wire_d=0
kids=[] 
wire_scan=[] 
noise_scan=[]
ID=''

def __WireScan_full_Analysis(BeamScanner, Noise = None, fsweep=None, doDf=False, fe='LFA', chains=None):
    """
    Tool to generate a standard set of plot to show current array performance. 

    @ param BeamScanner: Scan number of beam scanner
    @ type  BeamScanner: int
    @ param Noise: Scan number of the noise trace used for sensitivity calculation. 
    If None: a clean section of the BeamScanner is used
    @ type  Noise: int 
    @ param fsweep: Scan number of fsweep, if None: system look for the most recent one
    @ type  fsweep: int
    @ param fe: 'LFA' or 'HFA', array to be analyzed
    @ type  fe: string
    @ param chains: chains to be used in analysis. if None: all of them are used.
    @ param chains: iterable
    
    It returns a table with the following information:
    position=[ kid#,[pos_x,pos_y],[size_x,size_y],[noise_x,noise_y],[peak_x,peak_y] ]
    """
    global freq
    global T_sky 
    global wire_d
    global kids
    global wire_scan 
    global noise_scan
    global x
    global y
    global ID
    global data
    #Check we have a valid wire scan
    if not(scanIsBeamscan(BeamScanner,fe)):
        print('Scan is not a wire scan, use findWirescanForScan() to find a valid ScanNr')
        return
        
    # Initialize fe dependent variables
    frontEnd=fe    
    if Noise != None:    ID=fe + ', Scan %i & Scan %i, '%(BeamScanner, Noise)
    else: ID=fe + ', Scan %i '%(BeamScanner)

    if fe=='LFA':
        wire_d=0.5  #Effective wire size (from lab measurement)
    elif fe=='HFA':
        wire_d=0.65
    else: 
        print('fe must be HFA or LFA')
        return
    if fsweep == None:
        fsweep=findFsweepForScan(BeamScanner,fe=frontEnd,sweep_type='cal')
        #fsweep = findFsweepFromObslog(BeamScanner)
        #print("Using fsweep %i for calibration"%(fsweep))
    #initialize chains in case that all of them are required
    if chains==None:
        _,chains,n_kid=getFebe(fe)
    else:
        _,_,n_kid=getFebe(fe)	
        
    #Get Beam scanner positions
    try:
        _, x, y = readBeamscanPositions(BeamScanner, fe=fe)
        warn_relative_positions = False
    except RuntimeError as e:
        warn_relative_positions = True
        warn(e.message)
        warn('Unable to get beamscanner absolute positions - Getting relative positions.')
        _, x, y = readBeamscanPositions(Beamscanner, use_relative_positions=True, fe=fe)
    
    #if specified:load noise scan
    if Noise != None: 
        if doDf:
            getDf(Noise, fsweep=fsweep, fe=fe)
        else:
            getPhi(Noise, fsweep=fsweep, fe=fe)
        flagEmpty()
        flagNaN()
        correlChains()
        data.zeroStartEnd()
        base(order=3)
        noise_scan=data.Data
    #Get Tsky from Ambient Temp PWV and elevation and C-cabin temperature
    T_sky = scanTsky(febe=fe)
    values=fetchMonPoint(data, subscan=1, MonPo=['T_CABIN_C'])[0]
    maskTC=np.where(values > 0)[0]
    if len(maskTC) > 0:
        T_cabin_C=np.nanmedian(values)+273.15
    else:
        T_cabin_C=18.0+273.15

    #Get Beam scanner data used during analysis
    if doDf:
        getDf(BeamScanner,fsweep=fsweep,fe=fe)
    else:
        getPhi(BeamScanner,fsweep=fsweep,fe=fe)
    flagEmpty()
    flagNaN()      
    data.zeroStart()
    flagFractionRms()
    ##The new correlChips function does not need to flag the wire
    #flagArray = __findWire()
    #data.FlagHandler.setOnMask(flagArray,iFlags=[8])
    #correlChains()
    correlChips()
    base(order=1)
    #data.unflag(flag=8)   
    
    #List of kids to be analyzed
    kids = data.BolometerArray.checkChanList([])
    #print('FRONTEND1: '+getFrontend(data.BolometerArray.FeBe))
    #data.zeroStartEnd()
    wire_scan=copy.deepcopy(data.Data)

    #Get the sampling frequency from Beam scanner. We assume that is the same for the noise scan
    t = data.getChanData('mjd', flag='None', getFlagged=0)
    t_int = np.median(t[1:] - t[:-1])
    freq=1/t_int
    elevation=np.nanmean(data.ScanParam.El)

    #Now start data processing
    positions,status=__find_Wire(kids,checkPos=True)
    if status == -1:
        print( 'Warning: Could not determine positions and beam width from Wire Scanner')
        print( 'Warning: This hints to a problem with the time syncronisation of the Wire Scanner')
        print( 'Warning: Contact Mirek and Nico to look at this problem')
        positions,_=__find_Wire(kids,checkPos=False)
        
        
    positions2=find_good(positions)
    #print('FRONTEND2: '+getFrontend(data.BolometerArray.FeBe))
    if Noise != None: 
        #replace noise estimation by one taken from noise scan
        positions4=get_Noise(positions2)
    else:
        positions4=positions2
        # get the noise trace from from beam_scanner trace. to check if this code is necessary!
        """N,M=len(wire_scan[0,:]),len(wire_scan[:,0])
        noise_scan=[]
        for kid in range(N):
            l=int(M/11)
            trace_aux=wire_scan[0:l,kid]
            aux=np.nanstd(trace_aux)
            noise_scan.append(trace_aux)
            if aux!= aux: 
                continue
            for j in range(1,11):
                trace_aux=wire_scan[j*l:(j+1)*l,kid]
                noise=np.nanstd(trace_aux)
                if noise != noise or noise>aux: continue
                else: 
                    aux=noise
                    noise_scan[kid]= trace_aux
        noise_scan=np.transpose(np.array(noise_scan))"""
    #print('FRONTEND3: '+getFrontend(data.BolometerArray.FeBe))
    positions5=__get_KidFreq(positions4,BeamScanner,fe=fe)
    #positions4=find_good(positions3)
    parameters2=[]
    #print('FRONTEND4: '+getFrontend(data.BolometerArray.FeBe))
    return positions5,wire_d,T_sky,T_cabin_C,status
    
def plot_NoiseSpectra(kids=[], chains=[], doReturn=False,factor=1, method='median',
                      doPlot=True, ax=None, color='red', label=None, log=True):   
    '''Function to plot the Noise spectra from current data structure
    kids= list of kids to be used to generate the noise spectra
    chains= list of chains used to generate the noise spectra, if both kids and chain are given it just consider chains
    method= 'mean', 'median' or 'weight-mean'
    factor= number of bin used in the smoothing of the spectra
    ax= figure where to plot, if None new figure is generated
    Figure controls: log, label, color
    doReturn= Trur for function to return relevant data frequency, Spectra, KID_list, NET_kids
   '''
    #get system info
    global data
    noise=data.Data
    fe=data.tone_dict_fe
    _,_,n_kid=getFebe(fe) 
    if ax==None and doPlot: 
        fig,ax=plt.subplots(1)
    #Get the sampling frequency
    t = data.getChanData('mjd', flag='None', getFlagged=0)
    t_int = np.nanmedian(t[1:] - t[:-1])
    freq=1/t_int
    N=len(noise[:,1])
    f_0=freq #sampling frequency
    df=f_0/N
    f=np.linspace(0,f_0,N)[0:int(N/2)]
    #get the list of kids to operate
    kid_list=[]
    if kids!= [] and chains==[]:
        kid_list=kids
    elif kids== [] and chains==[]:
        print("error: both kids and chains cannot be empty!")
        return
    elif kids!= [] and chains!=[]:
        print("error: both kids and chains cannot coexist, you must chose!")
        return
    elif chains=='all':
        _,chains,n_kid=getFebe(fe)
    else:
        if type(chains) is int: chains=[chains]
        for chain in chains: 
            kid_list += range((chain-1)*n_kid+1,chain*n_kid) 
    if (np.nansum(noise[:, kid_list])==0):
        print('INFO: No data for selected chains')
        return
    spectra=[]
    NET=[]
    KIDs=[]
    repaired=0
    #generate Kermnel for smoothing
    if factor >1:
        aux=np.ones(round(factor))
        aux=np.convolve(aux,aux)  #generate a triangular kernel
        aux=(aux-factor)*2/factor
        aux=np.exp(-aux**2)       #generate a gaussian kernel
        aux=aux/np.sum(aux)
        plt.plot(aux)
    for kid in kid_list:
        P_total=np.zeros(int(N/2))
        traces=0
        if kids ==[] and (kid not in data.tone_dict.keys()):
            print('f')
            continue
        dataset=noise[:, kid-1]
        #Check if data has Nans
        numberNan=np.sum(np.isnan(dataset))
        if numberNan==len(dataset): #nonexisting data:
            continue        
        if numberNan>0:#  #fixing trace:
            for index in np.argwhere(np.isnan(dataset)):
                dataset[index]= dataset[index-1]  
                repaired+=1
        #This is for the doreturn option
        NET.append(np.std(dataset)/np.sqrt(freq))
        KIDs.append(kid)
        #Now calcualte power spectra
        P=(np.abs(np.fft.ifft(dataset))**2)[0:int(N/2)] #Signal power per bin, 
        #ANG=np.angle(np.fft.ifft(dataset)[0:int(N/2)])
        P=P/df #Signal power per Hz
        #Smoothing must be done in the power spectra domain!
        if factor >1:
            f,P = np.convolve(f,aux), np.convolve(P,aux)
            f,P=f[factor+2:-factor*2],P[factor+2:-factor*2]
        #Now going to the NET domain in mK/sqrt(s)
        P=np.sqrt(P) 
        spectra.append(P)
        
    if repaired !=0:
        print("Warning: %i data points with Nans in data! Attemp to repair trace."%repaired)
    spectra=np.array(spectra) 
    if method == 'median':
        P_total=np.nanmedian(spectra, axis=0)
    elif method == 'mean':
        P_total=np.nanmean(spectra, axis=0)
    elif method =='weight-mean':
        w=1/(np.array(NET)**2)
        for i in range(len(w)):
            spectra[i]=spectra[i]*w[i]
        P_total=np.nanmean(spectra, axis=0)/np.mean(w)
    else:
        print("Error: method not supported")
        return
    
    if not doPlot and doReturn:
        return f,P_total,np.array(KIDs),np.array(NET)
    ax.plot(f,P_total,color=color,linewidth=1,label=label)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('PSD('+ data.Unit +')')
    ax.set(title='Scan %i, Noise spectra'%data.ScanParam.ScanNum)
    ax.legend()
    if log:
        ax.set_xscale('log')
        ax.set_yscale('log')  
    if label:
        ax.legend()
    if doReturn:
        return f,P_total,np.array(KIDs),np.array(NET)
    else:
        return
    
def __find_Wire(kids=[12],checkPos=True):
    """Function to analyze a wire scan. It returns an array with kid 
    information=[kid, position,sizes,maxima,noise]
    noise is retrieved from the part of the traces were the beam is not present
    """
    positions=[]
    scan=copy.deepcopy(wire_scan)
    N=len(scan[:,1])
    pos_failed=0
    pos_ok=0
    for kid in kids:        
        dataset=scan[:,kid-1]
        #check for nan in dataset, if more than 10% skip KID
        if np.sum(np.isnan(dataset))/len(dataset) > 0.1: continue
        #Estimation of the noise signal in beam scanner measurement
        section_noise=np.zeros(N)
        for i in range(20):
            noise=np.nanstd(dataset[int(i*N/20):int((i+1)*N/20)])
            section_noise[int(i*N/20):int((i+1)*N/20)]=noise
        xp=np.arange(0,len(scan),1.)
        
        #plt.clf()
        #plt.plot(xp,scan[:,kid-1])
        #mask=np.where(scan[:,kid-1] > 20*np.min(section_noise))[0]
        #plt.plot(xp[mask],scan[mask,kid-1],c='green')
        #plt.show()
        #raw_input()
        
        noise=np.nanmin(section_noise)/(freq**0.5) #Noise in Phi/Hz**0.5
        
        #Now estimate beam position and size
        dx=x[1:]-x[:-1]
        x1=np.where(dx>0.1)[0]
        data_X1=scan[x1, kid-1]
        X1=x[x1]
        Xpos1,Xsize1,Xpeak1=gauss_fit(X1, data_X1)
        
        x2=np.where(dx<-0.1)[0]
        data_X2=scan[x2, kid-1]
        X2=x[x2]
        Xpos2,Xsize2,Xpeak2=gauss_fit(X2, data_X2)
        
        if Xpos1 == -9999 and Xsize1 == -9999 and Xpeak1 == -9999 and Xsize2 != -9999:
            Xpos1  = Xpos2
            Xsize1 = Xsize2
            Xpeak1 = Xpeak2
        if Xpos2 == -9999 and Xsize2 == -9999 and Xpeak2 == -9999 and Xsize1 != -9999:
            Xpos2  = Xpos1
            Xsize2 = Xsize1
            Xpeak2 = Xpeak1
            
        dy=y[1:]-y[:-1]
        y1=np.where(dy>0.1)[0]
        data_Y1=scan[y1, kid-1]
        Y1=y[y1]    
        Ypos1,Ysize1,Ypeak1=gauss_fit(Y1, data_Y1)
        
        y2=np.where(dy<-0.1)[0]
        data_Y2=scan[y2, kid-1]
        Y2=y[y2]    
        Ypos2,Ysize2,Ypeak2=gauss_fit(Y2, data_Y2)

        if Ypos1 == -9999 and Ysize1 == -9999 and Ypeak1 == -9999 and Ysize2 != -9999:
            Ypos1  = Ypos2
            Ysize1 = Ysize2
            Ypeak1 = Ypeak2
        if Ypos2 == -9999 and Ysize2 == -9999 and Ypeak2 == -9999 and Ysize1 != -9999:
            Ypos2  = Ypos1
            Ysize2 = Ysize1
            Ypeak2 = Ypeak1
         
        if 150 > Xsize1 > 0 and 150 > Xsize2 > 0 and 150 > Ysize1 > 0 and 150 > Ysize2 > 0:
            ok = 1
        else:
            ok = 0

        if checkPos == True:
            if np.abs(Xpos1-Xpos2) < 15.0 and np.abs(Ypos1-Ypos2) < 15.0 and ok == 1: #test positions and beam width (default)
                pos1=np.mean([Xpos1,Xpos2])
                pos2=np.mean([Ypos1,Ypos2])
                size1=np.mean([Xsize1,Xsize2])
                size2=np.mean([Ysize1,Ysize2])
                peak1=np.mean([Xpeak1,Xpeak2])
                peak2=np.mean([Ypeak1,Ypeak2])
                kid_freq=-1
                kid_power = -1
                positions.append([kid,[pos1,pos2],[size1,size2],[noise,noise],[peak1,peak2],kid_freq,kid_power])
                pos_ok=pos_ok+1
            else:
                pos_failed=pos_failed+1
            
        else:
            if 0.7*Xpeak1 < Xpeak2 < 1.3*Xpeak1 and 0.7*Ypeak1 < Ypeak2 < 1.3*Ypeak1: #test peak intensity of wire signal in case postion info is corrupted
                pos1=-9999
                pos2=-9999
                size1=-1
                size2=-1
                peak1=np.mean([Xpeak1,Xpeak2])
                peak2=np.mean([Ypeak1,Ypeak2])
                kid_freq=-1
                kid_power = -1
                positions.append([kid,[pos1,pos2],[size1,size2],[noise,noise],[peak1,peak2],kid_freq,kid_power])
                pos_ok=pos_ok+1
            else:
                pos_failed=pos_failed+1
            
    if checkPos == True:
        print("-----------------------------------------------------------------")
        print("Found %i Kid positions on wire"%(pos_ok))
        print("Failed to determine Kid positions for %i Kids"%(pos_failed))
    else:
        print("-----------------------------------------------------------------")
        print("No position information from wire scanner used in reduction")
        print("Determined response for %i Kids from the wirescan"%(pos_ok))
        print("Failed to determined response for %i Kids"%(pos_failed))

        #frequency estimated from KID #. To be changed: frequency is avalaible in HK data
        #kid_freq=0
        #positions.append([kid,[pos1,pos2],[size1,size2],[noise,noise],[max(data_X),max(data_Y)],kid_freq])
    if len(positions)==1:
        plt.clf()
        plt.plot(dataset)

    if pos_ok < 150:
        status = -1
    else:
        status = 1
    
    return positions,status

def gauss_fit(x, y,checkRise=False):
    """Simple function to do a gauss fit of the wire signal"""
    # Fitting is done over a +-20mm section of the wire scan
    delta=20
    try:
        maxsimple=np.nanargmax(y) #position of maxima
    except:
        return -9999,-9999,-9999
    ## New code to find the maximum avoiding that on picks up a noise spike
    if checkRise==True:
        smoothed = smoothSignal(y, window_len=5, window='hanning')
        maxsimple=np.nanargmax(smoothed)
        nums=np.where(np.diff(smoothed) > 0)[0]
        rise=[0]
        for k, g in groupby(enumerate(nums), lambda (i, vx): i-vx):
            result = list(map(itemgetter(1), g))
            result = result+list([result[-1]+1])
            if len(list(result)) >= 3 and np.max(smoothed[list(result)]) > np.max(smoothed[rise]):
                rise = list(result)

        nums=np.where(np.diff(smoothed) < 0)[0]
        fall=[0]
        for k, g in groupby(enumerate(nums), lambda (i, vx): i-vx):
            result = map(itemgetter(1), g)
            if len(list(result)) >= 3 and np.max(smoothed[list(result)]) > np.max(smoothed[fall]):
                fall = list(result)
        
    
        if len(fall) > 1 and maxsimple in fall: 
            maxind = maxsimple
        elif len(rise) > 1 and maxsimple in rise: 
            maxind = maxsimple
        elif np.max(rise) <= maxsimple <= np.min(fall) and np.abs(np.min(fall)-np.max(rise)) < 10:
            maxind = maxsimple
        elif len(fall) > 1 and len(rise) > 1 and np.max(rise) == np.min(fall):
            maxind = np.max(rise)
        else:
            maxind = 0

        if maxind > 0:
            #subtract baseline
            fitx=np.concatenate((x[0:np.min(rise)],x[np.max(fall):len(x)]))
            fity=np.concatenate((y[0:np.min(rise)],y[np.max(fall):len(x)]))
            p = np.polyfit(fitx, fity, 1)
            y= y - (p[1]+p[0]*x)
            #print("ok")
            #xp=np.array(range(len(x)))
            #plt.clf()
            #xp=np.array(range(len(x)))
            #plt.plot(xp,y)
            #plt.scatter(xp,y)
            #plt.plot(xp[rise],y[rise],c='green')
            #plt.plot(xp[fall],y[fall],c='red')
            #plt.show()
            #raw_input()
    else:
        maxind = maxsimple

    if maxind > 0:
        pos=x[maxind]
        y=y[(x>(pos-delta)) & (x<(pos+delta))]
        y=y-np.min([y[0],y[-1]])
        x=x[(x>(pos-delta)) & (x<(pos+delta))]
        N=len(x)
        step=(np.nanmax(x)-np.nanmin(x))/N
        pos=np.nansum(x*y)/np.nansum(y)
        size=np.abs(np.nansum(y)*step/np.nanmax(y))
        peak=np.nanmax(y)
        
        #plt.clf()
        #plt.plot(x,y)
        #plt.show()
        #raw_input()
        return pos,size,peak
    else:
        #print("failed")
        #print(np.max(rise), maxsimple, np.min(fall), np.min(fall)-np.max(rise))
        #plt.clf()
        #xp=np.array(range(len(x)))
        #plt.plot(xp,y)
        #plt.scatter(xp,y)
        #plt.plot(xp[rise],y[rise],c='green')
        #plt.plot(xp[fall],y[fall],c='red')
        #plt.show()
        #raw_input()
        return -1, -1, -1


def __get_KidFreq(positions,scannerID,fe='LFA'):
    """Function to add KID resonance frequency and power information to the position table
    returns a new position table"""
    #scan=findFsweepForScan(scannerID,fe=fe,sweep_type='cal')
    #febe,_,n_kid=getFebe(fe)
    #getIQ(scan,fe=fe)
    #data.read(scan,febe=febe+'_HK')
    #getHK(scannerID)
    aux_positions=copy.deepcopy(positions)
    tone_list=data.tone_dict.keys()
    for pos in aux_positions:
        tone=pos[0]
        if tone in tone_list:
            pos[5]= data.tone_dict[tone]['freq']#data.Data[1,tone]#
            pos[6]= data.tone_dict[tone]['power']
    return aux_positions

def get_Noise(positions,fstart=3,fstop=15):   
    '''Function to replace the initial estimation of noise, take from the wire scan data, 
    by the noise measurement on a stable load. Noise is calculated in the fourier domain 
    as the average noise spectra above 3Hz. 
    Input is list of kid parameters obtained by '''
    aux_positions=copy.deepcopy(positions)
    noise=(noise_scan)
    N=len(noise[:,1])
    for pos in aux_positions:
        kid=pos[0]
        dataset=noise[:, kid-1]
        f_0=freq #sampling frequency
        df=f_0/N
        f=np.linspace(0,f_0,N)[0:int(N/2)]
        P=(np.abs(np.fft.ifft(dataset))**2)[0:int(N/2)] #Signal power per bin
        P=P/df #Signal power per Hz
        aux=np.sqrt(2)*np.mean(P[(f>fstart) & (f<fstop)])**0.5
        pos[3][0]=aux
        pos[3][1]=aux
    if len(positions)==1:
        plt.plot(f,P,'r')
        plt.plot(f[(f>fstart) & (f<fstop)],P[(f>fstart) & (f<fstop)])
        plt.yscale('log')
    return aux_positions
"""
def plot_NoiseSpectra(kids=[], chains=[], color='red',ax=None, label=None, factor=1):   
    '''Function to plot the Noise spectra from current data structure
    If ax is given: the plot is generated on the indicated axes
    If chains is given: plot the average kid of all kids in the specified chains
    If kids is given: plot the average spectra of KIDS in list
    color and label are passed to the plot.
    factor give user control over the smoothing factor (just average over number of bins)
   '''
    #get system info
    global data
    noise=data.Data
    fe=data.tone_dict_fe
    _,_,n_kid=getFebe(fe) 
    if ax==None: 
        fig,ax=plt.subplots(1)
    #Get the sampling frequency
    t = data.getChanData('mjd', flag='None', getFlagged=0)
    t_int = np.nanmedian(t[1:] - t[:-1])
    freq=1/t_int
    N=len(noise[:,1])
    f_0=freq #sampling frequency
    df=f_0/N
    f=np.linspace(0,f_0,N)[0:int(N/2)]
    #get the list of kids to operate
    if kids!= [] and chains==[]:
        kid_list=kids
    elif kids== [] and chains==[]:
        print("error: both kids and chains cannot be empty!")
        return
    elif kids!= [] and chains!=[]:
        print("error: both kids and chains cannot coexist, you must chose!")
        return
    elif chains=='all':
        _,chains,n_kid=getFebe(fe)
    else: # chains != []:
        for chain in chains: kid_list = range((chain-1)*n_kid+1,chain*n_kid) 
    spectra=[]
    net=[]
    repaired=0
    for kid in kid_list:
        P_total=np.zeros(int(N/2))
        traces=0
        if (np.nansum(noise[:, kid_list])==0):
            print('INFO: No data for selected chains')
            return
        if kid not in data.tone_dict.keys():
            continue
        dataset=noise[:, kid-1]
        #Check if data has Nans
        numberNan=np.sum(np.isnan(dataset))
        if numberNan>0: #fixing trace:
            for index in np.argwhere(np.isnan(dataset)):
                dataset[index]= dataset[index-1]  
                repaired+=1
        P=(np.abs(np.fft.ifft(dataset))**2)[0:int(N/2)] #Signal power per bin
        P=P/df #Signal power per Hz
        P=np.sqrt(P)
        spectra.append(P)
    if repaired !=0:
        print("Warning: %i data points with Nans in data! Attemp to repair trace."%repaired)
    spectra=np.array(spectra) 
    P_total=np.nanmedian(spectra, axis=0)
    factor=3
    aux=np.ones(factor)/factor
    f,P_total = np.convolve(f,aux), np.convolve(P_total,aux)
    ax.plot(f[:-factor],P_total[:-factor],color=color,linewidth=1,label=label)
    ax.set_yscale('log')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('PSD('+ data.Unit +')')
    ax.set_xscale('log')
    ax.legend(fontsize='small')
    return

"""
def inspect_Traces(kids=[],trace='wire'):
    """ function to inspect the 'noise' or 'wire' data for a specific KID.
    to be depracated"""
    print("function to be depracated")
    return
    plt.figure()
    plt.xlabel('samples')
    plt.ylabel('Phi [rad] decorrelated')
    for kid in kids:
        if trace =='noise':plt.plot(noise_scan[:,kid-1],label='Noise, kid%i'%kid)
        if trace =='wire' :plt.plot( wire_scan[:,kid-1],label='Wire,  kid%i'%kid)
    plt.legend()
    return

#def plot_Positions(wire_dict):
#    '''Function to plot the kid positions in array. '''
#    plt.figure()
#    kids_list=wire_dict.keys()
#    for k in kids_list:
#        kid=k
#        x,y        =wire_dict[k]['x0'],wire_dict[k]['y0'] 
#        plt.plot(x,y,'.r') #plot a red spot at pixel location
#        plt.annotate(kid,[x,y]) #annotate pixel numbers
#    plt.xlabel('X position [mm]')
#    plt.ylabel('Y position [mm]')
#    plt.grid()
#    plt.xlim(-300,300)
#    plt.ylim(-300,300)
#    plt.tight_layout()
#    plt.show()
#    return

def find_KID(wire_dict,kid):
    """    Return position x,y of a KID in the list
    """
    kids_list=wire_dict.keys()
    for k in kids_list:
        if k==kid:
            return wire_dict[k]['x0'],wire_dict[k]['y0']
    print('KID not in list')
    return 0,0

def find_XY(positions,xy):
    """kid=find_xy(positions,[x0,y0])
     Return KID at position x0,y0 in the list"""
    x0,y0=xy[0],xy[1]
    distance=300**2
    kid=[]
    kids_list=wire_dict.keys()
    for k in kids_list:
        dist=(x0-wire_dict[k]['x0'])**2+(y0-wire_dict[k]['y0'])**2   
        if dist < distance:
            distance, kid =dist,k
    return kid

def plot_correlParameters(wire_dict, parameter1='freq', parameter2='NET',chains=None):
    #clear current figure to avoid overplotting
    plt.clf()
    
    #get basic information from dictionary
    kids_list=wire_dict.keys()
    T_sky      =wire_dict[kids_list[0]]['T_sky']
    wire_d     =wire_dict[kids_list[0]]['wire_d']
    fe         =wire_dict[kids_list[0]]['fe'] 
    scan       =wire_dict[kids_list[0]]['scan']
    ID= fe + ', '+scan
    
    keys=wire_dict[kids_list[0]].keys()
    if parameter1 not in keys or parameter2 not in keys:
       print('Invalid parameter name, possible values are;')
       print(keys)
       return
    
    if chains==None:
        _,chains,n_kid=getFebe(fe)
    else:
        _,_,n_kid=getFebe(fe)
        	
    parameter3='freq'
    X,Y,Z=[],[],[]
    for k in kids_list:
        kid=k
        chain = kid/n_kid +1 
        if chain not in chains: continue
        x = wire_dict[k][parameter1]
        y = wire_dict[k][parameter2]
        if len(chains)==1: z = wire_dict[k][parameter3]
        else:              z= chain 
        X.append(x)
        Y.append(y)
        Z.append(z)
    plt.scatter(X,Y, c=Z, marker='o', s=15, cmap='jet',alpha=0.9, edgecolors='black')
    plt.xlabel(parameter1)
    plt.ylabel(parameter2)
    X=sorted(X)
    Y=sorted(Y)
    n=len(X)
    #print(n)
    X=X[int(0.05*n):int(0.95*n)]
    Y=Y[int(0.05*n):int(0.95*n)]
    #print(min(Y),max(Y))
    plt.xlim(np.mean(X) - 5*np.std(X),np.mean(X) + 5*np.std(X))
    plt.ylim(np.mean(Y) - 5*np.std(Y),np.mean(Y) + 5*np.std(Y))
    cbar=plt.colorbar()
    if len(chains)==1:
        plt.title(ID+ ', chain%i'%chains[0])
        cbar.set_label(parameter3)
    else:
        plt.title(ID+ 'selected chains')
        cbar.set_label('Chain number')
    plt.grid(color='gray', linestyle='-.',linewidth=1)
    return

def __WireScannerParameters():
    """Internal MARS function. It returns the wire scanner parameters as installed at APEX"""
    theta= -33.7 #Angle of coordinate system regarding antenna coordinate system
    X,Y= 30,47  #coordinates of scanner center to antenna optical axes
    return theta,X,Y
        
def plot_WireParameter(scan, parameter='NET',fe='LFA' ,  chains=None, scale=None, RCPxy=True):
    '''Function to plot the main parameters of a list of kids in positions.
    The index 1 or 2 indicate which direction is used, if no number we use the average
    parameter= 'ellipticity' , 'noise', 'SN', 'size' ,'signal','NET','max' ,'responsivity'
    using 1 or 2 will display the parameter in just 1 direction (x or y). example 'size1' or 'NET1' 
chains= iterable s is used to display results on a specific chain
scale=  can be 'auto' for 0-3*median(data), None for default, or iterable form minvalue and maxvalu
RCPxy= True: transform the measured x,y values into RCP standard for easy of comparison 
       False: Display measured XY positions
    '''
    wire_dict=reduceWireScan(scan,fe=fe,doReturn=True)
    #get basic information from dictionary
    kids_list=wire_dict.keys()
    T_sky      =wire_dict[kids_list[0]]['T_sky']
    T_cabin_C  =wire_dict[kids_list[0]]['T_cabin_C']
    wire_d     =wire_dict[kids_list[0]]['wire_d']
    fe         =wire_dict[kids_list[0]]['fe'] 
    scan       =wire_dict[kids_list[0]]['scan']
    ID= fe + ', '+scan
    
    if chains==None:
        _,chains,n_kid=getFebe(fe)
    else:
        _,_,n_kid=getFebe(fe)	
        
    if str.upper(fe) == 'LFA':
            nu = 345.0
    if str.upper(fe) == 'HFA':
            nu = 850.0

    if parameter=='ellipticity' : 
        min_value,max_value, unit = 1.0, 3.0, '' #default is 2, 3 for misaligned cases
    elif parameter=='noise' or parameter=='noise1' or parameter=='noise2': 
        if fe=='LFA': min_value,max_value, unit = 0 , 35, '$\Delta$ KHz / Hz**0.5'
        else: min_value,max_value, unit = 0 , 35, '$\Delta$ KHz / Hz**0.5' 
    elif parameter=='SN' or parameter=='SN1' or parameter=='SN2'         : 
        min_value,max_value, unit = 0, 5000, ''
    elif parameter=='size' or parameter=='size1' or parameter=='size2'   : 
        min_value,max_value, unit = 2, 24.0, 'mm' #default is 12mm, 24 for very misaligned cases
    elif parameter=='signal' or parameter=='signal1' or parameter=='signal2':    
        if fe=='LFA': min_value,max_value, unit = 0, 10, 'MHz*mm' 
        else: min_value,max_value, unit = 0, 4, 'MHz*mm'   
    elif parameter=='NET1' or parameter=='NET' or parameter=='NET2'      : 
        min_value,max_value, unit = 0, 6, 'mK * s**0.5 @ Tsky=%iK'%T_sky
    elif parameter=='max' or parameter=='max1' or parameter=='max2'      : 
        if fe=='LFA': min_value,max_value, unit = 0, 1, 'KHz'
        else: min_value,max_value, unit = 0, 0.3, 'KHz'
    elif parameter=='responsivity': 
        if fe=='LFA': min_value,max_value, unit = 0, 30, 'MHz/K'
        else: min_value,max_value, unit = 0, 100, 'MHz/K'
    elif parameter=='freq':
        min_value,max_value, unit = 3, 9, 'GHz'
    else: 
        print("Error: Parameter does not exist, read help for list of valid parameters")
        return
        

    fig,ax=plt.subplots(1,2,figsize=(14, 6),gridspec_kw={'width_ratios':[2,3]})
    fig.suptitle(ID + ', '+ parameter+' ['+unit+']' ,fontsize=16)
    ax[1].axis('square')
    values=[]
    X,Y,kids=[],[],[]
    for k in kids_list:
        kid=k
        chain = kid/n_kid +1 
        if chain not in chains: continue
        x,y        =wire_dict[k]['x0'],wire_dict[k]['y0'] #pos[1][0],pos[1][1]
        size1,size2=wire_dict[k]['sizeX'],wire_dict[k]['sizeY'] #pos[2][0],pos[2][1]
        n1,n2      =wire_dict[k]['noiseX'],wire_dict[k]['noiseY'] #pos[3][0],pos[3][1]
        sig1,sig2  =wire_dict[k]['signalX'],wire_dict[k]['signalY'] #pos[4][0],pos[4][1]
        freq       =wire_dict[k]['freq']/1E3 #pos[5]
        e=max(size1/size2,size2/size1)
        
        #Rotate x,y to show in coherent way with RCP based functions:
        theta, _, _ = __WireScannerParameters()
        theta=theta*np.pi/180
        x2=-x*np.cos(theta) - y*np.sin(theta)
        y2=-x*np.sin(theta) + y*np.cos(theta)
        x,y=x2,y2
        if parameter=='ellipticity'        : value = e
        elif parameter=='noise' : value = (n1+n2)*500
        elif parameter=='noise1': value = n1*1000
        elif parameter=='noise2': value = n2*1000
        
        elif parameter=='SN' : value = (sig1/n1+sig2/n2)/2
        elif parameter=='SN1': value = sig1/n1
        elif parameter=='SN2': value = sig2/n2

        elif parameter=='size' : value = (size1+size2)/2
        elif parameter=='size1': value = size1
        elif parameter=='size2': value = size2
        
        elif parameter=='signal' : value = (sig1*size1+sig2*size2)/2
        elif parameter=='signal1': value = sig1*size1
        elif parameter=='signal2': value = sig2*size2
        
        elif parameter=='max' : value=max(sig1,sig2)
        elif parameter=='max1': value = sig1
        elif parameter=='max2': value = sig2

        elif parameter=='NET' : 
            diff = B(nu,T_cabin_C)-B(nu,T_sky)
            dT=TRJ(nu,diff)
            value= (0.5*((dT)*1E3*n1/((sig1*size1)/wire_d) + (dT)*1E3*n2/((sig2*size2)/wire_d)))
        elif parameter=='NET1': 
            diff = B(nu,T_cabin_C)-B(nu,T_sky)
            dT=TRJ(nu,diff)
            value= (dT)*1E3*n1/((sig1*size1)/wire_d) 
        elif parameter=='NET2': 
            diff = B(nu,T_cabin_C)-B(nu,T_sky)
            dT=TRJ(nu,diff)
            value= (dT)*1E3*n2/((sig2*size2)/wire_d)
        elif parameter=='responsivity': 
            diff = B(nu,T_cabin_C)-B(nu,T_sky)
            dT=TRJ(nu,diff)
            value= (dT)/(((sig1)*size1)/wire_d)
            value=1/value
        elif parameter=='freq':
            value=freq
        else: value=0
        X.append(x)
        Y.append(y)
        kids.append(kid)
        values.append(value)
    if wire_d==0.65: markersize=10
    else:markersize=30

    if scale=='auto':
        min_value,max_value=0, 3*np.median(np.array(values))
    elif scale=='default' or scale==None:
        pass
    else:
        min_value,max_value=scale[0],scale[1] 
        
    z_plot=ax[1].scatter(X,Y, c=values,marker='o', s=markersize, cmap='jet',alpha=0.9, vmin=min_value, vmax=max_value, edgecolors='black')
    cbar = plt.colorbar(z_plot,ax=ax[1])
    ax[1].set(xlabel='X position [mm]', ylabel='Y position [mm]')
    ax[1].grid()
    if RCPxy:
       ax[1].set(xlim=[-330,280] , ylim=[-280,330])
       ax[1].annotate('2',(-250,250),color='r',va='center',ha='center')
       ax[1].annotate('1',(-250,-150),color='r',va='center',ha='center')
       ax[1].annotate('4',(150,-200),color='r',va='center',ha='center')
       ax[1].annotate('3',(250,250),color='r',va='center',ha='center')
    else:
       ax[1].set(xlim=[-280,280] , ylim=[-280,280])
       ax[1].annotate('3',(-250,250),color='r',va='center',ha='center')
       ax[1].annotate('4',(-250,-250),color='r',va='center',ha='center')
       ax[1].annotate('1',(250,-250),color='r',va='center',ha='center')
       ax[1].annotate('2',(250,250),color='r',va='center',ha='center')

    values=np.array(values)
    ax[1].set_title( 'N=%i'%len(values), fontsize=12)
    values2=values[(values>min_value)&(values<max_value)]
    counts, bins = np.histogram(values2,bins=40,range=(min_value,max_value))
    ax[0].hist(bins[:-1], bins, weights=counts)
    ax[0].set_title( 'N=%i, median=%.1f, mean=%.1f, std=%.1f'%(len(values2), np.median(values2), np.mean(values2), np.std(values2)), fontsize=12)
    ax[0].set(xlabel=parameter + ' ' + '['+unit+']')

    #if parameter[0:3]=='NET':
    #    Z=np.sort(values2)**2
    #    aux=[]
    #    for i in range(1,len(Z)):
    #        aux.append(np.sqrt(np.sum(Z[0:i]))/i)
    #    index=np.argmin(aux)
    #    N_eff=index
    #    NET_eff=aux=aux[index]*np.sqrt(index+1)
    #    ax[0].plot([NET_eff,NET_eff],[0,max(counts)],color='blue')
    #    ax[0].annotate('NET_eff=%.1f mK sqrt(s) \nN_eff=%i'%(NET_eff,N_eff),
    #                    (NET_eff+0.5,0.95*max(counts) ) )
        
    plt.savefig(parameter)
    plt.show()
    return
   
def print_Positions(positions,kids=[]): 
    print('Function to be deprecated, was used during development')
    return
    if True: 
        for pos in positions:
            if kids == []: 
                print('Kid=%i, pos=(%i,%i)[mm] , W=(%.1f,%.1f) [mm], Noise=%.3f,%.3f [deg/Hz**0.5], Max_signal=%.1f,%.1f [deg]'%(pos[0],pos[1][0],pos[1][1],pos[2][0],pos[2][1],pos[3][0],pos[3][1],pos[4][0],pos[4][1]))
            else:
                for kid in kids:
                    if kid==pos[0]: print('Kid=%i, pos=(%i,%i)[mm] , W=(%.1f,%.1f) [mm], Noise=%.3f,%.3f [deg/Hz**0.5], Max_signal=%.1f,%.1f [deg]'%(pos[0],pos[1][0],pos[1][1],pos[2][0],pos[2][1],pos[3][0],pos[3][1],pos[4][0],pos[4][1]))
            
def find_good(positions2):
    '''Function to filter kids'''
    positions3=[]
    for pos in positions2:
        kid=pos[0]
        x,y        =pos[1][0],pos[1][1]
        size1,size2=pos[2][0],pos[2][1]
        n1,n2      =pos[3][0],pos[3][1]
        sig1,sig2  =pos[4][0],pos[4][1]
        if n1==0 or n2==0: continue
        SN=min(np.abs(sig1/n1) ,np.abs(sig2/n2))
        if SN>20: positions3.append(pos)
    return positions3

def inspect_IQ(scan, fe='LFA', kids=[]):
    '''Function to inspect IQ data of a scan. plot all kids in list kids
    If scan is an fsweep use subscan 3 (Noise scan)'''
    
    print("Function to be deprecated, please inform if you need to use it")
    return
    
    # get IQ data of scan
    if scanIsFsweep(scan,fe=fe):getIQ(scan,sub=3,fe=fe)
    else: getIQ(scan,fe=fe)
    #Get timestamps
    t = data.getChanData('mjd', flag='None', getFlagged=0)
    fig,ax=plt.subplots(3,1,figsize=(9,12))
    ax[0].set_title('IQ data, Scan %i'%scan)
    maxim=0
    for kid in kids:
        trace=data.Data[:,kid-1]
	maxim=max(maxim,max(np.abs(trace)))
	ax[0].plot(np.real(trace),np.imag(trace),label=' %i'%kid)
	ax[1].plot(t,np.abs(trace))
        ax[2].plot(t,np.angle(trace))
    ax[0].set(xlim=[-maxim,maxim],ylim=[-maxim,maxim])
    #ax[0].axis('square')
    ax[0].set(xlabel='I [a.u.]',ylabel='Q [a.u.]')
    ax[1].set(xlabel='Time [s]', ylabel='Amplitude [a.u]')
    ax[2].set(xlabel='Time [s]', ylabel='Phase [rad]')
    if len(kids)<10: ax[0].legend()
    return fig,ax


def makeSkyRCP(scans=[], fe='LFA',frame='sky',distlim=None):
    """Function to generat an RCP file from beamscanner measurement"""
    
    if isinstance(scans,int):
        scans=[scans]
        
            
    wireDicts=[]
    for scan in scans:
        wireDicts.append(reduceWireScan(scan,fe=fe,doReturn=True))
    
    combinedDict=combineWireScans(wireDicts=wireDicts,fe=fe,param='signal',distlim=distlim)
    KidsWithConistentPos=[]
    posX=[]
    posY=[]
    gain=[]
    for kid in combinedDict:
        KidsWithConistentPos.append(kid)
        posX.append(combinedDict[kid]['x0'])
        posY.append(combinedDict[kid]['y0'])
        gain.append(combinedDict[kid]['signal'])

    gain=np.array(gain)/np.mean(np.array(gain))

    if str.upper(frame) == 'SKY':
        rcpposX,rcpposY=RCP_transform(np.array(posX),np.array(posY))
    else:
        rcpposX=np.array(posX)
        rcpposY=np.array(posY)
    
    nr = len(KidsWithConistentPos)
    #Initialize variables for RCP file    
    rcp_nr = np.zeros(nr,'int')
    rcp_x = np.zeros(nr,'f')
    rcp_y = np.zeros(nr,'f')
    rcp_g = np.zeros(nr,'f')
    for i,kid in enumerate(KidsWithConistentPos):
        rcp_nr[i] = kid
        rcp_x[i] = rcpposX[i]
        rcp_y[i] = rcpposY[i]
        rcp_g[i] = gain[i]
        
    #norm rcp gains
    rcp_g = rcp_g/np.nanmean(rcp_g)
    
    #Write file 
    ScanNr = scans[0]
    if str.upper(frame) == 'SKY':
        outname = "AMKID-"+str.upper(fe)+"-"+str(ScanNr)+"_WireScan_sky.rcp"
    else:
        outname = "AMKID-"+str.upper(fe)+"-"+str(ScanNr)+"_WireScan_cabin.rcp"
    f=open(outname,'w')
    f.write("# Feed   Gain  Gain  X [\"]    Y [\"]\n")
    for i in range(len(rcp_nr)):
        f.write("%5i %5.3f %5.3f %7.2f %7.2f\n"%(rcp_nr[i], rcp_g[i], rcp_g[i], rcp_x[i], rcp_y[i]))
    f.close()
    print("combind rcp writen to %s"%(outname)) 
    #Copy file to list of valid RCPs
    #command="cp "+str(outname)+" "+str(BoaConfig.rcpPath)
    #os.system(command)
    #updateRCP(outname)
    #plotArray()
    return    

def combineWireScans(wireDicts=[],fe='LFA',distlim=None):
    """Function to combine information from several wirescanner dictionaries. 
    Pixel information is just merged when they have spatial coherence (same KID 
    name on same XY position)
    """
    if len(wireDicts) == 0:
        print("no wireDicts given")
        return combinedDict

    if len(wireDicts) == 1:
        print("only a single wireDict given, nothing to combined")
        return wireDicts[0]
    
    if not distlim:
        if str.upper(fe) == 'LFA':
            distlim = 15. #distlim is in mm, conversion to arcsec ~ x2, so default limit for LFA is 30''
        if str.upper(fe) == 'HFA':
            distlim = 10. 
    combinedDict = {}


    KidsWithConistentPos=[]
    if len(wireDicts) > 1:
        for i in range(len(wireDicts)):
            refDict=wireDicts[i]
            for j in range(len(wireDicts)):
                if i != j:
                    compDict=wireDicts[j]
                    for kid in refDict.keys():
                        if kid in compDict.keys():
                            dist=np.sqrt((refDict[kid]['x0']-compDict[kid]['x0'])**2+(refDict[kid]['y0']-compDict[kid]['y0'])**2)
                            if dist <= distlim and kid is not KidsWithConistentPos:
                                KidsWithConistentPos.append(kid)
                                
    for kid in KidsWithConistentPos:
        keys=wireDicts[0][kid].keys()
        #This dictionary store a list for each parameter.
        mydict={}
        for key in keys:
            if key == 'scan' or key=='fe':
                mydict[key]=''
            else:
                mydict[key]=[]
        #fill the empty list with values from the different dictionaries
        for Dict in wireDicts:
            for key in keys:
                if key == 'scan': 
                    mydict[key]+= ('_'+ Dict[kid][key])
                elif key == 'fe': 
                    mydict[key]=Dict[kid][key]
                else:                
                    mydict[key].append(Dict[kid][key])
                
        #now calculate the average value of each key        
        for key in keys:
            if key == 'scan' or key == 'fe': 
                mydict[key]=mydict[key]
            else:
                mydict[key] = np.nanmean(np.array(mydict[key]))
        #just add the kid to the final dictionary if beam size make sense
        if mydict['sizeX'] < distlim*2 and mydict['sizeY'] < distlim*2:
            combinedDict[kid] = mydict
    
    return combinedDict

def RCP_transform(x,y,p=None):
    #p=[2.20,-2.19,33.1, 217.4,-535.3] #fit on LFA rcp 56907 and wire scanner scan 56908 from Nov 2023
    #p=[-2.10,2.16,33.4,  59.8, -95.6] #fit on LFA RCP AMKID-LFA-12Nov2023-geoCenter.rcp to wire scanner 580 Apr2024
    if not p:
        p=[-2.15,2.20,33.7,  -1.2,-108.1] #fit on LFA RCP 16492_LFA_Mars.rcp to AMKID-LFA-16490_WireScan_cabin.rcp
    # scalingX: arc/mm, scalingY: arc/mm, rotation: rot_ang, translation: X0, Y0
    x_rcp,y_rcp=myrotate(p,x,y)
    return x_rcp,y_rcp


def Fit_RCP_transform(scanner_rcp,sky_rcp,fe='LFA',doPlot=1,doReturn=False):
    if isinstance(sky_rcp,int):
        BeamMapDict=readBeamMapResult(scan=sky_rcp_scan,fe=fe)
    if isinstance(scanner_rcp,int):        
        ScannerDict=readScannerResult(scan=scanner_rcp,fe=fe)
    if isinstance(sky_rcp,str):
        kidNr_sky,_,gain, offX_sky, offY_sky = openRCP(rcpfile=sky_rcp)
        BeamMapDict={}
        for i,kid in enumerate(kidNr_sky):
            BeamMapDict[kid]={'xoff': offX_sky[i],'yoff': offY_sky[i],'gain': gain[i]}
    if isinstance(scanner_rcp,str):        
        kidNr_scanner,_,_, offX_scanner, offY_scanner = openRCP(rcpfile=scanner_rcp)
        ScannerDict={}
        for i,kid in enumerate(kidNr_scanner):
            ScannerDict[kid]={'x0': offX_scanner[i],'y0': offY_scanner[i]}

    if isinstance(sky_rcp,dict):
        BeamMapDict = sky_rcp
    if isinstance(scanner_rcp,dict):
        ScannerDict=scanner_rcp

    common = []
    pos_scanner = []
    pos_sky = []
    zparam = [] #only using in case doPlot = 1
    for kid in BeamMapDict.keys():
        if kid in ScannerDict.keys() and np.abs(ScannerDict[kid]['x0']-45) < 200 and np.abs(ScannerDict[kid]['y0']-12) < 200:
            common.append(kid)
            pos_scanner.append([ScannerDict[kid]['x0'],ScannerDict[kid]['y0']])
            pos_sky.append([BeamMapDict[kid]['xoff'],BeamMapDict[kid]['yoff']])
            if doPlot == 1:
                if isinstance(scanner_rcp,int):
                    zparam.append(1000*ScannerDict[kid]['response'])
                    wedge_label = r'response [mrad/K]'  
                    value_label = r'Sky RCP'             
                    pixel_unit = 'mrad/K'
                    mykey='responce'
                    factor = 1000.
                    limitsZ=[10,150]
                else:
                    zparam.append(BeamMapDict[kid]['gain'])
                    wedge_label = r'relative gain'  
                    value_label = r'Sky RCP'            
                    pixel_unit = ''
                    mykey='gain'
                    factor = 1.0
                    limitsZ=[0.2,1.5]
 
    pos_sky = np.array(pos_sky).T
    pos_scanner = np.array(pos_scanner).T
    Xscanner=pos_scanner[0]
    Yscanner=pos_scanner[1]
    Xsky = pos_sky[0]
    Ysky = pos_sky[1]
    #print np.min(Xscanner),np.min(Yscanner),np.min(Xsky), np.min(Ysky),len(Xscanner),len(Xsky)
    guess=[-2.0,2.0,30.0,0.0,0.0] 
    fitparams,success = optimize.leastsq(fit_transform, guess[:], args=(Xscanner,Yscanner,Xsky,Ysky), maxfev=50000, xtol=1e-10) #[0]
    
    xrot,yrot=myrotate(fitparams,Xscanner,Yscanner)  
    dist=np.sqrt((Xsky-xrot)**2+(Ysky-yrot)**2)
    badind=np.where(dist > 3* np.std(dist))[0]
    goodind=np.where(dist <= 3* np.std(dist))[0]
    print("")
    print("Nr of reliable Kids positions: %i"%(len(dist[goodind])))
    print("Nr of unreliable Kids positions: %i"%(len(dist[badind])))
    print("Position accuracy: %5.1f arcsec"%(np.std(dist[goodind])))
    print("")

    print("Fitted parameters:")
    print("X scale arc/mm: %5.2f"%(fitparams[0]))
    print("Y scale arc/mm: %5.2f"%(fitparams[1]))
    print("Rotation: %4.1f deg"%(fitparams[2]))
    print("Shift: X0: %6.1f Y0: %6.1f arcsec"%(fitparams[3],fitparams[4]))
    print ""
    print "Update set of parameters for function RCP_transform:"
    print("p=[%4.2f,%4.2f,%4.1f,%6.1f,%6.1f]"%(fitparams[0],fitparams[1],fitparams[2],fitparams[3],fitparams[4]))
    
    if doPlot == 1:
        if str.upper(fe) == 'LFA':
            beamshape=16.0
            scale = 0.5
        else:
            beamshape=7.5
            scale=1.0
        
        showArrayPlot(xrot, yrot, zparam, beamshape,scale,
    				fileName=None,
    				wedge_label=wedge_label,
    				value_label=value_label,
    				pixel_ids=common,
    				pixel_label=r'KID',pixel_unit=pixel_unit,febe=fe,limitsZ=limitsZ)
        
    if doReturn == True:
        result = [fitparams[0],fitparams[1],fitparams[2],fitparams[3],fitparams[4]]
        return result
        
def myrotate(p,Xin,Yin):
    xscale = p[0]
    yscale = p[1]
    theta=np.radians(p[2])
    x0 = p[3]
    y0 = p[4]
    xrot = xscale*(np.cos(theta)*Xin-np.sin(theta)*Yin)+x0
    yrot = yscale*(np.sin(theta)*Xin+np.cos(theta)*Yin)+y0
    return xrot,yrot


def fit_transform(p,Xscanner,Yscanner,Xsky,Ysky):
    xrot,yrot = myrotate(p,Xscanner,Yscanner)
    xdiff = xrot-Xsky
    ydiff = yrot-Ysky
    return np.array(list(xdiff)+list(ydiff))


def reduceWireScan(scan, noise=None, fsweep=None, doDf=True,fe='LFA',chains=None,overwrite=False,doReturn=False,doPlot=True, fig=None):
    """Script to process a WireScanner and save the results in a dictionary.
    Results are also stored as human legible *.csv file locaded in 'CalFiles/'. 
    
    @ param BeamScanner: Scan number of beam scanner
    @ type  BeamScanner: int
    @ param noise: Scan number of the noise trace used for sensitivity calculation. 
    If None: a clean section of the BeamScanner is used
    @ type  Noise: int 
    @ param fsweep: Scan number of fsweep, if None: system look for the most recent one
    @ type  fsweep: int
    @ param fe: 'LFA' or 'HFA', array to be analyzed
    @ type  fe: string
    @ param chains: chains to be used in analysis. if None: all of them are used.
    @ param chains: iterable
    @ param overwrite: if True force to overwrite existing file
                       if False do not reprocess data
    @ param overwrite: Boolean
    return: Dictionary
    """
    
    if not scanIsBeamscan(scan,fe=fe):
        print('ERROR:Scan is not a wirescan')
        return {}
    fe = str.upper(fe)
    outfile = 'CalFiles/beam_scanner_%i_%s.csv' % (scan, str.upper(fe))
    if os.path.exists(outfile) == False or overwrite:
        #if not noise:
        #    noise = scan  #if a noise scan is not specified use the same scan as noise.
        result,wire_d,Tsky,T_cabin_C,status=__WireScan_full_Analysis(BeamScanner=scan,Noise=noise,fsweep=fsweep, doDf=doDf, fe=fe,chains=chains)
        if isinstance(result,int):
            return
        responseDict={}
        #print(scan)
        if str.upper(fe) == 'LFA':
            nu = 345.0
            nominalsize1 = 11.7956 #in case the wire postions are wrong also the size may be corrupted
            nominalsize2 = 10.3439 #we then use theses sizes that have been determined from the wire scans in July2025
            ID='Scan %i, LFA'%scan
        elif str.upper(fe) == 'HFA':
            nu = 850.0
            nominalsize1 = 11.7956*350./870. #scaled from LFA
            nominalsize2 = 10.3439*350./870.
            ID='Scan %i, LFA'%scan
        else:
            print("Unknown frontend. allowed options are LFA or HFA")
            return responseDict

            
        for i in range(len(result)):
            kid = np.int(result[i][0])
        
            x,y = result[i][1][0],result[i][1][1]
            size1,size2=result[i][2][0],result[i][2][1]
            if size1 == -1:
                size1=nominalsize1
            if size2 == -1:
                size2=nominalsize2

            n1,n2 = result[i][3][0],result[i][3][1]  # VWO: these two are equal for wirescanners (minimum section noise)
            sig1,sig2  =result[i][4][0],result[i][4][1]
            KID_freq=result[i][5]
            KID_power=result[i][6]

            ell=max(size1/size2,size2/size1)
            noise = (n1+n2)/2. #in rad -> VWO: in KHz/sqrt(Hz) for wirescanner since doDf=True (default)
            SN = (sig1/n1+sig2/n2)/2.
            #signal = (sig1*size1+sig2*size2)/2. #Nico
            signal = (sig1+sig2)/2.

            diff = B(nu,T_cabin_C)-B(nu,Tsky)
            dT=TRJ(nu,diff)
            #NET =  (0.5*(dT*1E3*n1/((sig1*size1)/wire_d) + dT*1E3*n2/((sig2*size2)/wire_d))) #mK sqrt(s)
            # VWO: this computation introduced polarization-sensitive NETs by averaging NETs with different
            # wire signal responses (X-Y), but noise is not polarized

            r1 = sig1*size1/wire_d
            r2 = sig2*size2/wire_d
            r = (r1+r2)/2.
            phi2K= dT/r  # VWO: most of the time is dF2K since doDf=True (default). Change name everywhere to sig2K to generalize?

            NET = noise * phi2K * 1E3  # mK sqrt(s), VWO: equivalent to calculating noise in Kelvin-calibrated timelines
        
            ff = (size1/wire_d + size2/wire_d)/2.
            T_wire = dT/ff

            responseDict[kid] = {'T_sky':Tsky, 'T_cabin_C':T_cabin_C, 'T_wire': T_wire,'fe':fe, 'scan':str(scan), 'wire_d':wire_d, 'x0':x, 'y0':y, 'sizeX':size1, 'sizeY':size2, 'noiseX':n1, 'noiseY':n2, 'signalX':sig1,'signalY':sig2, 'freq':KID_freq, 'power':KID_power, 'ellipticity':ell, 'noise':noise, 'signal':signal, 'NET':NET, 'response':1./phi2K, 'phi2K': phi2K}
        
        if len(responseDict.keys()) > 0:
            writeScannerResult(responseDict)
        else:
            print('Warning: Wire scanner reduction failed - empty response dictionary')
        if status == -1:
            print('Warning: Could not determine positions and beam width from Wire Scanner')
    else:
        responseDict= readScannerResult(scan,fe)
        badposnr=0
        for kid in  responseDict.keys():
            if responseDict[kid]['x0'] == -9999:
                badposnr=badposnr+1
        if badposnr == 0:
            status = 1
        else:
            satus = -1
    if doReturn:
        return responseDict
    if doPlot:
        if len(responseDict.keys()) > 0:
            if status == 1:
                showWireScanResult(responseDict,fe=fe, fig=fig)
                print('Plotting positions and beam shape from Wire Scan ')
            else:
                showWireScanResult(responseDict,fe=fe,positions='b', fig=fig)
                
        

def writeScannerResult(responseDict):
    if os.path.exists("CalFiles") == False:
        os.makedirs("CalFiles")
    kid_list= responseDict.keys()
    
    fe=responseDict[kid_list[0]]['fe']
    scan=responseDict[kid_list[0]]['scan']
    outfile = 'CalFiles/beam_scanner_%s_%s.csv' % (scan, str.upper(fe))          
    with open(outfile,'w') as outFile:
        #copy common info from first KID
        kid=kid_list[0]
        outFile.write('#T_sky[K]=%6.2f, T_cabin_C[K]=%6.2f, T_wire[K]=%5.2f, fe=%s, scan=%s, wire_d[mm]=%f\n'%(responseDict[kid]['T_sky'], responseDict[kid]['T_cabin_C'], responseDict[kid]['T_wire'], responseDict[kid]['fe'], responseDict[kid]['scan'], responseDict[kid]['wire_d']))
                                                              
        outFile.write('#kidNr ,x0[mm], y0[mm] ,sizeX[mm], sizeY[mm], noiseX[rad/sqrt(Hz)], noiseY[rad/sqrt(Hz)], signalX[rad*mm], signalY[rad*mm], KID_freq[MHz], KID_power[dBm], ellipticity, noise[rad/sqrt(Hz)], signal[rad*mm], NET[mK sqrt(s)] ,response,phi2K[rad/K]\n')
        outFile.write('kid,x0,y0,sizeX,sizeY,noiseX,noiseY,signalX,signalY,freq,power,ellipticity,noise,signal,NET,response,phi2K\n')
        for kid in sorted(kid_list):
             outFile.write('%i,%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n' %(kid,
                                                              responseDict[kid]['x0'],
                                                              responseDict[kid]['y0'],
                                                              responseDict[kid]['sizeX'],
                                                              responseDict[kid]['sizeY'],
                                                              responseDict[kid]['noiseX'],
                                                              responseDict[kid]['noiseY'],
                                                              responseDict[kid]['signalX'],
                                                              responseDict[kid]['signalY'],
                                                              responseDict[kid]['freq'],
                                                              responseDict[kid]['power'],  
                                                              responseDict[kid]['ellipticity'],
                                                              responseDict[kid]['noise'],
                                                              responseDict[kid]['signal'],
                                                              responseDict[kid]['NET'],
                                                              responseDict[kid]['response'],
                                                              responseDict[kid]['phi2K']))
        outFile.close()
    print("Created new Wire scanner file: %s"%outfile) 
    print("Wire Scanner dictionary contains %i Kids"%(len(responseDict.keys())))

    return

def readScannerResult(scan,fe='LFA',account=None):
    if account:
        infile = 'CalFiles/%s/beam_scanner_%i_%s.csv' % (str.upper(account),scan, str.upper(fe))
    else:
        infile = 'CalFiles/beam_scanner_%i_%s.csv' % (scan, str.upper(fe))
    responseDict = {}
    dictKeys=[]
    if os.path.exists(infile) == True:
        with open(infile,'r') as inFile:
            reader = csv.reader(inFile)
            i=1
            for row in reader:
                if row[0][0:2]=='#T':
                    T_sky=float(row[0].split('=')[1])
                    T_cabin_C=float(row[1].split('=')[1])
                    T_wire=float(row[2].split('=')[1])
                    fe=row[3].split('=')[1]
                    scan=row[4].split('=')[1]
                    wire_d=float(row[5 ].split('=')[1])
                elif row[0][0]=='#':
                    continue
                else:  
                    if i == 1:
                        for j in range(len(row)):
                            if j > 0:
                                dictKeys.append(row[j])
                    else:
                        kid= np.int(row[0])
                        aux_dict = {}
                        for j in range(len(dictKeys)):
                            aux_dict[dictKeys[j]] = np.float(row[j+1])
                        if T_sky: aux_dict['T_sky']=T_sky
                        if T_cabin_C: aux_dict['T_cabin_C']=T_cabin_C
                        if fe   : aux_dict['fe']=fe
                        if scan : aux_dict['scan']=scan
                        if T_sky: aux_dict['wire_d']=wire_d
                        responseDict[kid] = aux_dict 
                    i=i+1
    else:
        messages.info('%s does not exist'% infile)
    return responseDict

def getScannerFile(scan='',fe='LFA'):
    if not scan:
        scan=data.ScanParam.ScanNum
    if not fe:
        if data.BolometerArray.FeBe == 'AMKID870-AMKID870BE':
            fe = 'LFA'
        if data.BolometerArray.FeBe == 'AMKID350-AMKID350BE':
            fe = 'HFA'
    search = "CalFiles/beam_scanner_*_%s.csv"%(str(fe))
    globlist=glob(search)
    mindiff = 1e6
    outfile = -1
    for file in globlist:
        calscan = np.int(file.split("_")[2])
        if scan >= calscan and scan-calscan < mindiff:
            outfile = file
            mindiff = scan-calscan
    if outfile == -1:
        messages.info('No Wire calibration found for %s'% fe)
    return outfile


def showWireScanResult(input=None,param='NET',fe='LFA',limitsZ=[],frame='sky',positions='W',doSave=0, fig=None):
    '''
    Function to plot the result from a wire scan
    By default the positions for each RX are shown in the Sky frame with 
    position from the last beam map

    @input: can be a scan number, a list of scan number or a wire scanner dictionary

    @param: Defines what is plotted color coded for each receiver.
            Implemented options are (case insensitive): 
            "NET" for the sensitivity in mK sqrt(s) (default)
            "Response" or "R" for response in radians/Kelvin or frequency_shift/Kelvin (depending on the reduction)
            "RCP" for the relative response 
            "Signal" or "S" for the signal in radians or frequency_shift (depending on the reduction)
            "Noise"  or "N" for the noise  in radians or frequency_shift (depending on the reduction)

    @positions: By default the position and beam shapes from the lastest wire scanner will be shown (positions="W"").
                One can use the positions and beam shapes from the latest beam map by using positions="B"
          
    @frame: defines if positions are plottend in the sky frame or the cabin frame
            allowed options are "SKY" (default, az-el) or "CABIN". "S" or "C" maybe used as well.
            The parameter is only effective when positions from the wirescanner are shown (positions="W")

    @limitsZ: limits for color range of the plotted quantity 
    @limitsZ type: list with min and max


    @doSave: save the plot to a pdf file or not
    @doSave type:  0 or 1 

    '''
    if input==None:
        print("no input provided: give a beam scanner dictionary, a scan number or a list of scan numbers")
        return
    else:
        if isinstance(input,int):
            scan=input
            wireDict=reduceWireScan(scan=scan,fe=fe,doPlot=False,doReturn=True)
        elif isinstance(input,list) and len(input) == 1:
            scans = input
            scan = scans[0]
            wireDict=reduceWireScan(scan=scan,fe=fe,doPlot=False,doReturn=True)
        elif isinstance(input,list) and len(input) > 1:
            scans = input
            scan=scans[0]
            wireDicts=[]
            for scan in scans:
                wireDicts.append(reduceWireScan(scan=scan,fe=fe,doPlot=False,doReturn=True))
            wireDict=combineWireScans(wireDicts=wireDicts,fe=fe)
        elif isinstance(input,dict):
            wireDict = input
            kid=wireDict.keys()[0]
            scan=int(wireDict[kid]['scan'])
            fe=str(wireDict[kid]['fe'])

    
    if str.upper(positions) == 'BEAMMAP' or str.upper(positions) == 'B': 
        beamfile=getBeamMapFile(scan=scan,fe=fe)
        if beamfile != -1:
            posDict=readBeamMapResult(infile=beamfile)
            info( 'Using positions and beam shape from  beam map %s'% beamfile)
        else:
            info( 'Warning: no %s beam map found, using positions from the wire scan'% fe)
            posDict = None
    else:
        posDict = None
    
    res=[]
    for kid in wireDict.keys():
        res.append(wireDict[kid]['response'])
    
    if np.median(res) < 0.1:
        unit = 'phi'
    else:
        unit = 'fshift'
    
    usedkids=[]
    offX=[]
    offY=[]
    fwhmX=[]
    fwhmY=[]
    ang=[]
    mypara=[]

    fe=str.upper(fe)
    factor = 1.

    title = None
    if str.upper(param) == 'R' or str.upper(param) == 'RESPONSE':
        mykey = 'response'
        value_label = r'Wire response'
        title="Wire response %s from scan %i"%(fe,scan)
        if unit == 'phi':
            wedge_label = r'response [mrad/K]'
            value_label = r'Wire response'
            pixel_unit = 'mrad/K'
            factor = 1000.
            if len(limitsZ) == 0:
                limitsZ=[10,150]
        else:
            wedge_label = r'response [kHz/K]'
            value_label = r'Wire response'
            pixel_unit = 'kHz/K'
            factor = 1.
            if len(limitsZ) == 0:
                limitsZ=[0,10]

        
    if str.upper(param) == 'NET':
        mykey = 'NET'
        wedge_label = r'NET [mK sqrt(s)]'
        value_label = r'NET'
        pixel_unit = 'mK sqrt(s)'
        title="%s NET from wire scan %i"%(fe,scan)
        factor = 1.
        if len(limitsZ) == 0:
            limitsZ=[0.3,4.]

    if str.upper(param) == 'RCP':
        mykey = 'signal'
        wedge_label = r'relative Gain'
        value_label = r'Gain'
        pixel_unit = ''
        title="%s RCP from wire scan %i"%(fe,scan)
        factor = 1.
        if len(limitsZ) == 0:
            limitsZ=[0.3,1.7]

    if str.upper(param) == 'SIGNAL' or str.upper(param) == 'S':
        mykey = 'signal'
        title="%s Signal from wire scan %i"%(fe,scan)
        value_label = r'Signal'
        if unit == 'phi':
            wedge_label = r'Signal [rad]'
            pixel_unit = 'rad'
            factor = 1.
            if len(limitsZ) == 0:
                limitsZ=[0.0,1.1]
        else:
            wedge_label = r'Signal [kHz]'
            pixel_unit = 'kHz'
            factor = 1.
            if len(limitsZ) == 0:
                limitsZ=[0.0,100]

    if str.upper(param) == 'NOISE' or str.upper(param) == 'N':
        mykey = 'noise'
        title="%s Noise from wire scan %i"%(fe,scan)
        value_label = r'Noise'
        if unit == 'phi':
            wedge_label = r'Noise [mrad x 10]'
            pixel_unit = 'mrad x 10'
            factor = 1000.*10.
            if len(limitsZ) == 0:
                limitsZ=[0.0,4.0]
        else:
            wedge_label = r'Noise [Hz]'
            pixel_unit = 'Hz'
            factor = 1000.
            if len(limitsZ) == 0:
                limitsZ=[0.0,50.]

    if posDict:
        for kid in posDict.keys():
            if kid in wireDict.keys():
                usedkids.append(kid)
                offX.append(posDict[kid]['xoff'])
                offY.append(posDict[kid]['yoff'])
                fwhmX.append(posDict[kid]['xfwhm'])
                fwhmY.append(posDict[kid]['yfwhm'])
                ang.append(posDict[kid]['tilt'])
                mypara.append(factor*wireDict[kid][mykey])
    else:
        for kid in wireDict.keys():
            usedkids.append(kid)
            if str.upper(frame) == 'CABIN' or str.upper(frame) == 'C':
                offX.append(wireDict[kid]['x0'])
                offY.append(wireDict[kid]['y0'])
            else:
                oX,oY=RCP_transform(np.array([wireDict[kid]['x0']]),np.array([wireDict[kid]['y0']]))
                offX.append(oX)
                offY.append(oY)
            fwhmX.append(wireDict[kid]['sizeX'])
            fwhmY.append(wireDict[kid]['sizeY'])
            ang.append(0.0)
            mypara.append(factor*wireDict[kid][mykey])

    if str.upper(param) == 'RCP':
        mypara = np.array(mypara)/np.mean(np.array(mypara))

    if len(limitsZ) == 0:
        limitsZ=[np.min(mypara),np.max(mypara)]
    
    beamshape=np.zeros((len(mypara),3),'f')
    beamshape[:,0]=np.array(fwhmX)
    beamshape[:,1]=np.array(fwhmY)
    beamshape[:,2] = np.array(ang)

    scale= 0.5
    
            
    
    showArrayPlot(offX, offY, mypara, beamshape,scale,
    				fileName=None,
    				wedge_label=wedge_label,
    				value_label=value_label,
    				pixel_ids=usedkids,
    				pixel_label=r'KID',pixel_unit=pixel_unit,febe=fe,limitsZ=limitsZ,title=title, fig=fig)


def writeWireScannerRCP(input=None,fe='LFA',outfile=None,skytransform=None):
    if isinstance(input,int):
        wD=readScannerResult(input,fe)
        scan=input
    if isinstance(input,dict):
        wD = input

    sky=""
    if skytransform:
        if isinstance(skytransform,int):
            p=None
        if isinstance(skytransform,list):
            p=skytransform

        for kid in wD.keys():              
            x=wD[kid]['x0']
            y=wD[kid]['y0']
            xr,yr=RCP_transform(x,y,p)
            wD[kid]['x0']=xr
            wD[kid]['y0']=yr
        sky="sky_"

    if not outfile:
        outfile='%i_%s_%swirescanner.rcp'%(scan,str.upper(fe),sky)

    responsenorm = 0
    for kid in wD.keys():
        responsenorm += wD[kid]['response']
    responsenorm = responsenorm/len(wD.keys())

    f=open(outfile,'w')
    for kid in wD.keys():
        #oX,oY=RCP_transform(np.array([wD[kid]['x0']]),np.array([wD[kid]['y0']]))
        gain=wD[kid]['response']/responsenorm
        f.write("%i %5.3f %5.3f %6.2f %6.2f\n"%(kid,gain,gain,wD[kid]['x0'],wD[kid]['y0']))
    f.close()
    
def CompareRCP_Wire(fe='LFA', masterFile=None,doPlot=False,doReturn=False,threshold=0.5):
    """+
    Compare current data RCP with last wirescanner.
    return list of kids with positions that do not match together
    parameters:
    threshold (float): if difference between RCP and Wire is larger than threshold*beamSize 
                       the KID is considered to have an inconsistent response.
    doPlot (boolean): If True a figure of  the offset between RCP and wire scanner is prepared.
    return:
    list of kids with inconsistent response
    """
    #get MASTER
    if not masterFile:
       Master_dict=__readMaster(fe=fe)
    else:
       Master_dict=__readMaster(masterFile,fe=fe)
    #get RCP
    scan=data.ScanParam.ScanNum
    if fe=='LFA':
       beamsize=17
    elif fe=='HFA':
       beamsize=8
    else:
        print('Error: fe must be HFA or LFA')
        return

    #get current RCP positions
    modifyRCP()
    usedBolos=data.BolometerArray.checkChanList([])
    usedBolos=usedBolos#[0:150]
    #get wire positions     
    wirescan=findWirescanForScan(scan,fe=fe)
    wire_dict=reduceWireScan(wirescan,fe=fe,doPlot=False,doReturn=True)
    kids=list(set(wire_dict.keys()) & set(Master_dict.keys()))
    counter=len(list(set(wire_dict.keys()) - set(Master_dict.keys())))
    if counter>0: print("Warning: %i KIDS not in MASTER"%counter)                            
    x_rcp , y_rcp = [],[]
    x_wire, y_wire= [],[]
    for kid in kids:
       x_wire.append(wire_dict[kid]['x0'])
       y_wire.append(wire_dict[kid]['y0'])      
       x_rcp.append(data.BolometerArray.Offsets[0,kid-1])
       y_rcp.append(data.BolometerArray.Offsets[1,kid-1])
       #if 'xoff' in Master_dict[kid].keys():
       #   x_rcp.append(Master_dict[kid]['xoff'])
       #   y_rcp.append(Master_dict[kid]['yoff'])
       #else:
       #   x_rcp.append(np.nan)
       #   y_rcp.append(np.nan)
     
    x_wire,y_wire=np.array(x_wire),np.array(y_wire)
    x_rcp ,y_rcp =np.array(x_rcp) ,np.array(y_rcp)
    
    x_rcp[x_rcp==-9999]  =np.nan
    y_rcp[y_rcp==-9999]  =np.nan
    x_wire[x_wire==-9999]=np.nan
    y_wire[y_wire==-9999]=np.nan
    
    par=[-33.7,2.2,30,47]
    #Note that 33.7 is the nominal angle of the wire scanner with antenna coorcdinate system
    #2.2 arcsec/mm is the nominal platescale of telescope
    #30,47 is the nominal offset of wire scanner to antenna coordinate system
    
    print(len(x_wire),len(x_rcp)) 
    def Transform(par):
       theta=par[0]*np.pi/180
       x2=-(x_wire-par[2])*np.cos(theta) - (y_wire-par[3])*np.sin(theta)
       y2=-(x_wire-par[2])*np.sin(theta) + (y_wire-par[3])*np.cos(theta)
       x3,y3=par[1]*x2,par[1]*y2
       d2=(x3-x_rcp)**2+(y3-y_rcp)**2
       aux=np.nanmedian(d2)
       return aux
    
    result = minimize(Transform,par,method='Powell')
    #print(result)
    par=result['x']
    
    theta=par[0]*np.pi/180
    x2=-(x_wire-par[2])*np.cos(theta) - (y_wire-par[3])*np.sin(theta)
    y2=-(x_wire-par[2])*np.sin(theta) + (y_wire-par[3])*np.cos(theta)
    x_wire,y_wire=par[1]*x2,par[1]*y2

    if doPlot:
       fig,ax=plt.subplots(1)
       ax.axis('square')             
       #ax.plot(x_wire,y_wire,'.',color='r')
       ax.plot(x_rcp,y_rcp,'.',color='b',label='RCP')
       ax.set(xlim=[-800,800],ylim=[-800.,800])
       ax.set(xlabel=r"$\Delta$ Az [arcsec]",ylabel=r"$\Delta$ El [arcsec]")
       ax.set(title='RCP vss wire, scan %i'%scan)
    flagged_kids=[]
    limit=threshold*beamsize**2
    for i in range(len(kids)):
       x=[x_rcp[i],x_wire[i]]
       y=[y_rcp[i],y_wire[i]]
       d2=(x[1]-x[0])**2 + (y[1]-y[0])**2
       if d2>limit:
          flagged_kids.append(kids[i])
       if doPlot: ax.plot(x,y,color='white',label=None)
    if doPlot: 
       ax.plot(x,y,color='white',label=r"$\Delta$ pos")
       ax.legend()
    print("Detected %i detectors with inconsistent positions between RCP and wirescanner %i"%(len(flagged_kids), wirescan))  
    if doReturn:    
        return flagged_kids
    else:
        return


def __findWire(windowtime=0.2, sig=5, expspikefree=85, ignoreblinds=True, full_output=False, doplots=False):
    '''
    CURRENTLY A COPY OF findspikes_IQBT WITH CHANGED DEFAULT PARAMS. WORKS GOOD!
    
    Finds spiked windowtime-long windows in BT-corrected I and Q data for each tone 
    and outputs a mask for data.Data of whether a timestamp should be used for each tone.
    Timestamps within a spiked window are all flagged for a given tone, even if not all
    timestamps are spiked. This is to ensure at least <windowtime> seconds of stable data.
    Whether all non-spiked windows have the same KID state is unknown: the algorithm can
    only detect jumps, not if the tone went back to the original KID and KID state for which
    it was calibrated.
    '''
    
    info('Beginning wire detection...')
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
        warn('Windowtime of %1.2f seconds yields too few windows. Changing windowtime to %1.2f seconds to ensure 50 windows...'%(windowtime, totaltime/nwindows))
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
    tone_floor_speedMEANs = []
    tone_floor_speedSTDs = []
    tone_thresholds_speed = []

    # inizialize data flagging array
    # if blindtones are ignored in process, flag is false for all blindtones so they will not be affected
    flagmask = np.zeros_like(data.Data, dtype=bool)

    info('Detecting wire...')
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
            tone_floor_speedMEANs.append(floor_speedMEAN)

            # compute the spike-free mean of speed STDs for tone (STD floor)
            floor_speedSTD = np.nanmean(windows_speed_std[valid, toneidx])
            tone_floor_speedSTDs.append(floor_speedSTD)

            # compute the mean + sig STD threshold
            threshold_speed = floor_speedMEAN + sig * floor_speedSTD
            tone_thresholds_speed.append(threshold_speed)

            # Tone has spike if ANY speed sample in window surpasses threshold
            hasspike = windows_speed_max[:, toneidx] >= threshold_speed

            # fill data mask array
            flagged_windows = np.where(hasspike)[0]
            for winidx in flagged_windows:
                t_start = windows_tstart[winidx]
                t_end = t_start + windowtime
                flagmask[:, toneidx] |= (time >= t_start) & (time < t_end)

            # fill window mask array
            windowflag[:, toneidx] = hasspike

        else:
            # Tone is only NaNs, then fill masks to leave un-edited (no spikes)
            tone_floor_speedMEANs.append(np.nan)
            tone_floor_speedSTDs.append(np.nan)
            tone_thresholds_speed.append(np.nan)
            flagmask[:, toneidx] = False
            windowflag[:, toneidx] = False
    
    # convert to arrays
    tone_floor_speedMEANs = np.array(tone_floor_speedMEANs)
    tone_floor_speedSTDs = np.array(tone_floor_speedSTDs)
    tone_thresholds_speed = np.array(tone_thresholds_speed)

    # final info
    spikedfraction_alltones = np.sum(windowflag, axis=0) / float(np.shape(windowflag)[0])
    spikedfraction_avg = np.nanmean(spikedfraction_alltones)
    spikedtime_avg = spikedfraction_avg * totaltime
    info('%.1f'%(spikedfraction_avg*100) + '% (avg.) of the KID timelines is the wire')

    if doplots:
        info('Initializing interactive plot...')

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
        fig, ax = plt.subplots(2, 2, figsize=(15, 8))

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

            # raw dF
            ax[0, 0].plot(time, data.Data.real[:, toneidx], label='Full dF')
            ax[0, 0].plot(time, despiked_Z_testtone.real, label='Wire-masked dF', zorder=1e9)
            ax[0, 0].set_xlabel('Time [s]')
            ax[0, 0].set_ylabel('dF [kHz]')
            ax[0, 0].legend(framealpha=1)

            # dF speed
            ax[0, 1].plot(auxtime, speed_testtone)
            ax[0, 1].plot(auxtime, despiked_speed_testtone)
            ax[0, 1].axhline(tone_floor_speedMEANs[toneidx], c='green', lw=2, label='Average spike-free tone speed')
            ax[0, 1].axhline(tone_thresholds_speed[toneidx], c='red', lw=2, label='Mean + %s Sigma'%sig)
            try:
                ax[0, 1].set_ylim(0, tone_thresholds_speed[toneidx]*1.5)
            except:
                pass
            ax[0, 1].set_xlabel('Time (s)')
            ax[0, 1].set_ylabel('dF/dt (kHz/s)')
            ax[0, 1].legend(framealpha=1)

            # dF speed window mean
            ax[1, 0].plot(windows_time, windows_speed_mean[:, testtone-1])
            ax[1, 0].plot(windows_time, despiked_windows_speed_mean_testtone)
            ax[1, 0].axhline(0)
            ax[1, 0].axhline(tone_floor_speedMEANs[toneidx], c='green', lw=2, label='Average spike-free tone speed')
            ax[1, 0].legend(framealpha=1)
            ax[1, 0].set_xlabel('Time (s)')
            ax[1, 0].set_ylabel('Windowed dF-speed mean (kHz/s)')

            # IQBT speed window mean
            ax[1, 1].plot(windows_time, windows_speed_std[:, testtone-1])
            ax[1, 1].plot(windows_time, despiked_windows_speed_std_testtone)
            ax[1, 1].axhline(0)
            ax[1, 1].axhline(tone_floor_speedSTDs[toneidx], c='red', lw=2, label='Average spike-free tone speed STD')
            ax[1, 1].legend(framealpha=1)
            ax[1, 1].set_xlabel('Time (s)')
            ax[1, 1].set_ylabel('Windowed dF-speed STD (kHz/s)')

            fig.suptitle('Wire detection for tone %i in scan %i'%(testtone, scannum))

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
                    testtone = int(userInput)
                except:
                    print('Tone number must be a number.')
                    continue

            if testtone not in usedtones:
                print('Tone %s not in used tones.'%testtone)

            ax[0, 0].cla()
            ax[0, 1].cla()
            ax[1, 0].cla()
            ax[1, 1].cla()

    if full_output:
        return flagmask, spikedfraction_alltones
    else:
        return flagmask
