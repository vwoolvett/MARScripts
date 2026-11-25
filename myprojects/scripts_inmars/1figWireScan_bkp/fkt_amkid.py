'''Description: ?
$LastChangedDate: 2026-09-07 08:52:06 +0000 (Mon, 07 Sep 2026) $
$Revision: 160 $
$LastChangedBy: nreyes 
$2024-09-04: showSig() was implemented using Pyplot and not PGplot.
$2024-10-10: showSig() handle IQ complex data, add units in plot

'''

from mars.fortran import fStat
from mars.fortran import fUtilities
import os, commands
from scipy import optimize
from scipy.ndimage import label
import numpy as np
from matplotlib import pylab as pl
import matplotlib as mpl
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle as rect
from glob import glob
import copy as copy
import warnings

execfile(os.getenv('BOA_HOME_AMKID') + '/rcp.py')

def get(scannr=0000, febe=None, doPlot=0, nrch=0, limitsY=[0,0], cal=0, sub=0, fe='LFA'):
    # fe was febe before, but as we standardize parameter names, we now use fe
    # febe is still accepted for backward compatibility
    # fe has precedence over febe in case both are given
    if febe is not None:
        fe = febe
        warnings.warn(DeprecationWarning("get(): the 'febe' parameter is depricated. Use 'fe' instead."))
        
    if fe.upper() == 'LFA':
        if sub == 0:
            read(str(scannr),febe='AMKID870-AMKID870BE')
        else:
            read(str(scannr),febe='AMKID870-AMKID870BE',subscans=[sub])
    else:
        if fe.upper() == 'HFA':
            read(str(scannr),febe='AMKID350-AMKID350BE')
        else:
            print "FEBE unknown- use HFE or LFE"
            return

    flagEmpty()
    checkNaN()
    #checkMJD()

    nr=len(data.ScanParam.MJD)
    nr = nr -1
    mjdref = (data.ScanParam.MJD[nr]-data.ScanParam.MJD[0])/2.+data.ScanParam.MJD[0]
    rcp=getMKIDsRCP(mjdref)
    print rcp
    #raw_input()
    updateRCP(rcp)
    if cal == 1:
        data.zeroStart()
        calibrateAMKID()
        #data.flagFractionRms(ratio=5)
        #correlChains(chanRef=-2,factor=0.95,nbloop=1)
#         data.flattenFreq(below=0.2,hiref=0.35)

    usedbolos = data.BolometerArray.checkChanList([])
    nrbolos = len(usedbolos)
    if nrbolos == 0:
        print "##############################################"
        print("No KIDs with data found in scan %i"%(scannr))
        print "##############################################"
        return
    else:
        print "##############################################"
        print("%i KIDs with data found "%(nrbolos))
        print "##############################################"

    if nrch == 0:
        if fe.upper() == 'HFA':
            nrchperplot = 300
        else:
            nrchperplot = 100
    else:
        nrchperplot = nrch
    nrplots = np.round(nrbolos/nrchperplot,0)+1
    nrplots = np.int(nrplots)
    k=0
    if nrbolos > 0 and doPlot== 1:
        for i in range(nrplots):
            if k+nrchperplot < nrbolos:
                end=k+nrchperplot
            else:
                end=nrbolos-1
            if limitsY[0] == 0 and limitsY[1] ==0:
                sig(range(usedbolos[k],usedbolos[end]))
            else:
                sig(range(usedbolos[k],usedbolos[end]),limitsY=limitsY)
            print("KIDS %i to %i"%(usedbolos[k],usedbolos[end]))
            k=k+nrchperplot
            raw_input()

def showSig_MPL(nrch=25, usedbolos=[], limitsX=None, limitsY=None, mydata=None, mydata2=None, fft=False, logX=False, logY=False, isSweep=False, chain=None, plotFlagged=True):
    """
    Show nrch signals simultaneously.

    @param nrch: Number of signals to show.
    @type nrch: int
    @param usedbolos: List of channels.
    @type usedbolos: list
    @param limitsX: Lower and upper x-limits.
    @type limitsX: list
    @param limitsY: Lower and upper y-limits.
    @type limitsY: list
    @param mydata: Data to display. If none is given, the global data will be used
               (Default).
    @type mydata: BoaDataEntity
    @param mydata2: Second data to display. fft and isSweep option not supported.
    @type mydata2: BoaDataEntity
    @param fft: Show the FFT of the signal (Default False)
    @type fft: boolean
    @param logX: FFT only: use log-scale for X.
    @type logX: boolean
    @param logY: FFT only: use log-scale for Y.
    @type logY: boolean
    @param isSweep: plot info from tonedictionary, not from data.Data
    @type isSweep: boolean
    @param chain: Chain to display. -1 = show all (Default).
    @type chain: int
    @param plotFlagged: plot flagged data or not
    @type plotFlagged: boolean
    @type: void
    """
    #chack parameter range:
    if nrch>100:
        print("Error: Cannot display more than 100 channels per plot")
        return
    #get data information
    global data    
    if mydata is None:
        mydata = data
    source=mydata.ScanParam.Object
    kids_per_chain = getFebe(getFrontend(mydata.BolometerArray.FeBe))[2]
    chains=getFebe(getFrontend(mydata.BolometerArray.FeBe))[1]
    febe = mydata.BolometerArray.FeBe
    fe= getFrontend(febe)
    scannr = mydata.ScanParam.ScanNum
    sub=mydata.ScanParam.SubscanNum
    if len(sub)>1: sub='%i-%i'%(sub[0],sub[-1])
    else: sub='%i'%sub[0]
    
    #check if data is real or complex
    real_data = True
    real_data2= True
    if  np.sum(np.iscomplex(mydata.Data)) != 0:
        real_data=False 
    if mydata2:
        if np.sum(np.iscomplex(mydata2.Data)) != 0:
           real_data2=False
        if real_data != real_data2 or fft:
           print("Warning: mydata2 is not compatible to be on same plot")
           mydata2=None 
    #check that the option isSweep can be used:
    if isSweep:
       try: 
          mydata.tone_dict   
       except:
          print("Warning: Data does not have tone dictionary => Skip option isSweep")
          isSweep=False
    #get user list of kids or chain
    if len(usedbolos) == 0:
        if isSweep:
            usedbolos = data.tone_dict.keys()
            if len(usedbolos)==0: 
                print("ERROR: Sweep option requires that sweep was reduced")
                return
        else:
            usedbolos = mydata.BolometerArray.checkChanList([])
    if chain == None:
        nrbolos = len(usedbolos)
    else:
        if chain not in chains:
            print("ERROR: chain %i not avalaible in data"%chain)
            return
        mylist = list(set(range(np.int((chain-1)*kids_per_chain+1), np.int(chain*kids_per_chain+1))).intersection(usedbolos))
        mylist.sort()
        usedbolos = mylist 
        nrbolos = len(usedbolos)
        
    if nrbolos == 0:
        print "##############################################"
        print("No KIDs with data found in scan %i." % (scannr))
        print "##############################################"
        return
    else:
        print "##############################################"
        print("%i KIDs with data found in scan %i."%(nrbolos, scannr))
        print "##############################################"
    
    #prepare plot layout    
    if nrch == 0:
            nrchperplot = 25
    N=int(np.ceil(nrch**0.5))
    fig,ax=plt.subplots(N,N,sharex='col',sharey='row') 
    if N==1: ax=np.array([[ax]])
    fig.set_size_inches(10,8)
    title='Scan:%i(%s) - %s - %s'%(scannr ,sub,fe, source)
    fig.text(0.4,0.98,title )
    if mydata2:
        scannr2=mydata2.ScanParam.ScanNum
        sub2=(mydata2.ScanParam.SubscanNum)
        if len(sub2)>1: sub2='%i-%i'%(sub2[0],sub2[-1])
        else: sub2='%i'%sub2[0]
        source2=mydata2.ScanParam.Object
        title2='Scan:%i(%s) - %s - %s'%(scannr2 ,sub2,fe, source2)
        fig.text(0.4,0.95,title2 ,color='yellow')  
    if isSweep and real_data:
        fig.text(0.5,0.05,'Delta Frequency [MHz]')
        fig.text(0.05,0.5,'Phi [rad]',rotation='vertical')
    elif fft and real_data:
        fig.text(0.5,0.05,'Frequency [Hz]')
        fig.text(0.05,0.5,'Intensity [?]',rotation='vertical')
    elif not real_data: 
        fig.text(0.5,0.05,'I - min(I) %s'%data.Unit)
        fig.text(0.05,0.5,'Q - min(Q) %s'%data.Unit,rotation='vertical')    
    else:
        fig.text(0.5,0.05,'MJD - MJD[0] [s]')
        fig.text(0.05,0.5,'Intensity [%s]'%data.Unit,rotation='vertical')

    fig.subplots_adjust(hspace=0,wspace=0)
            
    #fig.supxlabel('text X') just avalaible in Python 3.7
    #fig.supylabel('test Y') just avalaible in Python 3.7
    
    #print interactive message
    msg = '\n------------------------------------------------\n'
    msg += 'Plot next:                                <Enter>\n'
    msg += 'Plot previous:                          - <Enter>\n'
    msg += 'Goto KID Number:                 <number> <Enter>\n'
    msg += 'Exit script execution:                  q <Enter>\n'
    msg += '-------------------------------------------------\n'
    msg += 'Enter choice: '
    
    #start of plot function
    if  real_data: 
        mjd= mydata.ScanParam.MJD #np.array(mydata.getChanData('mjd',usedbolos[0]))
        time=(mjd-mjd[0])*24*3600
    if real_data and mydata2:
        mjd2= mydata2.ScanParam.MJD #np.array(mydata.getChanData('mjd',usedbolos[0]))
        time2=(mjd2-mjd2[0])*24*3600
    index=0
    while True:
        if index >= len(usedbolos): 
            print('Data inspection finished')
            break
        print('\nPlotting KIDs from %i to %i'%(usedbolos[index],usedbolos[index]+N**2) ) 
        #plot desired traces in subplots
        for i in range(N):
            for j in range(N):
                if index >= len(usedbolos): break
                kid=usedbolos[index]     
                dataind = mydata.BolometerArray.getChanIndex([kid])[0]
                trace=mydata.Data[:,dataind]
                if mydata2:
                   trace2=mydata2.Data[:,dataind]
                #If data is complex=> plot in complex plane
                #If data is real   => plot in time domain
                #If isSweep plot tonedictionary information and not data.Data
                if isSweep and real_data:
                    time=data.tone_dict[kid]['dfData']
                    trace=data.tone_dict[kid]['phData']     
                if  not real_data: 
                    a=np.nanmin(np.real(trace))
                    b=np.nanmin(np.imag(trace))
                    time =np.real(trace)-a#-np.real(np.nanmin(trace))
                    trace=np.imag(trace)-b#-np.imag(np.nanmin(trace))
                    if mydata2: 
                       time2 =np.real(trace2)-a#-np.real(np.nanmin(trace))
                       trace2=np.imag(trace2)-b#-np.imag(np.nanmin(trace))  
                    if isSweep and ('circleR' in data.tone_dict[kid].keys()):
                       Z_tune=data.tone_dict[kid]['z_tuningpoint']-complex(a,b)
                       Z_I= data.tone_dict[kid]['circleI']-a
                       Z_Q= data.tone_dict[kid]['circleQ']-b
                       R=data.tone_dict[kid]['circleR']
                       circle=plt.Circle((Z_I,Z_Q),R, color='grey', fill=False)
                       ax[i,j].add_patch(circle)
                       ax[i,j].plot([Z_I,np.real(Z_tune)],[Z_Q,np.imag(Z_tune)],color='grey')
                if fft:
                    plotFlagged==True  # VWO: do not touch 'time' and 'trace' vars
                    N_nans= np.sum(np.isnan(trace))
                    if N_nans !=0:
                        print("WARNING: %i nans in data for KID %i, attemp to fix it"%(N_nans,kid))
                        trace = trace[~np.isnan(trace)]
                        time  =  time[0:len(trace)]
                    trace=np.abs(np.fft.ifft(trace))[0:len(time)/2]
                    time=np.linspace(0,1.0/np.mean(time[1:]-time[:-1]),len(time))[0:len(time)/2]  

                # plot Y vs X
                if plotFlagged==True:
                    ax[i,j].plot(time,trace)
                    if mydata2: 
                        ax[i,j].plot(time2,trace2,color='yellow')
                    

                else:
                    isnotflagged = mydata.FlagHandler._aFlags == 0
                    flagmask = isnotflagged[:, dataind]
                    time_toplot = time[flagmask]
                    trace_toplot = trace[flagmask]
                    ax[i,j].plot(time_toplot,trace_toplot)
                    if mydata2:
                        isnotflagged2 = mydata2.FlagHandler._aFlags == 0
                        flagmask2 = isnotflagged2[:, dataind]
                        time2_toplot = time2[flagmask2]
                        trace2_toplot = trace2[flagmask2]
                        ax[i,j].plot(time2_toplot,trace2_toplot,color='yellow')
                    

                #ax[i,j].scatter(time,trace,s=10**2)
                ax[i,j].annotate('%i'%kid,(0.75,0.85),xycoords='axes fraction',color='red')
                ax[i,j].minorticks_on()
                if limitsX:
                    ax[i,j].set(xlim=limitsX)
                if limitsY:
                    ax[i,j].set(ylim=limitsY)
                if (fft or real_data) and logY:
                    ax[i,j].set(yscale='log')
                if (fft or real_data) and logX:
                    ax[i,j].set(xscale='log')
                index+=1             
        plt.show()  
        
        #get user input on how to continue
        userInput = raw_input(msg)
        userInput=str(userInput)
        if userInput == '':
            index = index
        elif userInput == '-':
            index = index-2*N**2
            if index <1: index=1 
        elif userInput == 'q':
            break
        else:    
            try:
                kid = int(userInput)
                index= list(usedbolos).index(kid)
            except:
                print('KID number must be a number and exist in array')
                index-=N**2
                
        #erase traces from subplots
        for i in range(N):
            for j in range(N):       
                ax[i,j].cla()
                
    #close figure at end of execution    
    plt.close(fig)
    return
        


def showSig_PGP(nrch=25, usedbolos=[], 
            limitsX=[0,0], limitsY=[0,0],
            mydata=None, mydata2=None, 
            fft=0, logX=1, logY=1, 
            chain=-1, interactive=True):
    """
    Function replaced by showSig(). Main difference is that showSig() now uses pyplot 
    instead of PG plot.
    
    Show nrch signals simultaneously.
    
    @param nrch: Number of signals to show.
    @type nrch: int
    @param usedbolos: List of channels.
    @type usedbolos: list
    @param limitsX: Lower and upper x-limits.
    @type limitsX: list
    @param limitsY: Lower and upper y-limits.
    @type limitsY: list
    @param mydata: Data to display. If none is given, the global data will be used 
                   (Default).
    @type mydata: BoaDataEntity
    @param mydata2: Second data to display.
    @type mydata2: BoaDataEntity
    @param fft: Show the FFT of the signal (Default False)
    @type fft: boolean
    @param logX: FFT only: use log-scale for X.
    @type logX: boolean
    @param logY: FFT only: use log-scale for Y.
    @type logY: boolean
    @param chain: Chain to display. -1 = show all (Default).
    @type chain: int
    @rtype: void
    
    """
    global data
    kids_per_chain = getFebe(getFrontend(data.BolometerArray.FeBe))[2]
    
    if mydata is None:
        mydata = data
    if mydata2 is None:
        mydata2 = data
        
    febe = mydata.BolometerArray.FeBe
    scannr = mydata.ScanParam.ScanNum
    if len(usedbolos) == 0:
        usedbolos = mydata.BolometerArray.checkChanList([])
        
    if chain == -1:
        nrbolos = len(usedbolos)
    else:
        print chain
        mylist = list(set(range(np.int((chain-1)*kids_per_chain+1), np.int(chain*kids_per_chain+1))).intersection(usedbolos))
        mylist.sort()
        usedbolos = mylist 
        nrbolos = len(usedbolos)
    
    if nrbolos == 0:
        print "##############################################"
        print("No KIDs with data found in scan %i." % (scannr))
        print "##############################################"
        return
    elif interactive==False:
        pass
    else:
        print "##############################################"
        print("%i KIDs with data found in scan %i."%(nrbolos, scannr))
        print "##############################################"
        
    if nrch == 0:
        if (febe == 'HFA') or (febe=='hfa'):
            nrchperplot = 25
        else:
            nrchperplot = 25
    else:
        nrchperplot = nrch
        
    nrplots = np.int(np.round(nrbolos/nrchperplot,0) + 1)
    
    k = 0
    if nrbolos > 0:
        for _ in range(nrplots):
            if k+nrchperplot < nrbolos:
                end=k+nrchperplot
            else:
                end=nrbolos-1
            
            if usedbolos[k] == usedbolos[end]:
                first=usedbolos[k]
                last=first+1
            else:
                first=usedbolos[k]
                last=usedbolos[end]

            if limitsX[0] == 0 and limitsX[1] ==0 and limitsY[0] == 0 and limitsY[1] ==0 :
                if fft == 0:
                    mydata.signal(range(first,last))
                    if mydata2 != data:
                        mydata2.signal(range(first,last),overplot=1,ci=2)
                else:
                    mydata.plotFFT(range(first,last),logY=logY,logX=logX)
                    if mydata2 != data:
                        mydata2.plotFFT(range(first,last),logY=logY,logX=0,ci=2,overplot=1)

            else:
                if limitsX[0] == 0 and limitsX[1] ==0:
                    if fft == 0:
                        mydata.signal(range(first,last),limitsY=limitsY)
                        if mydata2 != data:
                            mydata2.signal(range(first,last),limitsY=limitsY,ci=2,overplot=1)
                    else:
                        mydata.plotFFT(range(first,last),logY=logY,logX=logY,limitsY=limitsY)
                        if mydata2 != data:
                            mydata2.plotFFT(range(first,last),logY=logX,logX=logX,limitsY=limitsY,ci=2,overplot=1)

                if limitsY[0] == 0 and limitsY[1] ==0:
                    if fft == 0:
                        mydata.signal(range(first,last),limitsY=limitsY)
                        if mydata2 != data:
                            mydata2.signal(range(first,last),limitsX=limitsX,ci=2,overplot=1)
                    else:
                        mydata.plotFFT(range(first,last),logY=logY,logX=logX,limitsX=limitsX)
                        if mydata2 != data:
                            mydata2.plotFFT(range(first,last),logY=logY,logX=logX,limitsX=limitsX,ci=2,overplot=1)

                if (limitsX[0] != 0 or limitsX[1] !=0) and (limitsY[0] != 0 or limitsY[1] !=0):
                    if fft == 0:
                        mydata.signal(range(first,last),limitsX=limitsX,limitsY=limitsY)
                        if mydata2 != data:
                            mydata2.signal(range(first,last),limitsX=limitsX,limitsY=limitsY,ci=2,overplot=1)
                    else:
                        mydata.plotFFT(range(first,last),logY=logY,logX=logX,limitsX=limitsX,limitsY=limitsY)
                        if mydata2 != data:
                            mydata2.plotFFT(range(first,last),logY=logY,logX=logX,limitsX=limitsX,limitsY=limitsY,ci=2,overplot=1)

            if interactive==False:
               return
            print("KIDS %i to %i"%(first,last-1))
            k=k+nrchperplot
            raw_input('Press Ctrl + C to quit, or Enter to continue')


def OverPlot(chans=[],limitsX=[],limitsY=[],scaleY=1):
    for c in chans:
        mjd= np.array(data.getChanData('mjd',c))
        flux= np.array(data.getChanData('flux',c))
        if scaleY==1:
            if len(limitsX) == 2:
                refaver=np.mean(flux[(mjd>limitsX[0])*(mjd < limitsX[1])])
                refmin=np.min(flux[(mjd>limitsX[0])*( mjd < limitsX[1])])
                refmax=np.max(flux[(mjd>limitsX[0])*(mjd < limitsX[1])])
                px = mjd[(mjd>limitsX[0])*(mjd < limitsX[1])]
                py = flux[(mjd>limitsX[0])*(mjd < limitsX[1])]
            else:
                refaver=np.mean(flux)
                refmin=np.min(flux)
                refmax=np.max(flux)
                px = mjd
                py = flux
            #limitsY=[refmin,refmax]
        ind=chans.index(c)
        #for i in range(len(chans[ind:])):
        for i in range(len(chans)):    
            oc=chans[i]
            mjd2= np.array(data.getChanData('mjd',oc))
            flux2= np.array(data.getChanData('flux',oc))
            if scaleY==1:
                if len(limitsX) == 2:
                    aver2=np.mean(flux2[(mjd2>limitsX[0])*(mjd2 < limitsX[1])])
                    min2=np.min(flux2[(mjd2>limitsX[0])*(mjd2 < limitsX[1])])
                    max2=np.max(flux2[(mjd2>limitsX[0])*(mjd2 < limitsX[1])])
                    px2=mjd2[(mjd2>limitsX[0])*(mjd2 < limitsX[1])]
                    py2=flux2[(mjd2>limitsX[0])*(mjd2 < limitsX[1])]
                    py2 = (py2)*(refmax-refmin)/(max2-min2)
                    py2 = py2-np.mean(py2)+refaver
                else:
                    aver2=np.mean(flux2)
                    min2=np.min(flux2)
                    max2=np.max(flux2)
                    px2 = mjd2
                    py2 = flux2*(refmax-refmin)/(max2-min2)
                    py2 = py2-np.mean(py2)+refaver
            
            if len(limitsY) != 2 :
                plot(px,py,ci=1,style='b',labelX='MJD',labelY='Intensitiy')
            else:
                plot(px,py,ci=1,style='b',limitsY=limitsY,labelX='MJD',labelY='Intensitiy')
                
            if c != oc:                
                plot(px2,py2,ci=2,style='b',overplot=1)
                print "Kid %i on %i"%(oc,c)
                raw_input()
                
                
def showSig(nrch=25, usedbolos=[], 
            limitsX=[0,0], limitsY=[0,0],
            mydata=None, mydata2=None, 
            fft=0, logX=1, logY=1, 
            chain=-1, interactive=True):

    showSig_PGP(nrch=nrch, usedbolos=usedbolos, 
            limitsX=limitsX, limitsY=limitsY,
            mydata=mydata, mydata2=mydata2, 
            fft=fft, logX=logX, logY=logY, 
                chain=chain, interactive=interactive)

def computeBeamMap(tol=30.):
    bolos=data.BolometerArray.checkChanList([])
    index = data.BolometerArray.getChanIndex(bolos)

    orgref=data.BolometerArray.RefChannel
    refindex=data.BolometerArray.RefChannel-1
    refoffset=[data.BolometerArray.Offsets[0][refindex],data.BolometerArray.Offsets[1][refindex]]

    xoff=data.BolometerArray.Offsets[0][index]-refoffset[0]
    yoff=data.BolometerArray.Offsets[1][index]-refoffset[1]
    azmin,azmax=fStat.minmax(xoff)
    elmin,elmax=fStat.minmax(yoff)
    elmax=elmax+60.
    azmin=azmin-100
    center=[(azmax+azmin)/2.,(elmax+elmin)/2.]

    mindist=1e6
    centerbolo=-1

    for i in range(len(xoff)):
        dist=  np.sqrt((xoff[i]-center[0])**2+(yoff[i]-center[1])**2)
        if dist < mindist:
            mindist=dist
            centerbolo=bolos[i]

    data.BolometerArray.RefChannel=centerbolo
    refindex=data.BolometerArray.RefChannel-1
    refoffset=[data.BolometerArray.Offsets[0][refindex],data.BolometerArray.Offsets[1][refindex]]
    xoff=data.BolometerArray.Offsets[0][index]-refoffset[0]
    yoff=data.BolometerArray.Offsets[1][index]-refoffset[1]

    #now find horizontal angle for beam map OTF and OTF size
    bols1=np.extract((bolos <= 2*880),bolos)
    bols2=np.extract((bolos > 2*880),bolos)
    ind1=data.BolometerArray.getChanIndex(bols1)
    ind2=data.BolometerArray.getChanIndex(bols2)

    x1=fStat.f_mean(data.BolometerArray.Offsets[0][ind1]-refoffset[0])
    y1=fStat.f_mean(data.BolometerArray.Offsets[1][ind1]-refoffset[1])
    x2=fStat.f_mean(data.BolometerArray.Offsets[0][ind2]-refoffset[0])
    y2=fStat.f_mean(data.BolometerArray.Offsets[1][ind2]-refoffset[1])

    #print x1,y1,x2,y2
    oldref=data.BolometerArray.RefChannel
    data.BolometerArray.RefChannel=centerbolo
    angle=-1*np.arctan((y1-y2)/(x1-x2))
    xrot= xoff*np.cos(angle)-yoff*np.sin(angle)
    yrot= xoff*np.sin(angle)+yoff*np.cos(angle)
    mapxmin,mapxmax=fStat.minmax(xrot)
    mapymin,mapymax=fStat.minmax(yrot)
    xlen=(mapxmax-mapxmin)+tol
    ylen=(mapymax-mapymin)+tol
    timeLFA=(int(ylen/6.0)+1)*xlen/120.
    timeHFA=(int(ylen/2.5)+1)*xlen/120.


    xoff=np.array([mapxmin-tol,mapxmin-tol,mapxmax+tol,mapxmax+tol])
    yoff=np.array([mapymin-tol,mapymax+tol,mapymax+tol,mapymin-tol])
    
    xrot= xoff*np.cos(-1*angle)-yoff*np.sin(-1*angle)
    yrot= xoff*np.sin(-1*angle)+yoff*np.cos(-1*angle)
    
    xrot=xrot.tolist()
    yrot=yrot.tolist()
    xrot.append(xrot[0])
    yrot.append(yrot[0])
    
   
    plotArray(limitsX=[-700,700],limitsY=[-700,700])
    plot(xrot,yrot,style='l',overplot=1,ci=2)
    

    data.BolometerArray.RefChannel=orgref
  
    if centerbolo > 0:
        print "best center pixel: %i"%(centerbolo)
        print "beam map: otf(xlen=%5.1f,ylen=%5.1f,angle=%5.1f,system=\'ho\')"%(xlen,ylen,-1*angle*180./np.pi)
        print "observing time for LFA/HFA: %5.1f/%5.1f min"%(timeLFA/60.,timeHFA/60.)
    else:
        print "failed to determine best pixel for reasons beyond my control"
    
                 

def getSkyTemp():
    indc1 = data.ScanParam.SubscanIndex[0,0]
    indc2 = data.ScanParam.SubscanIndex[1,0]
    usedbolos = data.BolometerArray.checkChanList([])
    skytemp=np.zeros_like(usedbolos)
    for i in range(len(usedbolos)):
        flux= np.array(data.getChanData('flux',usedbolos[i]))
        skytemp[i]=np.mean(flux[indc1:indc2])
    
    
    febe=getFeBe()
    scan = data.ScanParam.ScanNum
    op('%i_%s_skyTemp.ps/CPS' % (scan,febe))
    plot(usedbolos,skytemp,labelX='Nr Kids',labelY='Tsky [K]',caption='Sky Temperature [K]',style='p')
    close()
    meanskytemp=fStat.f_mean(skytemp)
    medianskytemp=fStat.f_median(skytemp)

    print "Median/Mean Sky Temperature for %s: %f5.1/ %f5.1"%(febe,medianskytemp,meanskytemp)
    
     

def getFeBe():
    if data.BolometerArray.FeBe == 'AMKID870-AMKID870BE':
        febe='LFA'
    if data.BolometerArray.FeBe == 'AMKID350-AMKID350BE':
        febe='HFA'
    if febe != 'LFA' and febe != 'HFA':
        febe=data.BolometerArray.FeBe
    return febe
        
        
def showArrayPlot(az, el, values, beamshape, scale=0.5, fileName=None, 
					wedge_label=r'Sensitivity [ mJy sqrt(s) ]',
					value_label=r'Sensitivity',
					pixel_ids=[], pixel_label=r'Pixel',pixel_unit='mK sqrt(s)',febe='LFA',title=None,
                    limitsZ=[], highlightKids=[],addRCP=[], fig=None):
    '''
    Display and optionally save the array sensitivity.
    
    @param az: 1D array of azimuth positions.
    @type az: numpy.array
    @param el: 1D array of elevation positions.
    @type el: numpy.array
    @param sensitivity: 1D array of sensitivities - will be color coded.
    @type sensitivity: numpy.array
    @param beamshape: Either single value indicating same beamsize for all 
                     channels, or 1D array defining the circular beamsize for each channel
                     or a 3D array defining the ellipses for each channel
    @type beamsize: float or 1D or 3D numpy.array
    @paraum scale: value to scale the beam size for plotting - default is 0.5
    @tyape scale:  float
    @param fileName: (optional - default:'') When given, additionally, a file
                     will be generated.
    @type fileName: string
    @param addRCP: list=[Az,El,beamsize] with Az and El positions of an additional RCP to display and the beam size
    @type        : list
    @return: void
    
    '''
    if not len(az) == len(el) == len(values):
        msg = 'Az, El and sensitivity must have same length. '
        msg += '(Found |az|=%i, |el|=%i, |sensitivity|=%i)' % (len(az), len(el), len(values))
        raise RuntimeError(msg)
    
    if len(pixel_ids)==len(az):
        pix_ids = pixel_ids
    else:
        pix_ids = np.arange(len(az))

    beams=np.zeros((len(az),3),'f')
    if np.isscalar(beamshape):
        searchradius = beamshape/2. * 10
        beams[:,0]=beamshape
        beams[:,1]=beamshape
        beamell = 0
    elif len(np.shape(beamshape)) == 1:
        searchradius = np.max(beamshape)/2. * 10
        beams[:,0]=beamshape
        beams[:,1]=beamshape
        beamell = 0
    else:
        searchradius = np.max(beamshape[:,0])/2. * 10
        beams = beamshape
        beams[:,2]=-1*beams[:,2] #-90.0 #to match astronomical definition that north is angle = 0.0
        beamell = 1
    def _onClick(event):
        '''
        Helper function to handle the mouse click events and show the tooltip
        when the mouse is hovering a pixel when click occurs.
        
        @param event: Mouse event handler
        @type event: event
        @return: void
        
        '''
        if event.xdata != None and event.ydata != None: # mouse is inside the axes
            dAz = np.abs(event.xdata - az)
            dEl = np.abs(event.ydata - el)
            dist = dAz**2 + dEl**2
            try:
                i = np.where(dist==min(dist[dist<searchradius]))[0][0]
            except IndexError:
                txt.set_visible(False)
            except ValueError:
                txt.set_visible(False)
            else:
                txt.set_visible(True)
                chain=_getChain(pix_ids[i],fe=febe)
                if value_label == 'Gain':
                    if beamell == 0:
                        txt.set_text('%s: %i (Chain: %i)\n%s: %.3f %s'%(pixel_label, pix_ids[i],chain, value_label, values[i],pixel_unit))
                    else:
                        txt.set_text('%s: %i (Chain: %i)\n%s: %.1f %s\nbeam: %5.1f\'\' x %5.1f\'\''%(pixel_label, pix_ids[i],chain, value_label, values[i],pixel_unit,beams[i][0],beams[i][1]))
                else:
                    if beamell == 0:
                        txt.set_text('%s: %i (Chain: %i)\n%s: %.1f %s'%(pixel_label, pix_ids[i],chain, value_label, values[i],pixel_unit))
                    else:
                        txt.set_text('%s: %i (Chain: %i)\n%s: %.1f %s\nbeam: %5.1f\'\' x %5.1f\'\''%(pixel_label, pix_ids[i],chain, value_label, values[i],pixel_unit,beams[i][0],beams[i][1]))
                txt.set_position((event.xdata, event.ydata))
            
            fig.canvas.draw()
        
    # Set up the figure
    xdiff = max(az)-min(az)
    ydiff = max(el)-min(el)
    if fig==None:
        fig = pl.figure(figsize=(10, 10))
    else:
        fig.clf()
    
    # Set up the drawing axes
#     ax1 = fig.add_subplot(111, axisbg='black')
    ax1 = fig.add_subplot(111)
    ax1.set_xlim([min(az)-(max(az)-min(az))*0.1, max(az)+(max(az)-min(az))*0.1])
    ax1.set_ylim([min(el)-(max(el)-min(el))*0.1, max(el)+(max(el)-min(el))*0.1])
    ax1.set_xlabel(r'Azimuth Offset')
    ax1.set_ylabel(r'Elevation Offset')
    ax1.set_aspect(1.)
    
    # Set up the color bar
    # ...first the min/max values
    if len(limitsZ) == 2:
        vmin=limitsZ[0]
        vmax=limitsZ[1]
    else:
        valmin, valmax = np.min(values), np.max(values)
        meanval = np.mean(values)
        rmsval = np.std(values)
        vmin = meanval - 2*rmsval
        vmax = meanval + 2*rmsval
        
        #vmin, vmax = np.min(values), np.max(values)
        if len(limitsZ) != 2:
            if vmin < valmin:
                vmin = valmin
            if vmax > valmax:
                vmax = valmax
    
    # ...colorbar ticks
    tlist = []
    for i in range(11):
        value = vmin + (vmax-vmin) * 0.1 * i
        tlist.append(value)
    #print vmax,vmin
    # ...colorbar ticks format
    labelFormat = "%.0f"
    if vmax > 30:
        labelFormat = "%.0f"
    if vmax <= 30 and vmax > 1:
        labelFormat = "%.1f"
    if vmax <= 1:
        labelFormat = "%.2f"
    
    # ...and finally the colorbar itself.
    try:
        cmap = mpl.cm.spectral
    except:
        cmap = mpl.cm.plasma

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = pl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(values)
        
    cbar = fig.colorbar(sm, ticks=tlist, format=labelFormat, ax=ax1, fraction=0.06, pad=0.01, aspect=80)
    cbar.set_label(wedge_label)

    #plot Title:
    if title:
        titlestr = title
    else:
        if value_label != 'Gain':
            titlestr = str(value_label)+' '+str.upper(febe)
        else:
            titlestr = 'RCP '+str.upper(febe)
    
    ax1.set_title(titlestr,fontsize=16)

    # Plot the pixels
    for x, y, c, k, b in zip(az, el, values, pix_ids, beams):
        colorVal = sm.to_rgba(c)
        if k in highlightKids:
#             ax1.scatter(x,y, s=48., marker='x', c='#FFFF00', zorder=100, linewidth=1)
            ax1.scatter(x,y, s=50., marker='x', c='#FF6688', zorder=100, linewidth=2)
            ax1.scatter(x,y, s=55., marker='x', c='k', zorder=99, linewidth=4)
        #ax1.add_artist(Circle(xy=(x, y),
        #                      color=colorVal,
        #                      radius=beamradius*scale))
        ax1.add_artist(mpl.patches.Ellipse(xy=(x, y),width=b[0]*scale,height=b[1]*scale,angle=b[2],color=colorVal))

        
    #Plot the additional RCP if provided:
    if len(addRCP) > 0:
        for entry in addRCP:
            #print(entry)
            #for x, y, bmaj, bmin, bpa in zip(entry[0],entry[1],entry[2],entry[3],entry[4]):
            ax1.add_artist(mpl.patches.Ellipse(xy=(entry[0], entry[1]),width=entry[2]*scale,height=entry[3]*scale,angle=entry[4],fc=None,ec='white',fill=False))

    # Set up the tooltip
    txt = ax1.text(0., 0., '', color='black', 
                   bbox=dict(facecolor='#ffffff', edgecolor='black', pad=10.0, alpha=0.9))
    txt.set_visible(False)
    txtHint = ax1.text(min(az)-(max(az)-min(az))*0.1 - (min(az)-(max(az)-min(az))*0.1)*0.01,
                       min(el)-(max(el)-min(el))*0.1 - (min(el)-(max(el)-min(el))*0.1)*0.01,
                       'Click %s to show details'%pixel_label, color='#666666',
                       size='small')
    
    # Let the mouse hover...
#     fig.canvas.mpl_connect('motion_notify_event', _onMotion)
    fig.canvas.mpl_connect('button_press_event', _onClick)

    #if np.max(np.abs(az)) > 1000 or np.max(np.abs(el)) > 1000.:
    #    ax1.set(xlim=[-750,750])
    #    ax1.set(ylim=[-750,750])
    
    ax1.set(xlim=[-750,750])
    ax1.set(ylim=[-750,750])

    if value_label == 'Coupling':
        _,validchains,kidsperchain=getFebe(fe=febe)
        posX=[-670,-530,400,250]
        posY=[-370,500,350,-510]
        if str.upper(febe) == 'LFA':
            for i,chainNum in enumerate(validchains):
                kidRange=[(chainNum-1)*kidsperchain + 1, (chainNum)*kidsperchain]
                mask=np.where( (np.array(pixel_ids) >= kidRange[0]) *  (np.array(pixel_ids) <= kidRange[1]))[0]
                if len(mask) > 0:
                    coup=np.median(np.array(values)[mask])
                    mytest="Coupling LFA %i:"%(i+1)
                    ax1.text(posX[i],posY[i],mytest,ha='left')
                    mytest="%4.1f Jy/K "%(coup)
                    ax1.text(posX[i],posY[i]-35,mytest,ha='left')
        else:
            chips=[]
            for chainNum in validchains:
                refchains=[[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20]]
                for i,rch in enumerate(refchains):
                    if chainNum in rch:
                        if i+1 not in chips:
                            chips.append(i+1)
            
            
            for i,chip in enumerate(chips):
                kidRange=[(chip-1)*kidsperchain*5 + 1, chip*kidsperchain*5]
                mask=np.where( (np.array(pixel_ids) >= kidRange[0]) *  (np.array(pixel_ids) <= kidRange[1]))[0]
                if len(mask) > 0:
                    coup=np.median(np.array(values)[mask])
                    mytest="Coupling HFA %i:"%(i+1)
                    ax1.text(posX[i],posY[i],mytest,ha='left')
                    mytest="%4.1f Jy/K "%(coup)
                    ax1.text(posX[i],posY[i]-35,mytest,ha='left')

    if value_label == 'NET':
        _,validchains,kidsperchain=getFebe(fe=febe)
        posX=[-650,-530,400,250]
        posY=[-350,500,350,-510]
        if str.upper(febe) == 'LFA':
            for i,chainNum in enumerate(validchains):
                kidRange=[(chainNum-1)*kidsperchain + 1, (chainNum)*kidsperchain]
                mask=np.where( (np.array(pixel_ids) >= kidRange[0]) *  (np.array(pixel_ids) <= kidRange[1]))[0]
                if len(mask) > 0:
                    sen=np.median(np.array(values)[mask])
                    mytest="NET LFA %i:"%(i+1)
                    ax1.text(posX[i],posY[i],mytest,ha='left')
                    if sen >=10:
                        mytest="%4.1f mK sqrt(s)"%(sen)
                    else:
                        mytest="%3.1f mK sqrt(s)"%(sen)
                    ax1.text(posX[i],posY[i]-30,mytest,ha='left')
        else:
            chips=[]
            for chainNum in validchains:
                refchains=[[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20]]
                for i,rch in enumerate(refchains):
                    if chainNum in rch:
                        if i+1 not in chips:
                            chips.append(i+1)
            for i,chip in enumerate(chips):
                kidRange=[(chip-1)*kidsperchain*5 + 1, chip*kidsperchain*5]
                mask=np.where( (np.array(pixel_ids) >= kidRange[0]) *  (np.array(pixel_ids) <= kidRange[1]))[0]
                if len(mask) > 0:
                    sen=np.median(np.array(values)[mask])
                    mytest="NET HFA %i:"%(i+1)
                    ax1.text(posX[i],posY[i],mytest,ha='left')
                    if sen >=10:
                        mytest="%4.1f mK sqrt(s)"%(sen)
                    else:
                        mytest="%3.1f mK sqrt(s)"%(sen)
                    ax1.text(posX[i],posY[i]-30,mytest,ha='left')
    
    #RC=rect((-410-540/2.,75-540/2.),540,540,-33.7,fill=False)
    #ax1.add_artist(RC)
    pl.show()
    if fileName:
        fig.savefig(fileName)


def showRCP(rcpFilename='',fe='LFA',beamsize=0,limitsZ=[], highlightKids=[],addRCP=None,align=True):
    if rcpFilename != '':
        kidNr, gain, _, offX, offY = openRCP(rcpFilename)
    else:
        rcpFilename = getLastRCP()
        kidNr, gain, _, offX, offY = openRCP(rcpFilename)
        
        
    if np.max(kidNr) <= 3520:
        fe='LFA'
        if beamsize == 0:
            beamsize = 12.
    else:
        fe='HFA'
        if beamsize == 0:
            beamsize = 8.

    if addRCP:
        if isinstance(addRCP,str):
            addRCPlist=[addRCP]
        else:
            addRCPlist=addRCP

        secondRCP=[]
        if align == True:
            for rcp in addRCPlist:
                scanDict=alignRCPs(rcp,rcpFilename)
                addkidNr,addoffX,addoffY = getOffsetsFromRCPDict(scanDict)
                secondRCP.append([addoffX,addoffY,beamsize])
        else:
            for rcp in addRCPlist:
                addkidNr, _, _, addoffX, addoffY = openRCP(rcp)
                secondRCP.append([addoffX,addoffY,beamsize])
    else:
        secondRCP=[]
    
    
    showArrayPlot(offX, offY, gain, beamsize,
    				fileName=None,
    				wedge_label=r'Gain',
    				value_label=r'Gain',
    				pixel_ids=kidNr,
    				pixel_label=r'KID',pixel_unit='',febe=fe,limitsZ=limitsZ,highlightKids=highlightKids,addRCP=secondRCP)

def openCalFile(calfile=''):
    '''
    Open and read an Cal file.
    
    @param calfile: Cal filename (in working directory)
    @return: List of  kidNr,responce, Tmirror, Tlinear, Phlinear, Sensitivity 
    
    '''
    kidNr, responce, tmirror, tlinear, tsky, phlinear, sensitivity = [], [], [], [], [], [], []
    if calfile != '':
        try:
            f = file(calfile)
        except IOError:
            messages.error("Could not open file %s" % calfile)
            return
        
        lines = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for line in lines:
            if not line[0].startswith(('!', '#')):
                kidNr.append(np.int(line[0]))
                responce.append(np.float(line[1]))
                tmirror.append(np.float(line[2]))
                tsky.append(np.float(line[3]))
                tlinear.append(np.float(line[4]))
                phlinear.append(np.float(line[5]))
                sensitivity.append(np.float(line[6]))
     
        kidNr = np.array(kidNr)
        responce = np.array(responce)
        tmirror = np.array(tmirror)
        tsky = np.array(tsky)
        tlinear = np.array(tlinear)
        phlinear = np.array(phlinear)
        sensitivity = np.array(sensitivity)
        f.close()
    else:
        messages.error("No Cal file specified")
       
    return kidNr, responce, tmirror, tsky, tlinear, phlinear, sensitivity


def openTempCalFile(calfile=''):
    '''
    Open and read an Cal file.
    
    @param rcpfile: Cal filename (in working directory)
    @return: List of  kidNr,responce, Tmirror, Tlinear, Phlinear, Sensitivity 
    
    '''
    kidNr, tmirror, para= [], [], []
    if calfile != '':
        try:
            f = file(calfile)
        except IOError:
            messages.error("Could not open file %s" % calfile)
            return
        
        lines = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for line in lines:
            if not line[0].startswith(('!', '#')):
                kidNr.append(np.int(line[0]))
                tmirror.append(np.float(line[1]))
                para.append((line[2]))
                
     
        kidNr = np.array(kidNr)
        tmirror = np.array(tmirror)
        f.close()
    else:
        messages.error("No Cal file specified")
       
    return kidNr, tmirror, para



def openCalBeamMapFile(calfile=''):
    '''
    Open and read an Cal file.
    
    @param calfile: Cal filename (in working directory)
    @return: List of  kidNr,responce, Tmirror, Tlinear, Phlinear, Sensitivity 
    
    '''
    kidNr, responce, sensibest,sensifull,noisebest,noisefull = [], [], [], [], [], []
    if calfile != '':
        try:
            f = file(calfile)
        except IOError:
            messages.error("Could not open file %s" % calfile)
            return
        
        lines = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for line in lines:
            if not line[0].startswith(('!', '#')):
                kidNr.append(np.int(line[0]))
                responce.append(np.float(line[1]))
                sensibest.append(np.float(line[2]))
                sensifull.append(np.float(line[3]))
                noisebest.append(np.float(line[4]))
                noisefull.append(np.float(line[5]))
                
     
        kidNr = np.array(kidNr)
        responce = np.array(responce)
        sensibest = np.array(sensibest)
        sensifull = np.array(sensifull)
        noisebest = np.array(noisebest)
        noisefull = np.array(noisefull)
        f.close()
    else:
        messages.error("No Cal file specified")
       
    return kidNr, responce, sensibest,sensifull,noisebest,noisefull



def showSensitivity(calscan=0,febe='LFA',beamsize=12.):
    if calscan != 0:
        search = "*%i*_cal.dat"%(calscan)
        globlist=glob(search)
        print globlist,len(globlist)
        if len(globlist) > 1:
            "print filename ambigious: %s"%(globlist)
        elif len(globlist) ==  0:
            "No calfile with scannr %i found"%(calscan)
        else:
            calfile=globlist[0]
            print calfile
            kidNr, responce, tmirror, tlinear, phlinear, sensitivity = openCalFile(calfile)
           
            tmp=calfile.split('_')
            febe=tmp[2]
            rcp=getLastRCP(febe=febe)
            kidNrRCP, _, _, offX, offY = openRCP(rcp)
            kidNrRCP=list(kidNrRCP)

            offXcal=[]
            offYcal=[]
            sen=[]

            for i,k in enumerate(kidNr):
                try:
                    index=kidNrRCP.index(k)
                except:
                    pass
                else:
                    if sensitivity[i] < 30.0:
                        offXcal.append(offX[index])
                        offYcal.append(offY[index])
                        sen.append(sensitivity[i])
                
    
            showArrayPlot(offXcal, offYcal, sen, beamsize,
    				fileName=None,
    				wedge_label=r'Sensitivity [mK sqrt(s)]',
    				value_label=r'Sensitivity',
    				pixel_ids=kidNr,
    				pixel_label=r'KID')


    else:
        print "No calfile for scan %i found."%(calscan)
        
def showCal(calscan=0,febe='LFA',para='S',limitsZ=[],beamsize=0.):
    if beamsize == 0:
        if str.upper(febe) == 'LFA':
            beamsize=12.
            backend = 'AMKID870-AMKID870BE'
        else:
            beamsize=8.0
            backend = 'AMKID350-AMKID350BE'
    if calscan != 0:
        #backend=str.upper(febe)
        search = "*_%i_*%s_KidCal*.dat"%(calscan,backend)
        globlist=glob(search)
       
        if len(globlist) > 1:
            print "filename ambigious: %s"%(globlist)
        elif len(globlist) ==  0:
            print "No calfile with scannr %i found"%(calscan)
        else:     
            calfile=globlist[0]
            kidNr, responce, tmirror, tsky, tlinear, phlinear, sensitivity = openCalFile(calfile)
            
            tmp=calfile.split('_')
            myfebe=tmp[2]
            rcp=getLastRCP(febe=myfebe)
            kidNrRCP, _, _, offX, offY = openRCP(rcp)
            kidNrRCP=list(kidNrRCP)
            
            usedkids=[]
            offXcal=[]
            offYcal=[]
            param=[]

            for i,k in enumerate(kidNr):
                try:
                    index=kidNrRCP.index(k)
                except:
                    pass
                else:
                    if len(limitsZ)==2:
                        findmin=limitsZ[0]
                        findmax=limitsZ[1]
                    else:
                        findmin=-1e9
                        findmax= 1e9
                    if str.upper(para) == 'S':
                        value=round(sensitivity[i],1)
                        messagestr='sensitivity'
                        wedge_label = r'Sensitivity [mK sqrt(s)]'
                        value_label = r'Sensitivity'
                        pixel_unit = 'mK sqrt(s)'
                    if str.upper(para) == 'R':
                        value=round(1./(np.pi/180.*responce[i]),1)
                        messagestr='resonce'
                        wedge_label = r'Responce [deg/K]'
                        value_label = r'Responce'
                        pixel_unit = 'deg/K'
                    if str.upper(para) == 'TM':
                        value=np.round(tmirror[i])
                        messagestr='mirror temperature'
                        wedge_label = r'Mirror Temperature [K]'
                        value_label = r'Tmirror'
                        pixel_unit = 'K'
                    if str.upper(para) == 'TL':
                        value=np.round(tlinear[i],1)
                        messagestr='linear temperature'
                        wedge_label = r'Linear Temperature Range [K]'
                        value_label = r'linear Temp'
                        pixel_unit = 'K'
                    if str.upper(para) == 'TS':
                        value=np.round(tsky[i],1)
                        messagestr='sky temperature'
                        wedge_label = r'Sky Temperature [K]'
                        value_label = r'Tsky'
                        pixel_unit = 'K'
                    if str.upper(para) == 'PL':
                        value=np.round(phlinear[i],1)
                        messagestr='linear readout phase'
                        wedge_label = r'Linear Readout Phase Range [deg]'
                        value_label = r'Linear Phase'
                        pixel_unit = 'deg'
                    if  findmin < value < findmax:
                        usedkids.append(k)
                        offXcal.append(offX[index])
                        offYcal.append(offY[index])
                        param.append(value)
            if len(param) == 0:
                print "No values for %s found in %s within the specified limits"%(messagestr,calfile)
            else:
                showArrayPlot(offXcal, offYcal, param, beamsize,
    				fileName=None,
    				wedge_label=wedge_label,
    				value_label=value_label,
    				pixel_ids=usedkids,
    				pixel_label=r'KID',pixel_unit=pixel_unit,febe=febe)


    else:
        print "No calfile for scan %i found."%(calscan)
 
def showCalBeamMap(calscan=0,febe='LFA',para='S',limitsZ=[],beamsize=0.,Tcal=1):
    if beamsize == 0:
        if str.upper(febe) == 'LFA':
            beamsize=12.
            backend = 'AMKID870-AMKID870BE'
        else:
            beamsize=8.0
            backend = 'AMKID350-AMKID350BE'
    if calscan != 0:
        if Tcal != 1:
            search = "%i*_sensitivity.dat"%(calscan)
        else:
            search = "%i*_sensitivity-Tcal.dat"%(calscan)
        globlist=glob(search)
       
        if len(globlist) > 1:
            print "filename ambigious: %s"%(globlist)
        elif len(globlist) ==  0:
            print "No calfile with scannr %i found"%(calscan)
        else:     
            calfile=globlist[0]
            kidNr, responce, sensibest,sensifull,noisebest,noisefull = openCalBeamMapFile(calfile)
        
            Febe = 'LFA'
            rcp=getLastRCP(febe=febe)
            kidNrRCP, _, _, offX, offY = openRCP(rcp)
            kidNrRCP=list(kidNrRCP)
            
            usedkids=[]
            offXcal=[]
            offYcal=[]
            param=[]

            for i,k in enumerate(kidNr):
                try:
                    index=kidNrRCP.index(k)
                except:
                    pass
                else:
                    if len(limitsZ)==2:
                        findmin=limitsZ[0]
                        findmax=limitsZ[1]
                    else:
                        findmin=-1e9
                        findmax= 1e9
                    if str.upper(para) == 'S':
                        value=round(sensibest[i],1)
                        messagestr='sensitivity'
                        wedge_label = r'Sensitivity [mJy sqrt(s)]'
                        value_label = r'Sensitivity'
                        pixel_unit = 'mJy sqrt(s)'
                    if str.upper(para) == 'R':
                        if Tcal != 1:
                            value=round(1000*180./np.pi*1./responce[i],5)
                            messagestr='responce'
                            wedge_label = r'Responce [mdeg/Jy]'
                            value_label = r'Responce'
                            pixel_unit = 'mdeg/Jy'
                        else:
                            value=round(responce[i],1)
                            messagestr='coupling'
                            wedge_label = r'Coupling [Jy/K]'
                            value_label = r'Coupling'
                            pixel_unit = 'Jy/K'
                    if str.upper(para) == 'C':
                        value=round(responce[i],1)
                        messagestr='coupling'
                        wedge_label = r'Coupling [Jy/K]'
                        value_label = r'Coupling'
                        pixel_unit = 'Jy/K'

                    if str.upper(para) == 'SF':
                        value=np.round(sensifull[i])
                        messagestr='sensitivity'
                        wedge_label = r'Sensitivity [mJy sqrt(s)]'
                        value_label = r'Sensitivity'
                        pixel_unit = 'mJy sqrt(s)'
                    if str.upper(para) == 'N':
                        value=np.round(noisebest[i],1)
                        messagestr='noise'
                        wedge_label = r'Noise [Jy]'
                        value_label = r'noise'
                        pixel_unit = 'Jy'
                    if str.upper(para) == 'NF':
                        value=np.round(noisefull[i],1)
                        messagestr='noise'
                        wedge_label = r'Noise [mJy]'
                        value_label = r'noise'
                        pixel_unit = 'mJy'
                    
                    #if  findmin < value < findmax:
                    usedkids.append(k)
                    offXcal.append(offX[index])
                    offYcal.append(offY[index])
                    param.append(value)
            if len(param) == 0:
                print "No values for %s found in %s within the specified limits"%(messagestr,calfile)
            else:
                print "mean responce: ",np.median(responce)
                showArrayPlot(offXcal, offYcal, param, beamsize,
    				fileName=None,
    				wedge_label=wedge_label,
    				value_label=value_label,
    				pixel_ids=usedkids,
    				pixel_label=r'KID',pixel_unit=pixel_unit,febe=febe,limitsZ=limitsZ)


    else:
        print "No calfile for scan %i found."%(calscan)
 

def showAP(infile='',febe='LFA',para='minS',limitsZ=[],beamsize=0.):
    if beamsize == 0:
        if str.upper(febe) == 'LFA':
            beamsize=12.
            backend = 'AMKID870-AMKID870BE'
        else:
            beamsize=8.0
            backend = 'AMKID350-AMKID350BE'
    if infile != '':
        kidNr,jyperk, sensibest,sensimean,stdsensi,djyperk = openCalBeamMapFile(infile)


        s20,s80,speak,ppeak,sdist,pdist=getSensitivityBreakdown(sensibest)
        
        ch1,ch2,ch3,ch4=getKidsperChain(kidNr)
        chains=[ch1,ch2,ch3,ch4]
        print "##########################################"
        print "        Mean best 20%     Mean best 80% "  
        print "    ALL:     %6.1f          %6.1f      mJy sqrt(s) (%i KIDS) "%(s20,s80,len(kidNr))
        for i in range(len(chains)):
            sen = []
            for j in range(len(chains[i])):
                ind=list(kidNr).index(chains[i][j])
                sen.append(sensibest[ind])
        
            s20,s80,speak,ppeak,sdist,pdist=getSensitivityBreakdown(sen)
            print "Chain %i:     %6.1f          %6.1f      mJy sqrt(s) (%i KIDS)"%(i+1,s20,s80,len(sen))
        print "##########################################"
        if len(kidNr) > 0:
            #print para
            Febe = 'LFA'
            rcp=getLastRCP(febe=febe)
            kidNrRCP, _, _, offX, offY = openRCP(rcp)
            kidNrRCP=list(kidNrRCP)
            
            usedkids=[]
            offXcal=[]
            offYcal=[]
            param=[]

            for i,k in enumerate(kidNr):
                try:
                    index=kidNrRCP.index(k)
                except:
                    pass
                else:
                    if len(limitsZ)==2:
                        findmin=limitsZ[0]
                        findmax=limitsZ[1]
                    else:
                        findmin=-1e9
                        findmax= 1e9
                    if str.upper(para) == 'MINS':
                        value=round(sensibest[i],1)
                        messagestr='best sensitivity'
                        wedge_label = r'Sensitivity [mJy sqrt(s)]'
                        value_label = r'best Sensitivity'
                        pixel_unit = 'mJy sqrt(s)'
                
                    
                    if str.upper(para) == 'C':
                        value=round(jyperk[i],1)
                        messagestr='mean coupling'
                        wedge_label = r'Coupling [Jy/K]'
                        value_label = r'Coupling'
                        pixel_unit = 'Jy/K'

                    if str.upper(para) == 'STDC' or str.upper(para) == 'RMSC':
                        value=round(djyperk[i],1)
                        messagestr='std coupling'
                        wedge_label = r'std Coupling [Jy/K]'
                        value_label = r'std Coupling'
                        pixel_unit = 'Jy/K'
                    if str.upper(para) == 'STDC' or str.upper(para) == 'RMSC':
                        value=round(djyperk[i],1)
                        messagestr='best coupling'
                        wedge_label = r'best Coupling [Jy/K]'
                        value_label = r'best Coupling'
                        pixel_unit = 'Jy/K'

                    if str.upper(para) == 'SM':
                        value=np.round(sensimean[i])
                        messagestr='mean sensitivity'
                        wedge_label = r'Sensitivity [mJy sqrt(s)]'
                        value_label = r'mean Sensitivity'
                        pixel_unit = 'mJy sqrt(s)'

                    if str.upper(para) == 'SSTD' or str.upper(para) == 'SRMS':
                        value=np.round(stdsensi[i])
                        messagestr='std sensitivity'
                        wedge_label = r'std Sensitivity [mJy sqrt(s)]'
                        value_label = r'std Sensitivity'
                        pixel_unit = 'mJy sqrt(s)'
                    

                                        
                    #if  findmin < value < findmax:
                    usedkids.append(k)
                    offXcal.append(offX[index])
                    offYcal.append(offY[index])
                    param.append(value)
            if len(param) == 0:
                print "No values for %s found in %s within the specified limits"%(messagestr,calfile)
            else:
                #print "mean responce: ",np.median(responce)
                showArrayPlot(offXcal, offYcal, param, beamsize,
    				fileName=None,
    				wedge_label=wedge_label,
    				value_label=value_label,
    				pixel_ids=usedkids,
    				pixel_label=r'KID',pixel_unit=pixel_unit,febe=febe,limitsZ=limitsZ)


    else:
        print "No file %s found."%(infile)
 

       
def showLN2Cal(calscan=0,febe='LFA',para='TM',angle=0.0,limitsZ=[],beamsize=0.,calfile=''):
    
    if beamsize == 0:
        if str.upper(febe) == 'LFA':
            backend = "AMKID870-AMKID870BE"
            beamsize=12.
        else:
            backend = "AMKID350-AMKID350BE"
            beamsize=5.0
    if calscan != 0 or calfile:
        if calfile:
            search = "%s"%(calfile)
        else:
            search = "*%i*%s*_TempCalspline.dat"%(calscan,backend)
        globlist=glob(search)
    
    
        if len(globlist) > 1:
            print "filename ambigious: %s"%(globlist)
        elif len(globlist) ==  0:
            print "No calfile with scannr %i found"%(calscan)
        else:     
            calfile=globlist[0]
            kidNr, tmirror, TempCalPara= openTempCalFile(calfile)
            tmp=calfile.split('_')
            myfebe=tmp[2]
            rcp=getLastRCP(febe=myfebe)
            kidNrRCP, _, _, offX, offY = openRCP(rcp)
            kidNrRCP=list(kidNrRCP)
            
            usedkids=[]
            offXcal=[]
            offYcal=[]
            param=[]

            for i,k in enumerate(kidNr):
                try:
                    index=kidNrRCP.index(k)
                except:
                    pass
                else:
                    if len(limitsZ)==2:
                        findmin=limitsZ[0]
                        findmax=limitsZ[1]
                    else:
                        findmin=-1e9
                        findmax= 1e9
                    if str.upper(para) == 'TM':
                        value=np.round(tmirror[i],1)
                        messagestr='mirror temperature'
                        wedge_label = r'Mirror Temperature [K]'
                        value_label = r'Tmirror'
                        pixel_unit = 'K'
                    if str.upper(para) == 'TA':
                        entries = str(TempCalPara[i]).replace('[','').replace(']','').split(',')
                        Thot = np.float(entries[3])
                        Tref = np.float(entries[4])
                        bestLN2angle = np.float(entries[5])
                        Tcold=((Tref-Thot*np.cos((bestLN2angle-np.float(entries[2]))*np.pi/180.)**2)/(1-np.cos((bestLN2angle-np.float(entries[2]))*np.pi/180.)**2)+np.float(entries[1])*(np.sin((angle-np.float(entries[0])-np.float(entries[2]))*np.pi/180.))-np.float(entries[1])*(np.sin((bestLN2angle-np.float(entries[0])-np.float(entries[2]))*np.pi/180.)))
                        temp = Thot*np.cos((angle-np.float(entries[2]))*np.pi/180.)**2+Tcold*(1.-np.cos((angle-np.float(entries[2]))*np.pi/180.)**2)
                        
                        value=np.round(temp,1)
                        messagestr='Pol Grid temperature'
                        wedge_label = r'Pol Grid Temperature [K]'
                        value_label = r'Tpol'
                        pixel_unit = 'K'
                    #if  findmin < value < findmax:
                    usedkids.append(k)
                    offXcal.append(offX[index])
                    offYcal.append(offY[index])
                    param.append(value)
            if len(param) == 0:
                print "No values for %s found in %s within the specified limits"%(messagestr,calfile)
            else:
                showArrayPlot(offXcal, offYcal, param, beamsize,
    				fileName=None,
    				wedge_label=wedge_label,
    				value_label=value_label,
    				pixel_ids=usedkids,
    				pixel_label=r'KID',pixel_unit=pixel_unit,febe=febe,limitsZ=limitsZ)


    else:
        print "No calfile for scan %i found."%(calscan)
        




def findPixelsInRadius(chanNr, radius=90):
    chanIndex = data.BolometerArray.getChanIndex(chanNr)
    az = data.BolometerArray.Offsets[0][chanIndex]
    el = data.BolometerArray.Offsets[1][chanIndex]
    usedChannels = data.BolometerArray.checkChanList([])
    allIndex = data.BolometerArray.getChanIndex(usedChannels)

    dAz = np.take(np.abs(az - data.BolometerArray.Offsets[0]),allIndex)
    dEl = np.take(np.abs(el - data.BolometerArray.Offsets[1]),allIndex)
    
    dist = np.sqrt(dAz**2+dEl**2) - radius
    mask = dist<radius
    
    nearbyChans = np.extract(mask, usedChannels).astype(int)

    return nearbyChans


def nrChans():
    usedbolos = data.BolometerArray.checkChanList([])
    return len(usedbolos)


def checkNaN():
    """
    Check if number of NaNs in timestream is >2% and pause script execution if so.
    @return: void
    
    """
    usedchans = data.BolometerArray.checkChanList([])
    dataY = np.asarray(data.getChanListData('flux', usedchans))
    
    nrnan = np.sum(np.isnan(dataY))
    
    nrdumps = dataY.size
    
    if float(nrnan)/nrdumps * 100. > 2.0:
        warn('%i NaNs (%8.4f percent) found on time stream!' % (nrnan, float(nrnan)/nrdumps*100.))
        raw_input('Press <enter> to continue:')


def checkMJD():
    mytime=1000*3600*24*(data.ScanParam.MJD-data.ScanParam.MJD[0])
    nrdumps=len(mytime)
    dt=[]
    for i in range(1,nrdumps):
        dt.append(mytime[i]-mytime[i-1])
    dt=np.array(dt)
    meansampling=np.median(dt)
    if np.max(dt) > 1.1*meansampling:
        nrofevents=len(np.extract(np.where(abs(np.array(dt))>1.1*meansampling),dt))               
        print "%i corrupted MJD entries found- check load on AMKID backend system"%(nrofevents)
        print "mean sampling rate is %5.1f ms (%5.1f Hz)"%(meansampling,1000./meansampling) 
        plot(dt,labelX="MJD index",labelY="dt [ms]")
        raw_input()
    else:
        if 1000./meansampling < 20:
            print "WARNING: mean sampling rate is %5.1f ms (%5.1f Hz)"%(meansampling,1000./meansampling)
            raw_input()
        else:
           print "mean sampling rate is %5.1f ms (%5.1f Hz)"%(meansampling,1000./meansampling) 

def writeFits3(data,outfile='boaMap.fits',overwrite=0,limitsX=[],limitsY=[],intensityUnit="Jy/beam"):
        """
        DES: store the current map (2D array with WCS info) to a FITS file
        INP: (str)   outfile: output file name (default boaMap.fits)
             (bool) overwrite: overwrite existing file -
                              default = 0: do not overwrite existing file
             (f list) limitsX/Y: optional map limits (in world coordinates)
             (string) intensityUnit: optional unit of the intensity (default: "Jy/beam")
        """
        #from boa           import BoaDataAnalyser, BoaMessageHandler, BoaFits, BoaFlagHandler
        #from boa           import BoaConfig, BoaError
        #from boa.Bogli     import Plot, MultiPlot
        #from boa.fortran   import fUtilities, fMap, fStat
        #from boa.Utilities import Timing, compressNan, tolist_boa, inPolygon, outPolygon
        #from boa.Utilities import modelBaseEllipticalGaussian,fitBaseEllipticalGaussian
        #import Numeric

        if os.path.exists(outfile):
            if not overwrite:
                print('File %s exists' % outfile)
                return
        
        try:
            dataset = BoaFits.createDataset("!" + outfile)
        except:
            print('Could not open file %s' % (outfile))
            return

        
        localMap = copy.deepcopy(data)
        localMap.Data=Numeric.array(localMap.Data)
        
        try:
            #write FLux plane
            localMap._Image__writeImage(dataset, "Intensity", intensityUnit=intensityUnit)

            #compute RMS and SNR plane
            rmsMap = copy.deepcopy(localMap)
            snrMap = copy.deepcopy(localMap)
            snrMap.Data *= sqrt(snrMap.Weight)
            a = snrMap.computeRms(rmsKappa=3.0)
            scale=snrMap.RmsBeam
            snrMap.Data /= array(scale,'f')

            rmsMap.Data =  (rmsMap.Data*0.0+1.0)/sqrt(rmsMap.Weight)
            rmsMap.Data *= array(scale,'f')
            
            
            rmsMap.Data=Numeric.array(rmsMap.Data)
            snrMap.Data=Numeric.array(snrMap.Data)
            #write RMS plane
            rmsMap._Image__writeImage(dataset, "Intensity", intensityUnit=intensityUnit)
            #write SNR plane
            snrMap._Image__writeImage(dataset, "Intensity", intensityUnit='SNR')
            dataset.close()
            
        except:
            print('Could not write data to file %s' % (outfile))
            return

        
        localMap = 0  # free memory
        snrMap = 0
        rmsMap = 0
        dataset = 0

def flagBeamShape(infile=None,scan=None,fe='LFA',maxE=1.3,minC=0.6):
    if not infile:
        infile=getBeamMapFile(scan=scan,fe=fe)
        print(infile)
    beamDict=readBeamMapDict(infile)
    if str.upper(fe) == 'LFA':
        couplingRef = 40.0
    if str.upper(fe) == 'HFA':
        couplingRef = 75.0

    badChans=[]
    for kid in beamDict.keys():
        fwhmmax = np.max([beamDict[kid]['xfwhm'],beamDict[kid]['yfwhm']])
        fwhmmin = np.min([beamDict[kid]['xfwhm'],beamDict[kid]['yfwhm']])
        if beamDict[kid]['JyperK'] > couplingRef/minC or fwhmmax/fwhmmin > maxE:
            badChans.append(kid)
    #flagC(badChans)
    return badChans

def checkLOFreqs(scan,fe='LFA'):
    getHK(scan=scan,fe=fe)
    _,chains,n_kid=getFebe(fe)
    kids=np.array(data.tone_dict.keys())
    for chain in chains:
        firstkid=(chain-1)*n_kid+1
        lastkid=firstkid+n_kid-1
        ind=np.where( (firstkid <= kids) * (kids <= lastkid) )[0]
        if len(ind) > 0:
            kid=kids[ind[0]]
            print("Chain %i: LO1: %9.4f LO2: %9.4f LO3: %9.4f MHz"%(chain,data.tone_dict[kid]['f_lo1'],data.tone_dict[kid]['f_lo2'],data.tone_dict[kid]['f_lo3']))
        else:
            print("Chain %i: no entries for this chain in HK Header"%(chain)) 
        

def optimizeP(wirescans=[],fe='LFA',tonelist=None):
    senDict={}
    for scan in wirescans:
        wsdict=readScannerResult(scan=scan,fe=fe)
        for kid in wsdict.keys():
            if kid not in senDict.keys():
                auxdict={}
                auxdict['NETs'] = [wsdict[kid]['NET']]
                auxdict['powers'] = [wsdict[kid]['power']]
                auxdict['freqs'] = [wsdict[kid]['freq']]
            else:
                auxdict=senDict[kid]
                auxdict['NETs'].append(wsdict[kid]['NET'])
                auxdict['powers'].append(wsdict[kid]['power'])
                auxdict['freqs'].append(wsdict[kid]['freq'])
            
            senDict[kid]=auxdict


    tonedict=readToneList(tonelist)

    for kid in tonedict.keys():
        if kid in senDict.keys():
            nets=np.array(senDict[kid]['NETs'])
            powers=np.array(senDict[kid]['powers'])
            freqs=np.array(senDict[kid]['freqs'])
            ind=np.where(nets == np.min(nets))[0][0]
            bestpower=powers[ind]
            bestfreq=freqs[ind]
            tonedict[kid]['power'] = bestpower
            tonedict[kid]['frequency'] = bestfreq
    
    filename = "ToneFiles/"+tonelist.split(".MSK")[0]+"_optimized.MSK"
    lines_kid = ''
    lines_ref = ''
    
    nrkids=0
    nrblinds=0
    for kid in tonedict.keys():
        if tonedict[kid]['type'] == 'D':
            lines_kid += '%i, %9.4f, %6.2f, %s, %i\n' % (kid, tonedict[kid]['frequency'], tonedict[kid]['power'], 'D', tonedict[kid]['hwchain'])
            nrkids=nrkids+1
        if tonedict[kid]['type'] == 'B':
            lines_ref += '%i, %9.4f, %6.2f, %s, %i\n' % (kid, tonedict[kid]['frequency'], tonedict[kid]['power'], 'B', tonedict[kid]['hwchain'])
            nrblinds=nrblinds+1

    header =  "#description=%s\n" % filename
        
    header += "#%s tone list based on optimised power and frequencies for tonefile %s\n" % (fe, tonelist)
    header += "#input wire scans: %s\n" %(str(wirescans))
    header += "#Number of KIDs: %i\n" % lines_kid.count('\n')
    header += "#KidNr, tone frequency [MHz], tone power [dBm], tone type [D = Detector, B = Blind tone], Hardware chain number\n"

    with open(filename, 'w') as f:
        f.write(header)
        f.write(lines_kid)
        f.write(lines_ref)
        
    info("New tonelist written to '%s'." % filename)
    info("Tonelist contains %i KIDs and %i Bind tones" %(nrkids,nrblinds))


    #return senDict

def getScansFromCalFiles(scantype='scanner',above=None, below=None,fe='LFA',cal_dir=None,verbose=True):
   if str.lower(scantype) not in ['scanner','beamscan','fsweep','sweep','beammap']:
      print("Error - unsupported scantype: allowed options are 'scanner','beamscan','fsweep','sweep','beammap'")
      return
   if cal_dir is None:
       cal_dir = os.path.join(os.getcwd(),'CalFiles')
   
   if not above:
       above = 0
       
   if not below:
       below = 99999
   if str.lower(scantype) in ['scanner','beamscan']:
      mytype = "scanner"
   elif str.lower(scantype) in ['fsweep','sweep']:
      mytype = "fsweep"
   else:
      mytype = "beam_map"
   
   myend=str.upper(fe)+".csv"
   search = cal_dir+"/*"+str(mytype)+"*"+str(myend)
   print(search)
   cals=glob(search)
   scans = []
   for cal in cals:
      scannr = np.int(cal.split("/")[-1].split(myend)[0].split("_")[-2])
      if above <= scannr <= below:
         scans.append(scannr)
   print("found %i scans of type %s"%(len(scans),mytype))
   return scans

def getScansFromObslog(searchstr='',keyword='comment', output_params='number',above=None, below=None, obslog_dir=None,verbose=True):
    """
    Get the scans matching a string pattern from the obslog
    
    @param searchstr: string pattern
    @param keyword: string ; indicating in which keyword to find the searchstr (e.g. in 'comment')
    @param output_params: string or list of strings indicating which parameters to return
    @param below: format YYYY-MM-DD: date up to which the search is carried out  
    @param above: format YYYY-MM-DD: date from which onwards the search is carried out  
    @param scan_dir: Scan directory containing the obslogs (optional).
    @param verbose: logical ; if true function prints the matchs of the searchstr.
    @return: list with a list of output paramters
    @rtype: str
    
    """
    validkeys=['CA','IE','azimuth','command','comment','date','declination','dewPoint','elevation','febes','focTempCorr','geom','humidity','lines','lst','mode','numSubScans','number','observerID','offsets','operatorID','pressure','projectID','pwv','reference','rightAscension','source','status','strokeMode','switchMode','temperature','totalScanTime','type','utc','velocity','windDirection','windSpeed','xFocus','yFocus','zFocus']

    if keyword not in validkeys:
       print("Error: keyword %s for search string is not a valid input"%(keyword))
       print("Valid keys are", validkeys)
       return 

    if isinstance(output_params,str):
       output_params=[output_params]
    
    for param in output_params:
        if param not in validkeys:
           print("Error: output parameter %s is not a valid input"%(param))
           print("Valid keys are", validkeys)
           return 

    if obslog_dir is None:
       obslog_dir = BoaConfig.inDir.replace('rawdata','obslogs')
       if os.path.exists(obslog_dir) == False:
           obslog_dir=os.path.join(os.getcwd(),'obslogs')
           if os.path.exists(obslog_dir) == False:
               print("error: can not find obslog directory")
               print("set directory manually by using the obslog_dir=XX parameter")
               return

    #print(np.int(above.split("-")[0]),np.int(above.split("-")[1]),np.int(above.split("-")[2]))
    if above:
       #print(above)
       #print(np.int(above.split("-")[0]),np.int(above.split("-")[1]),np.int(above.split("-")[2]))
       date_start = datetime.datetime(np.int(above.split("-")[0]),np.int(above.split("-")[1]),np.int(above.split("-")[2]),0,0,0)
    else:
       date_start = datetime.datetime(2000,01,01,0,0,0)
    
    if below:
       date_stop = datetime.datetime(np.int(below.split("-")[0]),np.int(below.split("-")[1]),np.int(below.split("-")[2]),0,0,0)
    else:
       date_stop = datetime.datetime(2100,01,01,0,0,0)   
    
    ## find all obslogs
    if obslog_dir[-1] != "/":
        search = obslog_dir+"/*.obslog"
    else:
        search = obslog_dir+"*.obslog"
    
    obslogs=glob(search)
    mylogs=[]
    for obslog in obslogs:
       #print(obslog)
       year= obslog.split("/")[-1].split("-")[5]
       month=obslog.split("/")[-1].split("-")[6]
       day=  obslog.split("/")[-1].split("-")[7].split(".")[0]
       this_date = datetime.datetime(np.int(year),np.int(month),np.int(day),0,0,0) 
       timediff1 = this_date - date_start
       timediff2 = date_stop - this_date
       if timediff1.days >= 0 and timediff2.days >= 0:
          mylogs.append(obslog)
    
    
    result=[]
    for filename in mylogs:
       obslog=xmltree.parse(filename)
       for scan in obslog.findall('scan'):
          item = scan.find(keyword).text
          if str.upper(searchstr) in str.upper(item):
             values =[]
             for key in output_params:
                values.append(scan.find(key).text)
             if verbose == True:
                print "scan: ",scan.find('number').text,"match:", item, scan.find(key).text
             result.append(values)
    return result
    
def getPlanetScans(start=None, stop=None,sources=None,mytype='cal',checkFoc=True,focTime=[180.0,20.0]):
    if not sources:
        sources=['Mars','Uranus','Neptune']
    else:
        if isinstance(sources,str):
            sources=list(sources)
    myscans = []
    for source in sources:
        result=getScansFromObslog(searchstr=source,keyword='source',output_params=['number','totalScanTime','command','type','date','utc'],above=start, below=stop,verbose=False)
        
        for res in result:
            if str.upper(mytype) == "CAL":
                if np.float(res[1]) > 10. and (str.lower(res[2]).find("point") >= 0 or str.lower(res[2]).find("raster") >= 0) and  str.lower(res[3]).find("focus") == -1:
                    if checkFoc == False:
                        myscans.append(np.int(res[0]))

                    else:
                        refdate =  datetime.datetime(np.int(res[4].split("-")[0]),np.int(res[4].split("-")[1]),np.int(res[4].split("-")[2]),
                                                     np.int(res[5].split(":")[0]),np.int(res[5].split(":")[1]) ,np.int(res[5].split(":")[2])) 
                        focDoneBefore = False
                        focDoneAfter = False

                        ##check only obslogs a day before and a day after the scans
                        day = np.int(res[4].split("-")[2])-1
                        if day < 10:
                            myday = "0"+str(day)
                        else:
                            myday = str(day)
                        above = str(res[4].split("-")[0])+"-"+str(res[4].split("-")[1])+"-"+myday
                        
                        day = np.int(res[4].split("-")[2])+1
                        if day < 10:
                            myday = "0"+str(day)
                        else:
                            myday = str(day)
                        below = str(res[4].split("-")[0])+"-"+str(res[4].split("-")[1])+"-"+myday
                        ####

                        #print(above,below)
                        resfoc=getScansFromObslog(searchstr='focus',keyword='command',output_params=['date','utc','number'],above=above,below=below,verbose=False)
                        for resf in resfoc:
                            thisdate = datetime.datetime(np.int(resf[0].split("-")[0]),np.int(resf[0].split("-")[1]),np.int(resf[0].split("-")[2]),
                                                         np.int(resf[1].split(":")[0]),np.int(resf[1].split(":")[1]) ,np.int(resf[1].split(":")[2])) 
                            dt=(thisdate - refdate).total_seconds()

                            ## OK if a focus scan has been done before
                            if -1*focTime[0]*60.0 <= dt <= 0:
                                focDoneBefore = True

                            ## remove if a focus has been done just afterwards
                            if 0 < dt <= focTime[1]*60.0:
                                focDoneAfter = True
                            
                        
                        if focDoneBefore == True and focDoneAfter == False:
                            myscans.append(np.int(res[0]))
                        

            if str.upper(mytype) == "BEAMMAP":
                if np.float(res[1]) > 10. and str.lower(res[2]).find("otf") >= 0 and  str.lower(res[3]).find("focus") == -1:
                    if checkFoc == False:
                        myscans.append(np.int(res[0]))
                    else: #for beam maps we skip the focus check since this should have been done
                        myscans.append(np.int(res[0]))
    myscans.sort()
    return myscans

def getSeconadryScans(start=None, stop=None,sources=None,checkFoc=True,focTime=[180.0,20.0]):
    if not sources:
        sources = []
        for source in calibFluxes.keys():
            sources.append(source)
    else:
        if isinstance(sources,str):
            sources=list(sources)
    myscans = []
    for source in sources:
        result=getScansFromObslog(searchstr=source,keyword='source',output_params=['number','totalScanTime','command','type','date','utc'],above=start, below=stop,verbose=False)
        
        for res in result:
            if np.float(res[1]) > 10. and (str.lower(res[2]).find("point") >= 0 or str.lower(res[2]).find("raster") >= 0 or str.lower(res[2]).find("otf") >= 0) and  str.lower(res[3]).find("focus") == -1:
                if checkFoc == False:
                    myscans.append(np.int(res[0]))
                else:
                    refdate =  datetime.datetime(np.int(res[4].split("-")[0]),np.int(res[4].split("-")[1]),np.int(res[4].split("-")[2]),
                                                     np.int(res[5].split(":")[0]),np.int(res[5].split(":")[1]) ,np.int(res[5].split(":")[2])) 
                    focDoneBefore = False
                    focDoneAfter = False

                    ##check only obslogs a day before and a day after the scans
                    day = np.int(res[4].split("-")[2])-1
                    if day == 0:
                        day=1
                    if day < 10:
                        myday = "0"+str(day)
                    else:
                        myday = str(day)
                    above = str(res[4].split("-")[0])+"-"+str(res[4].split("-")[1])+"-"+myday
                    
                    day = np.int(res[4].split("-")[2])+1
                    if day > 31:
                        day = 31
                    if day < 10:
                        myday = "0"+str(day)
                    else:
                        myday = str(day)
                    below = str(res[4].split("-")[0])+"-"+str(res[4].split("-")[1])+"-"+myday
                    ####

                        
                    resfoc=getScansFromObslog(searchstr='focus',keyword='command',output_params=['date','utc','number'],above=above,below=below,verbose=False)
                    for resf in resfoc:
                        thisdate = datetime.datetime(np.int(resf[0].split("-")[0]),np.int(resf[0].split("-")[1]),np.int(resf[0].split("-")[2]),
                                                         np.int(resf[1].split(":")[0]),np.int(resf[1].split(":")[1]) ,np.int(resf[1].split(":")[2])) 
                        dt=(thisdate - refdate).total_seconds()

                        ## OK if a focus scan has been done before
                        if -1*focTime[0]*60.0 <= dt <= 0:
                            focDoneBefore = True

                        ## remove if a focus has been done just afterwards
                        if 0 < dt <= focTime[1]*60.0:
                            focDoneAfter = True
                            
                        
                    if focDoneBefore == True and focDoneAfter == False:
                        myscans.append(np.int(res[0]))
    myscans.sort()
    return myscans

def createSourceModel(m=None,mtype='SNR',sm=0.0,highcut=5.0,lowcut=1.5,clip=3.0):
    '''
    VWO: new implementation without the need to rescale the data. Now noise map is accurate.
    '''
    info('Creating source model...')
    if not m and len(data.Map.Data) == 0:
        print("No input map found")
        return
    if str.upper(mtype) not in ['SNR','FLUX','S','F']:
        print("unknown type. Expected options are  'SNR','FLUX','S','F'")
        return
    
    if not m:
        flux=copy.deepcopy(data.Map)
    else:
        flux=copy.deepcopy(m)
        
    if sm > 0.0:    
        if flux.WCS['CTYPE1'] == 'RA---GLS' or flux.WCS['CTYPE1'] == 'GLON-GLS':
            __smoothBy(flux, sm/3600.)
        else:
            __smoothBy(flux, sm)
    
    rms = copy.deepcopy(flux)
    rms.Data = np.where(rms.Weight > 0.0, 1.0 / np.sqrt(rms.Weight), np.NaN)
    snr = copy.deepcopy(flux)
    snr.Data = flux.Data / rms.Data
    
    if clip > 0.0:
        minrms = np.nanpercentile(rms.Data, 5)  # instead of min
        medrms = np.nanmedian(rms.Data)
        info('min r.m.s. = %.1f mJy/b | median r.m.s. = %.1f mJy/b'%(minrms*1000,medrms*1000))
        info('Clipping model at RMS > %.1f*minrms = %.1f mJy/b'%(clip, clip*minrms*1000))
        mask = np.where(rms.Data > clip*minrms)
        snr.Data[mask] = np.NaN
        flux.Data[mask] = np.NaN

    modelmask=np.zeros(snr.Data.shape)

    ahigh = np.zeros(snr.Data.shape)
    alow = np.zeros(snr.Data.shape)
    ahigh[snr.Data >= highcut] = 1
    alow[snr.Data >= lowcut] = 1


    labels_low = label(alow)[0]
        
    ind = np.unique(labels_low * ahigh)[1:]
    print("%i regions found with SNR > %5.1f"%(len(ind),highcut))
    for i in ind:
        modelmask[labels_low == i] = 1
        
    if str.upper(mtype) in ['SNR','S']:
        snr.Data = snr.Data*modelmask
        return snr
    else:
        flux.Data = flux.Data*modelmask
        return flux


        
def createSourceModel_old(m=None,mtype='SNR',sm=0.0,highcut=5.0,lowcut=1.5,clip=3.0):
    if not m and len(data.Map.Data) == 0:
        print("No input map found")
        return
    if str.upper(mtype) not in ['SNR','FLUX','S','F']:
        print("unknown type. Expected options are  'SNR','FLUX','S','F'")
        return
    
    if not m:
        flux=copy.deepcopy(data.Map)
    else:
        flux=copy.deepcopy(m)
    
    rms=copy.deepcopy(flux)
    snr = copy.deepcopy(flux)

    snr.Data*=np.sqrt(snr.Weight)

    
    if sm > 0.0:
        smflux=copy.deepcopy(flux)
        if flux.WCS['CTYPE1'] == 'RA---GLS' or flux.WCS['CTYPE1'] == 'GLON-GLS':
            snr.smoothBy(sm/3600.)
            smflux.smoothBy(sm/3600.)
        else:
            snr.smoothBy(sm)
            smflux.smoothBy(sm)
            
    tmp=copy.deepcopy(snr)
    tmp.iterativeSigmaClip(below=-3,above=3)
    tmp.computeRms()
    scale=tmp.RmsBeam

    snr.Data /= array(scale,'f')
    mask=np.where(rms.Weight > 0.0)
    rms.Data[mask] = np.sqrt(1./rms.Weight[mask])/scale
    mask=np.where(rms.Weight <= 0.0)
    rms.Data[mask] = np.NaN
    
    if clip > 0.0:
        mask=np.where(rms.Data > clip*np.nanmin(rms.Data))
        snr.Data[mask] = np.NaN
        flux.Data[mask] = np.NaN

    modelmask=np.zeros(snr.Data.shape)

    ahigh = np.zeros(snr.Data.shape)
    alow = np.zeros(snr.Data.shape)
    ahigh[snr.Data >= highcut] = 1
    alow[snr.Data >= lowcut] = 1


    labels_low = label(alow)[0]
        
    ind = np.unique(labels_low * ahigh)[1:]
    print("%i regions found with SNR > %5.1f"%(len(ind),highcut))
    for i in ind:
        modelmask[labels_low == i] = 1
        
    if str.upper(mtype) in ['SNR','S']:
        snr.Data = snr.Data*modelmask
        return snr
    else:
        flux.Data = flux.Data*modelmask
        return flux 
    

def reject_outliers(data, m=2.):
    if len(data) == 1 or len(np.unique(np.array(data))) == 1:
        return data
    if isinstance(data,list):
        data = np.array(data)
        values=data[np.abs(data-np.median(data)) < m * np.std(data)]
        return list(values)
    else:
        return data[np.abs(data-np.median(data)) < m * np.std(data)]
    
    
