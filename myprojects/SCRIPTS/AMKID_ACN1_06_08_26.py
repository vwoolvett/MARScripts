#!/usr/bin/env python
# APEX - Atacama Pathfinder EXperiment Project
#
# who       when        what
# --------  ----------  ----------------------------------------------
#aweiss     Sep-2013     created
#vwo        Jul-2026    add ignore numpy warnings
'''
Name: AMKID.py
This macro contains reduction functions (shortcuts) for MKID data.
Files are expected in $BOA_HOME_AMKID.

Functions currenty available:
    - redfoc(ScanNr,fe,tau)    : Focus reduction
    - redmfoc(ScanNr,fe,tau)  : Focus-mapping reduction
    - redpnt(ScanNr,fe,tau)   : Pointing reduction
    - redweak(ScanNr,fe,tau)  : Pointing on weak source
    - redcal(ScanNr,fe,tau)   : Map reduction for calibration/pointing
    - redmap(ScanNr,fe,tau)   : Map reduction for science targets
    - redbeam(ScanNr,fe)      : Reduce a beam map
    - redscans(ScanNrList,fe): Reduce and coadd science maps
'''
import os, sys, re
import signal
from time import time
import time as ti
import IPython
from subprocess32 import CalledProcessError

from amkid.mars_help import findAMKIDFunctions
import matplotlib as mpl
import numpy as np
import scipy as sp
from mars.fortran import fMap

# to ignore python warnings
import warnings


class TimeoutError(Exception):
    def __init__(self, value = "Timed Out"):
        self.value = value
    def __str__(self):
        return repr(self.value)

def timeout(seconds_before_timeout):
    def decorate(f):
        def handler(signum, frame):
            raise TimeoutError()
        def new_f(*args, **kwargs):
            old = signal.signal(signal.SIGALRM, handler)
            old_time_left = signal.alarm(seconds_before_timeout)
            if 0 < old_time_left < seconds_before_timeout: # never lengthen existing timer
                signal.alarm(old_time_left)
            start_time = time()
            try:
                result = f(*args, **kwargs)
            finally:
                if old_time_left > 0: # deduct f's run time from the saved timer
                    old_time_left -= time() - start_time
                signal.signal(signal.SIGALRM, old)
                signal.alarm(old_time_left)
            return result
        new_f.func_name = f.func_name
        return new_f
    return decorate


BoaConfig.rcpPath = os.getenv('BOA_LOCAL_RCP')

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'),'boaAMKID.py'))
except Exception, e:
    print e
    print "File boaAMKID.py not found"
    print "AMKID reduction functions will not work!!!"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'),'cabling.py'))
except Exception, e:
    print e
    print "File boaAMKID.py not found"
    print "AMKID reduction functions will not work!!!"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'),'calib.py'))
except Exception, e:
    print e
    print "File calib.py not found"
    print "AMKID reduction functions will not work!!!"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'secondary_fluxes.py'))
except Exception, e:
    print e
    print "File secondary_fluxes.py not found"
    print "AMKID reduction of calibrators cannot run"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'planet-flux.py'))
except Exception, e:
    print e
    print "File planet-flux.py  not found"
    print "AMKID reduction of primary calibrators (Planets) cannot run"


try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'focus.py'))
except Exception, e:
    print e
    print "File focus.py not found."
    print "No reduction of map-focus possible."

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'reduce_rcp_amkid.py'))
except:
    print "reduce_rcp_amkid.py not found"
    print "No reduction of beam maps possilbe"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'fit_rcp.py'))
except:
    print "fit_rcp.py not found"
    print "No fitting of chip models to the RCP possible"
    
try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'boacorrections.py'))
except Exception, e:
    print e
    print "boacorrections.py not found"
    print "Correct map smoothing and FITS writing will not work"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'fkt_amkid.py'))
except Exception, e:
    print e
    print "fkt_amkid.py not found"
    print "some functions may not work"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'Wire_Scan.py'))
except Exception, e:
    print e
    print "Wire_Scan.py not found"
    print "some functions may not work"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'arrayStatus.py'))
except Exception, e:
    print e
    print "arrayStatus.py not found"
    print "some functions may not work"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'MasterList.py'))
except Exception, e:
    print e
    print "MasterList.py not found"
    print "some functions may not work"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'kids_diag.py'))
except Exception, e:
    print e
    print "kids_diag.py not found"
    print "some functions may not work"

try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'spikedetect.py'))
except Exception, e:
    print e
    print "spikedetect.py not found"
    print "IQBT despiking will not work"

ip = IPython.get_ipython()
execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'amkid/mars_styling.py'), ip.ns_table['user_global'], ip.ns_table['user_global'])

"""
if '_boa_get' not in globals().keys():
    info('adding data.Unit to Boa get')
    global _boa_get
    global get
    _boa_get = get

    def get(*args, **kwargs):
        data.Unit = 'phi [rad]'
        _boa_get(*args, **kwargs)
"""
        
@timeout(5)    
def pathExists(p):
    return os.path.exists(p)
    
    
def getIndirLocations():
    """
    Get a list of all available BoaConfig.indir locations.
    
    @return: List of all available indir locations {locationName: location_path}.
    @rtype: dict
    
    """
    indir_locations = {'testing': '/home/amkid/sim_rawdata_paruma',
                       'mpifrdata': '/home/amkid/rawdata_mpifr',
                       'apexdata': '/apexdata/rawdata/T-0112.F-9992A-2023'
                       }

    # we use our pathExists function here, which times out after 5 seconds in case
    # a network resource is not available
    available_indir_locations = []
    for l, p in indir_locations.items():
        info('Adding %s...' % p)
        from subprocess32 import check_call, TimeoutExpired
        try:
            # test if path exists and is reachable with a timeout of 3 seconds
            # we use check_call (instead of os.path.exists), as a timout is needed, 
            # as NFS mounts might not be reachable
            check_call(['test', '-d', p], timeout=3)
        except (TimeoutExpired, CalledProcessError) as e:
            warn('Archive "%s" not reachable!' % p)
            continue
        else:
            info('Archive "%s" added.' % p)
            available_indir_locations.append((l, p))
                                             
    available_indir_locations = dict(available_indir_locations)
            
    return available_indir_locations


def setInDirLocation(location='apexdata'):
    """
    Set the BoaConfig.inDir from a set of given locations.
    
    @param location: Location name. One of 'testing' (default) or 'archive'
    @type location: str
    @return: void
    
    """
    info('Setting indir location to "%s"...' % location)
    locations = getIndirLocations()
    if location.lower() in locations.keys():
        indir(locations[location.lower()])
    else:
        raise RuntimeError("Unknown indir location '%s'!\nAllowed locations:%s" % (location, str(locations.keys())))
    info('Indir set.')

#BONN    
# Set up input directory:
#indir(os.path.join(os.getenv('HOME'),'rawdata'))
#setInDirLocation('apexdata')
if os.getenv('HOME').split("/")[2] == 'amkid':
#    #indir('/home/amkid/sim_rawdata_paruma')
#    #indir('/data2/rawdata_tmp/T-0112.F-9992A-2023')
#    #indir('/homes/t-0115.f-9992a-2025/rawdata')
    indir('/apex-archive/rawdata/T-0117.F-9992A-2026')
else:
    myaccount=str.upper(os.getenv('USER'))
    myindir="/apex-archive/rawdata/"+myaccount
    indir(myindir)

#    #indir(os.path.join(os.getenv('HOME'),'rawdata')) #APEX
#    #indir('/display_data/rawdata/T-0115.F-9992A-2025') #current storage as of Nov 2025
#    indir('/homes/t-0115.f-9992a-2025/rawdata')  #OLD

## APEX
#if os.getenv('HOME').split("/")[2] == 'amkid':
#    indir('/display_data/rawdata/T-0117.F-9992A-2026/')
#else:
#    myaccount=str.upper(os.getenv('USER'))
#    myindir="/apexdata/rawdata/"+myaccount
#    indir(myindir)



try:
    execfile(os.path.join(os.getenv('BOA_HOME_AMKID'), 'analysis.py'))
except Exception, e:
    print e
    print "analysis.py not found"
    print "some functions may not work"
    
########################################################################

print "------------------------------------------------"
print " "
print " \033[1;32mQuick reduction functions for AMKID:\033[0m"
print " - \033[1mredfoc\033[0m(ScanNr,fe): Focus"
print " - \033[1mredmfoc\033[0m(ScanNr,fe,tau): Focus (mapping mode)"
print " - \033[1mredpnt\033[0m(ScanNr,fe,tau,sm): Pointing map"
print " - \033[1mredweak\033[0m(ScanNr,fe,tau,sm): Pointing on weak source"
print " - \033[1mredcal\033[0m(ScanNr,fe,tau): Map for calibration/pointing"
#print " - redmap(ScanNr,tau,ra1,ra2,dec1,dec2)  : Map on target source"
print " - \033[1mredscans\033[0m([ [Scan1,Scan2,...],calscans,calfile,taufile,modelfile,oldmap):"
print "       Reduces multiple scans and coadds the result"
print " - \033[1mredbeam\033[0m(ScanNr,fe): Beam map reduction"
#print " - \033[1mredOTFbeam\033[0m([Scan1,Scan2..],doPlot=0/1):"
#print "       Create beam map from a list of OTF scans covering different parts of the "
#print "       array."
#print " - \033[1mredKidCal\033[0m(ScanNr,febe,doPlot=0/1): Kid calibration reduction"
#print " - \033[1mredPolGridTempCal\033[0m(ScanNr,febe,doPlot=0/1,c1=X,c2=Y):"
#print "       Pol Grid Temperture calibration reduction"
print
print("To find AMKIDs related functions use \033[1mfindAMKIDFunctions\033[0m().")
print("To get detailed help For any function, type 'functionName?'.")

#print " For the reduction of data from other frontends, use the appropriate"
#print " of the following macros:"
#print " * SABOCA:  \33[1;34mexecfile('%s/saboca/saboca.boa')\33[0m" %generalBoaDir 
#print " * ARTEMIS: \33[1;34mexecfile('%s/artemis/artemis.boa')\33[0m" %generalBoaDir 
#print " * CHAMP+:  \33[1;34mexecfile('%s/champ/champ.py')\33[0m" %generalBoaDir 


#########################################################################################
#                             FUNCTIONS START HERE
#########################################################################################


def redfoc(ScanNr=0,fe='LFA'):
    '''
    Reduces a focus measurement with AMKID.
    Optional parameter: Scan number, otherwise it takes the current scan.
    '''
    if (ScanNr != 0):
        if (str.upper(fe) == 'LFA') or (str.upper(fe)=='HFA' ):
            try:
                getPhi(ScanNr,fe=fe)
                tst=None
            except:
                tst=1
        
        if  fe == 'Laboca' or fe == 'LABOCA' or fe == 'laboca':
            tst = read(str(ScanNr))
        if tst:
            print "Scan %s not readable"%(str(ScanNr))
            return
        execfile(os.getenv('BOA_HOME_AMKID')+'/reduce_focus.py')
    else:
        print "No scan defined."


def redmfoc(ScanNr=0,fe='LFA',tau=0.0,subscans=[],processData=True,fsweep=None):
    '''
    Reduces a focus observed in spiral-mapping mode with LABOCA.
    Required parameter: Scan number
    Optional parameters: zenith opacity, list of subscans
    '''
    if not ScanNr:
        print "No scan defined."
        return
    tau0 = tau
    redfocmap(scan=ScanNr,fe=fe,tau=tau,subscans=subscans,displayMaps=True,displayFit=True,processData=processData,fsweep=fsweep)
    #if (tau0 == 0.0):
    #    print 'No opacity value passed, using ATM model: tau_z = %3.3f'%(tau)	
    
    
def redpnt(ScanNr=0,fsweep=None,fe='LFA',chains=[],tau=0.0,sm=0.0,size=150,oversamp=5.,doPlot=1,doPrint=True,doDf=True,doFlagSpeed=True,model=False,doReturn=False,myrcp=None):
    '''
    Reduces a LABOCA map on a pointing source.
    Optional parameters: Scan number, zenith opacity
    '''
    if (ScanNr != 0):
        if (str.upper(fe) == 'LFA') or (str.upper(fe)=='HFA' ):
            try:
                if doDf == False:
                    if fsweep:
                        getPhi(ScanNr,fsweep,fe=fe)
                    else:
                        getPhi(ScanNr,fe=fe)
                else:
                    if fsweep:
                        good,medium,bad=getDf(ScanNr,fsweep,fe=fe,returnPhaseSetting=True)
                    else:
                        good,medium,bad=getDf(ScanNr,fe=fe,returnPhaseSetting=True)
                tst=None
            except:
                tst=1
        #if (febe == 'LFA') or (febe=='lfa'):
        #    tst=read(str(ScanNr),febe='AMKID870-AMKID870BE')
        #if (febe == 'HFA') or (febe=='hfa'):
        #    tst=read(str(ScanNr),febe='AMKID350-AMKID350BE')
        if  fe == 'Laboca' or fe == 'LABOCA' or fe == 'laboca':
            tst = read(str(ScanNr))
        
        if tst:
            print "Scan %s not readable"%(str(ScanNr))
            return
        if doFlagSpeed:
            data.flagSpeed (below=20.)
        data.flagSpeed(above=400.)
        data.flagAccel(above=800.)
        

	tau0 = tau
        if (tau0 == 0.0):
	    tau = scanTau()
	    print 'No opacity value passed, using ATM model: tau_z = %3.3f'%(tau)

        if fe == 'Laboca' or fe== 'LABOCA' or fe == 'laboca':
            CntstoV(data)
            data.Data *= np.array(VtoJy,'f')

        if (str.upper(fe) == 'LFA') or (str.upper(fe) == 'HFA'):
            flagEmpty()
            #checkMJD()
            calibrateAMKID(fe=fe)
            if len(chains) !=0:
                if str.upper(fe) == 'HFA':
                    flagC(range(1,20*800+1))
                    for c in chains:
                        unflagC(range((c-1)*800+1,c*800+1))
                if str.upper(fe) == 'LFA':
                    flagC(range(1,4*880+1))
                    for c in chains:
                        unflagC(range((c-1)*880+1,c*880+1))
                flagEmpty()
            flagNaN()

        febe=data.BolometerArray.FeBe
        ## Apply Tau corr
        data.correctOpacity(tau)

        nr=len(data.ScanParam.MJD)
        nr = nr -1
        mjdref = (data.ScanParam.MJD[nr]-data.ScanParam.MJD[0])/2.+data.ScanParam.MJD[0]
 
        febe=data.BolometerArray.FeBe

        if febe == 'LABOCA-ABBA':
            execfile(os.path.join(os.getenv('BOA_HOME_LABOCA'), 'cabling.py'))
            execfile(os.path.join(os.getenv('BOA_HOME_LABOCA'), 'Laboca-RCPs.py'))
            BoaConfig.rcpPath = os.getenv('BOA_HOME_RCP')
            flagC(resistor)
            flagC(sealed_may07)
            cro=getLabocaCross(mjdref)
            rcp=getLabocaRCP(mjdref)
        else:
            BoaConfig.rcpPath = os.getenv('BOA_LOCAL_RCP')
            if not myrcp:
                rcp=getMKIDsRCP(mjdref)
            else: 
                try:
                    rcp=myrcp
                except:
                    print("RCP %s"/"%s not found"%(BoaConfig.rcpPath,rcp))
                    return
            #cro=getAMKIDCross(mjdref)

        
        messages.info('Updating RCP to %s'% rcp)
        try:
            updateRCP(rcp)
        except:
            updateRCP_Local(rcp)
        flagRCP(rcp)
        #flat()
        
        data.zeroStart()

        source_name = data.ScanParam.Object
        if str.upper(source_name) in ['JUPITER','VENUS','SATURN']:
            data.flagPosition(radius=60,flag=8)
        else:
            data.flagPosition(radius=30,flag=8)
        if febe == 'LABOCA-ABBA':
            medianNoiseRemoval(chanRef=-1,factor=0.9,nbloop=3)              
        if febe == 'AMKID870-AMKID870BE':
            #correlChains(chanRef=-2,factor=0.95,nbloop=3)
            correlChips(chanRef=-2,factor=0.95,nbloop=3)
        if febe == 'AMKID350-AMKID350BE':
            #try:
            #    bad=flagBeamShape('CalFiles/beam_map_25099_HFA_combined.csv',maxE=1.5,minC=0.3)
            #    flagC(bad)
            #except:
            #    pass
            #if len(chains) == 0:
            #    mychains=range(1,6)
            #    correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
            #    mychains=range(6,11)
            #    correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
            #    mychains=range(11,16)
            #    correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
            #    mychains=range(16,21)
            #    correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
            #else:
            #    correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=chains)
            correlChips(chanRef=-2,factor=0.95,nbloop=3)

        #data.flagFractionRms(ratio=5)
        #despike()
        #data.flattenFreq(below=0.2,hiref=0.3)
        base(order=1)
        weight()
        unflag(flag=8)

        if size == 150 and str.upper(febe) == 'AMKID350-AMKID350BE':
            size = 100.
        
        if size > 0:
            mapping(oversamp=oversamp,sizeX=[-1*size,size],sizeY=[-1*size,size],noPlot=1)
        else:
            mapping(oversamp=oversamp,noPlot=1)

        if sm >0:
            smoothBy(sm)
        elif sm == 0 and str.upper(fe) == 'LFA':
            smoothBy(5.0)
        elif sm == 0 and str.upper(fe) == 'HFA':
            smoothBy(2.5)
        else:
            pass
        
        try:
            if size > 0:
                data.solvePointingOnMap(plot=doPlot,radius=size/2.)
            else:
                data.solvePointingOnMap(plot=doPlot,radius=1000)
            result = data.PointingResult
            daz=result['gauss_x_offset']['value']
            delev=result['gauss_y_offset']['value']
            if doPrint == True:
                print("############################################")
                print("Phase setting")
                print("Good: %5.1f - Medium: %5.1f - Bad: %5.1f"%(good,medium,bad))
                print("############################################")
                print("# pcorr %5.1f, %5.1f "%(daz,delev))
                print("############################################")
            if doReturn == True:
                return daz,delev,tau,good,medium,bad

        except:
            if doPlot == 1:
                display(aspect=1)
            print "Failed to fit source"
            if doReturn == True:
                return -1,-1,tau,good,medium,bad
        
        if (tau0 == 0.0):
            print 'No opacity value passed, using ATM model: tau_z = %3.3f'%(tau)
            if doReturn == True:
                return daz,delev,tau,good,medium,bad
    else:
        print "No scan defined."


def redweak(ScanNr=0,fsweep=None,fe='LFA',chains=[],tau=0.0,sm=0.0,size=150,oversamp=5.,doPlot=1,model=None,subtract=False,extremeFilter=False,flagJumps=False,writeSummary=False,skipMapping=False, debug=False):
    '''
    Reduces a AMKID  map on a weak pointing source.
    Optional parameters: Scan number, zenith opacity
    '''
    
    if (ScanNr != 0):
        if (str.upper(fe) == 'LFA') or (str.upper(fe)=='HFA' ):
            try:
                if fsweep:
                    good,medium,bad=getDf(ScanNr,fsweep,fe=fe,flagJumps=flagJumps,returnPhaseSetting=True)
                else:
                    good,medium,bad=getDf(ScanNr,fe=fe,flagJumps=flagJumps,returnPhaseSetting=True)
                    
                tst=None
            except:
                tst=1
        #if (febe == 'LFA') or (febe=='lfa'):
        #    tst=read(str(ScanNr),febe='AMKID870-AMKID870BE')
        #if (febe == 'HFA') or (febe=='hfa'):
        #    tst=read(str(ScanNr),febe='AMKID350-AMKID350BE')
        if  fe == 'Laboca' or fe == 'LABOCA' or fe == 'laboca':
            tst = read(str(ScanNr))
        
        if tst:
            raise RuntimeError("Scan %s not readable"%(str(ScanNr)))
            return

        #flag high an low speed in case telescope was scaning
        azvel=data.ScanParam.get('azspeed')
        elvel=data.ScanParam.get('elspeed')
        if np.nanmax(np.sqrt(azvel**2+elvel**2)) > 20.:
            data.flagSpeed (below=20.)
            data.flagSpeed(above=400.)
            data.flagAccel(above=800.)

        if fe == 'Laboca' or fe== 'LABOCA' or fe == 'laboca':
            CntstoV(data)
            data.Data *= np.array(VtoJy,'f')

        if (str.upper(fe) == 'LFA') or (str.upper(fe) == 'HFA'):
            flagEmpty()
            #checkMJD()
            calibrateAMKID(fe=fe)
            if len(chains) !=0:
                if str.upper(fe) == 'HFA':
                    flagC(range(1,20*800+1))
                    for c in chains:
                        unflagC(range((c-1)*800+1,c*800+1))
                if str.upper(fe) == 'LFA':
                    flagC(range(1,4*880+1))
                    for c in chains:
                        unflagC(range((c-1)*880+1,c*880+1))
                flagEmpty()
            flagNaN()

        febe=data.BolometerArray.FeBe
        

        nr=len(data.ScanParam.MJD)
        nr = nr -1
        mjdref = (data.ScanParam.MJD[nr]-data.ScanParam.MJD[0])/2.+data.ScanParam.MJD[0]
 

        if febe == 'LABOCA-ABBA':
            execfile(os.path.join(os.getenv('BOA_HOME_LABOCA'), 'cabling.py'))
            execfile(os.path.join(os.getenv('BOA_HOME_LABOCA'), 'Laboca-RCPs.py'))
            BoaConfig.rcpPath = os.getenv('BOA_HOME_RCP')
            flagC(resistor)
            flagC(sealed_may07)
            cro=getLabocaCross(mjdref)
            rcp=getLabocaRCP(mjdref)
        else:
            BoaConfig.rcpPath = os.getenv('BOA_LOCAL_RCP')
            rcp=getMKIDsRCP(mjdref)
            #cro=getAMKIDCross(mjdref)

        messages.info('Updating RCP to %s'% rcp)
        updateRCP(rcp)
        flagRCP(rcp)
        
        data.zeroStart()

        messages.info('Retrieving atmospheric parameters...')
        scanel  = fStat.f_mean(data.ScanParam.El)
        tau0 = tau
        if (tau0 == 0.0):
	    tau = scanTau()
            print('No opacity value passed, using ATM model: tau_z = %3.3f'%(tau))
        else:
            print('Using opacity tau_z = %3.3f provided by user'%(tau))
        taucorr = np.exp(tau/np.sin(scanel * np.pi / 180.))
        tsky=scanTsky(febe=febe)

        if debug==True:
            return

        # Model flagging or subtraction before correlchains:
        if model:
            if subtract == False:
                print("flagging source Model") 
                data.flagSource(threshold=0.1,model=model)
            else:
                print("subtracting source Model")
                modeln=copy.deepcopy(model)
                modeln.Data = -1*model.Data/taucorr
                addSource_fixed(data,model=modeln)
                
        #if str.upper(fe) == 'LFA':
        #    data.flagPosition(radius=30,flag=8)

        if febe == 'LABOCA-ABBA':
                medianNoiseRemoval(chanRef=-1,factor=0.9,nbloop=3)
        if febe == 'AMKID870-AMKID870BE':
                correlChains(chanRef=-2,factor=0.95,nbloop=3)
        if febe == 'AMKID350-AMKID350BE':
            #bad=flagBeamShape('CalFiles/beam_map_25099_HFA_combined.csv',maxE=1.5,minC=0.3)
            #flagC(bad)
            if len(chains) == 0:
                mychains=range(1,6)
                correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
                mychains=range(6,11)
                correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
                mychains=range(11,16)
                correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
                mychains=range(16,21)
                correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
            else:
                correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=chains)

        base(order=3)
        #data.flagFractionRms(above=5,below=3)
        # VWO this does not flag >5*median and <median/3, it flags
        # >10*median and <median/10. Fix:
        data.flagFractionRms(ratio=5, above=True, below=False)
        data.flagFractionRms(ratio=3, above=False, below=True)
        #despike()


        rms=np.zeros(4,'f')
        _,_,n_kids=getFebe(fe)
        usedbolos = data.BolometerArray.checkChanList([])

        # VWO: info so we don't think it's stuck
        messages.info('Flattening low frequencies in timelines...')

        if str.upper(fe) == 'LFA':
            if model and subtract == True:
                if extremeFilter == True:
                    #data.blankFreq(below=3.0)
                    data.flattenFreq(below=3.0,hiref=5.0)
                else:
                    data.flattenFreq(below=1.0,hiref=2.0)
                #print("now despiking")
                nrspikes=data.despike()
            elif model:
                #print("now despiking")
                nrspikes=data.despike()
                data.flattenFreq(below=0.5,hiref=0.7)
            else:
                data.flattenFreq(below=0.5,hiref=0.7)
                nrspikes = 0
            base(order=1)
            nch = 1
            
            
        if str.upper(fe) == 'HFA':
            data.flattenFreq(below=0.2,hiref=0.35)
            #data.blankFreq(below=1.5)
            base(order=1)
            nch = 5
            nrspikes=0

        data._DataAna__statistics() 

        for i in range(1,5):
            start = (i-1)*n_kids*nch+1
            stop = i*n_kids*nch
            mask=np.where((usedbolos >=  start) * (usedbolos <= stop))[0]
            if len(mask) > 0:
                rms[i-1]=np.nanmedian(data.getChanListData('rms',usedbolos[mask]))
            
       
        dt=np.nanmedian(data.ScanParam.get('deltat'))
        sensi = rms * np.sqrt(dt) * 1000
       
        if model and subtract == True:
            print("adding source Model back to data") 
            modelp = copy.deepcopy(modeln)
            modelp.Data = -1*modeln.Data
            data.addSource(model=modelp)
            #addSource_fixed(data,model=modelp)
        messages.info('Projecting fluxes outside the atmosphere...')
        if (tau0 == 0.0):
            print 'No opacity value passed, using ATM model: tau_z = %3.3f'%(tau)
        data.correctOpacity(tau)
        data._DataAna__statistics() 
        #data.slidingWeight()
        weight()
        #return
        unflag(flag=8)

        if writeSummary:
            outname="%s-%s-%i_summary.txt"%(fe,data.ScanParam.Object,data.ScanParam.ScanNum)
            f=open(outname,'w')
            f.write("#scanel,taucorr,tsky,tint[s],sensi[0],sensi[1],sensi[2],sensi[3],goodPhase[%], mediumPhase[%], badPhase[%], nrspikes\n")
            f.write("%4.1f,%5.2f,%5.1f,%7.1f,%4.1f,%4.1f,%4.1f,%4.1f,%4.1f,%4.1f,%4.1f,%i\n"%(scanel,taucorr,tsky,np.max(data.ScanParam.get('MJD')),sensi[0],sensi[1],sensi[2],sensi[3],good,medium,bad,nrspikes))
            f.close()

        # ===========================
        # HERE REDUCTION IS FINISHED,
        # NEXT COMES POINTING INFO
        # ONLY NEEDED FOR WEAK POINTING
        # SOURCES
        # ==========================        
            
        # Non-pointing sources should end here, mapping is done in reduction script
        # in the intended coordinate system
        if skipMapping == True:
            print("############################################")
            print("Phase setting (Tsky = %.1f K):"%(tsky))
            print("Good: %5.1f - Medium: %5.1f - Bad: %5.1f"%(good,medium,bad))
            print("############################################")
            for i in range(4):
                if sensi[i] > 0:
                    print("# Sensitivity %s chip %i: %6.1f mJy sqrt(s)"%(fe,i+1,sensi[i]))
            print("############################################")
            return

        # And for pointing sources continue...
        if size == 150 and str.upper(febe) == 'AMKID350-AMKID350BE':
            size = 100.
        
        if size > 0:
            mapping(oversamp=oversamp,sizeX=[-1*size,size],sizeY=[-1*size,size],noPlot=1)
        else:
            mapping(oversamp=oversamp,noPlot=1)
        
        
        if sm >0:
            smoothBy(sm)
        elif sm == 0 and str.upper(fe) == 'LFA':
            smoothBy(10.0)
        elif sm == 0 and str.upper(fe) == 'HFA':
            smoothBy(4.5)
        else:
            pass
        
        
        try:
            if size > 0:
                data.solvePointingOnMap(plot=doPlot,radius=30)
            else:
                data.solvePointingOnMap(plot=doPlot,radius=300)
            result = data.PointingResult
            daz=result['gauss_x_offset']['value']
            delev=result['gauss_y_offset']['value']
            print("############################################")
            print("Phase setting (Tsky = %.1f K):"%(tsky))
            print("Good: %5.1f - Medium: %5.1f - Bad: %5.1f"%(good,medium,bad))
            print("############################################")
            for i in range(4):
                if sensi[i] > 0:
                    print("# Sensitivity %s chip %i: %6.1f mJy sqrt(s)"%(fe,i+1,sensi[i]))
            print("############################################")
            print("# pcorr %5.1f, %5.1f "%(daz,delev))
            print("###########################################")
            

        except:
            display(aspect=1)
            print("############################################")
            print("Phase setting (Tsky = %.1f K):"%(tsky))
            print("Good: %5.1f - Medium: %5.1f - Bad: %5.1f"%(good,medium,bad))
            print("############################################")
            for i in range(4):
                if sensi[i] > 0:
                    print("# Sensitivity %s chip %i: %6.1f mJy sqrt(s)"%(fe,i+1,sensi[i]))
            print("############################################")
            print("Failed to fit source")
            print("############################################")
    else:
        print "No scan defined."
    #t2=time()
    #runtime=t2-t1
    #print("Run Time: %5.1f sec"%(runtime))



def _auxsaveData(fileName='MarsData', overwrite=False):
    """
    FUNCTION IMPLEMENTED TO SAVE DATA BEFORE MEDIAN SUBTRACTION
    ON REDUCTION OF SCIENCE DATA. 
    MAY NOT WORK FOR OTHER PURPOSES!

    DES: save current data object to pickle and NumPy files
    INP: (string) fileName: basename of the output files
    """

    ipy = IPython.get_ipython()
    ipy_globals = ipy.ns_table['user_global']

    # Pickle data object with metadata, but without large arrays
    # (see BoaDataEntity::DataEntity:__getstate__() method).
    pckFileName = os.path.join(BoaConfig.outDir, '%s.pck' % (fileName))
    if os.path.exists(pckFileName) and overwrite==False:
        ans = raw_input('\nFile %s exists. Overwrite y|N? ' % (pckFileName))
        if ans.lower() != 'y':
            return

    try:
        print('         Saving dictionary to pickle...')
        with open(pckFileName, 'wb+') as fd:
            cPickle.dump(ipy_globals['data'], fd, 2)
    except:
        raise RuntimeError("Could not open file\n%s" % (pckFileName))
        #return

    # Write large arrays via np.save() because it is too
    # large to be pickled and the size on disk would increase
    # about 2.5 times compared to the intrinsic size.
    try:
        npyFileName = os.path.join(BoaConfig.outDir, '%s_Data.npy' % (fileName))
        print('         Saving NumPy files...')
        with open(npyFileName, 'w+') as fd:
            np.save(fd, ipy_globals['data'].Data, allow_pickle=False)

        npyFileName = os.path.join(BoaConfig.outDir, '%s_DataWeights.npy' % (fileName))
        with open(npyFileName, 'w+') as fd:
            # try and save disk space
            if np.all(ipy_globals['data'].DataWeights==1):
                np.save(fd, np.array([1]), allow_pickle=False)
            else:
                np.save(fd, ipy_globals['data'].DataWeights, allow_pickle=False)
    except:
        raise RuntimeError("Could not open file\n%s" % (npyFileName))
        #return

    print("         Current data successfully written to:\n         %s" % (fileName))



def _auxloadData(fileName='MarsData'):
    """
    FUNCTION IMPLEMENTED TO LOAD DATA BEFORE MEDIAN SUBTRACTION
    FOR REDUCTION OF SCIENCE DATA IN ITERATIONS 2 AND 3.
    MAY NOT WORK FOR OTHER PURPOSES!

    DES: load current data object from pickle and NumPy files
    INP: (string) fileName: basename of the output files
    """

    ipy = IPython.get_ipython()
    ipy_globals = ipy.ns_table['user_global']
    #if 'data' in ipy_globals.keys():
    #    # del(ipy_globals['data'])
    #    # VWO: don't delete, otherwise some functions don't work

    pckFileName = os.path.join(BoaConfig.outDir, '%s.pck' % (fileName))
    print('         Loading dictionary from pickle...')
    try:
        # VWO: don't replace 'data', just its dictionary keys.
        with open(pckFileName, 'rb') as fd:
            mynewdata = cPickle.load(fd)
        for dictkey in mynewdata.__dict__.keys():
            ipy_globals['data'].__dict__[dictkey] = mynewdata.__dict__[dictkey]
    except:
        raise RuntimeError("Could not open file\n%s" % (pckFileName))
        #return

    # Add data back since it doesn't exist in the pickle.
    try:
        npyFileName = os.path.join(BoaConfig.outDir, '%s_Data.npy' % (fileName))
        print('         Loading NumPy files...')
        ipy_globals['data'].Data = np.load(npyFileName)

        npyFileName = os.path.join(BoaConfig.outDir, '%s_DataWeights.npy' % (fileName))
        loaded = np.load(npyFileName)
        if len(loaded)==1 and loaded[0]==1:
            ipy_globals['data'].DataWeights = np.ones_like(ipy_globals['data'].Data)
        else:
            ipy_globals['data'].DataWeights = loaded
    except:
        raise RuntimeError("Could not open file\n%s" % (npyFileName))


def redscience(ScanNr=None, fsweep=None, fe='LFA', chains=[], src=None, tau=0.0, model=None, subtract=False, extremeFilter=False, correctbeam=True, flagJumps=False, writeSummary=False, debug=False):
    '''
    VWO: Redweak implementation specific for science reduction scripts.
    Reduces a AMKID data down to Jy/beam units projected outside the atmosphere. Mapping should be handled externally with intended coord system.
    Accepts scan number for full reduction. If model==None, iteration is 1 and will dump full data structure before correlated noise removal
    for future iterations. If model!=None, parameter "src" (for source) must be provided to retrieve data before correlated noise removal.
    Returns nothing; the global "data" variable is manipulated.
    '''
    assert ScanNr!=None and type(ScanNr)==int, 'ScanNr must be an existing scan number in the current input directory.'
    assert str.upper(fe)=='LFA' or str.upper(fe)=='HFA', 'This function is implemented for AMKID only: fe="LFA" or fe="HFA".'
    if model!=None:
        assert src!=None, 'For iterations that are not 1, a valid (reduced in iter 1) source name must be provided in argumnent "src".'

    dataObjname = '%s-%s'%(str.upper(fe), src)
    if flagJumps==True:
        dataObjname += '-flagJumps'
    dataObjname += '-%s'%ScanNr
    
    # Create folders if first run
    if os.path.exists('ReducedFiles') == False:
        os.makedirs('ReducedFiles')
    if os.path.exists('ReducedFiles/DataObjects/') == False:
        os.makedirs('ReducedFiles/DataObjects/')
        
    fulldataObjname = 'ReducedFiles/DataObjects/%s'%(dataObjname)
    dataiter1exists = len(glob(fulldataObjname+'*'))==3  # Data, dataWeight and pickle
    
    # If calibrated data from iteration 1 does not exist, then reduce from scratch
    if dataiter1exists==False:
        messages.info('Reducing data from start (IQ)...')
        try:
            if fsweep!=None:
                good,medium,bad=getDf(ScanNr,fsweep,fe=fe,flagJumps=flagJumps,returnPhaseSetting=True)
            else:
                good,medium,bad=getDf(ScanNr,fe=fe,flagJumps=flagJumps,returnPhaseSetting=True)
        except:
            print("Scan %s not readable"%(str(ScanNr)))
            return

        # Save phase setting
        data.ScanParam.PhaseSetting = [good, medium, bad]

        #flag high an low speed in case telescope was scaning
        azvel=data.ScanParam.get('azspeed')
        elvel=data.ScanParam.get('elspeed')
        if np.nanmax(np.sqrt(azvel**2+elvel**2)) > 20.:
            data.flagSpeed (below=20.)
            data.flagSpeed(above=400.)
            data.flagAccel(above=800.)

        if str.upper(fe) == 'LFA' or str.upper(fe) == 'HFA':
            flagEmpty()
            calibrateAMKID(fe=fe)  # Jy/beam after this
            if len(chains) !=0:
                if str.upper(fe) == 'HFA':
                    flagC(range(1,20*800+1))
                    for c in chains:
                        unflagC(range((c-1)*800+1,c*800+1))
                if str.upper(fe) == 'LFA':
                    flagC(range(1,4*880+1))
                    for c in chains:
                        unflagC(range((c-1)*880+1,c*880+1))
                flagEmpty()
            flagNaN()

        # Update RCP
        # AMKID-only
        febe = data.BolometerArray.FeBe
        mjdref = (data.ScanParam.MJD[-1]+data.ScanParam.MJD[0])/2.  # mid-scan
        BoaConfig.rcpPath = os.getenv('BOA_LOCAL_RCP')
        rcp=getMKIDsRCP(mjdref)
        messages.info('Updating RCP to %s'% rcp)
        updateRCP(rcp)  # this affects data object, not variables in MARS
        flagRCP(rcp)  # so they will be saved as well.. no?
        
        messages.info('Applying data zero start...')
        data.zeroStart()
        messages.info('Retrieving atmospheric parameters...')
        scanel  = fStat.f_mean(data.ScanParam.El)
        tau0 = tau
        if (tau0 == 0.0):
	    tau = scanTau()
            print('No opacity value passed, using ATM model: tau_z = %3.3f'%(tau))
        else:
            print('Using opacity tau_z = %3.3f provided by user'%(tau))
        taucorr = np.exp(tau/np.sin(scanel * np.pi / 180.))
        tsky=scanTsky(febe=febe)

        # Save these atmospheric parameters in data.ScanParam.ATMparams, bc they take
        # way too long to retrieve.
        #                           elev,   givenTau, retrievedTau, taucorr, tsky
        data.ScanParam.ATMparams = [scanel, tau0,     tau,          taucorr, tsky]

        # and dump data structure so we don't have to do EVERYTHING again for iter 2 and 3,
        # only model flagging/substraction and correlated noise removal.
        messages.info('Saving data structure for next iterations...')
        if not os.path.exists('ReducedFiles/DataObjects/'):
            os.makedirs('ReducedFiles/DataObjects/')
        dataObjname = '%s-%s'%(str.upper(fe), data.ScanParam.Object)
        if flagJumps==True:
            dataObjname += '-flagJumps'
        dataObjname += '-%s'%ScanNr
        fulldataObjname = 'ReducedFiles/DataObjects/%s'%(dataObjname)
        _auxsaveData(fulldataObjname, overwrite=True)

    # Otherwise, we should be able to retrieve the data structure before correlated
    # noise removal saved in iter 1.
    # The src parameter is needed for this.
    else:
        messages.info('Retrieving saved data:')
        print('         Scan %i on %s...'%(ScanNr, src))
        try:
            initData()
            dataObjname = '%s-%s'%(str.upper(fe), src)
            if flagJumps==True:
                dataObjname += '-flagJumps'
            dataObjname += '-%s'%ScanNr
            fulldataObjname = 'ReducedFiles/DataObjects/%s'%(dataObjname)
            _auxloadData(fulldataObjname)
        except:
            raise RuntimeError('Scan %i (source %s) does not have a saved iteration 1 reduction before correlated noise removal.'%(ScanNr, src)+\
                               '\nTry running iteration 1 through all scans first.')

        # redefine these from saved values in iter 1
        good, medium, bad = data.ScanParam.PhaseSetting
        scanel, tau0, tau, taucorr, tsky = data.ScanParam.ATMparams

        # Update RCP
        # AMKID-only
        febe = data.BolometerArray.FeBe
        mjdref = (data.ScanParam.MJD[-1]+data.ScanParam.MJD[0])/2.  # mid-scan
        

    if model!=None:
        if subtract == False:
            # this uses flag 8
            print("-----------------------------------------")
            print(">   Flagging source Model (S/N ratio)    ") 
            data.flagSource(threshold=0.1,model=model)
        else:
            print("-----------------------------------------")
            print(">    Subtracting source Model (Flux)     ")
            modeln=copy.deepcopy(model)
            modeln.Data = -1*model.Data/taucorr
            addSource_fixed(data,model=modeln)

    # debug check data is identical when loaded
    if debug==True:
        return

    print("-----------------------------------------")
    # Regardless of model or not and subtract or not, now we do the rest
    messages.info('Performing median subtraction...')
    if febe == 'AMKID870-AMKID870BE':
        correlChains(chanRef=-2,factor=0.95,nbloop=3)
    if febe == 'AMKID350-AMKID350BE':
        if len(chains) == 0:
            mychains=range(1,6)
            correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
            mychains=range(6,11)
            correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
            mychains=range(11,16)
            correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
            mychains=range(16,21)
            correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=mychains)
        else:
            correlChains(chanRef=-2,factor=0.95,nbloop=2,chains=chains)

    messages.info('Removing baseline...')
    # For some reason base() does not work after loading the data in iters 2 and 3...
    # Now it does, it was because loadData() deletes the data variable.
    # Now using _auxloadData() instead.
    base(order=3)
    #data.flagFractionRms(above=5,below=3)
    # VWO this does not flag >5*median and <median/3, it flags
    # >10*median and <median/10. Fix:
    data.flagFractionRms(ratio=5, above=True, below=False)
    data.flagFractionRms(ratio=3, above=False, below=True)
    _,_,n_kids=getFebe(fe)
    usedbolos = data.BolometerArray.checkChanList([])

    # Flattenfreq
    messages.info('Flattening low frequencies in timelines...')
    if str.upper(fe) == 'LFA':
        if model and subtract == True:
            if extremeFilter == True:
                data.flattenFreq(below=3.0,hiref=5.0)
            else:
                data.flattenFreq(below=1.0,hiref=2.0)
            nrspikes=data.despike()
        elif model:
            nrspikes=data.despike()
            data.flattenFreq(below=0.5,hiref=0.7)
        else:
            data.flattenFreq(below=0.5,hiref=0.7)
            nrspikes = 0
        base(order=1)
        #data.polynomialBaseline(order=1)
        nch = 1
            
    if str.upper(fe) == 'HFA':
        data.flattenFreq(below=0.2,hiref=0.35)
        base(order=1)
        #data.polynomialBaseline(order=1)
        nch = 5
        nrspikes=0
    
    # Statistics before projecting outside the atmosphere
    messages.info('Computing statistics at telescope...')
    data._DataAna__statistics()
    rms = np.zeros(4,'f')
    for i in range(1,5):
        start = (i-1)*n_kids*nch+1
        stop = i*n_kids*nch
        mask=np.where((usedbolos >=  start) * (usedbolos <= stop))[0]
        if len(mask) > 0:
            rms[i-1]=np.nanmedian(data.getChanListData('rms',usedbolos[mask]))
    dt = np.nanmedian(data.ScanParam.get('deltat'))
    sensi = rms * np.sqrt(dt) * 1000
    # Save for tests and in case one wants to see them again
    data.ScanParam.ChipSens = sensi  # chip 1, 2, 3, 4
    
    # Model back if subtracted
    if model and subtract == True:
        print("-----------------------------------------")
        print("> Adding source Model (Flux) back to data")
        modelp = copy.deepcopy(modeln)
        modelp.Data = -1*modeln.Data
        data.addSource(model=modelp)
    print("-----------------------------------------")
    
    # Project outside atmosphere
    messages.info('Projecting fluxes outside the atmosphere...')
    if (tau0 == 0.0):
        print 'No opacity value passed, using ATM model: tau_z = %3.3f'%(tau)
    data.correctOpacity(tau)
    
    # Statistics outside the atmosphere
    messages.info('Computing statistics outside the ATM...')
    data._DataAna__statistics()
    weight()
    unflag(flag=8)
    
    messages.info('Correcting AMKID beam:')
    # update beam if required
    if correctbeam:
        if str.upper(fe)=='LFA':
            oldbeam = data.BolometerArray.BeamSize
            newbeam = 16.7  # median 16.668" @June2026
        else:
            oldbeam = data.BolometerArray.BeamSize
            newbeam = data.BolometerArray.BeamSize  # no beam measurement
        data.BolometerArray.BeamSize = newbeam
        print('         Beam: %.3f" -> %.1f"'%(oldbeam, newbeam))

    # Summary if needed
    if writeSummary:
        outname="%s-%s-%i_summary.txt"%(fe,data.ScanParam.Object,data.ScanParam.ScanNum)
        f=open(outname,'w')
        f.write("#scanel,taucorr,tsky,tint[s],sensi[0],sensi[1],sensi[2],sensi[3],goodPhase[%], mediumPhase[%], badPhase[%], nrspikes\n")
        f.write("%4.1f,%5.2f,%5.1f,%7.1f,%4.1f,%4.1f,%4.1f,%4.1f,%4.1f,%4.1f,%4.1f,%i\n"%(scanel,taucorr,tsky,np.max(data.ScanParam.get('MJD')),sensi[0],sensi[1],sensi[2],sensi[3],good,medium,bad,nrspikes))
        f.close()

    # ===========================
    # HERE REDUCTION IS FINISHED
    # PRINT RESULT INFO
    # ===========================
    print("############################################")
    print("Phase setting (Tsky = %.1f K):"%(tsky))
    print("Good: %5.1f - Medium: %5.1f - Bad: %5.1f"%(good,medium,bad))
    print("############################################")
    for i in range(4):
        if sensi[i] > 0:
            print("# Sensitivity %s chip %i: %6.1f mJy sqrt(s)"%(fe,i+1,sensi[i]))
    print("############################################")
    


def redcal(ScanNr=0,fsweep=None,fe='LFA',tau=0.0,sm=0.0,size=150,oversamp=5.,doPlot=1,calscan='',doFlagSpeed=True,outfile=None,rcp=None):
    daz,delev,tau,good,medium,bad=redpnt(ScanNr=ScanNr,fsweep=fsweep,fe=fe,tau=tau,sm=sm,size=size,oversamp=oversamp,doPlot=doPlot,doPrint=False,doDf=True,doFlagSpeed=doFlagSpeed,doReturn=True,myrcp=rcp)

    chan0 = data.BolometerArray.checkChanList ([])[0]
    el = np.nanmedian(data.getChanData('el',chan0))
    taucorr = np.exp(tau/np.sin(el * np.pi / 180.))
    obsfluxmap = data.PointingResult['gauss_peak']['value']
    expect_flux = 0.

    source_name = data.ScanParam.Object
    if str.upper(fe) == 'LFA':
        beam = 16.0
        freq = 345.
        JypK = 40.0

    if str.upper(fe) == 'HFA':
        beam = 8.0
        freq = 810
        JypK=75.0
        

    if source_name in ['Uranus','Neptune','Mars','Saturn','Jupiter','Venus']:
        astrotime,astrodate=getAstroDate(data)
        expect_flux = PlanetFlux(source_name,astrotime,astrodate,beam,freq) 
    else:
        if str.upper(fe) == 'LFA':
            if calibFluxes.has_key(string.upper(source_name)):
                expect_flux = calibFluxes[string.upper(source_name)]
        if str.upper(fe) == 'HFA':
            if calibFluxes350.has_key(string.upper(source_name)):
                expect_flux = calibFluxes350[string.upper(source_name)]

    if expect_flux > 0:
        calcorr = obsfluxmap/expect_flux
        obsflux = obsfluxmap
        percent = 100.0*obsflux/expect_flux

        
        rms=np.nanmedian(data.getChanListData('rms'))*expect_flux/obsfluxmap
        dt=np.nanmedian(data.ScanParam.get('deltat'))


        print("############################################")
        print("Phase setting")
        print("Good: %5.1f - Medium: %5.1f - Bad: %5.1f"%(good,medium,bad))
        print "-------------------------------------------------------------"
        print "FLUX %s:  %7.2f Jy [expected: %7.2f Jy, %6.2f percent]"%(source_name,obsflux,expect_flux,percent)
        print "COUPLING: %5.1f Jy/K"%(JypK*expect_flux/obsfluxmap)
        print "SENSITIVITY: %5.1f mJy sqrt(s) [%5.1f mK sqrt(s)] "%(rms * np.sqrt(dt) * 1000. / taucorr, rms * np.sqrt(dt) * 1000. / taucorr * obsfluxmap/expect_flux/JypK)
        print "-------------------------------------------------------------"
        print "pcorr %5.1f, %5.1f "%(daz,delev)
        print "-------------------------------------------------------------"
        if outfile:
            mjdref=np.nanmean(data.ScanParam.MJD)
            scandate=data.ScanParam.DateObs
            f=open(outfile,'a')
            f.write('%i %s %s %18.12f %4.1f %5.3f %5.3f\n' %(ScanNr,source_name,scandate,mjdref,el,calcorr,taucorr))
            f.close()


def redcal_old(ScanNr=0,fe='LFA',showSen=1,tau=0.0,taufile='',temptau=0,outfile='',Tamb=273.0,doPlot=1,calscan='',doSkyFF=1):
    '''
    Reduces a MKIDs map of primary and secondary calibrator.
    Optional parameters: Scan number, zenith opacity
    '''
    if (ScanNr != 0):
        if (str.upper(fe) == 'LFA') or (str.upper(fe)=='HFA' ):
            try:
                getPhi(ScanNr,fe=fe)
                tst=None
            except:
                tst=1
            flagEmpty()
        
        if  fe == 'Laboca' or fe== 'LABOCA' or fe == 'laboca':
            tst = read(str(ScanNr))

        if tst:
            print "Scan %s not readable"%(str(ScanNr))
            return
        useRadiometer = 0

	if taufile:
            nr=len(data.ScanParam.MJD)
            nr = nr -1
            mjdref = (data.ScanParam.MJD[nr]-data.ScanParam.MJD[0])/2.+data.ScanParam.MJD[0]
            tau = getTau(mjdref,'linear',taufile)
        else:
            if (tau == 0.0):
                tau = scanTau()
                if temptau == 0:
                    useRadiometer = 1
                    print 'No opacity value passed, using ATM model: tau_z = %3.3f'%(tau)
                
	
        

        if fe == 'Laboca' or fe== 'LABOCA' or fe == 'laboca':
            CntstoV(data)
            data.Data *= np.array(VtoJy,'f')
            data.correctOpacity(tau)

        if (str.upper(fe) == 'LFA') or (str.upper(fe) == 'HFA'):
            checkMJD()
            #if calscan:
            #    search = "*_%i_*%s_KidCalspline.dat"%(calscan,str(data.BolometerArray.FeBe))
            #    #print search
            #    globlist=glob(search)
            #    if len(globlist) == 1:
            #        calfile=globlist[0]
            #    else:
            #        print "No calibration file for scan %i found - check scan nr or reduce kid caibration first"%(scan)
            #        return
            #    applyTempCal(calfile,doSkyFF=doSkyFF)
            #    #flagFractionRms()
            #    dataTsky= data._DataAna__computeMedianSignal()
            #    obsTsky = fStat.f_median(np.array(dataTsky))
            #    if (febe == 'LFA') or (febe=='lfa'):
            #        obsTsky = obsTsky-(1-Feff_LFA)*Tamb
            #    if (febe == 'HFA') or (febe=='hfa'):
            #        obsTsky = obsTsky-(1-Feff_HFA)*Tamb
            #    if temptau == 1:
            #        taucorrtemp=applyTempTauCorr()
            #        print 'ATM model: tau_z = %3.3f'%(tau)
            #    #applySkyFlat()
            #    applyJyperK()
            #else:
            calibrateAMKID(fe=fe)
                
            chan0 = data.BolometerArray.checkChanList ([])[0]
            el = data.getChanData('el',chan0)
            med_el = np.median(el)
            tau_los = tau/np.sin(med_el * np.pi / 180.)
            Modeltsky=Tamb*(1-np.exp(-1*tau_los))
            if temptau == 0:
                data.correctOpacity(tau)
                taucorr = np.exp(tau/np.sin(med_el * np.pi / 180.))
            
            

        nr=len(data.ScanParam.MJD)
        nr = nr -1
        mjdref = (data.ScanParam.MJD[nr]-data.ScanParam.MJD[0])/2.+data.ScanParam.MJD[0]

        febe=data.BolometerArray.FeBe

        if febe == 'LABOCA-ABBA':       
            execfile(os.path.join(os.getenv('BOA_HOME_LABOCA'), 'cabling.py'))
            execfile(os.path.join(os.getenv('BOA_HOME_LABOCA'), 'Laboca-RCPs.py'))
            BoaConfig.rcpPath = os.getenv('BOA_HOME_RCP')
            flagC(resistor)
            flagC(sealed_may07)
            cro=getLabocaCross(mjdref)
            rcp=getLabocaRCP(mjdref)
        else:
            BoaConfig.rcpPath = os.getenv('BOA_LOCAL_RCP')
            rcp=getMKIDsRCP(mjdref)
            

        print rcp 
        updateRCP(rcp)
        flagRCP(rcp)
        if febe == 'LABOCA-ABBA':       
            flat()
        
        data.zeroStart()
    
        #data.flagSpeed (below=30.)
        #data.flagSpeed(above=500.)
        #data.flagAccel(above=800.)

        source_name = data.ScanParam.Object

        data.flagPosition(radius=30,flag=8)
        if source_name in ['Saturn','Jupiter']:
            data.flagPosition(radius=50,flag=8)

        
        if source_name not in ['Saturn','Jupiter','Venus']:
            data.flagFractionRms(ratio=5)
            if febe == 'LABOCA-ABBA':
                medianNoiseRemoval(chanRef=-1,factor=0.9,nbloop=3)
            #if febe == 'AMKID870-AMKID870BE' or febe == 'AMKID350-AMKID350BE':
            #    correlChains(chanRef=-2,factor=0.95,nbloop=3)
                
        if febe == 'LABOCA-ABBA':
            medianNoiseRemoval(chanRef=-1,factor=0.9,nbloop=3)
        if febe == 'AMKID870-AMKID870BE' or febe == 'AMKID350-AMKID350BE':
                correlChains(chanRef=-2,factor=0.95,nbloop=3)
                #if  febe == 'LFA'  or  febe=='lfa' or febe == 'AMKID870-AMKID870BE':
                #    chainIndices = range(4)
                #if  febe == 'HFA'  or  febe=='hfa' or febe == 'AMKID350-AMKID350BE':
                #    chainIndices = range(20)
                #for chainIndex in chainIndices:
	 #		correlGroups(chainIndex, chanRef=-2, factor=0.97, nbloop=3,numRepeats=1, groupingThreshold=0.7)

        data.flagFractionRms(ratio=5)
        despike(below=-5,above=5)
        if source_name == 'HLTAU' or source_name == 'VYCma' or source_name == 'V883-ORI':
            if febe == 'LABOCA-ABBA':
                correlbox(data,factor=0.8,nbloop=2)
                correlgroup(data,factor=0.9,nbloop=2)
            
        

        base(order=1)
        data.flagFractionRms(ratio=5)
       

        if source_name not in ['Saturn','Jupiter','Venus']:
            ## FOR P-MODEL ONLY
            #data.flattenFreq(below=0.2,hiref=0.35)
            data.computeWeight()

        data._DataAna__statistics()    
        unflag(flag=8)
        if doPlot == 1:
            mapping(sizeX=[-150,150],sizeY=[-150,150],oversamp=5,noPlot=0,aspect=1)
        else:
            mapping(sizeX=[-150,150],sizeY=[-150,150],oversamp=5,noPlot=1,aspect=1)
        print "Solving for pointing on map..."
        try:
            solvePointingOnMap(plot=0,radius=30)
            obsfluxmap = data.PointingResult['gauss_peak']['value']
            daz=data.PointingResult['gauss_x_offset']['value']
            delev=data.PointingResult['gauss_y_offset']['value']
        except:
            return
            

        expect_flux = 0.

        if source_name in ['Uranus','Neptune','Mars','Saturn','Jupiter','Venus']:
            astrotime,astrodate=getAstroDate(data)
            if string.find(data.BolometerArray.FeBe,'LABOCA') >= 0:
                beam = 19.2
                freq = 345
            if febe == 'AMKID870-AMKID870BE':
                beam = 19.2
                freq = 345
            
	    if febe == 'AMKID350-AMKID350BE':
                beam = 7.8
                freq = 810	
        
            print source_name,astrotime,astrodate,str(beam),str(freq)
            expect_flux = PlanetFlux(source_name,astrotime,astrodate,beam,freq)
        else:
            if febe == 'AMKID870-AMKID870BE':
                if calibFluxes.has_key(string.upper(source_name)):
                    expect_flux = calibFluxes[string.upper(source_name)]
            if febe == 'AMKID350-AMKID350BE':
                if calibFluxes350.has_key(string.upper(source_name)):
                    expect_flux = calibFluxes350[string.upper(source_name)]

        if expect_flux > 0:
            calcorr = obsfluxmap/expect_flux
            obsflux = obsfluxmap


            percent = 100.0*obsflux/expect_flux


            if showSen == 1:
                data.flagPosition(radius=30,flag=8)
                if  str.upper(fe) == 'LFA'  or febe == 'AMKID870-AMKID870BE':
                    data.flagPosition(radius=30,flag=8)
                    chainIndices = range(4)
                    
                if  str.upper(fe) == 'HFA' or febe == 'AMKID350-AMKID350BE':
                    data.flagPosition(radius=15,flag=8)
                #    chainIndices = range(4)
                #for chainIndex in chainIndices:
	 #		correlGroups(chainIndex, chanRef=-2, factor=0.97, nbloop=3,numRepeats=1, groupingThreshold=0.7)
                data.flattenFreq(below=0.5,hiref=0.7)
                despike()
                data._DataAna__statistics()
                if source_name not in ['Saturn','Jupiter','Venus']:
                    data.computeWeight()
                unflag(flag=8)
                if temptau == 1:
                    sensi,wsensi,allsensi,bolos = calcsensitivity(unit='mJy sqrt(s)',scaling=1./(calcorr*np.median(taucorrtemp)))
                else:
                    sensi,wsensi,allsensi,bolos = calcsensitivity(unit='mJy sqrt(s)',scaling=1./(calcorr*np.median(taucorr)))
            #print calcorr,np.median(taucorrtemp)
            if calfile:
                print "-------------------------------------------------------------"
                print "Calibration based on Pol-Grid scan %s"%(calfile)
            print "-------------------------------------------------------------"
            print "FLUX %s:  %7.2f [expected: %7.2f, %6.2f percent]"%(source_name,obsflux,expect_flux,percent)
            print "-------------------------------------------------------------"
            print "pcorr %5.1f, %5.1f "%(daz,delev)
            print "-------------------------------------------------------------"

            if useRadiometer == 1:
                print 'No opacity value passed, using ATM model: tau_z = %3.3f'%(tau)

            if calfile:
                print 'Sky Temperture: %5.1f/%5.1f K (observed/model)'%(obsTsky,modelTsky)
                if fe == 'AMKID870-AMKID870BE':
                    effJyperK=JyperK_LFA/percent*100.
                if fe == 'AMKID350-AMKID350BE':
                   effJyperK=JyperK_LFA/percent*100.
                print 'System Gain: %5.1f Jy/K'%(effJyperK)
            #if showSen == 1:
            #    print "Sensitivity: %6.1f mJy sqrt(s)"%(calcorr*sensi/np.median(taucorrtemp))
            print "-------------------------------------------------------------"

            scandate = data.ScanParam.DateObs
            chan0 = data.BolometerArray.checkChanList ([])[0]
            el = data.getChanData('el',chan0)
            scanel  = fStat.f_mean(el)
            taucorr = np.exp(tau/np.sin(scanel * np.pi / 180.))
            
            if calcorr < 0.3:
                writecal = 0
            else:
                writecal = 1


            if outfile and writecal == 1:
            	out=file(outfile,'a')
                out.write('%i  %s %18.12f %5.3f %5.3f\n' %(ScanNr,scandate,mjdref,calcorr,taucorr))
                out.close()

                of2='long_'+outfile
                out=file(of2,'a')
                out.write('%i %s %s %18.12f %4.1f %5.3f %5.3f %5.3f\n' %(ScanNr,source_name,scandate,mjdref,scanel,tau,calcorr,taucorr))
                out.close()

            
                print '#######################################'
                print "Sensitivity: %6.1f mJy sqrt(s)"%(calcorr*sensi/np.median(taucorrtemp))
                print '#######################################'
        else:
            print "Calibrator not found in secondary_fluxes.py"
        
    else:
        print "No scan defined."


def redcalbeammap(ScanNr=0,calfile=0,febe='LFA',doDespike=0,mode='normal',fsweep=0):
    '''
    Reduces a MKIDs map of primary or secondary calibrator
    and derives RCP, JyperK and Sensitivity for each detector

    Optional parameters: Scan number, PolGird cal scan number, FeBe:
    redcalbeammap(ScanNr=0,calfile=0,febe='LFA')
    '''
    print "START"
    if str.upper(mode) != 'IQ':
    	if (ScanNr != 0):
    	    if (febe == 'LFA') or (febe=='lfa'):
    	        tst=read(str(ScanNr),febe='AMKID870-AMKID870BE')
    	        updateScanParam(data)
    	        flagEmpty()
    	    if (febe == 'HFA') or (febe=='hfa'):
    	        tst=read(str(ScanNr),febe='AMKID350-AMKID350BE')
    	        updateScanParam(data)
    	        flagEmpty()
    	    
    	    
    	    if tst:
    	        print "Scan %s not readable"%(str(ScanNr))
    	        return
    	    
    	
    	    scan = data.ScanParam.ScanNum
    	    mjdref = (data.ScanParam.MJD[-1]+data.ScanParam.MJD[0])/2.  
    	    rcp=getMKIDsRCP(mjdref)        
    	    updateRCP(rcp)
    	
    	    #### APPLY POLGRID TEMP CAL AND SKY FLAT FIELD
    	    if calfile:
    	        if str(calfile)[len(str(calfile))-4:len(str(calfile))] != '.dat':
    	            if str.upper(febe) == 'LFA':
    	                febename = 'AMKID870-AMKID870BE'
    	            else:
    	                febename = 'AMKID350-AMKID350BE'
    	            calnr = np.int(calfile)
    	            search = '*%i_%s_KidCalspline.dat' % (calnr,febename)
    	            #print search
    	            if len(glob(search)) > 0:
    	                calfile = glob(search)[-1]
    	            else:
    	                print "############################################################################"
    	                print "PolGrid calibration file %s not found - reduce PolGrid calibration first"%(calfile)
    	                print "############################################################################"
    	                return
    	
    	        Tamb = data.ScanParam.T_amb
    	        applyTempCal(calfile)
    	        #applySkyFlat()
    	        dataTsky= data._DataAna__computeMedianSignal()
    	        obsTsky = fStat.f_median(np.array(dataTsky))
    	        if (febe == 'LFA') or (febe=='lfa'):
    	            obsTsky = obsTsky-(1-Feff_LFA)*Tamb
    	        if (febe == 'HFA') or (febe=='hfa'):
    	            obsTsky = obsTsky-(1-Feff_HFA)*Tamb
    	    else:
    	        data.Data *= np.array(-1,'f')
    else:
        CalIQ(ScanNr,fsweep)
        data.BolometerArray.FeBe = 'AMKID870-AMKID870BE'

    flagEmpty()
    
    scan = data.ScanParam.ScanNum
    mjdref = (data.ScanParam.MJD[-1]+data.ScanParam.MJD[0])/2.  
    rcp=getMKIDsRCP(mjdref)        
    updateRCP(rcp)
        
    ######## GET OPACITY CORRECTION FROM RADIOMETER
    tau = scanTau(pwv=data.ScanParam.PWV)
    print 'ATM model opacity: tau_z = %3.3f'%(tau)
    chan0 = data.BolometerArray.checkChanList ([])[0]
    el = data.getChanData('el',chan0)
    scanel = np.median(el)
    taucorr = np.exp(tau/np.sin(scanel * np.pi / 180.))

    ### GET PLANET FLUX FROM ASTRO
    source = data.ScanParam.Object
    if source in ['Uranus','Neptune','Mars','Saturn','Jupiter','Venus']:
        astrotime,astrodate=getAstroDate(data)
        if (febe == 'LFA') or (febe=='lfa'):
            beam = 19.2
            freq = 345.
        if (febe == 'HFA') or (febe=='hfa'):
            beam = 7.8
            freq = 810.	
    
        print source ,astrotime,astrodate,str(beam),str(freq)
        expect_flux = PlanetFlux(source ,astrotime,astrodate,beam,freq)
    else:
        if calibFluxes.has_key(string.upper(source)):
            expect_flux = calibFluxes[string.upper(source)]  


    if expect_flux > 0:
        print "######################################"
        print "Expected flux for %s: %5.1f Jy"%(source,expect_flux)
        print "######################################"
    else:
        print "######################################"
        print "Expected flux of source undetermined - terminating reduction"
        print "######################################"
        return

    
    
    ############################
    
    ## PROCESS DATA
    data.zeroStart()
                    
    
    #data.flagSpeed (below=30.)
    #data.flagSpeed(above=500.)
    #data.flagAccel(above=800.)
    if source in ['Venus','Jupiter','jupiter','venus','saturn','Saturn','Mars','mars']:        
        data.flagPosition(radius=60,flag=8)
    else:
        data.flagPosition(radius=30,flag=8)
    base(order=0)
    #data.flagFractionRms(ratio=5)
    correlChains(chanRef=-2,factor=0.95,nbloop=3)
    if  febe == 'LFA'  or  febe=='lfa' or febe == 'AMKID870-AMKID870BE': 
        chainIndices = range(4)
    if  febe == 'HFA'  or  febe=='hfa' or febe == 'AMKID350-AMKID350BE':
        chainIndices = range(20)
    #for chainIndex in chainIndices:
    #    correlGroups(chainIndex, chanRef=-2, factor=0.97, nbloop=2,numRepeats=2, groupingThreshold=0.7)
    data.flattenFreq(below=0.5,hiref=0.7)
    #blankFreq(below=0.88,above=0.85)
    #blankFreq(below=5.30,above=5.26)
    #blankFreq(below=9.97,above=9.93)
    #blankFreq(below=9.77,above=9.74)
    #blankFreq(below=9.23,above=9.21)
    #blankFreq(below=8.75,above=8.72)
    #blankFreq(below=8.21,above=8.18)
    #blankFreq(below=7.19,above=7.16)
    #blankFreq(below=3.08,above=3.06)
    #blankFreq(below=2.06,above=2.04)
    
    base(order=1)
    weight()
    unflag(flag=8)
    
    ##### GET CAL CORR FOR EACH CHANNEL BASED ON FIT ON MAP
    mydata=copy.deepcopy(data)
    bolos = mydata.BolometerArray.checkChanList([])
    ChanIndex = mydata.BolometerArray.getChanIndex(bolos)
    for c in ChanIndex:
        if mydata.BolometerArray.Offsets[0][c] < -8000:
            mydata.BolometerArray.Offsets[0][c] = 0.0
        if mydata.BolometerArray.Offsets[1][c] < -8000:
            mydata.BolometerArray.Offsets[1][c] = 0.0
    if source in ['Venus','Jupiter','jupiter','venus','saturn','Saturn','Mars','mars']:
        radius_large=1250
    else:
        radius_large=100
    fit_beam_per_channel(mydata,bolos,radius_large=radius_large,noMap=False,max_peak=1e9)
    fitpeaktemp=[]
    fitchan=[]
    for i in range(len(mydata.arrayParamOffsets)):
        fitpeaktemp.append(mydata.arrayParamOffsets[i]['result']['gauss_peak']['value'])
        fitchan.append(mydata.arrayParamOffsets[i]['channel'])
            
    JyperK = np.array((expect_flux/taucorr)/np.array(fitpeaktemp))

    outrcp =  str(scan)+'_'+str.upper(febe)+'_'+str(source)+'.rcp'
    mydata.updateArrayParameters(outrcp)
    command = 'cp '+str(outrcp)+' '+str(BoaConfig.rcpPath)
    os.system(command)

    updateRCP(outrcp)
    flagRCP(outrcp)


    ## create rcp plot- to be put into a seperate fkt
    chanList= mydata.BolometerArray.checkChanList([])
    array_filename = '%i_%s_%s_array.ps'%(scan,febe,source)
    op('%i_%s_%s_array.ps/CPS'%(scan,febe,source))
    DeviceHandler.DevList
    plotRCP(bolos,ScanNr,num=0,data=mydata)
    Plot.nextpage()
    plotArray(num=0,ci=1,limitsX=[-650,650],limitsY=[-650,650])
    Plot.nextpage()
    
    
    if doDespike == 1:
        if source in ['Venus','Jupiter','jupiter','venus','saturn','Saturn','Mars','mars']:        
            data.flagPosition(radius=60,flag=8)
        despike(below=-5,above=5)
        base(order=1)
        weight()
        unflag(flag=8)

    mapping(sizeX=[-200,200],sizeY=[-200,200],oversamp=10,noPlot=0,aspect=1)
    Plot.nextpage()

    data.solvePointingOnMap(radius=60.0,plot=0)
    result = data.PointingResult
    norm = result['gauss_peak']['value']
    plYmax = 0.2*norm
    plYmin = -0.05*norm
    mapping(sizeX=[-600,600],sizeY=[-600,600],oversamp=10,aspect=1,limitsZ=[plYmin,plYmax])
    close()


    ## APPLY FITTED JY/K TO DATA (= to flatfielding on observed Mars Temperature)

    flagC(range(1,4000) )
    data.unflagChannels(fitchan)
    chanIndices = data.BolometerArray.getChanIndex(fitchan)
    data.Data[:,chanIndices] *= JyperK
    
    mapping(sizeX=[-200,200],sizeY=[-200,200],oversamp=7,noPlot=0,aspect=1)
    
    print "Solving for pointing on map..."
    try:
         data.solvePointingOnMap(plot=0,radius=30)
         obsfluxmap = data.PointingResult['gauss_peak']['value']
    except:
        return

    percent = 100.0*obsfluxmap/(expect_flux/taucorr)
    if 95.0 < percent < 105.0: 
        data.Data[:,chanIndices] /= np.array(percent/100.,'f')
        data.solvePointingOnMap(plot=0,radius=30)
        obsfluxmap = data.PointingResult['gauss_peak']['value']
        percent = 100.0*obsfluxmap/(expect_flux/taucorr)
        JyperK = JyperK/(percent/100.)
        data.Data[:,chanIndices] /= percent/100.


    print "-------------------------------------------------------------"
    print "FLUX %s:  %7.2f [expected: %7.2f, %6.2f percent]"%(source,obsfluxmap,expect_flux/taucorr,percent)
    print "-------------------------------------------------------------"

    ## COMPUTE SLIDING RMS

    nbInteg = 500
    chanList = data.BolometerArray.checkChanList([])
    noise = []
    for c in chanList:
        index = data.BolometerArray.getChanIndex(c)[0]
        rms = fStat.slidingrms(data.Data[:,index],
                               data.FlagHandler._aFlags[:,index],
                               nbInteg)
        noise.append(np.min(rms))
       
    calibrated_noise = np.array(noise)

    data.flagPosition(radius=60,flag=8)
    base(order=1)
    weight()
    calibrated_noisefull = np.array((data.getChanListData('rms')),'f')
    unflag(flag=8)

    ###  COMPUTE SENSITIVIES
    
    dt = np.median((data.ScanParam.get('deltat')))
    sensi = np.array(calibrated_noise*np.sqrt(dt)*1000.0)
    sensifull = np.array(calibrated_noisefull*np.sqrt(dt)*1000.0)

    array_sensitivity = sum(sensi/(calibrated_noise**2))/sum(1./(calibrated_noise**2))
    s20,s80,speak,ppeak,sdist,pdist=getSensitivityBreakdown(sensi)

    ## DO REST 

    if (febe == 'LFA') or (febe=='lfa'):
        nrchains=4
        refnr=880
        chains=LFAChains
    if (febe == 'HFA') or (febe=='hfa'):
        nrchains=20
        refnr=1080
        chains=HFAChains

    unit = 'mJy sqrt(s)'
    chaincount =0
    WeightedChanSen=[]
    for chain in chains:
           inChain=list(set(chain) & set(fitchan))
           chaincount = chaincount+1
           if len(inChain) > 0:
               chainsen = []
               chainrms = []
               for num in inChain:
                   index=list(fitchan).index(num)
                   chainsen.append(list(sensi)[index])
                   chainrms.append(list(calibrated_noise)[index])
               meanchainsen = sum(np.array(chainsen)/(np.array(chainrms)**2))/sum(1./(np.array(chainrms)**2))
               WeightedChanSen.append(meanchainsen)
               print "Sensitivity of Chain %i: %5.1f %s"%(chaincount,meanchainsen,unit)
           else:
               meanchainsen = -1
               WeightedChanSen.append(meanchainsen)
               print "No KIDS in Chain %i:found in this scan"%(chaincount)
               print "####################################"

    chain_sensitivity = WeightedChanSen
    scan = data.ScanParam.ScanNum
    if calfile:
        outfile = str(scan)+'_'+str(source)+'_sensitivity-Tcal-summary.txt'
    else:
        outfile = str(scan)+'_'+str(source)+'_sensitivity-summary.txt'
    output=file(outfile,'w')
    output.write("\n")
    output.write("#######################        Array Sensitivity        ################################\n")
    output.write("              noise weighted mean | mean best 20% | mean best 80% \n")
    output.write("                [mJy sqrt(s)]        [mJy sqrt(s)]    [mJy sqrt(s)] \n\n")
    output.write(" Full array:       %5.2f              %5.2f            %5.2f     (%4i KIDS)\n\n"%(array_sensitivity,s20,s80,len(sensi)))


    for i in range(nrchains):
        chainsen=np.array(np.extract( (np.array(fitchan)>i*refnr) * (np.array(fitchan)<=(i+1)*refnr) ,sensi))
        if len(chainsen) > 0:
            s20,s80,speak,ppeak,sdist,pdist=getSensitivityBreakdown(chainsen)      
            output.write("   Chain %2i:       %5.2f              %5.2f            %5.2f     (%4i KIDS)\n"%(i+1,chain_sensitivity[i],s20,s80,len(chainsen)))
        else:
            output.write("   Chain %2i:                                                      no Kids  \n"%(i+1))
    output.close()
    print "Sensitivity summary writen to file %s"%(outfile)
    if calfile:
       outfile = str(scan)+'_'+str(source)+'_sensitivity-Tcal.dat'
    else:
        outfile = str(scan)+'_'+str(source)+'_sensitivity.dat'

    output=file(outfile,'w')
    output.write("#######################        Pixel Sensitivities        ################################\n")
    if calfile:
        output.write("## KID    coupling     sensitivity (best)  sensitivity (full)  rms (best)  rms (full) \n")
        output.write("##         [Jy/K]          [mJy sqrt(s)]     [mJy sqrt(s)]       [Jy]      [Jy] \n")
    else:
        output.write("## KID    coupling     sensitivity (best)  sensitivity (full)  rms (best)  rms (full) \n")
        output.write("##        [Jy/rad]          [mJy sqrt(s)]     [mJy sqrt(s)]       [Jy]      [Jy] \n")
    
    for i in range(len(fitchan)):
        output.write("%5i       %5.1f              %5.2f          %5.2f           %5.2f      %5.2f\n"%(fitchan[i],JyperK[i],sensi[i],sensifull[i],calibrated_noise[i],calibrated_noisefull[i]))
    output.close()



    #SENSI PLOTS 
    psfile = str(scan)+'_'+str(source)+'_sensitivity-all.ps/CPS'
    op(psfile)
    plot(np.array(fitchan),np.array(sensi),labelX='Kid Nr',labelY='Sensitivity [mJy sqrt(s)]',limitsY=[70,600])
    close()
  
    psfile = str(scan)+'_'+str(source)+'_sensitivity-ch1.ps/CPS'
    op(psfile)
    plot(np.array(fitchan),np.array(sensi),labelX='Kid Nr',labelY='Sensitivity [mJy sqrt(s)]',limitsX=[0,880],limitsY=[70,600])
    close()

    psfile = str(scan)+'_'+str(source)+'_sensitivity-ch2.ps/CPS'
    op(psfile)
    plot(np.array(fitchan),np.array(sensi),labelX='Kid Nr',labelY='Sensitivity [mJy sqrt(s)]',limitsX=[881,1760],limitsY=[70,600])
    close()


    psfile = str(scan)+'_'+str(source)+'_sensitivity-ch3.ps/CPS'
    op(psfile)
    plot(np.array(fitchan),np.array(sensi),labelX='Kid Nr',labelY='Sensitivity [mJy sqrt(s)]',limitsX=[1761,2640],limitsY=[70,600])
    close()

    psfile = str(scan)+'_'+str(source)+'_sensitivity-ch4.ps/CPS'
    op(psfile)
    plot(np.array(fitchan),np.array(sensi),labelX='Kid Nr',labelY='Sensitivity [mJy sqrt(s)]',limitsX=[2641,3520],limitsY=[70,600])
    close()


    #else:
    #    print "No scan defined."


def test(ScanNr=0,calfile=0,febe='LFA'):
    if (ScanNr != 0):
        if (febe == 'LFA') or (febe=='lfa'):
            tst=read(str(ScanNr),febe='AMKID870-AMKID870BE')
            updateScanParam(data)
            flagEmpty()
        if (febe == 'HFA') or (febe=='hfa'):
            tst=read(str(ScanNr),febe='AMKID350-AMKID350BE')
            updateScanParam(data)
            flagEmpty()
        
        
        if tst:
            print "Scan %s not readable"%(str(ScanNr))
            return
        

        scan = data.ScanParam.ScanNum
        mjdref = (data.ScanParam.MJD[-1]+data.ScanParam.MJD[0])/2.  
        rcp=getMKIDsRCP(mjdref)        
        updateRCP(rcp)

    if calfile:
            if str(calfile)[len(str(calfile))-4:len(str(calfile))] != '.dat':
                if str.upper(febe) == 'LFA':
                    febename = 'AMKID870-AMKID870BE'
                else:
                    febename = 'AMKID350-AMKID350BE'
                calnr = np.int(calfile)
                search = '*%i_%s_KidCalspline.dat' % (calnr,febename)
                #print search
                if len(glob(search)) > 0:
                    calfile = glob(search)[-1]
                else:
                    print "############################################################################"
                    print "PolGrid calibration file %s not found - reduce PolGrid calibration first"%(calfile)
                    print "############################################################################"
                    return

            Tamb = data.ScanParam.T_amb
            applyTempCal(calfile)
            #applySkyFlat()
            dataTsky= data._DataAna__computeMedianSignal()
            obsTsky = fStat.f_median(np.array(dataTsky))
            if (febe == 'LFA') or (febe=='lfa'):
                obsTsky = obsTsky-(1-Feff_LFA)*Tamb
            if (febe == 'HFA') or (febe=='hfa'):
                obsTsky = obsTsky-(1-Feff_HFA)*Tamb
    else:
        data.Data *= np.array(-1,'f')
    print calfile


def redmap(ScanNr=0,oversamp=3,system='ho'):
    '''
    Reduces a LABOCA map on a target source.
    Optional parameters: Scan number, zenith opacity
                         ra1, ra2, dec1, dec2 = abs. map limits in degree
    e.g. redmap(1234,tau=0.3,ra1=190.4,ra2=190.2,dec1=-63.1,dec2=-62.9)
         redmap(1234)
    '''
    if (ScanNr != 0):
        tst = read(str(ScanNr))
        if tst:
            print "Scan %s not readable"%(str(ScanNr))
            return
        tau = scanTau()
        print 'No tau value passed, using ATM model: tau_z = %3.3f'%(tau)
        data.correctOpacity(tau)
        
    else:
        print "No scan defined."


def readscan(infile='dummy'):
    f=file(infile)
    t=cPickle.load(f)
    f.close()
    return t

    
def redpoints(ScanList=[]):
    '''
    Reduces multiple LABOCA pointings on a target source 
    Optional parameters: Scan numbers (list), zenith opacity single value or list
    '''
    if (len(ScanList) != 0):
        mapList=[]
	w1=[]
	w2=[]
	pk=[]
        for i in range(len(ScanList)):
	
	    scan=ScanList[i]
	    redcal(scan)
	    mapping(oversamp=5,limitsZ=[-1,5],system='HO',sizeX=[-75,75],sizeY=[-75,75])
	    solvePointingOnMap(plot=1)
	    result = data.PointingResult
	    w1.append(result['gauss_x_fwhm']['value'])
	    w2.append(result['gauss_y_fwhm']['value'])
	    pk.append(result['gauss_peak']['value'])
        
	print w1 
	print w2
	print pk
		       
    else:
        print "No scans defined."


def calib(ScanNr=0,tau=0.0):
    if (ScanNr != 0):
        tst = read(str(ScanNr))
        if tst:
            print "Scan %s not readable"%(str(ScanNr))
            return
        if (tau == 0.0):
	    tau = scanTau()
            print 'No opacity value passed, using ATM model: tau_z = %3.3f'%(tau)

        CntstoV(data)
        data.Data *= np.array(VtoJy,'f')

        nr=len(data.ScanParam.MJD)
        nr = nr -1
        mjdref = (data.ScanParam.MJD[nr]-data.ScanParam.MJD[0])/2.+data.ScanParam.MJD[0]

        rcp=getLabocaRCP(mjdref)
        updateRCP(rcp)


        data.zeroStart()
        correctBlind(data)

        flagRCP(rcp)
        cross=getLabocaCross(mjdref)
        flagC(cross)

        flat()
        data.correctOpacity(tau)
        sig()





def redweakfkt(ScanNr=0,taufile='',calfile=''):
    '''
    Reduces a LABOCA map on a pointing source.
    Optional parameters: Scan number, zenith opacity
    '''
    if (ScanNr != 0):
        tst = read(str(ScanNr))
        if tst:
            print "Scan %s not readable"%(str(ScanNr))
            return -1,-1,-1,-1
   
        CntstoV(data)
        data.Data *= np.array(VtoJy,'f')

        nr=len(data.ScanParam.MJD)
        nr = nr -1
        mjdref = (data.ScanParam.MJD[nr]-data.ScanParam.MJD[0])/2.+data.ScanParam.MJD[0]

        rcp=getLabocaRCP(mjdref)
        updateRCP(rcp)
        print rcp
        
        calcorr = 0.0
        if calfile:
           calcorr = getCalCorr(mjdref,'linear',calfile)
        else:
           calcorr = 1.0
        if calcorr > 0.0:
            data.Data /= np.array((calcorr),'f')
        else:
	    print "No calbration correction found - assuming 1.0"
            calcorr = 1.0	
            

        data.zeroStart()
        correctBlind(data)

        flagRCP(rcp)
        cross=getLabocaCross(mjdref)
        print cross
        flagC(cross)

        flat()

        data.flagSpeed (below=30.)
        data.flagSpeed(above=500.)
        data.flagAccel(above=800.)

        data.flagPosition(radius=30,flag=8)

        data.flagFractionRms(ratio=5)
        despike_nr = 0
        medianNoiseRemoval(chanRef=-1,factor=0.8,nbloop=5)
        toFlag=jumps()
        jumps_nr = len(toFlag)
        data.flagFractionRms(ratio=3)
        medianNoiseRemoval(chanRef=-1,factor=0.8,nbloop=5)
        spikenr = despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr
        correlbox(data,factor=0.8,nbloop=5)
        data.flagFractionRms(ratio=4)
        spikenr =despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr
        correlgroup(data,factor=0.8,nbloop=5)
        spikenr =despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr
        base(order=1)
        data.flagFractionRms(ratio=4)
        #data.flattenFreq(below=0.5,hiref=0.7)
	data.blankFreq(below=0.5)
        spikenr =despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr
        data.slidingWeight(nbInteg=50)
        unflag(flag=8)
        spikenr =despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr

        sensi = calcsensitivity()
        print '#######################################'
        print "Sensitivity: %6.1f mJy sqrt(s)"%(sensi)
        print '#######################################'

        tau=0.0
        if taufile:
           tau = getTau(mjdref,'linear',taufile)
        else:
           tau = scanTau()
        if tau > 0.0:
            data.correctOpacity(tau)
            data._DataAna__statistics()
	    data.computeWeight()
            scanel  = fStat.f_mean(data.ScanParam.El)
            taucorr = exp(tau/sin(scanel * pi / 180.))
            return sensi,taucorr,calcorr,jumps_nr,despike_nr,tau
        else:
            print "No tau correction found!"
            return -1,-1,-1,-1,-1,-1
            

    else:
        print "No scan defined."
        return -1,-1,-1,-1,-1,-1


def redbrightfkt(ScanNr=0,taufile='',calfile=''):
    '''
    Reduces a LABOCA map on a pointing source.
    Optional parameters: Scan number, zenith opacity
    '''
    if (ScanNr != 0):
        tst = read(str(ScanNr))
        if tst:
            print "Scan %s not readable"%(str(ScanNr))
            return -1,-1,-1,-1
   
        CntstoV(data)
        data.Data *= np.array(VtoJy,'f')

        nr=len(data.ScanParam.MJD)
        nr = nr -1
        mjdref = (data.ScanParam.MJD[nr]-data.ScanParam.MJD[0])/2.+data.ScanParam.MJD[0]

        rcp=getLabocaRCP(mjdref)
        updateRCP(rcp)
        print rcp
        
        calcorr = 0.0
        if calfile:
           calcorr = getCalCorr(mjdref,'linear',calfile)
        else:
           calcorr = 1.0
        if calcorr > 0.0:
            data.Data /= np.array((calcorr),'f')
        else:
	    print "No calbration correction found - assuming 1.0"
            calcorr = 1.0	
            

        data.zeroStart()
        correctBlind(data)

        flagRCP(rcp)
        cross=getLabocaCross(mjdref)
        print cross
        flagC(cross)

        flat()

        data.flagSpeed (below=30.)
        data.flagSpeed(above=500.)
        data.flagAccel(above=800.)

        data.flagPosition(radius=30,flag=8)
        data.flagFractionRms(ratio=5)

        despike_nr = 0

        spikenr =despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr
        medianNoiseRemoval(chanRef=-1,factor=0.8,nbloop=5)
        toFlag=jumps()
        jumps_nr = len(toFlag)
        data.flagFractionRms(ratio=5)
        correlbox(data,factor=0.8,nbloop=5)
        correlgroup(data,factor=0.8,nbloop=5)
        data.computeWeight()
        unflag(flag=8)

        sensi = calcsensitivity()
        print '#######################################'
        print "Sensitivity: %6.1f mJy sqrt(s)"%(sensi)
        print '#######################################'

        tau=0.0
        if taufile:
           tau = getTau(mjdref,'linear',taufile)
        else:
           tau = scanTau()
        if tau > 0.0:
            data.correctOpacity(tau)
            data._DataAna__statistics()
	    data.computeWeight()
            scanel  = fStat.f_mean(data.ScanParam.El)
            taucorr = exp(tau/sin(scanel * pi / 180.))
            return sensi,taucorr,calcorr,jumps_nr,despike_nr,tau
        else:
            print "No tau correction found!"
            return -1,-1,-1,-1,-1,-1
            

    else:
        print "No scan defined."
        return -1,-1,-1,-1,-1,-1



def redweakmodelfkt(ScanNr=0,taufile='',calfile='',modelfile='',sub=0,cliplevel=3.0):
    '''
    Reduces a LABOCA map on a pointing source.
    Optional parameters: Scan number, zenith opacity
    '''
    if (ScanNr != 0):
        tst = read(str(ScanNr))
        if tst:
            print "Scan %s not readable"%(str(ScanNr))
            return -1,-1,-1,-1,-1,-1
   
        CntstoV(data)
        data.Data *= np.array(VtoJy,'f')

        nr=len(data.ScanParam.MJD)
        nr = nr -1
        mjdref = (data.ScanParam.MJD[nr]-data.ScanParam.MJD[0])/2.+data.ScanParam.MJD[0]

        rcp=getLabocaRCP(mjdref)
        updateRCP(rcp)
        print rcp
        
        calcorr = 0.0
        if calfile:
           calcorr = getCalCorr(mjdref,'linear',calfile)
        else:
           calcorr = 1.0
        if calcorr > 0.0:
            data.Data /= np.array((calcorr),'f')
        else:
            print "No calbration correction found"
            return -1,-1,-1,-1,-1,-1

        data.zeroStart()
        correctBlind(data)

        flagRCP(rcp)
        cross=getLabocaCross(mjdref)
        print cross
        flagC(cross)

        flat()

        data.flagSpeed (below=30.)
        data.flagSpeed(above=500.)
        data.flagAccel(above=800.)
	
	jumps_nr = 0
	despike_nr = 0

        model=restoreFile(modelfile)
        modeln=copy.deepcopy(model)
        model.Data*=sqrt(model.Weight)
        model.smoothBy(19.0/3600.0)
        model.computeRms()
        scale=model.RmsBeam
        model.Data /= np.array(scale,'f')
        if sub ==0:
           data.flagSource(threshold=cliplevel,model=model)
        else:
           nbX,nbY = shape(modeln.Data)
           for i in range(nbX):
              for j in range(nbY):
                 if is_nan(modeln.Data[i,j]) or model.Data[i,j]<1.5:
                        modeln.Data[i,j] = 0.
           tau = 0.0
           if taufile:
              tau = getTau(mjdref,'linear',taufile)
           else:
              tau = scanTau()
           if tau > 0.0:
              scanel  = fStat.f_mean(data.ScanParam.El)
              taucorr = exp(tau/sin(scanel * pi / 180.))
              modeln.Data = -1*modeln.Data/taucorr
              modelp=copy.deepcopy(modeln)
              modelp.Data = -1*modelp.Data
              data.addSource(model=modeln)

        data.flagFractionRms(ratio=5)
        despike_nr = 0
        medianNoiseRemoval(chanRef=-1,factor=0.8,nbloop=5)
        toFlag=jumps()
        jumps_nr = len(toFlag)
        data.flagFractionRms(ratio=3)
        medianNoiseRemoval(chanRef=-1,factor=0.8,nbloop=5)
        spikenr = despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr
        correlbox(data,factor=0.8,nbloop=5)
        data.flagFractionRms(ratio=4)
        spikenr =despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr
        correlgroup(data,factor=0.8,nbloop=5)
        spikenr =despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr
        base(order=1)
        data.flagFractionRms(ratio=4)
        #data.flattenFreq(below=0.5,hiref=0.7)
	data.blankFreq(below=0.5)
        spikenr =despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr
        data.slidingWeight(nbInteg=50)
        if sub == 0:
           unflag(flag=8)
        else:
           data.addSource(model=modelp)
        spikenr =despike(below=-5,above=5)
        despike_nr = despike_nr+spikenr
        sensi = calcsensitivity()
        print '#######################################'
        print "Sensitivity: %6.1f mJy sqrt(s)"%(sensi)
        print '#######################################'

        tau=0.0
        if taufile:
           tau = getTau(mjdref,'linear',taufile)
        else:
           tau = scanTau()

        if tau > 0.0:
            data.correctOpacity(tau)
            data._DataAna__statistics()
	    data.computeWeight()
            scanel  = fStat.f_mean(data.ScanParam.El)
            taucorr = exp(tau/sin(scanel * pi / 180.))
            return sensi,taucorr,calcorr,jumps_nr,despike_nr,tau
        else:
            print "No tau correction found!"
            return -1,-1,-1,-1,-1,-1

            

    else:
        print "No scan defined."
        return -1,-1,-1,-1,-1,-1


def rednoise(ScanNr=0,febe='LFA',calfile='',doPlot=0):
    if ScanNr > 0:
        if (febe == 'LFA') or (febe=='lfa'):
            tst=read(str(ScanNr),febe='AMKID870-AMKID870BE')
            flagEmpty()
        if (febe == 'HFA') or (febe=='hfa'):
            tst=read(str(ScanNr),febe='AMKID350-AMKID350BE')
            flagEmpty()

    if calfile:
        applyTempCal(calfile)

    data.zeroStart()
    flagFractionRms()
    correlChains(chanRef=-2,factor=0.95,nbloop=3)
    if  febe == 'LFA'  or  febe=='lfa' or febe == 'AMKID870-AMKID870BE': 
        chainIndices = range(4)
    if  febe == 'HFA'  or  febe=='hfa' or febe == 'AMKID350-AMKID350BE':
        chainIndices = range(20)
    for chainIndex in chainIndices:
	 correlGroups(chainIndex, chanRef=-2, factor=0.97, nbloop=3,numRepeats=3, groupingThreshold=0.7)
    base(order=1)
    #despike()
    if doPlot == 1:
        showSig()


def redOTFbeam(scanlist=[],calscan=0,febe='LFA',factor=0.75,doPlot=1,oversamp=5,sizeX=[-300,300],sizeY=[-300,300],radius=180.):
    centerind = -1
    ell=[]
    ellaver=[]
    arrayParamOffsets=[]

    if len(scanlist) == 1:
        scan = str(scanlist[0])
    else:
        scan=str(scanlist[0])+'-'+str(scanlist[-1])
   
    
    for i in range(len(scanlist)):
        good = []
        AP=[]
        if i == 4:
            centerind = len(ell)-1
        if str.upper(febe) == 'LFA':
            if calscan > 0:
                redcal(scanlist[i],calscan,doPlot=0,showSen=0)
            else:
                _,_,_,_,_,_=redpnt(scanlist[i],doPlot=0)    
        else:
            get(scanlist[i], fe='HFA', cal=1)
            data.zeroStart()
            data.flattenFreq()
        

        if calscan > 0:
            source_name = data.ScanParam.Object
            if i == 0:
                logfile =  str(scan)+'_coupling_'+str(source_name)+'.log' 
                outlog=open(logfile,'w')

            if source_name in ['Uranus','Neptune','Mars','Saturn','Jupiter','Venus']:
                astrotime,astrodate=getAstroDate(data)
                beam = 19.2
                freq = 345
                expect_flux = PlanetFlux(source_name,astrotime,astrodate,beam,freq)
            else:
                if calibFluxes.has_key(string.upper(source_name)):
                    expect_flux = calibFluxes[string.upper(source_name)]
                else:
                    expect_flux = 0
        else:
            expect_flux = 0

        xoffscan=data.ScanParam.get('AzimuthOffset')
        yoffscan=data.ScanParam.get('ElevationOffset')
        xoffrcp,yoffrcp=data.BolometerArray.Offsets
        refX=data.BolometerArray.RefOffX
        refY=data.BolometerArray.RefOffY
        rcpchans=data.BolometerArray.UsedChannels

        xoffrcp = xoffrcp #-refX
        xoffrcp = xoffrcp #-refY


        xoffscan= xoffscan+refX
        yoffscan= yoffscan+refY

        cchans=[]
        for k in rcpchans:
            mindist=np.min(np.sqrt((xoffscan-xoffrcp[k-1])**2+(yoffscan-yoffrcp[k-1])**2))
            if mindist < 30.0:
                cchans.append(k)
        centerchans = []
        usedbolos = data.BolometerArray.checkChanList([])
        for c in cchans:
            if c in usedbolos:
                centerchans.append(c)

        
        flux=np.array(data.getChanListData('flux',centerchans))
        maxflux=np.max(flux,axis=1)

        minflux = np.median(maxflux) #/1.5

        for c in centerchans:
            #mapping(c,oversamp=oversamp,sizeX=sizeX,sizeY=sizeY,noPlot=1)
            mapping(c,oversamp=oversamp,noPlot=1)
            try:
                solvePointingOnMap(plot=doPlot,radius=radius)
                result = data.PointingResult
                result['gain'] = 1.0
                fwhmx=result['gauss_x_fwhm']['value']
                fwhmy=result['gauss_y_fwhm']['value']
                ang=result['gauss_tilt']['value']
                peak=result['gauss_peak']['value']
                aratio = np.max([fwhmx,fwhmy])/ np.min([fwhmx,fwhmy])
                if fwhmx > 15. and fwhmy > 15. and peak > minflux and aratio < 3.0:
                    ell.append([c,fwhmx,fwhmy,ang])
                    arrayParamOffsets.append({'channel' : c,\
                                       'result'  : result})
                    AP.append({'channel' : c,\
                                       'result'  : result})
                    good.append(c)
                    if doPlot == 1:
                        raw_input()
            except:
                print ''

        
        outrcp = 'tmp.rcp'
        data.arrayParamOffsets = copy.deepcopy(AP)
        data.updateArrayParameters(outrcp)  
        command = 'cp '+str(outrcp)+' '+str(BoaConfig.rcpPath)
        os.system(command)
        updateRCP(outrcp)

        data.flagPosition(radius=30.,flag=8)
        base(order=1)
        data._DataAna__statistics() 
        weight()
        unflag(flag=8)

        ch=getKidsperChain(good)
        print len(good)
        print len(ch[0]), len(ch[1]), len(ch[2]), len(ch[3])
        for j in range(4):
            if len(ch[j]) > 4 and len(good) > 0:
                mapping(ch[j],oversamp=8,sizeX=sizeX,sizeY=sizeY,noPlot=0)
                chnr=j+1
                if i < 4:
                    name=str(scanlist[i])+'_CHAIN'+str(chnr)+'_center.data'
                else:
                    name=str(scanlist[i])+'_CHAIN'+str(chnr)+'_inner.data'
                data.Map.dumpMap(name)
                try:
                    solvePointingOnMap(plot=doPlot,radius=radius)
                    result = data.PointingResult
                    result['gain'] = 1.0
                    fwhmx=result['gauss_x_fwhm']['value']
                    fwhmy=result['gauss_y_fwhm']['value']
                    ang=result['gauss_tilt']['value']
                    ellaver.append([ch[j],fwhmx,fwhmy,ang])
                    if expect_flux > 0:
                        obsflux = result['gauss_peak']['value']
                        percent = 100.0*obsflux/expect_flux
                        outlog.write('Chain%i: observed: %6.1f Jy expected %6.1f Jy; %4.1f percent (coupling: %5.1f Jy/K)\n' %(j+1,obsflux,expect_flux,percent,40*expect_flux/obsflux))
                except:
                    print 'fit failed - not in ellaver'
                #print ellaver[-1][0]
                #print len(ellaver[-1][0])
                if doPlot == 1:
                    raw_input()



    ####
    data.arrayParamOffsets = copy.deepcopy(arrayParamOffsets)
    source = data.ScanParam.Object
    
    if data.BolometerArray.FeBe == 'AMKID870-AMKID870BE':
        febe = 'LFA'
    else:
        febe = 'HFA'
    
    if len(scanlist) == 1:
        scan = str(scanlist[0])
    else:
        scan=str(scanlist[0])+'-'+str(scanlist[-1])
    outrcp =  str(scan)+'_'+str(febe)+'_'+str(source)+'.rcp'    
    data.updateArrayParameters(outrcp)  
    command = 'cp '+str(outrcp)+' '+str(BoaConfig.rcpPath)
    os.system(command)

    psname=outrcp.replace('.rcp','')+'_array.ps/CPS'
    op(psname)

    cap = "MKIDS %s fitted beam - %s %s"%(data.BolometerArray.FeBe,data.ScanParam.Object,data.ScanParam.DateObs)

    sX=[-600,600]
    sY=[-600,600]


    pgsci(2) 
    BogliConfig.xyouttext['color'] = 1
    BogliConfig.xyouttext['charheight'] = 0.5
    kids, _, _, offX, offY = openRCP(outrcp)


    ### PLOT INDIVIDUAL BEAMS
    Plot.plot([0],[0],limitsX=sX,limitsY=sY,nodata=1,aspect=1,labelX='Az offset ["]',labelY='El offset ["]',caption=cap)
    

    refX = 0.0
    refY = 0.0

    for i in range(len(ell)):
        if ell[i][0] in kids:
            ind=list(kids).index(ell[i][0])
            x=offX[ind]-refX
            y=offY[ind]-refY
            wx=ell[i][1]*factor
            wy=ell[i][2]*factor
            tilt=ell[i][3]
            if centerind > -1 and i > centerind:
                Forms.ellipse(x,y,wx,wy, tilt*np.pi/180.,overplot=1,ci=2)
            else:
                Forms.ellipse(x,y,wx,wy, tilt*np.pi/180.,overplot=1,ci=1)
    
    print 'Fit result for beam shape writen to file %s'%(psname)
    
    ### PLOT AVERAGE BEAMS
    factor = 2.0
    Plot.nextpage()
    plot([0],[0],limitsX=sizeX,limitsY=sizeY,nodata=1,aspect=1,\
         labelX='Az offset ["]',labelY='El offset ["]',caption=cap)

    for i in range(len(ellaver)):
        myX=[]
        myY=[]
        KidsPos = list(set(ellaver[i][0]).intersection(kids))
        for j in range(len(KidsPos)):
            ind= list(kids).index(KidsPos[j])
            myX.append(offX[ind])
            myY.append(offY[ind])

        meanX=np.mean(np.array(myX))
        meanY=np.mean(np.array(myY))
        
        wx=ellaver[i][1]*factor
        wy=ellaver[i][2]*factor
        tilt=ellaver[i][3]
        if ellaver[i][1] < 100.:
            cutstrx= 4
        else:
            cutstrx= 5

        if ellaver[i][2] < 100.:
            cutstry= 4
        else:
            cutstry= 5

        text=str(str(ellaver[i][1])[0:cutstrx])+"x"+str(str(ellaver[i][2])[0:cutstry])
        Forms.ellipse(meanX,meanY,wx,wy, tilt*np.pi/180.,overplot=1)
        Plot.xyout(meanX, meanY+45., text)
    close()
    #return ellaver

    if calscan > 0:
        outlog.close()

def redscans(ScanList=[[]],calscan='',smooth=12.0,calfile='',taufile='',modelfile='',oldmap='',weak=1,covclip=0.15,writeScans=1,sub=0,cliplevel=3.0,overwrite=0,oversamp=3,febe='LFA',aspect=0):
	'''
	Reduces multiple LABOCA map on a target source and coadds the data.
	tau 1 --> [scan1,scan1N] first group
	tau 2 --> [scan2,scan2N] second group 
	.
	.
	.
	tau N --> [scanN,scanNN] N group
	
	The format to execute is :
	
	redSscans(ScanList=[ [scan1,scan1N],[scan2,scan2N] ],tau=[tau 1,tau 2],smooth=18.0, addModel=0,covclip=0.25)
	
	addModel: Remove data model from the data, it can use after the first primary reduction 
	
	Map size is determined by the first scan in the list.
	
	'''
	
	ms,count = 0,0
	SCAN_LIST = []
	CAL_LIST = []
	tt=0

	# Checking the scanReduced directory was already created 
	try:
		os.system('mkdir scanReduced')
	except:
		pass

        for s in (ScanList):
		 if len(s)>1:
                   for k in (s):
                       SCAN_LIST.append(k)
		 else:
		   SCAN_LIST.append(s[0])
		 tt+=1

	### GET SOURCE NAME FROM FIRST SCAN ###

        test=read(SCAN_LIST[0],febe='AMKID870-AMKID870BE')
	if test != -1:
		source=data.ScanParam.Object
                raref=data.ScanParam.RA0
                decref=data.ScanParam.Dec0
	else:
		print "Can not read first scan - quit and check Scan numbers"
		raw_input()
        
	#### Determine Map size in RA DEC #########

	# If we have an previously reduced map from another run, get RA-DEC limits from this map
        if oldmap:
		print "Size from old map"
		ms=restoreFile(oldmap)
		ra1,ra2,dec1,dec2,oversamp=getMapLimits(SCAN_LIST[0],oldmap)
	else:
        # else  if we want to reduce all scans in the list, get RA-DEC limits from the first scan
		if overwrite==1:
			print "Size from scan - overwrite"
			ra1,ra2,dec1,dec2=computeScanRaDecLimits(SCAN_LIST[0],oversamp=oversamp)
		else:
	           #check if there is already a map that has been reduced
		   com='ls scanReduced/'+source+'_'+str.upper(febe)+'*.data'
		   result = commands.getoutput(com)
		   if result.count('No such file or directory') == 0:
			   print "Size from reduced map"
			   DATA_ALREADY_REDUCED=result.split()
			   ra1,ra2,dec1,dec2,oversamp=getMapLimits(SCAN_LIST[0],DATA_ALREADY_REDUCED[0])
			   

		   else:
			   print "Size from Scan - no map yet"
			   ra1,ra2,dec1,dec2=computeScanRaDecLimits(SCAN_LIST[0],oversamp=oversamp)  
		   
		 
	

	### CREATE LOG FILE IF IT DOES NOT EXIST ###
        if modelfile:
            logname=source+'_'+str.upper(febe)+'_model.log'
        else:
            logname=source+'_'+str.upper(febe)+'.log'
	print logname
	com='ls '+logname
	result = commands.getoutput(com)
	if result.count('No such file or directory') == 1:
		output=open(logname,'w')
		output.close()
		
	if overwrite == 1:
		output=open(logname,'w')
		output.close()
	
	 

	### Loop over all Scans ###	

	for i in range(len(SCAN_LIST)):
            try: 
		  test = read(SCAN_LIST[i],febe='AMKID870-AMKID870BE')
		  if test != -1:
			  if overwrite==1:
				  doReduce = 1
			  else:   
                                  if modelfile:
                                      com='ls scanReduced/'+source+'_'+str.upper(febe)+'_'+str(SCAN_LIST[i])+'_model.data'
                                  else:
                                      com='ls scanReduced/'+source+'_'+str.upper(febe)+'_'+str(SCAN_LIST[i])+'.data'
				  result = commands.getoutput(com)
				  if result.count('No such file or directory') == 1:
					doReduce = 1
				  else:
					doReduce = 0
                                        if modelfile:
                                            file='scanReduced/'+source+'_'+str.upper(febe)+'_'+str(SCAN_LIST[i])+'_model.data'
                                        else:
                                            file='scanReduced/'+source+'_'+str.upper(febe)+'_'+str(SCAN_LIST[i])+'.data'
					m=restoreFile(file)
					if ms:
						ms = mapsumfast([ms,m])
					else:
						ms = copy.deepcopy(m)
		  
		  
			  
			  
			  if doReduce ==1:  
				  

                                  sensitivity,taucorr,calcorr,jumps_nr,despike_nr,tau=redmapfkt(SCAN_LIST[i],calscan=calscan,febe=febe,taufile=taufile,calfile=calfile,weak=weak,modelfile=modelfile,sub=sub,cliplevel=cliplevel)
				  
				  

				  scandate=data.ScanParam.DateObs
				  scannumber=data.ScanParam.ScanNum
				  scanel  = fStat.f_mean(data.ScanParam.El)
				  
				  
				  output=open(logname,'a')
				  output.write('%i  %s %5.1f %5.3f %5.3f %5.1f %5.3f %i %i\n' %(scannumber,scandate,sensitivity,taucorr,calcorr,scanel,tau,jumps_nr,despike_nr))
				  output.close()
				  
				  data.doMap(system='EQ',sizeX =[ra1,ra2],sizeY=[dec1,dec2], noPlot=1, oversamp=oversamp)
				  if writeScans==1:
                                          if modelfile:
                                              mapname='scanReduced/'+str(source)+'_'+str.upper(febe)+'_'+str(scannumber)+'_model.data'
                                          else:
                                              mapname='scanReduced/'+str(source)+'_'+str.upper(febe)+'_'+str(scannumber)+'.data'
					  data.Map.dumpMap(mapname)
		 
				  
				  if ms:
					  ms = mapsumfast([ms,data.Map])
				  else:
					  ms = copy.deepcopy(data.Map)
					  
  
		          #### PLOT COADDED DATA ####
                          if doReduce ==1:  
                              mp=copy.deepcopy(ms)
                              m1,m2 = fStat.minmax(mp.Coverage)
                              test = mp.Coverage > covclip*m2
                              mp.Data = np.where(test,mp.Data,float('nan'))

                              Plot.panels(2,2)
                              Plot.nextpage()
                              mp.display(noerase=1,caption="RAW PLOT",aspect=aspect)
			      

                              sm = copy.deepcopy(mp)
                              sm.smoothBy(smooth/3600.0)
                              Plot.nextpage()
                              sm.display(noerase=1, caption = "SMOOTH PLOT",aspect=aspect)
			      

                              snr = copy.deepcopy(sm)
                              snr.Data *= np.sqrt(snr.Weight)
                              a = snr.computeRms()
                              scale=snr.RmsBeam
                              snr.Data /= np.array(scale,'f')
			      Plot.nextpage()
			      snrmin,snrmax=fStat.minmax(snr.Data)
			      if snrmax < 20:
				      snr.display(noerase=1, caption = "S/N PLOT",aspect=aspect)
			      else:
				      snr.display(noerase=1, caption = "S/N PLOT",limitsZ=[-5,20],aspect=aspect)
                              
                              
                          #-------------------     
		    
            except Exception,e:
                print e
                pass
         
	#COVERAGE CLIPPING
	
        if ms != 0:
            m1,m2 = fStat.minmax(ms.Coverage)
            test = ms.Coverage > covclip*m2
            ms.Data = np.where(test,ms.Data,float('nan'))

            rmsm=copy.deepcopy(ms)
            m1,m2 = fStat.minmax(rmsm.Coverage)
            test = rmsm.Coverage > 0.7*m2
            rmsm.Data = np.where(test,rmsm.Data,float('nan'))
        
        
            Plot.panels(2,2)
            Plot.nextpage()
            source=data.ScanParam.Object
            smoothed_string = str(int(smooth))
            if modelfile:
               smoothed_string = str(int(smooth))
               outname_raw     = source+'_'+str.upper(febe)+'_model.fits'
               outname_smooth  = source+'_'+str.upper(febe)+'_model_sm'+smoothed_string+'.fits'
               outname2_raw    = source+'_'+str.upper(febe)+'_model.data'
               outname2_smooth = source+'_'+str.upper(febe)+'_model_sm'+smoothed_string+'.data'
            else:
               outname_raw     = source+'_'+str.upper(febe)+'.fits'
               outname_smooth  = source+'_'+str.upper(febe)+'_sm'+smoothed_string+'.fits'
               outname2_raw    = source+'_'+str.upper(febe)+'.data'
               outname2_smooth = source+'_'+str.upper(febe)+'_sm'+smoothed_string+'.data'
               

	#RAW DATA MAP 
            ms.display(noerase=1,caption="RAW PLOT",aspect=aspect)
            ms.dumpMap(outname2_raw)
            writeFits2(ms,outname_raw,overwrite=1)
            print '''___________________________________________________'''
            print "Intensity raw map",
            rmsm.computeRms()
            rms_val = rmsm.RmsBeam
            sm = copy.deepcopy(ms)
            sm.smoothBy(smooth/3600.0)
            rmssm=copy.deepcopy(rmsm)
            rmssm.smoothBy(smooth/3600.0)
            Plot.nextpage()
	#SMOOTH MAP DISPLAY
            sm.dumpMap(outname2_smooth)
            writeFits2(sm,outname_smooth,overwrite=1)

            sm.display(noerase=1, caption = "SMOOTH PLOT",aspect=aspect)
            print "Intensity smooth map",
            rmssm.computeRms()
            rms_sm=rmssm.RmsBeam

            resra = -1
            resfwhm1 = 0.0 
            resfwhm2 = 0.0
            try:
                resra,resdec,flux,dflflux,resfwhm1,resfwhm2=snrfit(mm=sm,smooth=0.)
            except:
                flux = -1
            ### get integrated flux on smoothed map if source is extended###

            if (resfwhm1>1.1*sm.BeamSize*3600) or (resfwhm2>1.1*sm.BeamSize*3600):
                  intflux,dintflux=IntFlux(sm,3.0)
            else:
                  intflux = -1
            #print '''	___________________________________________________	'''

            snr = copy.deepcopy(sm)
            snr.Data *= np.sqrt(snr.Weight)
            a = snr.computeRms()
            scale=snr.RmsBeam
            snr.Data /= np.array(scale,'f')
            rms=copy.deepcopy(sm)
            rms.Data =  (rms.Data*0.0+1.0)/np.sqrt(rms.Weight)
            rms.Data *= np.array(scale,'f')
            
            coords=rms.physicalCoordinates()
            rmsmin=1.0e5
            distance=1e5
            nbX,nbY=np.shape(rms.Data)
            for i in range(nbX):
                for j in range(nbY):
                    if not_nan(rms.Data[i,j]):
                        if rms.Data[i,j] < rmsmin:
                            rmsmin=rms.Data[i,j]
                            bestra=i
                            bestdec=j
                        disttest=np.sqrt((coords[0][i][j]-raref)**2+(coords[1][i][j]-decref)**2)
                        if disttest < distance:
                            distance=disttest
                            rapoint=i
                            decpoint=j
                        
            rmscenter = rms.Data[rapoint,decpoint]
    
            Plot.nextpage()
	#S/N plot
	    snrmin,snrmax=fStat.minmax(snr.Data)
	    if snrmax < 20:
		    snr.display(noerase=1, caption = "S/N PLOT",aspect=aspect)
	    else:
		    snr.display(noerase=1, caption = "S/N PLOT",limitsZ=[-5,20],aspect=aspect)
            Plot.nextpage()
	#Plot.pgpap(5.0,0.75)
	#Plot.pgsvp(0.15,0.8,0.2,0.8)
            Plot.pgswin(0,300,0,30)
            Plot.pgbox('ABC',0,0,'ABC',0,0)
            Plot.pgsch(1.5)
            Plot.pgtext(10,28,str(source))
            Plot.pgtext(10,25,'mean RMS raw map : ')
            Plot.pgtext(10,23,'mean RMS smooth map : ')
            Plot.pgtext(10,21,'best RMS in smoothed map: ')
            Plot.pgtext(10,19,'RMS at pointing center: ')
            if rms_val < 1.0:
                Plot.pgtext(180,25,str(round(1000*rms_val,1))+' [mJy/beam]')
                Plot.pgtext(180,23,str(round(1000*rms_sm,1))+' [mJy/beam]')
                Plot.pgtext(180,21,str(round(1000*rmsmin,1))+' [mJy/beam]')
                Plot.pgtext(180,19,str(round(1000*rmscenter,1))+' [mJy/beam]')
            else:
                Plot.pgtext(180,25,str(round(rms_val,2))+' [Jy/beam]')
                Plot.pgtext(180,23,str(round(rms_sm,2))+' [Jy/beam]')
                Plot.pgtext(180,21,str(round(rmsmin,2))+' [Jy/beam]')
                Plot.pgtext(180,19,str(round(rmscenter,2))+' [Jy/beam]')
            if flux > 1.0:
                Plot.pgtext(20,15,'Peak Flux: '+str(round(flux,2))+' +/- '+str(round(dflflux,2))+' [Jy/beam]')
            if flux > 0.0 and flux < 1.0:
                Plot.pgtext(20,15,'Peak Flux: '+str(round(1000*flux,1))+' +/- '+str(round(1000*dflflux,1))+' [mJy/beam]')
            if intflux > 1.0:
                Plot.pgtext(20,13,'Integrated Flux: '+str(round(intflux,2))+' +/- '+str(round(dintflux,2))+' [Jy]')
            if intflux > 0.0 and intflux < 1.0:
                Plot.pgtext(20,13,'Integrated Flux: '+str(round(1000*intflux,1))+' +/- '+str(round(1000*dintflux,1))+' [mJy]')
            if resra != -1:
                Plot.pgtext(20,11,'Position: '+str(resra)+' '+str(resdec))
                
	#Plot.pgtext(150,25,str(round(1000*rms_val,4))+' [mJy/beam]')
	#Plot.pgtext(150,20,str(round(1000*rms_sm,4))+' [mJy/beam]')
                print '''to reload the image use the command data = newRestoreData("source_name.data")'''
	Plot.panels(1,1)
        ### CREATE EPS IMAGE ###
        if modelfile:
            epsname=source+'_'+str.upper(febe)+'redscans_model.ps/CPS'
        else:
            epsname=source+'_'+str.upper(febe)+'redscans.ps/CPS'
        op(epsname)
        Plot.panels(2,2)
        Plot.nextpage()
        ms.display(noerase=1,caption="RAW PLOT")
        Plot.nextpage()
        sm.display(noerase=1, caption = "SMOOTH PLOT")
        Plot.nextpage()
        snr.display(noerase=1, caption = "S/N PLOT")
        Plot.nextpage()
        Plot.pgswin(0,300,0,30)
        Plot.pgbox('ABC',0,0,'ABC',0,0)
        Plot.pgsch(1.5)
        Plot.pgtext(10,28,str(source))
        Plot.pgtext(10,25,'mean RMS raw map : ')
        Plot.pgtext(10,23,'mean RMS smooth map : ')
        Plot.pgtext(10,21,'best RMS in smoothed map: ')
        Plot.pgtext(10,19,'RMS at pointing center: ')
        if rms_val < 1.0:
                Plot.pgtext(180,25,str(round(1000*rms_val,1))+' [mJy/beam]')
                Plot.pgtext(180,23,str(round(1000*rms_sm,1))+' [mJy/beam]')
                Plot.pgtext(180,21,str(round(1000*rmsmin,1))+' [mJy/beam]')
                Plot.pgtext(180,19,str(round(1000*rmscenter,1))+' [mJy/beam]')
        else:
                Plot.pgtext(180,25,str(round(rms_val,2))+' [Jy/beam]')
                Plot.pgtext(180,23,str(round(rms_sm,2))+' [Jy/beam]')
                Plot.pgtext(180,21,str(round(rmsmin,2))+' [Jy/beam]')
                Plot.pgtext(180,19,str(round(rmscenter,2))+' [Jy/beam]')
            
        if flux > 1.0:
                Plot.pgtext(20,15,'Peak Flux: '+str(round(flux,2))+' +/- '+str(round(dflflux,2))+' [Jy/beam]')
        if flux > 0.0 and flux < 1.0:
                Plot.pgtext(20,15,'Peak Flux: '+str(round(1000*flux,1))+' +/- '+str(round(1000*dflflux,1))+' [mJy/beam]')
        if intflux > 1.0:
                Plot.pgtext(20,13,'Integrated Flux: '+str(round(intflux,2))+' +/- '+str(round(dintflux,2))+' [Jy]')
        if intflux > 0.0 and intflux < 1.0:
                Plot.pgtext(20,13,'Integrated Flux: '+str(round(1000*intflux,1))+' +/- '+str(round(1000*dintflux,1))+' [mJy]')
        if resra != -1:
                Plot.pgtext(20,11,'Position: '+str(resra)+' '+str(resdec))
        close()
        ### WRITE ASCII RESULT FILE
        if modelfile:
            asciiname=source+'_'+str.upper(febe)+'_model_result.txt'
        else:
            asciiname=source+'_'+str.upper(febe)+'_result.txt'
        
        output=open(asciiname,'w')
        output.close()
        output=open(asciiname,'a')
        if rms_val < 1.0:
            output.write('mean RMS in unsmoothed map: %5.1f mJy/b\n' %(1000.*rms_val))
            output.write('mean RMS in smoothed map  : %5.1f mJy/b\n' %(1000.*rms_sm))
            output.write('best RMS in smoothed map  : %5.1f mJy/b\n' %(1000.*rmsmin))
            output.write('RMS at pointing center    : %5.1f mJy/b\n\n' %(1000.*rmscenter))
            if flux != -1:
                output.write('Peak Flux: %5.1f +/- %5.1f mJy/b\n' %(1000*flux,1000*dflflux))
            else:
                output.write('No source detected\n')
                
        else:
            output.write('mean RMS in unsmoothed map: %5.2f Jy/b\n' %(rms_val))
            output.write('mean RMS in smoothed map  : %5.2f Jy/b\n' %(rms_sm))
            output.write('best RMS in smoothed map  : %5.2f Jy/b\n' %(rmsmin))
            output.write('RMS at pointing center    : %5.2f Jy/b\n\n' %(rmscenter))
            if flux != -1:
                output.write('Peak Flux: %5.2f +/- %5.2f Jy/b\n' %(flux,dflflux))
            else:
                output.write('No source detected\n')
        if resra != -1:
                output.write('Position: %s %s\n'%(resra,resdec))
        output.close()    



def redmapfkt(ScanNr=0,calscan='',febe='LFA',taufile='',calfile='',weak=1,modelfile='',sub=0,cliplevel=3.0):
    '''
    Reduces a LABOCA map on a pointing source.
    Optional parameters: Scan number, zenith opacity
    '''
    if (ScanNr != 0):
        if (febe == 'LFA') or (febe=='lfa'):
            tst=read(str(ScanNr),febe='AMKID870-AMKID870BE')
        if (febe == 'HFA') or (febe=='hfa'):
            tst=read(str(ScanNr),febe='AMKID350-AMKID350BE')
        if  febe == 'Laboca' or febe== 'LABOCA' or febe == 'laboca':
            tst = read(str(ScanNr))
        
        if tst:
            print "Scan %s not readable"%(str(ScanNr))
            return -1,-1,-1,-1
   
        if febe == 'Laboca' or febe== 'LABOCA' or febe == 'laboca':
            CntstoV(data)
            data.Data *= np.array(VtoJy,'f')

        if (febe == 'LFA') or (febe=='lfa') or (febe == 'HFA') or (febe=='hfa'):
            flagEmpty()	
	    if calscan:
		    search = "*%i*%s_KidCalspline.dat"%(calscan,str(data.BolometerArray.FeBe))
		    globlist=glob(search)
		    if len(globlist) == 1:
			    tempcalfile=globlist[0]
		    else:
			    print "No calibration file for scan %i found - check scan nr or reduce kid caibration first"%(scan)
			    return
		    applyTempCal(tempcalfile,doSkyFF=1)
		    applyJyperK()
	    else:
		    calibrateAMKID()


        nr=len(data.ScanParam.MJD)
        nr = nr -1
        mjdref = (data.ScanParam.MJD[nr]-data.ScanParam.MJD[0])/2.+data.ScanParam.MJD[0]

	rcp=getMKIDsRCP(mjdref)
        print rcp 
        updateRCP(rcp)
        flagRCP(rcp)
        
	
        calcorr = 0.0
	
        if calfile:
           calcorr = getCalCorr(mjdref,'linear',calfile)
        else:
           calcorr = 1.0

        if calcorr > 0.0 and calcorr != 1.0:
            data.Data /= np.array((calcorr),'f')
        else:
	    print "No calbration correction found - assuming 1.0"
            calcorr = 1.0	
	
	
        data.zeroStart()
    
	

        febe = data.BolometerArray.FeBe
        despike_nr = 0
        jumps_nr = 0

        data.flagSpeed (below=30.)
        data.flagSpeed(above=500.)
        data.flagAccel(above=800.)



        if modelfile:
            model=restoreFile(modelfile)
            modeln=copy.deepcopy(model)
            model.Data*=np.sqrt(model.Weight)
            model.smoothBy(7.0/3600.0)
	    model.computeRms()
            scale=model.RmsBeam
            model.Data /= np.array(scale,'f')
            if sub ==0:
                data.flagSource(threshold=cliplevel,model=model)
            else:
                nbX,nbY = np.shape(modeln.Data)
                for i in range(nbX):
                    for j in range(nbY):
                        if is_nan(modeln.Data[i,j]) or model.Data[i,j]<cliplevel:
                            modeln.Data[i,j] = 0.
                tau = 0.0
                if taufile:
                    tau = getTau(mjdref,'linear',taufile)
                else:
                    tau = scanTau()
                if tau > 0.0:
                    scanel  = fStat.f_mean(data.ScanParam.El)
                    taucorr = np.exp(tau/np.sin(scanel * np.pi / 180.))
                    modeln.Data = -1*modeln.Data/taucorr
                    modelp=copy.deepcopy(modeln)
                    modelp.Data = -1*modelp.Data
                    data.addSource(model=modeln)
                else:
                    print "No tau correction found - Model subtraction not possible!!!"
                    return -1,-1,-1,-1,-1,-1
        else:
            data.flagPosition(radius=60,flag=8)
	#decorrelation per chain
	correlChains(chanRef=-2,factor=0.95,nbloop=3)
	data.flagFractionRms(ratio=5)
	if  febe == 'LFA'  or  febe=='lfa' or febe == 'AMKID870-AMKID870BE':
		#correlChains(chanRef=-2,factor=0.95,nbloop=3)
		chainIndices = range(4)
		#decorrelation on goups in chains
		for chainIndex in chainIndices:
			correlGroups(chainIndex, chanRef=-2, factor=0.97, nbloop=3,numRepeats=1, groupingThreshold=0.6)
	if  febe == 'HFA'  or  febe=='hfa' or febe == 'AMKID350-AMKID350BE':
		correlChains(chanRef=-2,factor=0.95,nbloop=2)
		#decorrelation on goups in chains
		#chainIndices = range(20)
		#for chainIndex in chainIndices:
		#	correlGroups(chainIndex, chanRef=-2, factor=0.97, nbloop=3,numRepeats=1, groupingThreshold=0.6)

		
	    
        base(order=1)

        if weak == 1:
            data.flattenFreq(below=0.5,hiref=0.7)
            weight()
            if modelfile:
                if sub == 0:
                    unflag(flag=8)
                else:
                    data.blankFreq(below=0.8)
                    data.addSource(model=modelp)
            else:
                unflag(flag=8)
            #spikenr =despike(below=-5,above=5)
            #despike_nr = despike_nr+spikenr
	    sensi,wsensi,allsensi,bolos =calcsensitivity()
        else:
            if modelfile:
                if sub == 0:
			data.flattenFreq(below=0.5,hiref=0.7)
			weight()
			sensi,wsensi,allsensi,bolos =calcsensitivity()
			unflag(flag=8)
                else:  
                    #despike(below=-5,above=5)
		    for i in 20:
			    correlGroups(i,nbloops=2,numRepeats=1)
                    data.flattenFreq(below=0.5,hiref=0.7)
		    weight()
		    sensi,wsensi,allsensi,bolos =calcsensitivity()
                    data.addSource(model=modelp)
		    print "MODEL SUB WEAK=0"
            else:
                #data.flattenFreq(below=0.5,hiref=0.7)
		weight()
		print "now sensitivity calc"
		sensi,wsensi,allsensi,bolos =calcsensitivity()
                unflag(flag=8)
            

        print '#######################################'
        print "Sensitivity: %6.1f mJy sqrt(s)"%(sensi)
        print '#######################################'

        tau=0.0
        if taufile:
           tau = getTau(mjdref,'linear',taufile)
        else:
           tau = scanTau()
        if tau > 0.0:
            data.correctOpacity(tau)
            data._DataAna__statistics()
	    data.computeWeight()
            scanel  = fStat.f_mean(data.ScanParam.El)
            taucorr = np.exp(tau/np.sin(scanel * np.pi / 180.))
            return sensi,taucorr,calcorr,jumps_nr,despike_nr,tau
        else:
            print "No tau correction found!"
            return -1,-1,-1,-1,-1,-1
            

    else:
        print "No scan defined."
        return -1,-1,-1,-1,-1,-1

    

def addSource_fixed(self, model, chanList=[], factor=1.):
        """
        DES: add data to time stream according to a model map
        INP: (i list) chanList: the list of channels to work with
             (f)     factor: multiply by this factor (default 1)
             (Image object) model: the input model map (with WCS)
                            (default: use current data.Map)
        """

        if not model:
            model = self.Map

        if not np.any(model.Data):
            self.MessHand.error("no map computed yet, and no model provided")
            return

        chanList = self.BolometerArray.checkChanList(chanList)
        if len(chanList)<1: 
            self.MessHand.error("no valid channel")
            return
        chanListIndexes = self.BolometerArray.getChanIndex(chanList)

        if str.find(model.WCS['CTYPE1'], 'GLON') > -1:
            if not len(self.ScanParam.GalAngle):
                if not len(self.ScanParam.GLon):
                    self.ScanParam.computeGal()
                self.ScanParam.computeGalAngle()
            rotAngles = np.array(self.ScanParam.ParAngle) + np.array(self.ScanParam.GalAngle)
            XYOffsets = np.array([self.ScanParam.get('Glon', flag='None'), \
                               self.ScanParam.get('Glat', flag='None')])
        else:
            rotAngles = np.array(self.ScanParam.ParAngle)
            XYOffsets = np.array([self.ScanParam.get('RA', flag='None'), \
                               self.ScanParam.get('Dec', flag='None')])

        chanListAzEl = np.array(self.BolometerArray.UsedChannels)-1
        OffsetsAzEl=np.array((np.take(self.BolometerArray.Offsets[0, :], chanListAzEl), \
                           np.take(self.BolometerArray.Offsets[1, :], chanListAzEl)))
        
        refChOffsets=np.array((self.BolometerArray.RefOffX, self.BolometerArray.RefOffY), 'f')
        AXIS1 = np.array([model.WCS['NAXIS1'], model.WCS['CRPIX1'],
                       model.WCS['CDELT1'], model.WCS['CRVAL1'], 1.])
        AXIS2 = np.array([model.WCS['NAXIS2'], model.WCS['CRPIX2'],
                       model.WCS['CDELT2'], model.WCS['CRVAL2'], 1.])

        # get the new data + factor x model array
        tmp = fMap.addsource(chanListIndexes, self.Data, model.Data, \
                             XYOffsets, OffsetsAzEl, rotAngles, refChOffsets, \
                             AXIS1, AXIS2, factor)
        # replace self.Data with updated one
        self.Data = copy.copy(tmp)
        self._DataAna__resetStatistics()
        tmp = 0  # free memory       

def fullMap(oversamp=4,system='ho',limitsZ=[],sm=0):
    if sm == 0:
        sm = data.BolometerArray.BeamSize/3.

    if str.upper(system) == 'HO':
        aspect=1
    else:
        sm=sm/3600.
        aspect=0
    
        
    mapping(oversamp=oversamp,system=system,aspect=aspect,noPlot=1)
    data.Map.smoothBy(sm)
    if len(limitsZ) == 0:
        limitsZ=getZlimitsMap()
    caption="%s %s Scan %i"%(data.ScanParam.Object,str(data.BolometerArray.FeBe).split("-")[0],data.ScanParam.ScanNum)
    display(aspect=aspect,limitsZ=limitsZ,caption=caption)


def getZlimitsMap(frac=0.3):
    mask=np.where(data.Map.Coverage.real > np.nanmax(data.Map.Coverage.real)*frac)
    zmin=np.nanmin(data.Map.Data[mask])
    zmax=np.nanmax(data.Map.Data[mask])
    return [zmin,zmax]

def getZlimits(values,plow=0.05,phigh=0.95):
    
    if len(np.shape(values)) > 1:
        values=np.ravel(values)
    values=values[~np.isnan(values)]
    counts,bins=np.histogram(values,bins=100,density=True)
    counts=counts/np.sum(counts)
    zmin=-9999.
    zmax=-9999.
    for i in range(len(counts)):
        if np.sum(counts[0:i]) >= plow and zmin == -9999:
            zmin = bins[i]
        if np.sum(counts[0:i]) >= phigh and zmax == -9999:
            zmax = bins[i+1]
    limitsZ=[zmin,zmax]
    return limitsZ

# Ignore python warnings to avoid terminal clogging
warnings.filterwarnings("ignore")
