# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================
# Currently LFA-only

# --- Source and map parameters ---
source  = 'G345'            # As in observing logs
fe      = 'LFA'             # Frontend, either 'LFA' or 'HFA'
system  = 'GAL'             # Coordinate system for map, 'EQ', 'GAL' or 'HO' (default)
center  = [345.3, 1.7]      # Center of map in CHOSEN COORDINATES in deg
sizex   = 2.2               # Size of map in deg for X direction
sizey   = 2.2               # Size of map in deg for Y direction
padding = 0.6               # Padding around the map in deg for grid (default ~2x array)
doPlot  = True              # Display maps at each scan. If False, only final
                            # coadded map per iteration will be displayed.

# ----- Reduction parameters -----
# SUGGESTED: run all scans with niters=1, figure out bad scans using showMaps.py and then
# run with niters=2 or 3 ignoring bad scans
writeSummary = True         # Write summary of reductions or not
niters       = 3            # Number of iterations to run, 1 to 3 (recommended: 2 + PLANCK data)
clip         = -1           # Sigma clipping level (-1 or >=1.5) on noise map: masked where 
                            # noisemap > clip * mediannoise, else no clipping
flagJumps    = True         # Flag jumps/spikes in the data:
                            # recommended to set to True for LFA
smoothby_arcsec = 8.        # By how much to smooth final iteration maps. Default 8. arcsec
correctbeam     = True      # Whether to correct beam bookkeeping in final iteration maps

# ----- Scans ------
# If scans is empty, automatically retrieves all scans of the source
# specified above from the obslogs directory below
scans = []
obslogsdir = '~/obslogs'  # at MPIfR: '/apex-archive/obslogs/M-PROJECT.CODE-IN-CAPS/obslogs'

# Manually exclude bad scans if needed            
badscans = [32439, 33340, 33568, 34066, 34685, 34693, 34950]  
# VWO:
# 25088 -> half the map is missing (likely cancelled scan). Still good.
# 30659 -> despiking worked, but left NaN patches. Still good.
# 31953 -> could not even reduce, not readable, I and Q have different lengths - NOW WE CAN REDUCE AND IT'S FINE!
# 32439 -> absolutely no signal in map. Apparently camera got warm: IQBTs look like Fsweeps.
# 33340 -> still has spikes on top of map after despiking
# 33568 -> same as 32439
# 34066 -> despiking worked well, but big patches of the map are NaNs. Up to you to use - I removed it.
# 34685 -> despiking worked well, but big patches of the map are NaNs. Up to you to use - I removed it.
# 34693 -> same as 32439
# 34950 -> same as 32439

# ==============================
# ===== END OF USER INUPUT =====
# ==============================









# ===== REDUCTION CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
# NOTE: VWO: BoA smoothBy smooths weights with the kernel, but weights are 1 / rms^2 which is
# a non-linear scale. All smoothing should be done as:
# Sky_smoothed = Kernel * Sky
# Weights_smoothed = 1 / Variance_smoothed | with Variance_smoothed = Kernel^2 * Variance
# Coverage_smoothed = Kernel * Coverage

# NOTE 2: smoothBy does handle the sky convolution correctly, it's just the weights that are not
# correct after convolution

import warnings
import copy as copy
import BoaMapping as BOAMAP
from mars.fortran import fMap

# define the good functions :)
def findSciTargetScans(source, obslogsdir, verbose=True):
    scanlist = []
    files = os.listdir(obslogsdir)
    c=0
    for file in files:
        fullfilename = obslogsdir + file if obslogsdir[-1]=='/' else obslogsdir + '/' + file
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
                    if key == 'Scan':
                        scan_int = int(line[4:-6])
                        scan = (line[4:-6].ljust(6) + ' | ')       # 0
                    if key == 'Source':
                        src = (line[4:-6].ljust(12) + ' | ')       # 1  
                    if key == 'Scan type':
                        scantype = (line[4:-6].ljust(12) + ' | ')  # 2
                    if key == 'Observ. mode':
                        mode = (line[4:-6].ljust(12) + ' | ')      # 3
                    if key == 'Scan duration':
                        duration = (line[4:-6].ljust(12) + ' | ')  # 4
                    if key == 'Scan status':
                        status = (line[4:-6].ljust(12) + ' | ')    # 5
                    if key == 'Comment':
                        comment = (line[4:-6])             # last                        
                        message += scan + src + scantype + mode + duration + status

                start = False

                if source in src:
                    if  '-999' not in duration:
                        if 'MAP' in scantype and 'OTF' in mode and 'OK' in status and 'warm' not in str.lower(comment):
                            message += 'SCAN CONSIDERED'.ljust(15) + ' | ' + comment
                            scanlist.append(scan_int)
                        else:
                            message += 'SCAN DISCARDED'.ljust(15) + ' | ' + comment
                    else:
                        message += 'SCAN ONGOING'.ljust(15) + ' | ' + comment
                    if verbose:
                        print(message)

        if c==0 and len(keys)!=0 and verbose:
            print('============')
            print('OBSLOG KEYS:')
            print('============')
            print(keys)
            c+=1
    scanlist.sort()
    info("Number of 'MAP' scans on science target %s: %i"%(source, len(scanlist)))
    return scanlist



def auxsmoothby(m, Size):
    '''
    BoA-like smoothing but with correct variance propagation.

    - Data: convolved with K (same as BoA)
    - Weight: propagated via variance (K^2)
    - Coverage: convolved with K (same as BoA)
    '''
    # Build kernel (not normalized)
    pixsize = abs(m.WCS['CDELT2'])
    K0 = BOAMAP.Kernel(pixsize, Size).Data.astype(float)

    # Normalize kernel
    K = K0 / np.sum(K0)

    # Create elementwise-squared kernel for variance
    K2 = K**2

    # Smooth INTENSITY (same as BoA)
    #   I' = K * I     =     (K0/sum(K0_i)) * I
    # and ksmooth does
    #   I' = (K * I) / sum(K_i), but since sum(K_i)=1
    # then ksmooth does effectively
    #   I' = K * I, all good
    I1 = fMap.ksmooth(m.Data, K)

    # Correct variance propagation for weights:
    #   V' = K2 * V     =     (K0/sum(K0_i))^2 * V
    # but ksmooth does
    #   V' = (K2 * V) / sum(K2_i), and now sum(K2_i)!=1
    # then ksmooth does effectively
    #   V' = K2/sum(K2_i) * V
    # so an additional multiplication by sum(K2_i) is needed
    # to get back from ksmooth:
    #   V' = K2/sum(K2_i) * V * sum(K2_i) = K2 * V
    V0 = np.where(m.Weight > 0.0, 1.0 / m.Weight, np.NaN)
    V1 = fMap.ksmooth(V0, K2) * np.sum(K2)

    # Smooth COVERAGE (same as BoA)
    C1 = fMap.ksmooth(m.Coverage, K)
    
    # new scale per beam for Jy/beam units
    newbeam = np.sqrt(m.BeamSize**2 + Size**2)
    scale = (newbeam**2 / m.BeamSize**2)
    I1 *= scale  # now in Jy/newbeam
    V1 *= scale**2  # now in Jy^2/newbeam^2

    # create new weight map
    W1 = np.where(V1 > 0.0, 1.0 / V1, 0.0)

    # Update map with correct Jy/beam scale
    m.Data = I1
    m.Weight = W1
    m.Coverage = C1
    m.BeamSize = newbeam



def auxwriteFits(data=None,outfile='boaMap.fits',overwrite=0,limitsX=[],limitsY=[],intensityUnit="Jy/beam",clip=-1):
    """
    DES: store the current map (2D array with WCS info) to a FITS file
    INP: (str)   outfile: output file name (default boaMap.fits)
         (bool) overwrite: overwrite existing file -
                          default = 0: do not overwrite existing file
         (f list) limitsX/Y: optional map limits (in world coordinates)
         (string) intensityUnit: optional unit of the intensity (default: "Jy/beam")
    """
    from mars import BoaFits

    if os.path.exists(outfile):
        if not overwrite:
            print('File %s exists' % outfile)
            return
    info('Exporting map to fits file:')
    print(outfile)
    if not data:
        data = data.Map
    try:
        dataset = BoaFits.createDataset("!" + outfile)
    except Exception, data:
        print('Could not open file %s: %s' % (outfile, data))
        return

        
    localMap = copy.deepcopy(data)
        
    try:
        # RMS map creation
        rmsMap = copy.deepcopy(localMap)  # Signal
        rmsMap.Data = np.where(rmsMap.Weight > 0.0, 1.0 / np.sqrt(rmsMap.Weight), np.NaN)  # Noise = 1/sqrt(weight)

        # SNR map creation
        snrMap = copy.deepcopy(localMap)  # Signal
        snrMap.Data = np.where(snrMap.Weight > 0.0, snrMap.Data * np.sqrt(snrMap.Weight), np.NaN)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

        if clip > 0:
            info('Clipping map to %.1f*medianRMS (inner contour on display)...'%clip)
            mediannoise = np.nanmedian(rmsMap.Data)
            mask = np.where(rmsMap.Data > clip * mediannoise)
            localMap.Data[mask] = np.NaN
            rmsMap.Data[mask] = np.NaN
            snrMap.Data[mask] = np.NaN
            del mask  # free memory
 
        #write FLux plane                                                            
        localMap._Image__writeImage(dataset, "Intensity", intensityUnit=intensityUnit)
        #write RMS plane
        rmsMap._Image__writeImage(dataset, "Intensity", intensityUnit=intensityUnit+" (RMS)")
        #write SNR plane
        snrMap._Image__writeImage(dataset, "Intensity", intensityUnit='SNR')
        dataset.close()
            
    except Exception, data:
        try:
            dataset.close()
        except:
            pass
        print('Could not write data to file %s: %s' % (outfile, data))
        return
    
    del localMap  # free memory
    del snrMap  # free memory
    del rmsMap  # free memory
    del dataset  # free memory



# variable checks
if fe not in ['LFA', 'HFA']:
    raise ValueError("fe must be either 'LFA' or 'HFA'.")
if system not in ['EQ', 'GAL', 'HO']:
    raise ValueError("system must be either 'EQ', 'GAL', or 'HO'.")
if niters < 1 or niters > 3:
    raise ValueError("niters must be 1, 2, or 3.")
if clip < 1.5 and clip!=-1:
    raise ValueError("clip must be -1 (no clipping) or >= 1.5.")
if sizex + 2*padding > 360 or sizey + 2*padding > 180:
    raise ValueError("Your map is bigger than the sky...")

# Create map bounds
info('Creating map boundaries...')
biggerX = center[0] + sizex/2 + padding
smallerX = center[0] - sizex/2 - padding
biggerY = center[1] + sizey/2 + padding
smallerY = center[1] - sizey/2 - padding

# These can't happen
if biggerY > 90:
    raise ValueError('STOPPING SCRIPT: The upper border of the map has Y coordinate > +90 degrees! are you sure this is intended?')
if smallerY < -90:
    raise ValueError('STOPPING SCRIPT: The lower border of the map has Y coordinate < -90 degrees! are you sure this is intended?')

# Check X reframing.
# Example with an X width of 10 deg:
# Case 1: left = 150, right = 140 is left untouched
# Case 1: left = 185, right = 175 -> frame was 0:360, now left = -175, right = 175
# Case 3: left = 200, right = 190 -> frame was 0:360, now left = -160, right = -170
# Case 4 : same as before but one of the boundaries ended up < -180: add 360
sysreframe = False
if biggerX > 180 and system!='EQ':
    biggerX -= 360
    sysreframe = True
if biggerX < -180 and system!='EQ':
    biggerX += 360
    sysreframe = True
if smallerX > 180 and system!='EQ':
    smallerX -=360
    sysreframe = True
if smallerX < -180 and system!='EQ':
    smallerX +=360
    sysreframe = True

# information
if sysreframe:
    info('Map X boundaries were wrapped into the range [-180, 180] deg')

# Define boundary list for functions
ysize = [smallerY, biggerY]
# For EQ or GAL biggerX is to the left because X angle
# follows right-hand rule with thumb pointing to EQ or GAL north pole
xsize = [biggerX, smallerX]

# For HO smallerX is to the left because X angle
# follows left-hand rule with thumb pointing to zenith (eastward in ground)
if system =='HO':
    xsize = [smallerX, biggerX]

# find project home folder based on where MARS loaded and re-define obslogsdir
if obslogsdir == '~/obslogs':
    currdir = os.getcwd()
    splitted = currdir.split('/')
    projectidx = None
    for i in range(len(splitted)):
        # project code is separated once with dot and thrice with dash
        if len(splitted[i].split('.')) == 2 and len(splitted[i].split('-')) == 4:
            projectidx = i
    if projectidx != None:
        obslogsdir = '/homes/%s/obslogs'%splitted[projectidx]
    else:
        raise ValueError("STOPPING SCRIPT: Project code could not be extracted from: %s"%currdir)

if len(scans) == 0 and not os.path.exists(obslogsdir):
    raise ValueError('STOPPING SCRIPT: Either enter scans or an existing obslogs directory...')

# find scans if not provided
if len(scans) == 0 and os.path.exists(obslogsdir):
    info('Retrieving source scan numbers from ObsLogs...')
    scans = findSciTargetScans(source=source, obslogsdir=obslogsdir)
    if len(scans) == 0:
        raise ValueError('No scans of source %s found in ObsLogs directory: %s!'%(source, obslogsdir))

# sort scans
scans.sort()

# Remove bad scans from the list of scans to be reduced
for badscan in badscans:
    if badscan in scans:
        scans.remove(badscan)
        info('Scan %i in bad scans manually removed'%badscan)

# Check removing bads did not leave scans empty
if len(scans) == 0:
    raise ValueError('There are no good scans after removing bad scans list...')

# Define standardized "myname" variable for output files
myname = str(fe) + "-" + str(source) + "-" + str(system)
if flagJumps:
    myname += "-flagJumps"

# Set noPlot
if not doPlot:
    noPlot = True
else:
    noPlot = False

# Create directory for reduced files if it doesn't exist
if os.path.exists("ReducedFiles") == False:
    os.makedirs("ReducedFiles")

# Create directory for reduced files if it doesn't exist
if writeSummary and os.path.exists("Summaries") == False:
    os.makedirs("Summaries")

# smoothby to deg
smoothby_deg = smoothby_arcsec / 3600.

print('')
print('''\
=====================
Reduction parameters:
=====================
Source:             %s
Frontend:           %s
Coordinate system:  %s
Map center:         %s, %s deg
Map size (x,y):     %s, %s deg
Padding:            %s deg
Map Boundaries:     %s, %s deg in x; %s, %s deg in y
Iterations:         %i
Sigmaclip level:    %s
Flag jumps:         %s
Smoothing:          %s arcsec
Number of scans     %s'''%(source, fe, system, center[0], center[1], sizex, sizey, padding,
                           xsize[0], xsize[1], ysize[0], ysize[1], niters, clip, flagJumps,
                           smoothby_arcsec, len(scans)))

# ===========================
# Beginning of reduction loop
# ===========================
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for iter in range(1, niters+1):
        print('')
        print("####################################################################")
        print("####################### Iteration %i starting #######################"%(iter))
        print("####################################################################")

        if iter == 1:
            # First iteration -- no model
            mymodel = None
            subtract = False
        else:
            # retrieve last iteration map
            mymodel = "ReducedFiles/" + str(myname) + "-coadded-flux-iter" + str(iter-1) + ".data"
            coadded = restoreFile(mymodel)

            if iter == 2:
                subtract = False
                mymodel = createSourceModel(coadded, highcut=5.5, lowcut=2.5, sm=0., mtype='snr', clip=3)
            
            if iter == 3:
                subtract = True
                mymodel = createSourceModel(coadded, highcut=5.5, lowcut=2.5, sm=0., mtype='flux', clip=3)

            del coadded  # free memory

  
        # Initialize co-added map
        ms = None

        for i, scan in enumerate(scans):
            scanname = "ReducedFiles/"+str(myname)+"-"+str(scan)+"-iter"+str(iter)+".data"
            globlist = glob(scanname)

            # Initialize map for this scan
            m = None

            # Check if reduction does not exist
            if len(globlist) ==  0:
                print('')
                print('')
                print('')
                info('Reducing scan %s (iteration %i)...'%(scan, iter))

                # Reduce it
                redweak(scan,fe=fe,size=-1,model=mymodel,subtract=subtract,doPlot=doPlot,extremeFilter=False,writeSummary=writeSummary,flagJumps=flagJumps)
                # NOTE: redweak's summary is everything about the timelines, nothing about map.
                # NOTE 2: redweak then runs mapping in horizontal coords, forces a 10" (LFA) or 4.5"(HFA) smoothing
                # and tries to solve for pointing corrections on smoothed map. Then prints timeline sensitivity
                # and pointing corrections in smoothed maps. This is fine.

                # Immediately rename summary and move to new folder
                if writeSummary:
                    # VWO: made it iteration-specific
                    origname = "%s-%s-%i_summary.txt"%(fe, data.ScanParam.Object, data.ScanParam.ScanNum)
                    newdir = "Summaries/"
                    newname = myname + "-" + str(scan) + "-iter" + str(iter) + "_summary.txt"
                    outname = newdir + newname
                    os.rename(origname, outname)

                # Flagging example to flag a certain time range in a map (seconds from the beining of the scan) 
                #if scan == 22919: 
                #    flagMJD(above=1430, below=1600,flag=2)

                # Flagging example to flag a certain tone/KID in a scan
                #if scan == 28517:
                #    flagC(3353, flag=2)

                # Create map in chosen system and chosen box
                # where pixsize = BEAM_FWHM / oversamp
                mapping(oversamp=4,system=system,sizeX=xsize,sizeY=ysize,limitsZ=[-0.8,1.5],noPlot=noPlot)
                # NOTE: this has a smooth parameter, but is default 0
                # NOTE 2: data.Map.BeamSize is taken from data.BolometerArray.BeamSize
                # NOTE 3: data.BolometerArray.BeamSize is just 1.22 * lambda / D * 180/pi, not from beammap!
            
                # Save unsmoothed map, "native" resolution (m.BeamSize = data.BolometerArray.BeamSize)
                data.Map.dumpMap(scanname)

                # Assign BoA map to variable m
                m = restoreFile(scanname)

                # For this scan, add non-noisy area and median noise info to summary
                if writeSummary:
                    info('Smoothing copy of map for summary at final resolution...')
                    # copy map and smooth with same kernel as final file
                    m_smooth = copy.deepcopy(m)
                    auxsmoothby(m_smooth, smoothby_deg)

                    # Create smoothed noise map
                    rmsArray = np.where(m_smooth.Weight > 0.0, 1.0 / np.sqrt(m_smooth.Weight), np.NaN)

                    # Statistics and write
                    minnoise = np.nanmin(rmsArray)
                    mask = (rmsArray > 5*minnoise)
                    rmsArray[mask] = np.NaN
                    pixelsize = np.abs(m.WCS['CDELT2'])
                    nrpix = np.sum(~np.isnan(rmsArray))
                    area = nrpix * pixelsize**2
                    noise = np.nanmedian(rmsArray)
                    f = open(outname,'r')
                    lines = f.readlines()
                    f.close()
                    myline = lines[0].replace("\n","")
                    myline = myline+",{:.3f},{:.4f}\n".format(area,noise)
                    f = open(outname,'w')
                    f.write(myline)
                    f.close()

                    del m_smooth  # free memory
                    del rmsArray  # free memory
                    del mask  # free memory
            
            else:
                # Retrieve BoA map
                info('Reduction for scan %i (iteration %i) found. Loading...'%(scan, iter))
                m = restoreFile(scanname)

            if np.all(np.isnan(m.Data)):
                raise ValueError("Scan %i produced an all-NaN map. This almost always indicates "%scan+\
                                 "incorrect map bounds or coordinate system. Aborting reduction.")

            info('Coadding...')
            if ms and m:
                ms = mapsumfast([ms,m])  
        
            elif not ms:
                ms = copy.deepcopy(m)

            if doPlot:
                # SNR map creation
                snrMap = copy.deepcopy(ms)  # Signal
                # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)
                snrMap.Data = np.where(snrMap.Weight > 0.0, snrMap.Data * np.sqrt(snrMap.Weight), np.NaN)
                # plotting
                snrMap.display(aspect=1,limitsZ=[-3, +10])
                del snrMap  # free memory

            del m  # free memory

            # Space between co-adding scans
            print('')

        # ==========================================================
        # ITERATION COMPLETE, NO SMOOTHING AT ALL UP TO HERE IN "ms"
        # ==========================================================
        del mymodel  # free memory

        # Now create final iter maps and FITS.
        # First, smooth co-added if required:
        if smoothby_deg > 0.0:
            info('Smoothing co-added map for iteration %i by %.1f"...'%(iter, smoothby_arcsec))
            nativebeam = ms.BeamSize
            auxsmoothby(ms, smoothby_deg)
            newbeam = ms.BeamSize
            print('Original beam: %.3f"     New beam: %.3f"'%(nativebeam*3600, newbeam*3600))

        # RMS map creation
        rmsMap = copy.deepcopy(ms)  # Signal
        rmsMap.Data = np.where(rmsMap.Weight > 0.0, 1.0 / np.sqrt(rmsMap.Weight), np.NaN)  # Noise = 1/sqrt(weight)

        # SNR map creation
        snrMap = copy.deepcopy(ms)  # Signal
        snrMap.Data = np.where(snrMap.Weight > 0.0, snrMap.Data * np.sqrt(snrMap.Weight), np.NaN)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

        # Compute statistics, let auxwriteFits handle clipping
        mediannoise = np.nanmedian(rmsMap.Data)  # on full map
        if clip != -1:
            minnoise = np.nanmin(rmsMap.Data[rmsMap.Data<clip*mediannoise])
            meannoise = np.nanmean(rmsMap.Data[rmsMap.Data<clip*mediannoise])
        else:
            minnoise = np.nanmin(rmsMap.Data[rmsMap.Data<2*mediannoise])
            meannoise = np.nanmean(rmsMap.Data[rmsMap.Data<2*mediannoise])

        # plotting
        minsnr = min(-5, np.nanpercentile(snrMap.Data[snrMap.Data<0], 90))
        maxsnr = max(5, np.nanpercentile(snrMap.Data[snrMap.Data>0], 90))
        maxabs = max(abs(minsnr), abs(maxsnr))
        snrMap.display(aspect=1,limitsZ=[-maxabs, maxabs])
        if clip != -1:
            rmsMap.display(aspect=1,limitsZ=[0, clip*mediannoise],doContour=1,levels=[clip*mediannoise],overplot=1)
        else:
            rmsMap.display(aspect=1,limitsZ=[0, 2*mediannoise],doContour=1,levels=[2*mediannoise],overplot=1)

        # Save full-iteration map (will be smoothed if smooth > 0.0)
        outname = "ReducedFiles/"+str(myname)+"-coadded-flux-iter"+str(iter)+".data"  # goes into ReducedFiles dir
        ms.dumpMap(outname)

        print("################### Iteration %i finished ###################"%(iter))
        print("minimum noise: %5.1f mJy/b, mean noise: %5.1f mJy/b"%(1000*minnoise,1000*meannoise))
        print("############################################################")

        outname = str(myname)+"-coadded-iter"+str(iter)+".fits" # Goes into current dir.
        auxwriteFits(ms, outfile=outname, overwrite=1, clip=clip)
        
        del ms  # free memory
        del rmsMap  # free memory
        del snrMap  # free memory


print('\n\n\n')
print('############################')
info('Reduction finished.')
print('############################')

# Beam corrections
if correctbeam:
    print('\n\n\n')
    execfile('correctAMKIDbeam.py')