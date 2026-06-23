# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================

# --- Source and map parameters ---
source  = 'OG259'           # As in observing logs
fe      = 'LFA'             # Frontend, either 'LFA' or 'HFA'
system  = 'GAL'             # Coordinate system for map, 'EQ', 'GAL' or 'HO'
center  = [259.3, -1.4]     # Center of map in CHOSEN COORDINATES in deg
sizex   = 1.5               # Size of map in deg for X direction
sizey   = 1.5               # Size of map in deg for Y direction
padding = 0.5               # Padding around the map in deg for grid (default ~ 2x array)
doPlot  = True              # Display maps at each scan. If False, only final
                            # coadded map per iteration will be displayed.

# ----- Reduction parameters -----
# SUGGESTED: run with niters=1, figure out bad scans using showMaps.py and then
# run with niters=2 ignoring bad scans
writeSummary = True         # Write summary of reductions or not
niters       = 1            # Number of iterations to run, 1 to 3 (recommended: 2 + PLANCK data)
clip         = 5.           # Sigma clipping level for masking high noise pixels
flagJumps    = True         # Flag jumps/spikes in the data:
                            # recommended to set to True for 'weak' sources in LFA
smoothby_arcsec = 8.        # Default 8. arcsec
correctbeam     = True      # Whether to correct beam bookkeeping in final iteration maps

# ----- Scans ------
# If scans is empty, automatically retrieves all scans of the source
# specified above from the obslogs directory below
scans = []  # 27974 is first
obslogsdir = '~/obslogs'  # at MPIfR: '/apex-archive/obslogs/M-PROJECT.CODE-IN-CAPS/obslogs'

# Manually exclude bad scans if needed            
badscans = [27979, 27991, 28217, 28498, 33735] # 33735 not readable 

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
def findSciTargetScans(source, obslogsdir, verbose=False):
    scanlist = []
    files = os.listdir(obslogsdir)
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
                        scan=int(line[4:-6])
                        message+=(line[4:-6].ljust(6) + ' | ')
                    if key == 'Source':
                        message+=(line[4:-6].ljust(12) + ' | ')       
                    if key == 'Scan type':
                        message+=(line[4:-6].ljust(12) + ' | ')        
                    if key == 'Scan status':
                        message+=(line[4:-6].ljust(12) + ' | ') 
                start = False

                if source in message.split('|')[1]:
                    if 'MAP' in message.split('|')[2] and 'OK' in message.split('|')[3]:
                        message += 'SCAN CONSIDERED'
                        scanlist.append(scan)
                    else:
                        message += 'SCAN DISCARDED'
                    if verbose:
                        print(message)
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
            mediannoise = np.nanmedian(rmsMap.Data)
            mask = np.where(rmsMap.Data > clip * mediannoise)
            localMap.Data[mask] = np.NaN
            rmsMap.Data[mask] = np.NaN
            snrMap.Data[mask] = np.NaN
 
        #write FLux plane                                                            
        localMap._Image__writeImage(dataset, "Intensity", intensityUnit=intensityUnit)
        #write RMS plane
        rmsMap._Image__writeImage(dataset, "Intensity", intensityUnit=intensityUnit+" (RMS)")
        #write SNR plane
        snrMap._Image__writeImage(dataset, "Intensity", intensityUnit='SNR')
        dataset.close()
            
    except Exception, data:
        print('Could not write data to file %s: %s' % (outfile, data))
        return
    
    # free memory
    localMap = 0
    snrMap = 0
    rmsMap = 0
    dataset = 0



# variable checks
if fe not in ['LFA', 'HFA']:
    raise ValueError("fe must be either 'LFA' or 'HFA'.")
if system not in ['EQ', 'GAL', 'HO']:
    raise ValueError("system must be either 'EQ', 'GAL', or 'HO'.")
if niters < 1 or niters > 3:
    raise ValueError("niters must be 1, 2, or 3.")
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
if biggerX > 180:
    biggerX -= 360
    sysreframe = True
if biggerX < -180:
    biggerX += 360
    sysreframe = True
if smallerX > 180:
    smallerX -=360
    sysreframe = True
if smallerX < -180:
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

  
        # Initialize co-added map
        ms = None

        for i,scan in enumerate(scans):
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
                if scan == 28517:
                    flagC(3353, flag=2)

                # Create map in chosen system and chosen box
                # where pixsize = BEAM_FWHM / oversamp
                mapping(oversamp=4,system=system,sizeX=xsize,sizeY=ysize,limitsZ=[-0.8,1.5],noPlot=noPlot)
                # NOTE: this has a smooth parameter, but is default 0
                # NOTE 2: data.Map.BeamSize is taken from data.BolometerArray.BeamSize
                # NOTE 3: data.BolometerArray.BeamSize is just 1.22 * lambda / D * 180/pi, not from beammap!
            
                # Save unsmoothed map, "native" resolution (m.BeamSize = data.BolometerArray.BeamSize)
                data.Map.dumpMap(scanname)

                # Assing BoA map to variable m
                m = restoreFile(scanname)

                # Smooth if necessary with the correct handling of Weight map using our new function.
                if smoothby_deg > 0.0:
                    info('Smoothing map before co-adding...')
                    nativebeam = m.BeamSize
                    auxsmoothby(m, smoothby_deg)
                    newbeam = m.BeamSize
                    print('Unsmoothed beam: %.3f "'%(nativebeam*3600))
                    print('Smoothing by:    %.3f "'%(smoothby_deg*3600))
                    print('New beam:        %.3f "'%(newbeam*3600))

                # For this scan, add non-noisy area and median noise info to summary
                if writeSummary:
                    # Create smoothed noise map
                    rmsArray = np.where(m.Weight > 0.0, 1.0 / np.sqrt(m.Weight), np.NaN)

                    # Statistics and write
                    minnoise = np.nanmin(rmsArray)
                    mask = (rmsArray > 5*minnoise)
                    rmsArray[mask] = np.NaN
                    pixelsize = np.abs(m.WCS['CDELT2'])  # taken from smoothed map
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
            
            else:
                # Retrieve BoA map
                info('Reduction for scan %i (iteration %i) found. Loading...'%(scan, iter))
                m = restoreFile(scanname)

                # Smooth if necessary
                if smoothby_deg > 0.0:
                    info('Smoothing map before co-adding...')
                    nativebeam = m.BeamSize
                    auxsmoothby(m, smoothby_deg)
                    newbeam = m.BeamSize
                    print('Unsmoothed beam: %.3f "'%(nativebeam*3600))
                    print('Smoothing by:    %.3f "'%(smoothby_deg*3600))
                    print('New beam:        %.3f "'%(newbeam*3600))


            if np.all(np.isnan(m.Data)):
                warn('Map data is all NaNs! Ensure that the map bounds are correct.')
                break

            info('Coadding...')
            print('')
            if ms and m:
                # both are smoothed
                ms = mapsumfast([ms,m])  
        
            elif not ms:
                # is already smoothed
                ms = copy.deepcopy(m)

            if doPlot:
                # SNR map creation
                snrMap = copy.deepcopy(ms)  # Signal
                snrMap.Data = np.where(snrMap.Weight > 0.0, snrMap.Data * np.sqrt(snrMap.Weight), np.NaN )  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

                # plotting
                snrMap.display(aspect=1,limitsZ=[-4, 12])

        # Iteration complete, now create final coadded maps and display
        # final SNR map + noise contours with optional clipping of high noise pixels.
        # RMS map creation
        rmsMap = copy.deepcopy(ms)  # Signal
        rmsMap.Data = np.where(rmsMap.Weight > 0.0, 1.0 / np.sqrt(rmsMap.Weight), np.NaN)  # Noise = 1/sqrt(weight)

        # SNR map creation
        snrMap = copy.deepcopy(ms)  # Signal
        snrMap.Data = np.where(snrMap.Weight > 0.0, snrMap.Data * np.sqrt(snrMap.Weight), np.NaN)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

        # clipping high noise pixels if clip > 0
        mediannoise = np.nanmedian(rmsMap.Data)
    
        if clip > 0:
            mask = np.where(rmsMap.Data > clip * mediannoise)
            ms.Data[mask] = np.NaN
            snrMap.Data[mask] = np.NaN
            rmsMap.Data[mask] = np.NaN

        minnoise = np.nanmin(rmsMap.Data[rmsMap.Data<1.5*mediannoise])
        meannoise = np.nanmean(rmsMap.Data[rmsMap.Data<1.5*mediannoise])

        # plotting (these are already smoothed if used)
        snrMap.display(aspect=1,limitsZ=[-4,12])
        rmsMap.display(aspect=1,limitsZ=[0, 1.5*mediannoise],doContour=1,levels=[1.5*mediannoise],overplot=1)

        # Save smoothed (if used) full-iteration map
        outname = "ReducedFiles/"+str(myname)+"-coadded-flux-iter"+str(iter)+".data"  # goes into ReducedFiles dir
        ms.dumpMap(outname)

        print("################### Iteration %i finished ###################"%(iter))
        print("minimum noise: %5.1f mJy/b, mean noise: %5.1f mJy/b"%(1000*minnoise,1000*meannoise))
        print("############################################################")

        outname = str(myname)+"-coadded-iter"+str(iter)+".fits" # Goes into current dir.
        auxwriteFits(ms, outfile=outname, overwrite=1)

print('\n\n\n')
print('############################')
info('Reduction finished.')
print('############################')

# Beam corrections
if correctbeam:
    print('\n\n\n')
    execfile('correctAMKIDbeam.py')