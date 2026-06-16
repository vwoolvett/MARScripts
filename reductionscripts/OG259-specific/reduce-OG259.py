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
writeSummary = True         # Write summary of reductions or not
niters       = 1            # Number of iterations to run, 1 to 3 (recommended: 2 + PLANCK data)
clip         = 5.           # Sigma clipping level for masking high noise pixels
flagJumps    = True         # Flag jumps/spikes in the data:
                            # recommended to set to True for 'weak' sources in LFA
smoothby_arcsec = 8.        # Default 8. arcsec

# ----- Scans ------
# If scans is empty, automatically retrieves all scans of the source
# specified above from the obslogs directory below
scans = []
obslogsdir = '~/obslogs'  # at MPIfR: '/apex-archive/obslogs/M-PROJECT.CODE-IN-CAPS/obslogs'

# Manually exclude bad scans if needed            
badscans = [] #[27991, 27979, 28217, 28498]

# ==============================
# ===== END OF USER INUPUT =====
# ==============================














# ===== REDUCTION CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
# NOTE: VWO: BoA smoothBy smooths weights with the kernel, but weights are 1 / rms^2 which is
# a non-linear scale. All smoothing should be done as:
# Sky_smoothed = Kernel * Sky
# Weights_smoothed = 1 / Variance_smoothed | with Variance_smoothed = Kernel^2 * Variance
# Coverage_smoothed = Kernel * Coverage

# NOTE 2: redweak still uses smoothBy! removed smoothing in redweak
import warnings
import copy as copy
import BoaMapping as BOAMAP
from mars.fortran import fMap

# variable checks
if fe not in ['LFA', 'HFA']:
    raise ValueError("fe must be either 'LFA' or 'HFA'.")
if system not in ['EQ', 'GAL', 'HO']:
    raise ValueError("system must be either 'EQ', 'GAL', or 'HO'.")
if niters < 1 or niters > 3:
    raise ValueError("niters must be 1, 2, or 3.")

# define the good functions :)
def findSciTargetScans(source, obslogsdir):
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
                    if key == 'Scan status':
                        message+=(line[4:-6].ljust(12) + ' | ')
                    if key == 'Scan type':
                        message+=(line[4:-6].ljust(12) + ' | ')
                    #if key=='Comment':
                    #    message+=(line[4:-6].ljust(20))              
                start = False

                if source in message:
                    if 'OK' in message and 'MAP' in message:
                        message += 'SCAN CONSIDERED'
                        print(message)
                        scanlist.append(scan)
                    else:
                        message += 'SCAN DISCARDED'
                        print(message)
    # If nothing was found, break script
    if len(scanlist) == 0:
        raise ValueError('No scans of source %s found in ObsLogs directory: %s!'%(source, obslogsdir))
    print('')
    info("Number of 'MAP' scans on science target %s: %i"%(source, len(scanlist)))
    print('')
    return scanlist



def auxsmoothby(m, Size):
    '''
    BoA-like smoothing but with correct variance propagation.

    - Data: convolved with K
    - Weight: propagated via variance (K^2)
    - Coverage: convolved with K (same as BoA)
    '''
    # Build kernel
    pixsize = abs(m.WCS['CDELT2'])
    K = BOAMAP.Kernel(pixsize, Size).Data.astype(float)
    K_norm = K / np.sum(K)

    # Smooth INTENSITY (same as BoA)
    I1 = fMap.ksmooth(m.Data, K_norm)

    # Correct variance propagation for weights
    #    V' = K^2 * V
    V0 = np.where(m.Weight > 0.0, 1.0 / m.Weight, np.NaN)
    V1 = fMap.ksmooth(V0, K_norm**2)
    W1 = np.where(V1 > 0.0, 1.0 / V1, 0.0)

    # Smooth COVERAGE (same as BoA)
    C1 = fMap.ksmooth(m.Coverage, K_norm)
    
    # new scale per beam for Jy/beam units
    newbeam = np.sqrt(m.BeamSize**2 + Size**2)
    scale = (newbeam**2 / m.BeamSize**2)

    # Convert weight map from pixel-based to beam-based
    beam_factor2 = np.sum(K_norm**2)
    W1 /= beam_factor2

    # Update map with correct Jy/beam scale
    m.Data = I1 * scale
    m.Weight = W1 / scale**2
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

        
        localMap = 0  # free memory
        snrMap = 0
        rmsMap = 0
        dataset = 0



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
    
# sort scans
scans.sort()

# Remove bad scans from the list of scans to be reduced
for badscan in badscans:
    if badscan in scans:
        scans.remove(badscan)

# Define standardized "myname" variable for output files
myname = str(fe) + "-" + str(source) + "-" + str(system)
if flagJumps:
    myname += "-flagJumps"

# Create map bounds, for EQ or GAL first
ysize = [center[1] - sizey/2 - padding, center[1] + sizey/2 + padding]
xsize = [center[0] + sizex/2 + padding, center[0] - sizex/2 - padding]  # Bigger number first
if system =='HO':
    # Invert X boundaries back to normal
    xsize = [xsize[1], xsize[0]]

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
Smoothing:          %s arcsec'''%(source, fe, system, center[0], center[1], sizex, sizey, padding,
     xsize[0], xsize[1], ysize[0], ysize[1], niters, clip, flagJumps,
     smoothby_arcsec))

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

            # Check if reduction exists
            if len(globlist) ==  0:
                print('')
                print('')
                print('')
                info('Reducing scan %s (iteration %i)...'%(scan, iter))

                # Reduce it
                redweak(scan,fe=fe,size=-1,model=mymodel,subtract=subtract,doPlot=doPlot,extremeFilter=False,writeSummary=writeSummary,flagJumps=flagJumps)
            
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

                # Create map
                mapping(oversamp=4,system=system,sizeX=xsize,sizeY=ysize,noPlot=noPlot)
            
                # Save unsmoothed map
                data.Map.dumpMap(scanname)

                # Create BoA map
                m = restoreFile(scanname)

                # Smooth if necessary
                if smoothby_deg > 0.0:
                    auxsmoothby(m, smoothby_deg)

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
                info('Reduction for scan %i (iteration %i) found'%(scan, iter))
                m = restoreFile(scanname)

                # Smooth if necessary
                if smoothby_deg > 0.0:
                    auxsmoothby(m, smoothby_deg)


            if np.all(np.isnan(m.Data)):
                warn('Map data is all NaNs! Ensure that the map bounds are correct.')
                break

            info('Coadding...')
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
                snrMap.display(aspect=1,limitsZ=[-4,12])

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