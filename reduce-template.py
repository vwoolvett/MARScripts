import copy as copy

# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================

# --- Source and map parameters ---
source  = 'Name'                # As in observing logs
fe      = 'LFA'                 # Frontend, either 'LFA' or 'HFA'
system  = 'EQ'                  # Coordinate system for map, 'EQ', 'GAL' or 'HO'
center  = [0, 0]                # Center of map in CHOSEN coordinates in deg
sizex   = 1.0                   # Size of map in x direction in deg
sizey   = 1.0                   # Size of map in y direction in deg
padding = 0.5                  # Padding around the map in deg for grid, default is 
                                # about the width of the array.
doPlot  = True                  # Whether to display maps at each scan. If False, only final
                                # coadded map per iteration will be displayed.

# ----- Reduction parameters -----
writeSummary = True             # Whether to write a summary file for each scan with 
                                # noise and area information.
niters       = 2                # Number of iterations to run, 1 to 3 (recommended 2)
clip         = 5.               # Sigma clipping level for masking high noise pixels in
                                # the final coadded map.
flagJumps    = False            # Whether to flag jumps/spikes in the data, recommended
                                # to set to True for weak sources in LFA.

# ----- Scans ------
# If empty, automatically retrieves all scans of source from Obslogs
# NOTE: CURRENTLY NOT FUNCTIONAL, PLEASE MANUALLY INPUT SCAN NUMBERS
scans = []

# Manually exclude bad scans if needed            
badscans = []                   

# ==============================
# ===== END OF USER INUPUT =====
# ==============================









# ===== REDUCTION CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
# variable checks
if fe not in ['LFA', 'HFA']:
    raise ValueError("fe must be either 'LFA' or 'HFA'.")
if system not in ['EQ', 'GAL', 'HO']:
    raise ValueError("system must be either 'EQ', 'GAL', or 'HO'.")
if niters < 1 or niters > 3:
    raise ValueError("niters must be 1, 2, or 3.")

# find scans if not provided
#def findScans(source, fe):
    


if len(scans) == 0:
    scans = findScans(source, fe)

# Define myname variable
myname = str(fe) + "-" + str(source) + "-" + str(system)
if flagJumps:
    myname += "-flagJumps"

# Transform 0-360 deg to -180 to 180 if needed
if center[0] > 180:
    warn('Center X is greater than 180 deg, transforming to -180 to 180 deg...')
    center[0] -= 360

# map bounds in absolute EQ or GAL coordinates in deg
xsize = [center[0] + sizex/2 + padding, center[0] - sizex/2 - padding]
ysize = [center[1] - sizey/2 - padding, center[1] + sizey/2 + padding]

# Remove bad scans from the list of scans to be reduced
for badscan in badscans:
    if badscan in scans:
        scans.remove(badscan)


print(
'''
=====================
Reduction parameters:
=====================
Source:             %s
Frontend:           %s
Coordinate system:  %s
Map center:         %s, %s deg
Map size (x,y):     %s, %s deg
Map Boundaries:     %s, %s deg in x; %s, %s deg in y
Padding:            %s deg
Iterations:         %i
Sigmaclip level:    %s
Flag jumps:         %s

'''%(source, fe, system, center[0], center[1], sizex, sizey, 
     xsize[0], xsize[1], ysize[0], ysize[1], padding, niters, clip, flagJumps))

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

# Beginning of reduction loop
for iter in range(1, niters+1):
    print("####################################################################")
    print("####################### Iteration %i starting #######################"%(iter))
    print("####################################################################")

    if iter == 1:
        mymodel=None
        subtract = False
    else:
        mymodel=str(myname)+"-coadded-flux-iter"+str(iter-1)+".data"
        m=restoreFile(mymodel)
        if iter == 2:
            subtract = False
            mymodel=createSourceModel(m,highcut=5.5,lowcut=2.5,sm=8.,mtype='snr',clip=3)
            
        if iter == 3:
            subtract = True
            mymodel=createSourceModel(m,highcut=5.5,lowcut=2.5,sm=0.,mtype='flux',clip=3)
            
    
    ms=None

    for i,scan in enumerate(scans):
        scanname = "ReducedFiles/"+str(myname)+"-"+str(scan)+"-iter"+str(iter)+".data"
        globlist = glob(scanname)

        m = None
        if len(globlist) ==  0:
            info('Reducing scan %s (iteration %i)...'%(scan, iter))
            redweak(scan,fe=fe,size=-1,model=mymodel,subtract=subtract,doPlot=doPlot,extremeFilter=False,writeSummary=writeSummary,flagJumps=flagJumps)
            #if scan == 22919: #flagging example to flag a certain time range in a map (seconds from the beining of the scan) 
            #    flagMJD(above=1430,below=1600,flag=2)
            mapping(oversamp=4,system=system,sizeX=xsize,sizeY=ysize,noPlot=noPlot)
            data.Map.dumpMap(scanname)
            m = restoreFile(scanname)
            m.smoothBy(8./3600.)
            if writeSummary:
                rmsMap = copy.deepcopy(m)
                rmsMap.Data =  (rmsMap.Data*0.0+1.0)/np.sqrt(rmsMap.Weight)
                minnoise = np.nanmin(rmsMap.Data)
                mask = np.where(rmsMap.Data > 5*minnoise)
                rmsMap.Data[mask] = np.NaN
                pixelsize = np.abs(rmsMap.WCS['CDELT2'])
                nrpix = np.sum(~np.isnan(rmsMap.Data))
                area = nrpix*pixelsize**2
                noise = np.nanmedian(rmsMap.Data)
                outname = "%s-%s-%i_summary.txt"%(fe,data.ScanParam.Object,data.ScanParam.ScanNum)
                f=open(outname,'r')
                lines=f.readlines()
                f.close()
                myline=lines[0].replace("\n","")
                myline=myline+",{:.3f},{:.4f}\n".format(area,noise)
                f=open(outname,'w')
                f.write(myline)
                f.close()

                # VWO: make iteration-specific and move to
                # Summmaries directory to clean up current dir.
                newdir = "Summaries/"
                newname = myname+"-"+str(scan)+"-iter"+str(iter)+"_summary.txt"
                os.rename(outname, newdir+newname)
            
        else:
            info('Reduction for scan %i (iteration %i) found'%(scan, iter))
            m = restoreFile(scanname)
            m.smoothBy(8./3600.)

        if np.all(np.isnan(m.Data)):
            warn('Map data is all NaNs! Ensure that the map bounds are correct.')
            break
        
        if ms and m:
            ms = mapsumfast([ms,m])
        elif not ms:
            ms = copy.deepcopy(m)

        if doPlot:
            # SNR map creation
            snrMap = copy.deepcopy(ms)  # Signal
            snrMap.Data *= np.sqrt(snrMap.Weight)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

            # rescaling SNR map by beam size
            a = snrMap.computeRms()
            scale = snrMap.RmsBeam
            snrMap.Data /= np.array(scale,'f')

            # plotting
            snrMap.display(aspect=1,limitsZ=[-4,12])

    # Iteration complete, now create final coadded maps and display final SNR map with optional clipping of high noise pixels
    snrMap = copy.deepcopy(ms)  # Signal
    snrMap.Data *= np.sqrt(snrMap.Weight)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

    # rescaling SNR map by beam size
    # snrMap.iterativeSigmaClip(below=-4,above=4)        
    a = snrMap.computeRms()
    scale = snrMap.RmsBeam
    snrMap.Data /= np.array(scale,'f')

    # creating rms map
    rmsMap = copy.deepcopy(ms)  # Signal
    rmsMap.Data =  (rmsMap.Data*0.0+1.0)/np.sqrt(rmsMap.Weight) # Noise = 1/sqrt(weight)

    # Rescaling rms map by beam size
    rmsMap.Data *= np.array(scale,'f')

    # clipping high noise pixels if clip > 0
    minnoise = np.nanmin(rmsMap.Data)
    meannoise = np.nanmean(rmsMap.Data)
    mediannoise = np.nanmedian(rmsMap.Data)
    
    if clip > 0:
        mask=np.where(rmsMap.Data > clip * mediannoise)
        ms.Data[mask] = np.NaN
        snrMap.Data[mask] = np.NaN
        rmsMap.Data[mask] = np.NaN

    snrMap.display(aspect=0,limitsZ=[-4,12])
    rmsMap.display(aspect=0,limitsZ=[0, 2*mediannoise],doContour=1,levels=[mediannoise],overplot=1)

    outname = "ReducedFiles/"+str(myname)+"-coadded-flux-iter"+str(iter)+".data"  # goes into ReducedFiles dir
    ms.dumpMap(outname)

    print("################### Iteration %i finished ###################"%(iter))
    print("minimum noise: %5.1f mJy/b, mean noise: %5.1f mJy/b"%(1000*minnoise,1000*meannoise))
    print("############################################################")

    outname = str(myname)+"-coadded-iter"+str(iter)+".fits" # Goes into current dir.
    writeFits2(ms, outfile=outname, overwrite=1)