import copy as copy

# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================

# --- Source and map parameters ---
source  = 'OG259'                # As in APECS/ObsLogs
fe      = 'LFA'                 # Frontend, either 'LFA' or 'HFA'
system  = 'GAL'                  # Coordinate system for map, either 'EQ' or 'GAL'
center  = [259.3, -1.4]            # Center of map in CHOSEN absolute coordinates in deg
xsize   = 1.5                   # Size of map in x direction in DEG
ysize   = 1.5                   # Size of map in y direction in DEG
padding = 0.25                  # Padding around the map in DEG for grid, default is 
                                # about the width of the array.
doPlot  = False                 # Whether to display maps at each iteration

# ----- Reduction parameters -----
writeSummary = False             # Whether to write a summary file for each scan with 
                                # noise and area information.
niters       = 2                # Number of iterations to run, 1 to 3 (recommended 2)
clip         = 3.               # Sigma clipping level for masking high noise pixels in
                                # the final coadded map.
flagJumps    = False            # Whether to flag jumps/spikes in the data, recommended
                                # to set to True for weak sources in LFA.

# ----- Scans ------
scans    = [27974]#,27975,27979,27980,27990,28212,28213,28217,28218,28231,28232,28235,28493,28494,28498,28499,28516,28517,28775,28776]                   # 'Auto' or list of scans to reduce
                                # NOTE: If using 'Auto', make sure to set the correct
                                # source name and frontend above
badscans = []#[27979,28213,28217,28498]                   # Manually exclude bad scans if needed

# ==============================
# ===== END OF USER INUPUT =====
# ==============================









# ===== REDUCTION CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
# Define myname variable
myname = str(fe) + "-" + str(source) + "-" + str(system)
if flagJumps:
    myname += "-flagJumps"

# map bounds in absolute EQ or GAL coordinates in deg
xsize = [center[0] - xsize/2 - padding, center[0] + xsize/2 + padding]
ysize = [center[1] - ysize/2 - padding, center[1] + ysize/2 + padding]

# Remove bad scans from the list of scans to be reduced
for badscan in badscans:
    if badscan in scans:
        scans.remove(badscan)

# Set noPlot
if not doPlot:
    noPlot = True
else:
    noPlot = False

# Create directory for reduced files if it doesn't exist
if os.path.exists("ReducedFiles") == False:
    os.makedirs("ReducedFiles")

# Beginning of reduction loop
for iter in range(1,3):
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
        info('Processing scan %s...'%(scan))
        globlist=glob(scanname)

        m = None
        if len(globlist) ==  0:
            info('Reducing scan %s...'%(scan))
            redweak(scan,fe='LFA',size=-1,model=mymodel,subtract=subtract,doPlot=doPlot,extremeFilter=False,writeSummary=writeSummary,flagJumps=flagJumps)
            mapping(oversamp=4,system=system,sizeX=xsize,sizeY=ysize,noPlot=noPlot)
            data.Map.dumpMap(scanname)
            m=restoreFile(scanname)
            m.smoothBy(8./3600.)
            if writeSummary:
                rmsMap = copy.deepcopy(m)
                rmsMap.Data =  (rmsMap.Data*0.0+1.0)/np.sqrt(rmsMap.Weight)
                minnoise=np.nanmin(rmsMap.Data)
                mask=np.where(rmsMap.Data > 5*minnoise)
                rmsMap.Data[mask] = np.NaN
                pixelsize=np.abs(rmsMap.WCS['CDELT2'])
                nrpix = np.sum(~np.isnan(rmsMap.Data))
                area = nrpix*pixelsize**2
                noise=np.nanmedian(rmsMap.Data)
                outname="%s-%s-%i_summary.txt"%(fe,data.ScanParam.Object,data.ScanParam.ScanNum)
                f=open(outname,'r')
                lines=f.readlines()
                f.close()
                myline=lines[0].replace("\n","")
                myline=myline+",{:.3f},{:.4f}\n".format(area,noise)
                f=open(outname,'w')
                f.write(myline)
                f.close()
            
        else:
            info('Reduction found at:')
            print(scanname)
            m=restoreFile(scanname)
            m.smoothBy(8./3600.)
            

        if ms and m:
            ms = mapsumfast([ms,m])
        elif not ms:
            ms = copy.deepcopy(m)
    
        
    
        rmsMap = copy.deepcopy(ms)
        snrMap = copy.deepcopy(ms)
        

        snrMap.Data *= np.sqrt(snrMap.Weight)

        tmp=copy.deepcopy(snrMap)
        a = tmp.computeRms()
        scale=tmp.RmsBeam
        snrMap.Data /= np.array(scale,'f')

    
        if doPlot:
            snrMap.display(aspect=1,limitsZ=[-4,12])

    rmsMap = copy.deepcopy(ms)
    snrMap = copy.deepcopy(ms)
        

    snrMap.Data *= np.sqrt(snrMap.Weight)

    tmp=copy.deepcopy(snrMap)
    tmp.iterativeSigmaClip(below=-4,above=4)        
    a = tmp.computeRms()
    scale=tmp.RmsBeam
    snrMap.Data /= np.array(scale,'f')

    rmsMap.Data =  (rmsMap.Data*0.0+1.0)/np.sqrt(rmsMap.Weight)
    rmsMap.Data *= np.array(scale,'f')

    minnoise=np.nanmedian(rmsMap.Data)
    meannoise=np.nanmedian(rmsMap.Data)
    
    if clip > 0:
        mask=np.where(rmsMap.Data > 5*minnoise)
        ms.Data[mask] = np.NaN
        snrMap.Data[mask] = np.NaN
        rmsMap.Data[mask] = np.NaN

    snrMap.display(aspect=0,limitsZ=[-4,12])
    outname=str(myname)+"-coadded-flux-iter"+str(iter)+".data"
    ms.dumpMap(outname)


    print("###################iteration %i ##########################"%(iter))
    print("minimum noise: %5.1f mJy/b, mean noise: %5.1f mJy/b"%(1000*minnoise,1000*meannoise))
    print("#####################################################")



    outname=str(myname)+"-coadded-iter"+str(iter)+".fits"
    writeFits2(ms,outfile=outname,overwrite=1)

