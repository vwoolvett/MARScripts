# =============================================================================
# ========================== BEGINNING OF USER INPUT ==========================
# =============================================================================
# Last edited by: VWO 10.08.2026

# ------------------------- OBSERVER or PI mode -------------------------------
observer    = True          # *NOTE1* True or False
projcode    = 'auto'        # Project code (automatic ONLY AT APEX)

# --- Source and map parameters ---
source      = 'SrcName'     # As in observing logs and source catalog
fe          = 'LFA'         # Frontend, either 'LFA' or 'HFA'
system      = 'EQ'          # Coordinate system for map, 'EQ', 'GAL' or 'HO'
center      = [0, 0]        # Center of map in DEG for **CHOSEN COORDINATES**
sizex       = 1.0           # Size of map in DEG for X direction
sizey       = 1.0           # Size of map in DEG for Y direction
padding     = 0.3           # Padding around the map in DEG for the grid
smoothing   = 'default'     # *NOTE2* By how much to smooth final iter. maps

# ------------------------- Reduction parameters ------------------------------
badscans        = []        # Manually exclude bad scans if needed
niters          = 1         # Number of iters., 1 to 3 (recomm.: 3 + PLANCK)
clip            = -1        # *NOTE3* Sigma clipping level (-1 or >=1.5)
flagJumps       = False     # Flag spikes: recomm. True if spikes are present
showsig         = True      # Show example signals after each calibration step
doPlot          = True      # Display coadded map after each scan or final only
writefits       = True      # Write FITS of final iteration maps.
correctbeam     = True      # Correct AMKID beam to nominal: recomm. True
writeSummary    = False     # Write summary of reductions (mostly debugging)

# ------------------------- Scans (automatic at APEX) -------------------------
scans       = []            # *NOTE4*
obslogsdir  = 'default'     # If not at MPIfR or APEX, input manually!
verbose     = False         # print ObsLog scan selection criteria (debugging)

# ------------------------- NOTES ---------------------------------------------
# *NOTE1*   If you are an observer, leave as True to assess AMKID performance
#           and/or calib at scan reduction. Prompts upon running script as 
#           observer give more information on what to do. PIs should set to 
#           observer = False and re-reduce Iteration 1 with an empty "badscans"
#           list variable in this reduction script; then assess with the script
#           "showMaps.py" which scans to actually discard.  For exceptional
#           cases, additional flagging is needed. Consult with Axel or Vicente
#           if this applies to you.

# *NOTE2*   The arrays have beam sizes of 16.7 (LFA) and 7.5 arcsec (HFA).
#           The default smoothings are 8.0 arcsec (LFA) and 3.6 arcsec (HFA).
#           If a target beam size is requested by the PI, consider:
#           smoothing^2 = targetbeam^2 - nativebeam^2

# *NOTE3*   The image is masked where noisemap > clip*mediannoise (clip>=1.5), 
#           or else (clip==-1) no clipping.

# *NOTE4*   If scans is empty, attempts to automatically retrieve all scans of
#           the source specified above from the specified obslogs directory.
#           If not reducing at APEX or MPIfR, you must manually input the
#           obslogs directory and likely set the (raw)data input directory via 
#           the "indir('directory')" function in mars.

# =============================================================================
# ============================= END OF USER INPUT =============================
# =============================================================================








# =============================================================================
# === REDUCTION CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING ====
# =============================================================================
def findSciTargetScans(source, obslogsdir, fe, verbose=False):
    assert fe=='LFA' or fe=='HFA', 'fe must be LFA or HFA'
    # no HFA-only mode, so:
    FeBedict = {'LFA': 'AMKID870-AMKID870BE', 'HFA':'AMKID350-AMKID350BE'}
    febe = FeBedict[fe]

    scanlist = []
    files = os.listdir(obslogsdir)
    c=0
    for file in files:
        fullfilename = obslogsdir + file if obslogsdir[-1]=='/'\
                                         else obslogsdir + '/' + file
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
                    if key == 'Scan':                               # 0
                        scan_int = int(line[4:-6])
                        scan = (line[4:-6].ljust(6) + ' | ')

                    if key == 'Source':                             # 1 
                        src = (line[4:-6].ljust(12) + ' | ') 

                    if key == 'Scan type':                          # 2
                        scantype = (line[4:-6].ljust(12) + ' | ')

                    if key == 'Observ. mode':                       # 3
                        mode = (line[4:-6].ljust(12) + ' | ')

                    if key == 'Frontend-backend':                   # 4
                        thisFeBe = line[4:-6]
                        if FeBedict['HFA'] in thisFeBe:
                            thisfe = ('LFA + HFA'.ljust(9) + ' | ')
                        elif FeBedict['LFA'] in thisFeBe:
                            thisfe = ('LFA'.ljust(9) + ' | ')
                        else:
                            thisfe = ('NOT AMKID!'.ljust(9) + ' | ')

                    if key == 'Command':
                        command = (line[4:-6].ljust(12) + ' | ')
                        # only add first 12 characters
                        command = command[:12] + ' | '              # 5

                    if key == 'Scan duration':
                        duration = (line[4:-6].ljust(8) + ' | ')   # 6

                    if key == 'Scan status':
                        status = (line[4:-6].ljust(12) + ' | ')     # 7

                    if key == 'Comment':                            # last
                        comment = (line[4:-6])[0:85]       
                        message += scan + src + scantype + mode + thisfe +\
                                   command + duration + status

                start = False

                if source in src:
                    if  '-999' not in duration:
                        #if 'MAP' in scantype and 'OTF' in mode\
                        #    and fe in thisfe and 'OK' in status:
                        #    message += 'SCAN CONSIDERED'.ljust(15) + ' | '
                        #    message += comment
                        #    scanlist.append(scan_int)
                        if 'calibrate(' not in command\
                            and 'beamscan(' not in command\
                            and 'go(' not in command\
                            and fe in thisfe\
                            and 'OK' in status:
                            message += 'SCAN CONSIDERED'.ljust(15) + ' | '
                            message += comment
                            scanlist.append(scan_int)
                        else:
                            message += 'SCAN DISCARDED'.ljust(15) + ' | '
                            message += comment
                    else:
                        message += 'SCAN ONGOING'.ljust(15) + ' | '
                        message += comment
                    if verbose:
                        print(message)

        if c==0 and len(keys)!=0 and verbose:
            print('============')
            print('OBSLOG KEYS:')
            print('============')
            print(keys)
            c+=1
    scanlist.sort()
    info("Number of valid scans on source %s (%s): %i"\
         %(source, fe, len(scanlist)))
    return scanlist

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

# -------------------------
# --- OBSERVER OVERRIDE ---
# -------------------------
if observer == True:
    doPlot = False          # Observer wants to check map of scan at reduction,
                            # not coadded until that scan.
                            # That is implemented separately below.
    writeSummary = False    # No need.
    niters = 1              # Source model is the most accurate when Iter1scans
                            # are complete. Also, saves time and disk space.
    clip = -1               # full map, no clipping
    writefits = False       # save time, don't clog directory

# Find project code if at APEX
if projcode == 'auto':
    curraccount = os.getenv('USER')
    # project code is separated once with dot and thrice with dash
    if len(curraccount.split('.')) == 2 and len(curraccount.split('-')) == 4:
        projcode = curraccount
        info('Project code extracted from current account: %s'\
             %(projcode))
        del curraccount
    else:
        raise ValueError("STOPPING SCRIPT: project code could not be" +\
                         " extracted from current account: %s."%(curraccount)+\
                         "\nIf you are at APEX, log in with a project"+\
                         " account and re-run script. If you are not at APEX"+\
                         ", you must manually set the project code variable")
else:
    # remove slash if present for some reason
    if projcode[-1] == '/':
        projcode = projcode[:-1]
    # project code is separated once with dot and thrice with dash
    if len(projcode.split('.')) != 2 or len(projcode.split('-')) != 4:
        raise ValueError("STOPPING SCRIPT: project code '%s' is not correct."\
                         %(projcode))

# create lowercase and CAPS version
projcode_low = str.lower(projcode)
projcode_caps = str.upper(projcode)

# find project obslogs folder and set indir if needed
if obslogsdir == 'default':
    obslogsdir = None
    APEX_obslogpath = '/homes/' + projcode_low + '/obslogs/'
    MPIfR_obslogpath = '/apex-archive/obslogs/' + projcode_caps + '/'
    if os.path.exists(APEX_obslogpath):
        # At APEX
        obslogsdir = APEX_obslogpath
        if BoaConfig.inDir != '/apexdata/rawdata/' + projcode_caps + '/':
            indir('/apexdata/rawdata/' + projcode_caps + '/')
    elif os.path.exists(MPIfR_obslogpath):
        # At MPIfR
        obslogsdir = MPIfR_obslogpath
        if BoaConfig.inDir != '/apex-archive/rawdata/' + projcode_caps + '/':
            indir('/apex-archive/rawdata/' + projcode_caps + '/')
        
    else:
        raise ValueError("STOPPING SCRIPT: project obslogs folder could not"\
                         " be extracted from project code: %s"\
                         %projcode + \
                         "\nPlease ensure project code is correct, or " +\
                         "manually set obslogsdir variable in " +\
                         "reduction script to the correct path.")

if len(scans) == 0 and not os.path.exists(obslogsdir):
    raise ValueError("STOPPING SCRIPT: Either enter scans or an existing "+\
                     "obslogs directory...")

# find scans if not provided
if len(scans) == 0 and os.path.exists(obslogsdir):
    info('Retrieving source scan numbers from ObsLogs...')
    scans = findSciTargetScans(source=source, obslogsdir=obslogsdir, fe=fe,
                               verbose=verbose)
    if len(scans) == 0:
        raise ValueError("No scans of source %s (%s) "%(source, fe) +\
                         "found in ObsLogs directory:\n%s"%(obslogsdir))

# sort scans
scans.sort()
badscans.sort()

# Remove bad scans from the list of scans to be reduced
for badscan in badscans:
    if badscan in scans:
        scans.remove(badscan)

# Check removing bads did not leave scans empty
if len(scans) == 0:
    raise ValueError('There are no good scans after removing bad scans list.')

# Create map bounds
info('Creating map boundaries...')
biggerX = center[0] + sizex/2 + padding
smallerX = center[0] - sizex/2 - padding
biggerY = center[1] + sizey/2 + padding
smallerY = center[1] - sizey/2 - padding

# These can't happen
if biggerY > 90:
    raise ValueError('STOPPING SCRIPT: The upper border of the map has Y '+\
                     'coordinate > +90 degrees! (comment this if intended)')
if smallerY < -90:
    raise ValueError('STOPPING SCRIPT: The lower border of the map has Y '+\
                     'coordinate < -90 degrees! (comment this if intended)')

# Check X reframing.
# Example with an X width of 10 deg:
# Case 1: left = 150, right = 140 is left untouched
# Case 2: left = 185, right = 175 
#           -> frame was 0:360, now left = -175, right = 175
# Case 3: left = 200, right = 190
#           -> frame was 0:360, now left = -160, right = -170
# Case 4 : same as before but one of the boundaries ended up < -180: add 360
sysreframe = False
if biggerX > 180 and system != 'EQ':
    biggerX -= 360
    sysreframe = True
if biggerX < -180 and system != 'EQ':
    biggerX += 360
    sysreframe = True
if smallerX > 180 and system != 'EQ':
    smallerX -=360
    sysreframe = True
if smallerX < -180 and system != 'EQ':
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

# Define standardized "myname" variable for output files
myname = str(fe) + "-" + str(source) + "-" + str(system)
if flagJumps:
    myname += "-flagJumps"

# Create directory for reduced files if it doesn't exist
if os.path.exists("ReducedFiles") == False:
    os.makedirs("ReducedFiles")

# Create directory for reduced files if it doesn't exist
if writeSummary and os.path.exists("Summaries") == False:
    os.makedirs("Summaries")

# Create directory for FITS files of final iteration maps
if writefits and os.path.exists("FITSfiles") == False:
    os.makedirs("FITSfiles")

if smoothing == 'default':
    # Default is 23% of beam in area
    if fe == 'LFA':
        smoothby_arcsec = 8.0  # sqrt(0.23)*16.7"
    else:
        smoothby_arcsec = 3.6  # sqrt(0.23)*7.5"
else:
    smoothby_arcsec = smoothing

# smoothby to deg
smoothby_deg = smoothby_arcsec / 3600.

# initialize MJD list for all scans
mymjdrefs = []


# SUMMARY PRINT
print('')
print('''\
=====================
Reduction parameters:
=====================
Observer:           %s
Project code:       %s
Source:             %s
Frontend:           %s
Coordinate system:  %s
Map center:         %.5f, %.5f deg
Map size (x,y):     %.3f, %.3f deg
Padding:            %.3f deg
Map Boundaries:     X: %.3f, %.3f deg | Y: %.3f, %.3f deg
Correct beam:       %s
Smoothing:          %s
Iterations:         %i
Sigmaclip level:    %s
Flag jumps:         %s
Number of scans     %s (valid)'''%(observer, projcode_caps, source, fe, system,
                           center[0], center[1], sizex, sizey, padding,
                           xsize[0], xsize[1], ysize[0], ysize[1],
                           correctbeam,
                           '%.1f arcsec (default)'%(smoothby_arcsec) \
                           if smoothing=='default' \
                           else '%.1f arcsec'%(smoothby_arcsec),
                           niters,
                           clip if clip != -1 else 'No clipping', flagJumps,
                           
                           len(scans)))
if len(badscans) > 0:
    info('Bad scans removed:')
    print('         %s'%badscans)


# ===========================
# Beginning of reduction loop
# ===========================
# impossible to keep PEP8 from here onward...
if True:  # Just to indent
    for iter in range(1, niters+1):
        print('')
        print("#####################################################################")
        print("####################### Iteration %i starting ########################"%(iter))
        print("#####################################################################")

        if iter == 1:
            # First iteration -- no model
            mymodel = None
            subtract = False
        else:
            # retrieve last iteration map
            mymodel = "ReducedFiles/" + str(myname) + "-coadded-flux-iter" + str(iter-1) + ".data"
            coadded = restoreFile(mymodel)

            # Only create source model if a scan in scans list is missing reduction in this iteration
            scanreduced = []
            for scan in scans:
                scanname = "ReducedFiles/"+str(myname)+"-"+str(scan)+"-iter"+str(iter)+".data"
                isreduced = True if len(glob(scanname))!=0 else False
                scanreduced.append(isreduced)

            scanreduced = np.array(scanreduced)

            if iter == 2:
                subtract = False
                if np.any(scanreduced==False):
                    mymodel = createSourceModel(coadded, highcut=5.5, lowcut=2.5, sm=0., mtype='snr', clip=3)
                else:
                    mymodel = None
            
            if iter == 3:
                subtract = True
                if np.any(scanreduced==False):
                    mymodel = createSourceModel(coadded, highcut=5.5, lowcut=2.5, sm=0., mtype='flux', clip=3)
                else:
                    mymodel = None

            del coadded  # free memory

  
        # Initialize co-added map
        ms = None
        tint = 0

        for i, scan in enumerate(scans):
            scanname = "ReducedFiles/"+str(myname)+"-"+str(scan)+"-iter"+str(iter)+".data"
            globlist = glob(scanname)

            # Initialize map for this scan
            m = None
            # avoid mapping() call for map N after CTRL+C during reduction->mapping gets previous data.Map (N-1)
            data.Data = None
            data.Map = None

            # Check if reduction does not exist
            if len(globlist) ==  0:
                print('')
                print('')
                info('Reducing scan %i (Iter%i | scan %i/%i)...'%(scan, iter, i+1, len(scans)))

                # Reduce it
                redscience(scan, fsweep=None, fe=fe, src=source, model=mymodel,
                           subtract=subtract, extremeFilter=False,
                           correctbeam=correctbeam, flagJumps=flagJumps,
                           writeSummary=writeSummary, showsig=showsig)

                # If we CTRL+C while in reduction, sometimes map is written and it is empty.
                # this is just a safe check to see if reduction finished, otherwise stop script.
                if data.Unit != 'Flux density [Jy/beam]':
                    raise RuntimeError('Stopping script: either CTRL+C was used or reduction failed.')
                if data.ScanParam.ScanNum != scan:
                    raise RuntimeError('Stopping script: either CTRL+C was used or reduction failed.')

                # Immediately rename summary if used and move to new folder
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
                mapping(oversamp=4, system=system, sizeX=xsize, sizeY=ysize,
                        limitsZ=[-0.8,1.5], noPlot=True)
                
                # Add MJD of middle and integration time to dumped map
                data.Map.MJDref = (data.ScanParam.MJD[-1] + data.ScanParam.MJD[0]) / 2  # MJD
                data.Map.Tint = np.sum(data.ScanParam.get('deltat'))  # seconds

                # Save unsmoothed map, "native" resolution (m.BeamSize = data.BolometerArray.BeamSize)
                data.Map.dumpMap(scanname)

                # Assign BoA map to variable m
                m = restoreFile(scanname)

                # For this scan, add non-noisy area and median noise info to summary
                if writeSummary:
                    info('Smoothing copy of map for summary at final resolution...')
                    # copy map and smooth with same kernel as final file
                    m_smooth = copy.deepcopy(m)
                    __smoothBy(m_smooth, smoothby_deg)

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

                #------------------------------------
                # --- OBSERVER PAUSE AT REDUCTION ---
                # -----------------------------------
                if observer == True:
                    # Show this scan's map, should be last finished scan, unsmoothed just as in reduction
                    rmsArray = np.where(m.Weight > 0.0, 1.0 / np.sqrt(m.Weight), np.NaN)
                    mediannoise = np.nanmedian(rmsArray)
                    meannoise = np.nanmean(rmsArray[rmsArray<2*mediannoise])  # no borders
                    del rmsArray  # free memory
                    caption = '%s - %s - Iter%i - Scan %i | Intensity (no smoothing): -3 to +10 sigma'%(source, fe, iter, scan)
                    m.display(aspect=1, limitsZ=[-3*meannoise, +10*meannoise], caption=caption)
                    print('')
                    msg  = "-------------------------------------------------------------------------\n"
                    msg += "-> Check map and answer whether it looks fine, regardless of sensitivity.\n"
                    msg += "   ** black patches (NaN values) are OK but please mention in obslogs ** \n"
                    msg += "                                                                         \n"
                    msg += "-> Consider phase setting (> ~80% good) and sensitivity (mJy sqrt(s))"+" \n"
                    msg += "   for re-calibration!                                                   \n"
                    msg += "-------------------------------------------------------------------------\n"
                    msg += "Map OK:                      <Enter>\n"
                    msg += "Map NOT OK:           no/n + <Enter>\n\n"
                    msg += "Observer input:"
                    obs_input = raw_input(msg)
                    obs_input = str(obs_input)

                    if str.upper(obs_input) in ['NO', 'N']:
                        print("------------------------------------------------------------------------")
                        info('Removing reduction of scan %i in file:'%(scan))
                        print('         %s'%scanname)
                        os.remove(scanname)
                        raise RuntimeError("Stopping script:"
                                           "\nMap of scan %i was reported as bad!"%scan +\
                                           "\n*** Remember to add this scan to 'badscans' list in"+\
                                           "\nreduction script before executing again ***'")
                    else:
                        print("------------------------------------------------------------------------")
                        info('Scan %i OK'%scan)
            
            else:
                # Retrieve BoA map
                info('Reduction for scan %i found (Iter%i | scan %i/%i). Loading...'%(scan, iter, i+1, len(scans)))
                m = restoreFile(scanname)

            if np.all(np.isnan(m.Data)):
                os.system('rm -f %s'%scanname)
                raise RuntimeError("Scan %i produced an all-NaN map. This almost always indicates "%scan+\
                                   "incorrect map bounds or coordinate system. Aborting reduction "+\
                                   "script. Please check your map bounds and coordinate system.")
            
            info('Coadding...')
            if ms and m:
                if np.shape(ms.Data)!=np.shape(m.Data):
                    raise RuntimeError("Coadded map and scan %i map have different grids. Cannot co-add!"%scan+\
                                       "\nDid you change map size or padding in reduction script?")
                ms = mapsumfast([ms,m])
        
            elif not ms:
                ms = copy.deepcopy(m)
            
            # Add integration time and delete map m of scan
            try:
                tint += m.Tint
                mymjdrefs.append(m.MJDref)
                mytints.append(m.Tint)
            except:
                pass
            del m  # free memory

            if doPlot:
                # SNR map creation
                snrMap = copy.deepcopy(ms)  # Signal
                # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)
                snrMap.Data = np.where(snrMap.Weight > 0.0, snrMap.Data * np.sqrt(snrMap.Weight), np.NaN)
                # plotting
                caption = '%s - %s - Iter%i - Coadded up to scan %i | SNR (no smoothing): -3 to +10'%(source, fe, iter, scan)
                snrMap.display(aspect=1,limitsZ=[-3, +10], caption=caption)
                del snrMap  # free memory


        # ==========================================================
        # ITERATION COMPLETE, NO SMOOTHING AT ALL UP TO HERE IN "ms"
        # ==========================================================
        del mymodel  # free memory

        # Now create final iter maps and FITS.
        # First, smooth co-added if required:
        if smoothby_deg > 0.0:
            print('')
            info('Smoothing co-added map for iteration %i by %.1f"...'%(iter, smoothby_arcsec))
            nativebeam = ms.BeamSize
            __smoothBy(ms, smoothby_deg)
            newbeam = ms.BeamSize
            print('         Original beam: %.3f"     New beam: %.3f"'%(nativebeam*3600, newbeam*3600))

        # RMS map creation
        rmsMap = copy.deepcopy(ms)  # Signal
        rmsMap.Data = np.where(rmsMap.Weight > 0.0, 1.0 / np.sqrt(rmsMap.Weight), np.NaN)  # Noise = 1/sqrt(weight)

        # SNR map creation
        snrMap = copy.deepcopy(ms)  # Signal
        snrMap.Data = np.where(snrMap.Weight > 0.0, snrMap.Data * np.sqrt(snrMap.Weight), np.NaN)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

        # Compute statistics, let __writeFits handle clipping
        messages.info('Computing aperture-based noise statistics...')
        # compute noise statistics in a circular aperture of radius 2 arcmin centered on map center
        radius_deg = 3.0 / 60.0  # 2 arcmin
        # create a mask for the circular aperture
        x_indices, y_indices = np.indices(ms.Data.shape)
        x_center = ms.WCS['CRPIX1']
        y_center = ms.WCS['CRPIX2']
        aperture_mask = (x_indices - x_center)**2 + (y_indices - y_center)**2 <= (radius_deg / abs(ms.WCS['CDELT1']))**2
        minnoise = np.nanmin(rmsMap.Data[aperture_mask])  # on aperture
        meannoise = np.nanmean(rmsMap.Data[aperture_mask])  # on aperture
        mediannoise = np.nanmedian(rmsMap.Data)  # on full map
        aperturegauss = 10*np.exp(-((x_indices - x_center)**2 + (y_indices - y_center)**2)/(2*(radius_deg/abs(ms.WCS['CDELT1']))**2))
        aperturegauss *= np.where(aperture_mask, 1., np.nan)  # cut it
        minap = np.nanmin(aperturegauss)
        # create an image for this aperture to display
        apertureMap = copy.deepcopy(rmsMap)
        apertureMap.Data = aperturegauss
        

        # plot SnR map
        caption = '%s - %s - Iter%i - Coadded up to scan %i | SNR (smoothed by %.1f"): -3 to +10 '%(source, fe, iter, scan, smoothby_arcsec)
        snrMap.display(aspect=1,limitsZ=[-3, 10], caption=caption)

        # plot noisemap contours
        if clip != -1:
            rmsMap.display(aspect=1,limitsZ=[0, clip*mediannoise],doContour=1,levels=[clip*mediannoise],overplot=1)
        else:
            # use 2*median noise to show "edges" of map, but not to clip
            rmsMap.display(aspect=1,limitsZ=[0, 2*mediannoise],doContour=1,levels=[2*mediannoise],overplot=1)

        # plot aperture map
        apertureMap.display(aspect=1,limitsZ=[0, minap],doContour=1,levels=[minap],overplot=1)#,colors=['cyan'])

        print('')
        print("####################### Iteration %i finished ########################"%(iter))
        print(" Time: %3.1f Hrs | min. noise: %3.1f mJy/b | mean noise: %3.1f mJy/b "%(tint/3600, 1000*minnoise,1000*meannoise))
        print("#####################################################################")

        # Save full-iteration map (will be smoothed if smooth > 0.0)
        outname = "ReducedFiles/"+str(myname)+"-coadded-flux-iter"+str(iter)+".data"  # goes into ReducedFiles dir
        ms.dumpMap(outname)

        # Save FITS file if requested
        if writefits:
            outname = str(myname)+"-coadded-iter"+str(iter)+".fits"
            outname = "FITSfiles/" + outname                         # goes into FITSfiles dir.
            __writeFits(ms, outfile=outname, overwrite=1, clip=clip)
        
        del ms  # free memory
        del rmsMap  # free memory
        del snrMap  # free memory
        del radius_deg, x_indices, y_indices, x_center, y_center  # free memory
        del aperture_mask, aperturegauss  # free memory
        del apertureMap  # free memory

if observer==False:
    print('')
    print("#####################################################################")
    print("                         Reduction finished                          ")
    print("#####################################################################")