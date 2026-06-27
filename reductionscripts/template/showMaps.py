# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================

# --- Source and map parameters ---
source    = 'SrcName'       # As in observing logs
fe        = 'LFA'           # Frontend, either 'LFA' or 'HFA'
system    = 'EQ'            # Coordinate system of reduced map, 'EQ', 'GAL' or 'HO'
iter      = 1               # Which iteration of the reduction to show (usual 1-2)
show      = 'snr'           # Show Signal-SKY (sig), Noise (rms), or SNR (snr)
flagJumps = False           # Whether the maps to show were de-jumped with
                            # 'flagJumps = True' at reduction
smoothby_arcsec = 0.        # Default 0 arcsec. Smoothing before showing takes more time.

# ----- Scans ------
# If scans is empty, automatically retrieves all scans of the source
# specified above from the obslogs directory below
scans = []
obslogsdir = '~/obslogs'  # at MPIfR: '/apex-archive/obslogs/M-PROJECT.CODE-IN-CAPS/obslogs'

# ==============================
# ===== END OF USER INUPUT =====
# ==============================
















# ===== BEGINNING OF CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
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



# variable checks
if fe not in ['LFA', 'HFA']:
    raise ValueError("fe must be either 'LFA' or 'HFA'.")
if system not in ['EQ', 'GAL', 'HO']:
    raise ValueError("system must be either 'EQ', 'GAL', or 'HO'.")
if iter < 1 or iter > 3:
    raise ValueError("iter must be 1, 2, or 3.")
if show not in ['sig', 'rms', 'snr']:
    raise ValueError("show must be 'sig', 'rms', or 'snr'.")

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

# Define myname variable
myname = str(fe) + "-" + str(source) + "-" + str(system)
if flagJumps:
    myname += "-flagJumps"

# smoothby to deg
smoothby_deg = smoothby_arcsec / 3600.

# display message
msg = '''\
Show next map:      Enter
Quit:               q
'''

# =============================
# Beginning of map display loop
# =============================
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for i,scan in enumerate(scans):
        scanname = "ReducedFiles/"+str(myname)+"-"+str(scan)+"-iter"+str(iter)+".data"
        print('')
        info('Retrieving reduction for scan %s (iter %i) ...'%(scan, iter))

        # check if reduction exists
        globlist = glob(scanname)
        if len(globlist) == 0:
            warn('File not found. Skipping...')
            print('')
            continue

        # retrieve unsmoothed map
        m = restoreFile(scanname)

        # smooth if needed
        if smoothby_arcsec > 0.0:
            info('File found, smoothing by %.1f arcseconds...'%smoothby_arcsec)
            auxsmoothby(m, smoothby_deg)
        
        else:
            info('File found')
    
        if show == 'sig':
            info('Displaying Signal map...')
            std = np.nanstd(m.Data)
            m.display(aspect=1,limitsZ=[-3*std,+10*std])

        elif show == 'rms':
            info('Displaying Noise map...')
            # RMS map creation
            rmsMap = copy.deepcopy(m)  # Signal
            rmsMap.Data = 1.0 / np.sqrt(rmsMap.Weight)  # Noise = 1/sqrt(weight)

            # median noise
            mediannoise = np.nanmedian(rmsMap.Data)

            # plotting
            rmsMap.display(aspect=1, limitsZ=[0, 2*mediannoise])
            
        else:
            info('Displaying Signal-to-Noise map...')
            # SNR map creation
            snrMap = copy.deepcopy(m)  # Signal
            snrMap.Data *= np.sqrt(snrMap.Weight)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

            # plotting
            snrMap.display(aspect=1,limitsZ=[-3, 10])

        usrinput = raw_input(msg)
        if str.upper(str(usrinput)) == 'Q':
            break