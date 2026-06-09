# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================

# --- Source and map parameters ---
source    = 'SrcName'       # As in observing logs
fe        = 'LFA'           # Frontend, either 'LFA' or 'HFA'
system    = 'EQ'            # Coordinate system of reduced map, 'EQ', 'GAL' or 'HO'
iter      = 2               # Which iteration of the reduction to show (usual 1-2)
show      = 'sig'           # Show Signal (sig), Noise (noise) or SNR (snr)
flagJumps = False           # Whether the maps to show were de-jumped with
                            # 'flagJumps = True' at reduction
smoothby_arcsec = 8.        # Default 8. arcsec. Use same smoothing as reduction.

# ----- Scans ------
# If empty, automatically retrieves all scans of source from Obslogs
# NOTE: CURRENTLY NOT FUNCTIONAL, PLEASE MANUALLY INPUT SCAN NUMBERS
scans = []



# ==============================
# ===== END OF USER INUPUT =====
# ==============================










# ===== BEGINNING OF CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
import copy as copy
import BoaMapping as BOAMAP
from mars.fortran import fMap

# Define myname variable
myname = str(fe) + "-" + str(source) + "-" + str(system)
if flagJumps:
    myname += "-flagJumps"

# display message
msg = '''\
Show next map:      Enter
Quit:               q
'''

# smoothby to deg
smoothby_deg = smoothby_arcsec / 3600.

# define the good functions :)
def auxsmoothby(m, Size=smoothby_deg):
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

with np.errstate(divide='ignore', invalid='ignore'):
    for i,scan in enumerate(scans):
        scanname = "ReducedFiles/"+str(myname)+"-"+str(scan)+"-iter"+str(iter)+".data"
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
            auxsmoothby(m, smoothby_deg)
    
        if show == 'sig':
            info('File found, displaying Signal map...')
            m.display(aspect=1,limitsZ=[-0.2,0.5])

        else:
            if show =='noise':
                info('File found, displaying Noise map...')
                # RMS map creation
                rmsMap = copy.deepcopy(m)  # Signal
                rmsMap.Data = 1.0 / np.sqrt(rmsMap.Weight)  # Noise = 1/sqrt(weight)

                # median noise
                mediannoise = np.nanmedian(rmsMap.Data)

                # plotting
                rmsMap.display(aspect=1, limitsZ=[0, 2*mediannoise])
            

            else:
                info('File found, displaying Signal-to-Noise map...')
                # SNR map creation
                snrMap = copy.deepcopy(m)  # Signal
                snrMap.Data *= np.sqrt(snrMap.Weight)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

                # plotting
                snrMap.display(aspect=1,limitsZ=[-4,12])

        usrinput = raw_input(msg)
        if str.upper(str(usrinput)) == 'Q':
            break