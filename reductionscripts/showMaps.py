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

# ----- Scans ------
# If empty, automatically retrieves all scans of source from Obslogs
# NOTE: CURRENTLY NOT FUNCTIONAL, PLEASE MANUALLY INPUT SCAN NUMBERS
scans = []

# ==============================
# ===== END OF USER INUPUT =====
# ==============================










# ===== BEGINNING OF CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
import copy as copy

# Define myname variable
myname = str(fe) + "-" + str(source) + "-" + str(system)
if flagJumps:
    myname += "-flagJumps"

for i,scan in enumerate(scans):
    scanname = "ReducedFiles/"+str(myname)+"-"+str(scan)+"-iter"+str(iter)+".data"
    info('Retrieving reduction for scan %s (iter %i) ...'%(scan, iter))
    globlist = glob(scanname)
    if len(globlist) == 0:
        warn('File not found. Skipping...')
        print('')
        continue

    # retrieve unsmoothed map
    m = restoreFile(scanname)
    
    if show == 'sig':
        info('File found, displaying Signal map...')
        m.smoothBy(8./3600.)
        m.display(aspect=1,limitsZ=[-0.2,0.5])

    else:
        # extract unsmoothed RMS, then smooth
        rmsMap = copy.deepcopy(m)
        rmsMap.Data = 1.0 / np.sqrt(rmsMap.Weight) # Noise = 1/sqrt(weight)
        
        m.smoothBy(8./3600.)
        rmsMap.smoothBy(8./3600)

        if show =='noise':
            info('File found, displaying Noise map...')
            # median noise
            mediannoise = np.nanmedian(rmsMap.Data)

            # plotting
            rmsMap.display(aspect=1, limitsZ=[0, 2*mediannoise])
            

        else:
            info('File found, displaying Signal-to-Noise map...')

            # creatie SNR map
            snrMap = copy.deepcopy(m)  # already smoothed Signal
            snrMap.Data /= rmsMap.Data  # divide by smoothed noise

            # plotting
            snrMap.display(aspect=1,limitsZ=[-4,12])

    raw_input()