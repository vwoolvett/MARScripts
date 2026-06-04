import copy as copy

# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================

# --- Source and map parameters ---
source    = 'Name'          # As in observing logs
fe        = 'LFA'           # Frontend, either 'LFA' or 'HFA'
system    = 'EQ'            # Coordinate system of reduced map, 'EQ', 'GAL' or 'HO'
iter      = 2               # Which iteration of the reduction to show (usual 1-3)
show      = 'sig'           # Show Signal (sig), Noise (rms) or SNR (snr)
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
    m = restoreFile(scanname)
    m.smoothBy(8./3600.)
    info('File found, displaying scan %s...'%(scan))

    if show == 'sig':
        m.display(aspect=1,limitsZ=[-0.2,0.5])

    else:
        # creating SNR map
        snrMap = copy.deepcopy(m)  # Signal
        snrMap.Data *= np.sqrt(snrMap.Weight)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

        # rescaling SNR map by beam size 
        a = snrMap.computeRms()
        scale = snrMap.RmsBeam
        snrMap.Data /= np.array(scale,'f')

        if show =='snr':
            # plotting
            snrMap.display(aspect=1,limitsZ=[-4,12])

        else:
            # creating rms map
            rmsMap = copy.deepcopy(m)  # Signal
            rmsMap.Data =  (rmsMap.Data*0.0+1.0)/np.sqrt(rmsMap.Weight) # Noise = 1/sqrt(weight)

            # Rescaling rms map by beam size
            rmsMap.Data *= np.array(scale,'f')

            # median noise
            mediannoise = np.nanmedian(rmsMap.Data)

            # plotting
            rmsMap.display(aspect=1, limitsZ=[0, 2*mediannoise])

    raw_input()