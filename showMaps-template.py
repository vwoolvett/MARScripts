import copy as copy

# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================

# --- Source and map parameters ---
source    = 'OG259'              # as in APECS/ObsLogs
fe        = 'LFA'               # frontend, either 'LFA' or 'HFA'
system    = 'GAL'                # coordinate system of reduced maps, either 'EQ' or 'GAL' is usual
iter      = 2                   # which iteration of the reduction to show
flagJumps = False               # whether the maps were deJumped with 'flagJumps = True' at reduction

# ----- Scans ------
scans = [27974]#,27975,27979,27980,27990,28212,28213,28217,28218,28231,28232,28235,28493,28494,28498,28499,28516,28517,28775,28776]                      # 'Auto' or list of scans to reduce
                                # NOTE: If using 'Auto', make sure to set the correct
                                # source name and frontend above

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
    info('Retrieving reduction for scan %s ...'%(scan))
    print('File: %s'%scanname)

    globlist = glob(scanname)
    if len(globlist) == 0:
        warn('File not found, check input params. Skipping...')
        continue
    m=restoreFile(scanname)
    m.smoothBy(8./3600.)
    info('Displaying scan %s...'%(scan, scanname))
    m.display(aspect=1,limitsZ=[-0.2,0.5])    
    raw_input()
    
