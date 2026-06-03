import copy as copy

# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================

# --- Source and map parameters ---
source    = 'OG259'             # as in APECS/ObsLogs
fe        = 'LFA'               # frontend, either 'LFA' or 'HFA'
system    = 'EQ'                # coordinate system of reduced maps, either 'EQ' or 'GAL' is usual
iter      = 1                   # which iteration of the reduction to show
flagJumps = False               # whether the maps were deJumped with 'flagJumps = True' at reduction

# ----- Scans ------
# 'Auto' or list of scans to reduce. If using 'Auto', make sure to set the correct
# source name and frontend above
scans = [27974]

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
        warn('File not found, check script params. Skipping...')
        continue
    m=restoreFile(scanname)
    m.smoothBy(8./3600.)
    info('File found, displaying scan %s...'%(scan))
    m.display(aspect=1,limitsZ=[-0.2,0.5])    
    raw_input()
    
