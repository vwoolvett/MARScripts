import copy as copy

# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================

# --- Source and map parameters ---
source    = 'Name'                # as in APECS/ObsLogs
fe        = 'LFA'                 # frontend, either 'LFA' or 'HFA'
system    = 'EQ'                  # coordinate system of reduced maps, either 'EQ' or 'GAL' is usual
iter      = 2                     # which iteration of the reduction to show
flagJumps = False                 # whether the maps were deJumped with 'flagJumps = True' at reduction

# ----- Scans ------
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
    print('Retrieving scan %s (file: %s)...'%(scan, scanname))

    globlist = glob(scanname)
    if len(globlist) == 0:
        print('File %s not found, skipping...'%scanname)
        continue
    m=restoreFile(scanname)
    m.smoothBy(8./3600.)
    print('Displaying scan %s (file: %s)...'%(scan, scanname))
    m.display(aspect=1,limitsZ=[-0.2,0.5])    
    raw_input()
    
