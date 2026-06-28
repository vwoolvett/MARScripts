# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================
# NOTE: script intended to be used after reducing all scans (iteration 1)
# to assess whether any scan should be added to or removed from
# badscans. It takes all variables but 'iter' and 'show' from reduction
# script, so it will not work properly if not ran after reduction
# in the same MARS session.

# --- Plotting parameters ---
iter      = 1               # Which iteration of the reduction to show (usual 1)
show      = 'sig'           # Show Signal-SKY (sig), Noise (rms), or SNR (snr)

# ==============================
# ===== END OF USER INUPUT =====
# ==============================
















# ===== BEGINNING OF CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
import warnings
import copy as copy

if iter < 1 or iter > 3:
    raise ValueError("iter must be 1, 2, or 3.")
if show not in ['sig', 'rms', 'snr']:
    raise ValueError("show must be 'sig', 'rms', or 'snr'.")

# NO SMOOTHING
smoothby_arcsec = 0.

if len(scans) == 0 and not os.path.exists(obslogsdir):
    raise ValueError('STOPPING SCRIPT: execute reduction script to pre-define required variables.')

# display message
msg = '''\
-------------------------------
Show next map:          <Enter>
Quit:               q + <Enter>
-------------------------------
user input:
'''

# =============================
# Beginning of map display loop
# =============================
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for i, scan in enumerate(scans):
        scanname = "ReducedFiles/"+str(myname)+"-"+str(scan)+"-iter"+str(iter)+".data"
        print('')
        info('Retrieving reduction for scan %s (iter %i) ...'%(scan, iter))

        # check if reduction exists
        globlist = glob(scanname)
        if len(globlist) == 0:
            warn('File not found. Skipping...')
            print('')
            continue
        
        try:
            # retrieve unsmoothed map
            m = restoreFile(scanname)
            info('File found')

        except:
            # failed
            warn('File found, but could not open, Skipping...')
            print('')
            continue
    
        if show == 'sig':
            info('Displaying Signal map...')
            rmsArray = np.where(m.Weight > 0.0, 1.0 / np.sqrt(m.Weight), np.NaN)
            mediannoise = np.nanmedian(rmsArray)
            meannoise = np.nanmean(rmsArray[rmsArray<2*mediannoise])  # no borders
            del rmsArray  # free memory
            caption = '%s - %s - Iter%i - Scan %i | Intensity (no smoothing): -3 to +10 sigma'%(source, fe, iter, scan)
            m.display(aspect=1, limitsZ=[-3*meannoise, +10*meannoise], caption=caption)
            del m  # free memory

        elif show == 'rms':
            info('Displaying Noise map...')
            rmsMap = copy.deepcopy(m)  # Signal
            rmsMap.Data = 1.0 / np.sqrt(rmsMap.Weight)  # Noise = 1/sqrt(weight)
            mediannoise = np.nanmedian(rmsMap.Data)
            caption = '%s - %s - Iter%i - Scan %i | RMS (no smoothing): 0 to 2 x median'%(source, fe, iter, scan)
            rmsMap.display(aspect=1, limitsZ=[0, 2*mediannoise], caption=caption)
            del m  # free memory
            del rmsMap  # free memory
            
        else:
            info('Displaying Signal-to-Noise map...')
            snrMap = copy.deepcopy(m)  # Signal
            snrMap.Data *= np.sqrt(snrMap.Weight)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)
            caption = '%s - %s - Iter%i - Scan %i | SNR (no smoothing): -3 to +10'%(source, fe, iter, scan)
            snrMap.display(aspect=1,limitsZ=[-3, 10], caption=caption)
            del m  # free memory
            del snrMap  # free memory

        usrinput = raw_input(msg)
        if str.upper(str(usrinput)) == 'Q':
            break