# =================================
# ==== BEGINNING OF USER INPUT ====
# =================================

# SCRIPT TO CORRECT THE BEAM SIZE AND RE-SCALING IN A CO-ADDED MAP
# AND EXPORT IT INTO "./BeamCorrected" directory as FITS

# --- Source and map parameters ---
source    = 'OG259'       # As in observing logs
fe        = 'LFA'           # Frontend, either 'LFA' or 'HFA'
system    = 'GAL'            # Coordinate system of reduced map, 'EQ', 'GAL' or 'HO'
iter      = 2               # Which iteration is the map to be corrected
flagJumps = True           # Whether the maps to show were de-jumped with
                            # 'flagJumps = True' at reduction
smoothby_arcsec = 8.        # By how much was the map smoothed?
AMKID_beamsize  = 17.       # Nominal native-res beam FWHM of AMKID

# ==============================
# ===== END OF USER INUPUT =====
# ==============================






# ===== BEGINNING OF CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
import warnings
import copy as copy
import BoaMapping as BOAMAP
from mars.fortran import fMap

# variable checks
if fe not in ['LFA', 'HFA']:
    raise ValueError("fe must be either 'LFA' or 'HFA'.")
if system not in ['EQ', 'GAL', 'HO']:
    raise ValueError("system must be either 'EQ', 'GAL', or 'HO'.")
if iter < 1 or iter > 3:
    raise ValueError("iter must be 1, 2, or 3.")

def auxwriteFits(data=None,outfile='boaMap.fits',overwrite=0,limitsX=[],limitsY=[],intensityUnit="Jy/beam",clip=-1):
        """
        DES: store the current map (2D array with WCS info) to a FITS file
        INP: (str)   outfile: output file name (default boaMap.fits)
             (bool) overwrite: overwrite existing file -
                              default = 0: do not overwrite existing file
             (f list) limitsX/Y: optional map limits (in world coordinates)
             (string) intensityUnit: optional unit of the intensity (default: "Jy/beam")
        """
        from mars import BoaFits

        if os.path.exists(outfile):
            if not overwrite:
                print('File %s exists' % outfile)
                return
        if not data:
            data = data.Map
        try:
            dataset = BoaFits.createDataset("!" + outfile)
        except Exception, data:
            print('Could not open file %s: %s' % (outfile, data))
            return

        
        localMap = copy.deepcopy(data)
        
        try:
            # RMS map creation
            rmsMap = copy.deepcopy(localMap)  # Signal
            rmsMap.Data = np.where(rmsMap.Weight > 0.0, 1.0 / np.sqrt(rmsMap.Weight), np.NaN)  # Noise = 1/sqrt(weight)

            # SNR map creation
            snrMap = copy.deepcopy(localMap)  # Signal
            snrMap.Data = np.where(snrMap.Weight > 0.0, snrMap.Data * np.sqrt(snrMap.Weight), np.NaN)  # SNR = signal * sqrt(weight) = signal / sqrt(noise^2)

            if clip > 0:
                mediannoise = np.nanmedian(rmsMap.Data)
                mask = np.where(rmsMap.Data > clip * mediannoise)
                localMap.Data[mask] = np.NaN
                rmsMap.Data[mask] = np.NaN
                snrMap.Data[mask] = np.NaN
 
            #write FLux plane                                                            
            localMap._Image__writeImage(dataset, "Intensity", intensityUnit=intensityUnit)
            #write RMS plane
            rmsMap._Image__writeImage(dataset, "Intensity", intensityUnit=intensityUnit+" (RMS)")
            #write SNR plane
            snrMap._Image__writeImage(dataset, "Intensity", intensityUnit='SNR')
            dataset.close()
            
        except Exception, data:
            print('Could not write data to file %s: %s' % (outfile, data))
            return

        
        localMap = 0  # free memory
        snrMap = 0
        rmsMap = 0
        dataset = 0

# Define myname variable
myname = str(fe) + "-" + str(source) + "-" + str(system)
if flagJumps:
    myname += "-flagJumps"

# Extract file name of corresponding map
mycoadded_fname = str(myname) + "-coadded-flux-iter" + str(iter-1) + ".data"
mycoadded_fullfname = "ReducedFiles/" + mycoadded_fname

# Try to retrieve it
try:
    ms = restoreFile(mycoadded_fullfname)
except:
    print('File %s not found!'%mycoadded_fullfname)

# Extract uncorrected, convolved beam FWHM. Might be actually convolved or not.
UNCORRECTED_CONVOLVED_FWHM = ms.BeamSize

# If map was smoothed
if smoothby_arcsec > 0.0:
    # Continue
    pass
else:
    # Will leave map untouched but change the written beam
    smoothby_arcsec = 0.0

# Recover native beam of map before smoothing (if any)
UNCORRECTED_NATIVE_FWHM = np.sqrt(UNCORRECTED_CONVOLVED_FWHM^2 - smoothby_arcsec^2)

# The correct native FWHM is AMKID's.
CORRECT_NATIVE_FWHM = AMKID_beamsize
# If UNCORRECTED_NATIVE_FWHM == CORRECT_NATIVE_FWHM, then everything was done correctly
# And nothing changes here onwards.

# Compute AMKID's convolved beam after smoothing (if any)
CORRECT_CONVOLVED_FWHM = np.sqrt(CORRECT_NATIVE_FWHM^2 + smoothby_arcsec^2)  # = AMKID's native if smooth=0

# Undo wrong convolution rescaling (if any)
# (native^2 + smoothing^2) / native^2 was multiplied to 
# convolved intensity map and its square to the variance = 1/Weight map.
uncorrected_scale = (UNCORRECTED_CONVOLVED_FWHM**2 / UNCORRECTED_NATIVE_FWHM**2)  # =1 if smooth=0
ms.Data /= uncorrected_scale  # undo uncorrected scale
ms.Weight /= (1./uncorrected_scale**2)  # undo uncorrected scale^2

# Redo correct convolution rescaling (if any)
correct_scale = (CORRECT_CONVOLVED_FWHM**2 / CORRECT_NATIVE_FWHM**2)  # =1 if smooth=0
ms.Data *= correct_scale  # redo correct scale
ms.Weight *= (1./correct_scale**2)  # redo correct scale^2

# =================================================================
# Now the map is in the correct AMKID^2 + SMOOTHING^2 Jy/beam units.
# =================================================================
# All that's left to do is correct the written beam size
ms.BeamSize = CORRECT_CONVOLVED_FWHM  # = AMKID's native if smooth=0

# now export to fits
if os.path.exists('BeamCorrected') == False:
    os.makedirs("BeamCorrected")

outname = 'BeamCorrected/' + str(myname)+"-coadded-iter"+str(iter)+"-beamCorrected.fits" # Goes into ./BeamCorrected directory.
auxwriteFits(ms, outfile=outname, overwrite=1)

# free memory
ms = None