# SCRIPT TO CORRECT THE BEAM SIZE AND RE-SCALING IN CO-ADDED MAPS
# AND EXPORT THEM AS FITS INTO "./BeamCorrected" directory
# INTENDED TO BE RAN AUTOMATICALLY AFTER REDUCTION SCRIPT
# NO USER INPUT NEEDED
# TODO: extract AMKID_beamsize from average of beammaps instead!

# ===== BEGINNING OF CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
import warnings
import copy as copy
import BoaMapping as BOAMAP
from mars.fortran import fMap

print('=======================================================================')
info('Beginning SKY, RMS, and SNR map corrections for iteration maps')
print('=======================================================================')

# Create dir if missing
if os.path.exists('BeamCorrected') == False:
    os.makedirs("BeamCorrected")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")

    # determine nominal AMKID beam size
    AMKID_beamsize  = 17.  # arcsec
    from_where = 'Nominal value'

    info('Attempting to extract AMKID beam size from beammap info...')
    # Try and read nominal beam size from merged beammap
    if True:
        # find beam_map reduced file
        beammap_fnames = []
        for filename in os.listdir('CalFiles'):
            if 'beam_map' in filename and 'merged' in filename:
                beammap_fnames.append(filename)
        if len(beammap_fnames) > 1:
            # use last?
            beammap_fname = beammap_fnames[-1]
        elif len(beammap_fnames) == 1:
            beammap_fname = beammap_fnames[0]
        else:
            raise ValueError # stop trying
        
        # read and extract average beam size
        beammap_fname_full = 'CalFiles/' + beammap_fname
        beamdict = readBeamMapDict(infile=beammap_fname_full, fe=fe)  # fe defined at reduction

        # median of geometric average FWHM of all kids
        # produces same beam area in case beam is ellyptical
        AMKID_beamsize = np.nanmedian(np.array([(beamdict[kid]['xfwhm']*beamdict[kid]['yfwhm'])**0.5 for kid in beamdict.keys()]))

        # succesfully got a beammap estimate
        from_where = 'Extracted from beam map: %s'%beammap_fname_full

        info('Success extracting beam size from:')
        print(beammap_fname_full)
        info('Beam size is %.3f "'%AMKID_beamsize)
        

    else:
        warn('Beam size extraction was not possible. Using nominal value of %.3f "'%AMKID_beamsize)
        pass

    for iter in range(1, niters+1):
        # Extract file name of corresponding map
        mycoadded_fname = str(myname) + "-coadded-flux-iter" + str(iter) + ".data"
        mycoadded_fullfname = "ReducedFiles/" + mycoadded_fname
        globlist = glob(mycoadded_fullfname)

        # Try to retrieve it
        if len(globlist) != 0:
            print('')
            info('Loading map in file:')
            print(mycoadded_fullfname)
            ms = restoreFile(mycoadded_fullfname)
        else:
            print('')
            warn('File %s not found!'%mycoadded_fullfname)
            print('')
            continue

        # Extract uncorrected, convolved beam FWHM. Might be actually convolved or not.
        UNCORRECTED_CONVOLVED_FWHM = ms.BeamSize

        # If map was not smoothed
        if smoothby_arcsec <= 0.0:
            # Will leave map untouched but change the written beam
            smoothby_arcsec = 0.0

        # define smoothing in deg
        smoothby_deg = smoothby_arcsec / 3600.

        # Recover native beam of map before smoothing (if any)
        UNCORRECTED_NATIVE_FWHM = np.sqrt(UNCORRECTED_CONVOLVED_FWHM**2 - smoothby_deg**2)

        # The correct native FWHM is AMKID's.
        CORRECT_NATIVE_FWHM = AMKID_beamsize / 3600.
        # If UNCORRECTED_NATIVE_FWHM == CORRECT_NATIVE_FWHM, then everything was done correctly
        # And nothing changes here onwards.

        # Compute AMKID's convolved beam after smoothing (if any)
        CORRECT_CONVOLVED_FWHM = np.sqrt(CORRECT_NATIVE_FWHM**2 + smoothby_deg**2)  # = AMKID's native if smooth=0

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
        fluxfactor = UNCORRECTED_CONVOLVED_FWHM/CORRECT_CONVOLVED_FWHM  # corrects header beam on integration
        fluxfactor *= correct_scale/uncorrected_scale  # corrects image
        percentofreal =  1./fluxfactor * 100.

        # create summary for first iteration, will always be there...
        if iter==1:
            correctionsummary = 'Summary of corrections (should be same for all iterations):'
            correctionsummary += '\n------------------------------------------------------------'
            correctionsummary += '\nBeam read from files:            %.3f "'%(UNCORRECTED_CONVOLVED_FWHM*3600.)
            correctionsummary += '\nSmoothing was:                   %.3f "'%(smoothby_deg*3600.)
            correctionsummary += '\nUnsmoothed beam:                 %.3f "'%(UNCORRECTED_NATIVE_FWHM*3600.)
            correctionsummary += '\nAMKID median beam:               %.3f "'%(CORRECT_NATIVE_FWHM*3600.) + ' (%s)'%from_where
            correctionsummary += '\nImage was rescaled by:           %.3fx'%(correct_scale/uncorrected_scale)
            correctionsummary += '\nCorrected beam after smoothing:  %.3f "'%(CORRECT_CONVOLVED_FWHM*3600.)
            correctionsummary += '\n------------------------------------------------------------'
            correctionsummary += '\nFluxes were %.1f'%(percentofreal) + r'% of the expected flux'
            correctionsummary += '\n------------------------------------------------------------'

        # now export to fits
        outname = 'BeamCorrected/' + str(myname)+"-coadded-iter"+str(iter)+"-beamCorrected.fits" # Goes into ./BeamCorrected directory.
        auxwriteFits(ms, outfile=outname, overwrite=1)
        info('Beam-corrected FITS written to:')
        print(outname)

        # free memory
        ms = None

print('')
print(correctionsummary)
info('Check "BeamCorrected/" directory for beam-corrected FITS!')