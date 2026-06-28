# SCRIPT TO CORRECT THE BEAM SIZE AND RE-SCALING IN CO-ADDED MAPS
# AND EXPORT THEM AS FITS INTO "./BeamCorrected" directory.
# INTENDED TO BE RAN AUTOMATICALLY AFTER REDUCTION SCRIPT
# NO USER INPUT NEEDED

# ===== BEGINNING OF CODE, DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE DOING =====
import warnings
import copy as copy
import BoaMapping as BOAMAP
from mars.fortran import fMap

print('=======================================================================')
info('Beginning SKY, RMS, and SNR map corrections for iteration maps')
print('=======================================================================')

# Create dir if missing
if os.path.exists("FITSfiles/BeamCorrected") == False:
    os.makedirs("FITSfiles/BeamCorrected")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")

    # determine nominal AMKID beam size
    if fe == 'LFA':
        AMKID_beamsize = 16.5  # arcsec
    elif fe == 'HFA':
        AMKID_beamsize = 7.63  # arcsec
    from_where = 'Nominal value'

    info('Attempting to extract AMKID %s beam size from beammap info...'%fe)
    # Try and read nominal beam size from merged beammap
    printwarn = False
    try:
        # find beam_map reduced file
        beammap_fnames = []
        for filename in os.listdir('CalFiles'):
            if 'beam_map' in filename and 'merged' in filename and fe in filename:
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
        print('         %s'%beammap_fname_full)
        info('Beam size is %.3f "'%AMKID_beamsize)
        

    except:
        printwarn = True
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
            print('         %s'%mycoadded_fullfname)
            ms = restoreFile(mycoadded_fullfname)
        else:
            print('')
            warn('File for Iteration %i not found!'%iter)
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
        # Redo correct convolution rescaling (if any)
        uncorrected_scale = (UNCORRECTED_CONVOLVED_FWHM**2 / UNCORRECTED_NATIVE_FWHM**2)  # =1 if smooth=0
        correct_scale = (CORRECT_CONVOLVED_FWHM**2 / CORRECT_NATIVE_FWHM**2)  # =1 if smooth=0
        imagecorrfactor = correct_scale / uncorrected_scale  # =1 if no smoothing, then only beam in header changes

        # apply factor
        ms.Data *= imagecorrfactor
        ms.Weight /= imagecorrfactor**2  # apply correction factor squared to variance map

        # =================================================================
        # Now the map is in the correct AMKID^2 + SMOOTHING^2 Jy/beam units.
        # =================================================================
        # All that's left to do is correct the written beam size
        ms.BeamSize = CORRECT_CONVOLVED_FWHM  # = AMKID's native if smooth=0

        # Compute changes to fluxes:
        # Flux of a source is sum_apperture(Jyb_i * pixarea / beamarea)
        # so  Fakeflux = sum_apperture(Jyb_original_i)  * pixarea / uncorr_convolved_beam
        # and Realflux = sum_apperture(Jyb_corrected_i) * pixarea / correct_convolved_beam
        # where Jyb_corrected_i = imagecorrfactor * Jyb_original_i
        # so fluxfactor = Realflux / Fakeflux = imagecorrfactor * uncorr_convolved_beam/correct_convolved_beam
        # or fluxfactor = Realflux / Fakeflux = imagecorrfactor * uncorr_convolved_FWHM^2/correct_convolved_FWHM^2
        fluxfactor = imagecorrfactor * UNCORRECTED_CONVOLVED_FWHM**2/CORRECT_CONVOLVED_FWHM**2

        # RealFlux = FakeFlux * fluxfactor -> FakeFlux = (1/fluxfactor) * RealFlux = fraction * RealFlux
        fractionofreal = 1./fluxfactor
        percentofreal =  fractionofreal * 100.

        # create summary for first iteration, will always be there...
        if iter==1:
            correctionsummary = 'Summary of corrections (should be same for all iterations):'
            correctionsummary += '\n------------------------------------------------------------'
            correctionsummary += '\nBeam read from files:                       %s "'%(str(np.round(UNCORRECTED_CONVOLVED_FWHM*3600., 3)).ljust(6))
            correctionsummary += '\nSmoothing was:                              %s "'%(str(np.round(smoothby_deg*3600., 3)).ljust(6))
            correctionsummary += '\nUnsmoothed beam:                            %s "'%(str(np.round(UNCORRECTED_NATIVE_FWHM*3600., 3)).ljust(6))
            correctionsummary += '\nAMKID median beam:                          %s "'%(str(np.round(CORRECT_NATIVE_FWHM*3600., 3)).ljust(6)) + ' (%s)'%from_where
            correctionsummary += '\nSky map (image) was rescaled by:            %sx '%(str(np.round(imagecorrfactor, 3)).ljust(5))
            correctionsummary += '\nVariance map (1 / Weight) was rescaled by:  %sx '%(str(np.round(imagecorrfactor**2, 3)).ljust(5))
            correctionsummary += '\nCorrected beam after smoothing:             %s "'%(str(np.round(CORRECT_CONVOLVED_FWHM*3600., 3)).ljust(6))
            correctionsummary += '\n------------------------------------------------------------'
            correctionsummary += '\nFluxes before correction were %.1f'%(percentofreal) + r'% of the expected flux'
            correctionsummary += '\n------------------------------------------------------------'

        # now export to fits
        outname = 'FITSfiles/BeamCorrected/' + str(myname)+"-coadded-iter"+str(iter)+"-beamCorrected.fits" # Goes into ./BeamCorrected directory.
        auxwriteFits(ms, outfile=outname, overwrite=1, clip=clip)
        info('Beam-corrected FITS written to:')
        print('         %s'%outname)

        # free memory
        ms = None

print('')
print(correctionsummary)
if printwarn ==True:
    warn('Beam size extraction from beam maps was not possible: the nominal value of %.3f " for %s was used. '%(AMKID_beamsize, fe)+\
         'Check if a merged beammap (beam_map_SCAN_%s_merged.csv) file exists in CalFiles/ directory.'%fe)
info('Check "BeamCorrected/" directory for beam-corrected FITS!')