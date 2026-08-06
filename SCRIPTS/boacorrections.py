def __smoothBy(m, Size):
    '''
    BoA-like smoothing but with correct variance propagation.

    - Data: convolved with K (same as BoA)
    - Weight: propagated via variance (K^2)
    - Coverage: convolved with K (same as BoA)
    '''
    # Build kernel (not normalized)
    pixsize = abs(m.WCS['CDELT2'])
    K0 = BoaMapping.Kernel(pixsize, Size).Data.astype(float)

    # Normalize kernel
    K = K0 / np.sum(K0)

    # Create elementwise-squared kernel for variance
    K2 = K**2

    # Smooth INTENSITY (same as BoA)
    #   I' = K * I     =     (K0/sum(K0_i)) * I
    # and ksmooth does
    #   I' = (K * I) / sum(K_i), but since sum(K_i)=1
    # then ksmooth does effectively
    #   I' = K * I, all good
    I1 = fMap.ksmooth(m.Data, K)

    # Correct variance propagation for weights:
    #   V' = K2 * V     =     (K0/sum(K0_i))^2 * V
    # but ksmooth does
    #   V' = (K2 * V) / sum(K2_i), and now sum(K2_i)!=1
    # then ksmooth does effectively
    #   V' = K2/sum(K2_i) * V
    # so an additional multiplication by sum(K2_i) is needed
    # to get back from ksmooth:
    #   V' = K2/sum(K2_i) * V * sum(K2_i) = K2 * V
    V0 = np.where(m.Weight > 0.0, 1.0 / m.Weight, np.NaN)
    V1 = fMap.ksmooth(V0, K2) * np.sum(K2)

    # Smooth COVERAGE (same as BoA)
    C1 = fMap.ksmooth(m.Coverage, K)
    
    # new scale per beam for Jy/beam units
    newbeam = np.sqrt(m.BeamSize**2 + Size**2)
    scale = (newbeam**2 / m.BeamSize**2)
    I1 *= scale  # now in Jy/newbeam
    V1 *= scale**2  # now in Jy^2/newbeam^2

    # create new weight map
    W1 = np.where(V1 > 0.0, 1.0 / V1, 0.0)

    # Update map with correct Jy/beam scale
    m.Data = I1
    m.Weight = W1
    m.Coverage = C1
    m.BeamSize = newbeam



def __writeFits(data=None,outfile='boaMap.fits',overwrite=0,limitsX=[],limitsY=[],intensityUnit="Jy/beam",clip=-1):
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
    info('Exporting map to fits file:')
    print('         %s'%outfile)
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
            info('Clipping map to %.1f*medianRMS (inner contour on display)...'%clip)
            mediannoise = np.nanmedian(rmsMap.Data)
            mask = np.where(rmsMap.Data > clip * mediannoise)
            localMap.Data[mask] = np.NaN
            rmsMap.Data[mask] = np.NaN
            snrMap.Data[mask] = np.NaN
            del mask  # free memory
 
        #write FLux plane                                                            
        localMap._Image__writeImage(dataset, "Intensity", intensityUnit=intensityUnit)
        #write RMS plane
        rmsMap._Image__writeImage(dataset, "Intensity", intensityUnit=intensityUnit+" (RMS)")
        #write SNR plane
        snrMap._Image__writeImage(dataset, "Intensity", intensityUnit='SNR')
        dataset.close()
            
    except Exception, data:
        try:
            dataset.close()
        except:
            pass
        print('Could not write data to file %s: %s' % (outfile, data))
        return
    
    del localMap  # free memory
    del snrMap  # free memory
    del rmsMap  # free memory
    del dataset  # free memory
