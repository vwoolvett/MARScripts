mymaps = {}
for i in range(1, niters+1):
    myname = 'ReducedFiles/%s-%s-%s'%(source, fe, system)
    if flagJumps:
        myname += '-flagJumps'
    myname+='-coadded-flux-iter%i.data'%i
    mymap = restoreFile(myname)
    mymaprms = copy.deepcopy(mymap)
    mymaprms.Data = np.where(mymaprms.Weight>0, 1.0/np.sqrt(mymaprms.Weight), np.NaN)
    mymapsnr  = copy.deepcopy(mymap)
    mymapsnr.Data = mymap.Data/mymaprms.Data
    mymapcov = copy.deepcopy(mymap)
    mymapcov.Data = mymapcov.Coverage

    mymaps['ITER%i'%i] = {'sig':mymap, 'rms':mymaprms, 'snr':mymapsnr, 'cov':mymapcov}
    del mymap, mymaprms, mymapsnr, mymapcov
