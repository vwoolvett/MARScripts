mymaps = {}
for i in range(1, niters+1):
    myname = 'ReducedFiles/%s-%s-%s'%(fe, source, system)
    if flagJumps:
        myname += '-flagJumps'
    myname+='-coadded-flux-iter%i.data'%i
    mymap = restoreFile(myname)
    mymaprms = copy.deepcopy(mymap)
    mymaprms.Data = np.where(mymaprms.Weight>0, 1.0/np.sqrt(mymaprms.Weight), 0)
    mymapcov = copy.deepcopy(mymap)
    mymapcov.Data = mymapcov.Coverage

    mymaps['ITER%i'%i] = {'sig':mymap, 'rms':mymaprms, 'cov':mymapcov}
    del mymap, mymaprms, mymapcov
