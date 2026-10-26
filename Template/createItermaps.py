# =================================
# ==== BEGINNING OF USER INFO =====
# =================================
# NOTE: script intended to be used after reducing all scans and can be used
# for any iteration of the reduction.
# It creates the "mymaps" variable in the MARS session, which is a dictionary
# containing all maps for all iterations. Purpose is whatever user intends.

# ---------------
# Usage examples:
# ---------------
# >>> mymaps['ITER1']['sig'] contains a BoA object: signal map of iteration 1

# display the signal map of iteration 1:
# >>> mymaps['ITER1']['sig'].display(aspect=1)

# display the noise map of iteration 2:
# >>> mymaps['ITER2']['rms'].display(aspect=1)

# display the SNR map of iteration 3:
# >>> mymaps['ITER3']['snr'].display(aspect=1)

# ==============================
# ===== END OF USER INFO =======
# ==============================

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
