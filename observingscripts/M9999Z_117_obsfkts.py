import random

# define map dictionary (edit for each project accordingly). Entries are:
# 'sourceName': [xsize (arcmin), ysize (arcmin), mapangle (degree)]
# if equatorial map, angle is 0, if galactic map contact PI for angle.

# ex: mymapdict = {'source1':[sizeRA, sizeDEC, 0], 'source2':[sizeRA, sizeDEC, 0]}
mymapdict = {}



if len(mymapdict)==0:
    raise RuntimeError('Source/map dictionary is empty! edit the obsfkts.apecs script of the project accordingly.')

def obs(mysource=mymapdict.keys()[0], n=1, dir='xy', tiltangle=15., doCals=True):
    '''
    mysource:       [str], Source name as in source catalog
    n:              [int], Number of loops (kid_center + scanner + X/Y/XY OTF).
                           Each loop is ~ 40 minutes for X+Y. Evaluate number of loops
                           Based on Tsky stability along path of the source.
                           Limited by max 105 minutes without tonelist calibration.
    dir:            [str], scanning direction 'xy', 'x', 'y', default 'xy'.
    tiltangle:    [float], every OTF will be tilted by a random angle from [-tiltangle, +tiltangle]
    doCals:        [bool], Do calibrations before OTFs. Default true. Set to False only
                           if, for example, the second OTF in an X+Y was cancelled, you
                           want to re-do the Y dir, it and you did a kid_center + scanner
                           less than ~40 min ago. 
    '''

    assert str.upper(dir) in ['X', 'Y', 'XY'], "Invalid dir argument. Either 'x', 'y' or 'xy' OTFs allowed."
    

    # Define map
    xlen = float(mymapdict[mysource][0])*60.       # Xsize (arcsec)
    ylen = float(mymapdict[mysource][1])*60.       # Ysize (arcsec)
    sourceang = float(mymapdict[mysource][2])      # Angle (degree)
    arraysize = 15.*60.                          # 15 arcmin

    # Check map is smaller than array
    minsize = min(xlen, ylen)
    if minsize < arraysize:
        msg = 'WARNING: one or more dimensions of the map are smaller than the AMKID array:'
        msg = '\nX length = %.3f arcmin, Y length = %.3f arcmin, AMKID size = %.3f arcmin'%(xlen/60., ylen/60., 15)
        msg +='\nRectangular OTFs might not be the optimal scanning pattern.'
        msg +='\n\nContinue anyways? (y/n):     '
        userInput = str(raw_input(msg))
        if str.upper(userInput) in ['NO', 'N']:
            print('obs() call aborted. Source was not changed!')
            return

    # add array size to get uniform sensitivity within requested map:
    # the entire requested map is visited by the entire array.
    xlen += arraysize
    ylen += arraysize

    # Define speed based on nyquist sampling: 2.5 samples per beam for HFA
    # and LFA will be better than nyquist (not a problem).
    f_samp = 25.                                 # 25 samples / sec nominal
    HFA_fwhm = 7.                                # arcsec, conservative nominal
    oversamp = 2.5
    # oversamp = f_samp / speed * HFA_fwhm [samples / beam]:
    speed = f_samp * HFA_fwhm / oversamp         # is around 70 arcsec / sec
    scanstep = speed * 1.0                       # time in OTF command set to 1.0 sec

    # Define perpendicular step:
    # each OTF should be no longer than ~ 1000 sec, otherwise even chaiten
    # won't have enough memory to reduce them in HFA (confirmed).
    maxScanTime = 1000.                          # sec
    minperpstep = xlen*ylen/(speed*maxScanTime)  # step if we want scan to last 1000 sec

    # Now, if the perpendicular step is greater
    # than the array size there will be missing patches.
    # Moreover if it is greater than 1/2 the array,
    # there will be lines that will only be visited "once".
    # This can happen if the map is too big and we want
    # the OTF to last less than 1000 seconds...
    if minperpstep > arraysize/2 * 0.9:
        perpstep = arraysize/2 * 0.9
        scantime = xlen*ylen/(speed*perpstep)
        msg = 'WARNING: the map is too big to be observed in less than %.1f seconds per OTF.'%maxScanTime
        msg +='\nRequested OTFs will last %.1f seconds and might not be reduceable for HFA.'%scantime
        msg +='\nLFA-only scans can last up to about 2000 seconds.'
        msg +='\n\nContinue anyways? (y/n):     '
        userInput = str(raw_input(msg))
        if str.upper(userInput) in ['NO', 'N']:
            print('obs() call aborted. Source was not changed!')
            return
    else:
        perpstep = minperpstep                   # all good, each OTF will last ~17 min
        scantime = xlen*ylen/(speed*perpstep)

    # Check not more than 90 min without doing pointing + calibrator
    # and checking if new tonelist (power and freq optim.) is needed
    # With nominal 1000s OTFs, this is max n=2 loops.
    otfsperloop = 2 if str.upper(dir)=='XY' else 1
    calibtime = 220.  # sec, 100 fsweep, 120 wirescan
    looptime = otfsperloop * scantime + (calibtime if doCals==True else 0.)  # sec

    if n * looptime / 60 > 90:
        print('Requested number of loops leaves A-MKID without')
        print('pointing + calibrator obs. and without checking')
        print('tone power/frequency optimization for more than')
        print('90 minutes. Reduce number of loops "n".')
        return
    

    # --------------------
    # START OBS - all good
    # --------------------
    source(mysource, cats='user')                  # Center and set source
    
    # loop of OTFs
    for i in range(n):
        print('Starting loop %i/%i...'%(i+1, n))
        # Do calibrations (this is a must before a long science observation!)
        if doCals==True:
            print('Executing kid_center() (Fsweep)')
            kid_center()
            print('Executing scanner() (WireScan)')
            scanner()

        print('Executing %s...'%('X+Y OTFs' if str.upper(dir)=='XY' else '%s OTF'%str.upper(dir)))
        continuous_data('on')
        refcenter()
        
        myang1 = random.uniform(-tiltangle, tiltangle) + sourceang
        if str.upper(dir) == 'X':
            otf(xlen=xlen, ylen=ylen, xstep=scanstep, ystep=perpstep, time=1.0, angle=myang1, direction='x', zigzag=1, size_unit='arcsec', system='EQ')

        elif str.upper(dir) == 'Y':
            otf(xlen=xlen, ylen=ylen, xstep=perpstep, ystep=scanstep, time=1.0, angle=myang1, direction='y', zigzag=1, size_unit='arcsec', system='EQ')

        else:
            otf(xlen=xlen, ylen=ylen, xstep=scanstep, ystep=perpstep, time=1.0 ,angle=myang1, direction='x', zigzag=1, size_unit='arcsec', system='EQ')
            myang2 = random.uniform(-tiltangle, tiltangle) + sourceang
            otf(xlen=xlen, ylen=ylen, xstep=perpstep, ystep=scanstep, time=1.0, angle=myang2, direction='y', zigzag=1, size_unit='arcsec', system='EQ')


