'''Functions to visualize calibrations over time and, e.g., fit F_KID (T_sky)'''

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.stats import sigma_clip


'''Functions for building the main dictionary/table'''
def find_sweep_wire(dir='CalFiles', startscan=0, endscan=int(1e9), fe = 'LFA', keep_digits=False, verbose=True):
    '''
    Returns a list of Fsweep scan numbers and a list of their corresponding WireScan scan numbers, ordered.
    NOTE: RECOMMENDED FOR LFA 2025: startscan=25000, endscan=60000, everything else default.

    @param dir:         directory of sweep and wire scan files. Default = 'CalFiles'.
    @type dir:          str
    @param startscan:   starting scan number, inclusive
    @type startscan:    int
    @param endscan:     ending scan number, inclusive
    @type endscan:      int
    @param fe:          frontend
    @type fe:           str
    @param keep_digits: digits the scan number must have, ordered
    @type keep_digits:  int or False
    @param verbose:     display information or not
    @type verbose:      bool

    '''

    all_entries = os.listdir(dir)

    # filter to get only the correct files
    filename_list = [f for f in all_entries if os.path.isfile(os.path.join(dir, f))]  # only files
    filename_list = [f for f in filename_list if (('fsweep' in f) or ('beam_scanner' in f))]  # only fsweeps or wirescans
    filename_list = [f for f in filename_list if (f[-3:]=='csv') and (fe in f)]  # only csv files with the correct frontend

    scan_list = []
    type_list = []

    for filename in filename_list:
        scan = int(filename.split('_')[-2])

        if scan < startscan or scan > endscan:
            continue

        elif keep_digits==False:
            if 'beam_scanner' in str(filename):
                scan_type = 'wire'
                scan_list.append(scan)
                type_list.append(scan_type)

            elif 'fsweep' in str(filename):
                scan_type = 'sweep'
                scan_list.append(scan)
                type_list.append(scan_type)

            else:
                print('Scan #{} is not an fsweep or beam_scanner scan.'.format(scan))
                break

        else:
            assert type(keep_digits)==int, 'keep_digits must be False or type int'
            if str(keep_digits) in str(scan):
                if 'beam_scanner' in str(filename):
                    scan_type = 'wire'
                    scan_list.append(scan)
                    type_list.append(scan_type)

                elif 'fsweep' in str(filename):
                    scan_type = 'sweep'
                    scan_list.append(scan)
                    type_list.append(scan_type)

                else:
                    print('Scan #{} is not an fsweep or beam_scanner scan.'.format(scan))
                    break
            else:
                if verbose==True:
                    print('Scan #{} ignored: "{}" not in scan number'.format(scan, keep_digits))
                    pass
    
    assert len(scan_list)!=0, 'No scans found in the current directory'

    # sorting the lists
    # Zip the numbers and labels together
    zipped = list(zip(scan_list, type_list))

    # Sort by the number (first item in each tuple)
    zipped_sorted = sorted(zipped)

    # Unzip the sorted pairs back into two lists
    scans_sorted, types_sorted = zip(*zipped_sorted)

    # Convert to lists (optional)
    scan_list = list(scans_sorted)
    type_list = list(types_sorted)

    # detecting pairs
    sweep_scan_list = []
    wire_scan_list = []

    new_scan_list = []
    new_type_list = []

    current = 0
    while current <= len(scan_list)-2:  # can't compare last scan with next ones
        # if the current scan is sweep and the next scan is sweep, this sweep
        # is on masterlist, not used for calibration
        if type_list[current]=='sweep' and type_list[current+1]=='sweep':
            # skip it
            current +=1

        # if current is wire and next is wire, skip
        elif type_list[current]=='wire' and type_list[current+1]=='wire':
            current += 1
            
        # If current is wire and next is sweep, not related
        elif type_list[current]=='wire' and type_list[current+1]=='sweep':
            current +=1
        
        # if the current scan is sweep and the next scan is wire, both are a pair
        # Only if one is right after the other
        elif type_list[current]=='sweep' and type_list[current+1]=='wire' and (scan_list[current+1] == scan_list[current] + 1):
            # add them to the lists
            sweep_scan_list.append(scan_list[current])
            wire_scan_list.append(scan_list[current+1])
            new_scan_list.append([scan_list[current], scan_list[current+1]])
            new_type_list.append([type_list[current], type_list[current+1]])

            # and skip the next
            try:
                dummy = scan_list[current+2]
                current += 2
            except:
                if verbose==True:
                    print('\nScan #{} is the last sweep in directory \"{}\"'.format(scan_list[current], dir))
                break

        # If it's a wire after a sweep but they're not immediate one after the other, skip            
        elif type_list[current]=='sweep' and type_list[current+1]=='wire' and (scan_list[current+1] != scan_list[current] + 1):
            current +=1
        
    if verbose==True:
        print('\n=============================')
        print('========== SUMMARY ==========')
        print('=============================\n')
        print('Number of Fsweep/WireScan scans ({}) found in directory "{}": {}'.format(fe, dir, len(scan_list)))
        print('\nPercentage of sweep+wire scans in directory "{}" used (paired): {:.3f}% ({} out of {})'.format(
            dir, 100.0*(len(sweep_scan_list)+len(wire_scan_list))/len(scan_list),
            (len(sweep_scan_list)+len(wire_scan_list)), len(scan_list)))

    return sweep_scan_list, wire_scan_list


def empty_table():
    '''
    A substitution for Pandas dataframe rows
    '''

    columns=['sweep_scan', 'wire_scan', 'tone_id', 'chain', 'freq', 'responsivity', 'overdriven',
             'freq_max', 'df_max', 'max_response', 'T_sky', 'tone_power', 'Ndf', 'dF2K', 'NET', 'x0', 'y0']
    return {column:[] for column in columns}


def dict_SingleCalib(sweep_scan, wire_scan, fe='LFA'):
    '''
    Builds a dictionary with many KID attributes from the provided sweep and wire scans (one pair).

    @param sweep_scan:  Fsweep scan number
    @type sweep_scan:   int
    @param wire_scan:   Wire scanner scan number
    @type wire_scan:    int
    @param fe:          Frontend, 'LFA' or 'HFA'
    @type fe:           str

    '''

    table = empty_table()

    # full dictionaries, ignoring blindtones in sweep
    full_dict_sweep = readFsweepResult(sweep_scan, fe=fe)
    full_dict_sweep = {kid:full_dict_sweep[kid] for kid in full_dict_sweep.keys() if full_dict_sweep[kid]['type']=='KID'}  # filter blindtones
    full_dict_wire = readScannerResult(wire_scan, fe=fe)

    kids_list = [kid_idx for kid_idx in full_dict_wire.keys()]

    # for every kid:
    for kid_key in kids_list:
        status_sweep = True
        status_wire = True

        try:
            dict_wire = full_dict_wire[kid_key]
        except:
            # no wire data for this KID
            status_wire = False

        try:
            dict_sweep = full_dict_sweep[kid_key]
            # if bad circle fit in this KID (no wire scan counterpart):
            if isinstance(dict_sweep['phData'], int):
                status_sweep = False
        except:
            # no sweep data for this KID
            status_sweep = False

        if not (status_sweep and status_wire):
            continue

        # ----- SWEEP -----
        chain = int(dict_sweep['chain'])
        freq = dict_sweep['freq']
        overdriven = dict_sweep['overdriven']
        df = dict_sweep['df']
        dfs_arr = dict_sweep['dfData']
        phi_arr = dict_sweep['phData']
        responsivity = dict_sweep['responsivity']
        T_sky = dict_sweep['T_sky']
        power = dict_sweep['power']

        dphi_df_arr = (phi_arr[1:] - phi_arr[:-1]) / df
        aux_dfs_arr = (dfs_arr[1:] + dfs_arr[:-1]) / 2

        dphi_df_max = np.max(dphi_df_arr)
        idxmax = np.argmax(dphi_df_arr)
        df_max = aux_dfs_arr[idxmax]
        freq_max = freq + df_max

        dphi_df_arr = np.where(dphi_df_arr > 0, dphi_df_arr, 0)

        # ----- WIRE -----
        Ndf = dict_wire['noise']
        dF2K = dict_wire['phi2K']
        NET = Ndf * dF2K * 1E3
        x0 = dict_wire['x0']
        y0 = dict_wire['y0']

        # ----- APPEND -----
        table['sweep_scan'].append(sweep_scan)
        table['wire_scan'].append(wire_scan)
        table['tone_id'].append(kid_key)
        table['chain'].append(chain)
        table['freq'].append(freq)
        table['responsivity'].append(responsivity)
        table['overdriven'].append(overdriven)
        table['freq_max'].append(freq_max)
        table['df_max'].append(df_max)
        table['max_response'].append(dphi_df_max)
        table['T_sky'].append(T_sky)
        table['tone_power'].append(power)
        table['Ndf'].append(Ndf)
        table['dF2K'].append(dF2K)
        table['NET'].append(NET)
        table['x0'].append(x0)
        table['y0'].append(y0)

    # convert lists into numpy arrays (important!)
    for key in table:
        table[key] = np.array(table[key])

    return table


def concat_tables(table1, table2):
    '''
    Concatenates two tables of the above type
    '''

    if table1 is None:
        return table2

    out = {}
    for key in table1:
        out[key] = np.concatenate([table1[key], table2[key]])
    return out


def build_calDict(sweep_scan=None, wire_scan=None, badscans=[], fe='LFA', sortby=None, verbose=True):
    '''
    Builds a dictionary with many KID attributes from the provided sweep and wire scan lists.

    @param sweep_scan:  Fsweep scan number (or list).
                        Default None: all found in CalFiles by find_sweep_wire() for LFA.
    @type sweep_scan:   None/int/list
    @param wire_scan:   Wire scanner scan number (or list).
                        Default None: all found in CalFiles by find_sweep_wire() for LFA.
    @type wire_scan:    None/int/list
    @param badscans:    Bad scans. E.g., a different masterlist was used.
    @type badscans:     list
    @param fe:          Frontend, 'LFA' or 'HFA'
    @type fe:           str
    @param sortby:      sort by what. Default is first by tone_id and then by sweep_scan.
                        Ex: sortby = ['tone_id', 'sweep_scan']
    @type sortby:       list.
    @param verbose:     whether to display progress in building the dictionary or not.
                        Recommended to keep as True; this function takes long as it must
                        read a lot of reduction files.
    @type verbose:      bool

    '''

    if isinstance(sweep_scan, int) and isinstance(wire_scan, int):
        return dict_SingleCalib(sweep_scan, wire_scan, fe)
    
    elif isinstance(sweep_scan, list) and isinstance(wire_scan, list) and len(sweep_scan)==len(wire_scan):
        pass
    
    elif sweep_scan==None and wire_scan==None:
        info('Using all {} Fsweeps and WireScanners in directory "./CalFiles"'.format(fe))
        sweep_scan, wire_scan = find_sweep_wire()

    else:
        warn('Cannot specify only Fsweep scan or only WireScan scan. Specify both in lists or none.')
        return

    table = None
    total = len(sweep_scan)

    for idx, (sweep, wire) in enumerate(zip(sweep_scan, wire_scan)):
        if sweep in badscans or wire in badscans:
            if verbose:
                warn('Skipping bad scan pair: Fsweep={} | WireScan={}'.format(sweep, wire))
            continue

        # ---- PROGRESS ----
        if verbose:
            info('Processing batch {}/{}: Fsweep={} | WireScan={}'.format(idx+1, total, sweep, wire))

        new_table = dict_SingleCalib(sweep, wire, fe)
        table = concat_tables(table, new_table)

    
    # ----- FILTERING ENTRIES WITH NaNs OR WITHOUT TONE POWER  -----
    mask = np.ones_like(table['sweep_scan'], dtype=bool)  # start with all True
    for key in table:
        if np.issubdtype(table[key].dtype, np.number):
            mask &= ~np.isnan(table[key])  # keep only non-NaN values

    mask &= (table['tone_power'] != -1)
    table = {key: table[key][mask] for key in table}


    # ----- SORTING -----
    if sortby is None:
        sortby = ['tone_id', 'sweep_scan']

    if isinstance(sortby, str):
        sortby = [sortby]

    for col in sortby:
        if col not in table:
            raise ValueError("Column '{}' not found".format(col))
    
    # lexsort expects reversed order
    sort_keys = [table[col] for col in reversed(sortby)]
    idx = np.lexsort(sort_keys)


    for key in table:
        table[key] = table[key][idx]

    return table



'''Main function for writing/loading the main calibration dictionary/table'''
def calDict(sweep_scan=None, wire_scan=None, badscans=[], fe='LFA', mode='load', sortby=None, verbose=True):
    '''
    Builds a calibration dictionary with many KID attributes from the provided sweep and wire scan lists.
    In the default mode="load", it loads an existing file "./CalFiles/calib_dict_FRONTEND.npz".
    If the calibration dictionary file does not exist in "./CalFiles", you should build it and write it
    using mode="write" and choosing the desired parameters.
    See function find_sweep_wire to find suitable sweep-wire pairs.
    NOTE: RECOMMENDED find_sweep_wire parameters FOR LFA 2025: startscan=25000, endscan=60000, everything else default.

    @param sweep_scan:  Fsweep scan number (or list).
                        Default None: all found in CalFiles by find_sweep_wire() for LFA.
    @type sweep_scan:   None/int/list
    @param wire_scan:   Wire scanner scan number (or list).
                        Default None: all found in CalFiles by find_sweep_wire() for LFA.
    @type wire_scan:    None/int/list
    @param badscans:    Bad scans. E.g., a different masterlist was used.
    @type badscans:     list
    @param fe:          Frontend, 'LFA' or 'HFA'
    @type fe:           str
    @param mode:        whether to try loading or to write (build) the table from scratch.
                        Writing overwrites the saved file (if existing) in the process.
                        If mode='load' and a master file is found, all other parameters are ignored.
    @type mode:         str, 'load' or 'write'.
    @param sortby:      sort by what. Default is first by tone_id and then by sweep_scan.
                        Ex: sortby = ['tone_id', 'sweep_scan']
    @type sortby:       list.
    @param verbose:     whether to display progress in building the dictionary or not.
                        Recommended to keep as True; this function takes long as it must
                        read a lot of reduction files.
    @type verbose:      bool

    '''

    assert mode in ['load', 'write'], 'Parameter {} must be either "load" or "write"!'.format(mode)
    assert str.upper(fe) in ['LFA', 'HFA'],'Parameter {} must be either "LFA" or "HFA"!'.format(fe)

    if mode=='write':
        info('Building and saving calibration dictionary/table from provided scans and parameters')
        table = build_calDict(sweep_scan=sweep_scan, wire_scan=wire_scan, badscans=badscans, fe=fe, sortby=sortby, verbose=verbose)
        np.savez_compressed("CalFiles/calib_dict_{}.npz".format(fe), **table)
        info('Calibration dictionary file written to "./CalFiles/calib_dict_{}.npz"'.format(fe))
        return table
    
    else:
        if verbose:
            info('Loading calibration dictionary from "./CalFiles/calib_dict_{}.npz"'.format(fe))
        try:
            data = np.load("CalFiles/calib_dict_{}.npz".format(fe), allow_pickle=True)
            table_loaded = {key: data[key] for key in data}
            return table_loaded
        except:
            warn('Failed to load calibration dictionary for {}, file "./CalFiles/calib_dict_{}.npz" not found.'.format(fe, fe))
            warn('Build the calibration dictionary with function calDict(mode="write") and according to desired parameters!')
            return None



'''Functions for visualizing the table'''
def plot_CalibKids(calibdict, xaxis='T_sky', yaxis='freq_max', color_by='NET', startkid=None, endkid=None, interactive=True):
    '''
    Plots parameter yaxis vs xaxis from calibrations in calibdict (see function calDict) interactively
    for all kids.

    @param calibdict:       Calibration dictionary, built or loaded from function calDict.
    @type calibdict:        dictionary of type attribute:array
    @param xaxis:           attribute from calibdict to plot in the X axis. Attribites include tone_id, freq, freq_max, T_sky, etc.
    @type xaxis:            str
    @param yaxis:           attribute from calibdict to plot in the Y axis. Attribites include tone_id, freq, freq_max, T_sky, etc.
    @type yaxis:            str
    @param color_by:         attribute from calibdict to color the points. Attribites include tone_id, freq, freq_max, T_sky, etc.
    @type color_by:          str
    @param startkid:        initial kid to start plotting. Default None starts from KID 1
    @type startkid:         int
    @param endkid:          final kid to plot. Default None ends in last KID with data in calibdict.
    @type endkid:           int
    @param interactive:     whether to plot interactivelt or not. Value False only plots startkid.
    @type interactive:      bool 

    '''

    assert xaxis in calibdict.keys() and yaxis in calibdict.keys(),\
           'Arguments "xaxis" and "yaxis" must be keys of the calibration dictionary. Possible values are {}'.format(calibdict.keys())
    
    kidlist = np.sort(np.unique(calibdict['tone_id']))
    if startkid==None:
        startkid = kidlist[0]

    if endkid==None:
        endkid = kidlist[-1]

    currentkididx = np.where(kidlist==startkid)[0][0]
    currentkid = int(kidlist[currentkididx])

    #Message to interact with user
    msg = '\n------------------------------------------------\n'
    msg += 'Plot next KID:                            <Enter>\n'
    msg += 'Plot previous KID:                      - <Enter>\n'
    msg += 'Go to KID Number:                <number> <Enter>\n'
    msg += 'Exit script execution:                  q <Enter>\n'
    msg += '-------------------------------------------------\n'
    msg += 'Enter choice: '

    fig, ax = plt.subplots(1)
    fig.set_size_inches(10, 8)
    while True:
        if currentkid>=startkid and currentkid<=endkid:

            mask = (calibdict['tone_id']==currentkid)  # select kid
            dict_kid = {key: calibdict[key][mask] for key in calibdict}

            vmin, vmax = np.percentile(dict_kid[color_by], 1), np.percentile(dict_kid[color_by], 95)
            sc = ax.scatter(dict_kid[xaxis], dict_kid[yaxis], marker='o', s=8, c=dict_kid[color_by], vmin=vmin, vmax=vmax, cmap='brg')
            cbar = fig.colorbar(sc, ax=ax, orientation='horizontal')
            cbar.set_label(color_by)
            ax.set_xlabel(xaxis)
            ax.set_ylabel(yaxis)
            ax.set_title('%s vs %s for KID %i'%(yaxis, xaxis, currentkid))
            ax.ticklabel_format(axis='both', style='plain', useOffset=False)
            plt.show()
            
            if interactive==False:
                return

            inp = str(raw_input(msg))
            if inp == '':
                currentkididx += 1
                currentkid = int(kidlist[currentkididx])
            elif inp =='-':
                currentkididx -=1
                currentkid = int(kidlist[currentkididx])
            elif inp =='q':
                plt.close(fig)
                break
            
            else:
                try:
                    if int(inp) in kidlist:
                        currentkid = int(inp)
                        currentkididx = np.where(kidlist==currentkid)[0][0]
                    else:
                        print('\n=================================================================================')
                        print('Value "{}" entered not in kid list.'.format(inp))
                        print('=================================================================================')
                except:
                    print('\n=================================================================================')
                    print('Value "{}" entered is not valid.'.format(inp))
                    print('=================================================================================')
                
        else:
            print('\n=================================================================================')
            print('Reached maximum or minimum kid. Going back to startkid: {}'.format(startkid))
            print('=================================================================================')
            currentkid = startkid
            currentkididx = np.where(kidlist==currentkid)[0][0]

        ax.cla()
        cbar.remove()



def plot_NET_Tsky(calibdict, mode='median'):
    '''
    Plots the median, mode, or mean NET of the array as a function of T_sky based on the data in calibdict.
    @param calibdict:       Calibration dictionary, built or loaded from function calDict.
    @param mode:            The mode of the NET to plot ('median', 'mode', or 'mean').
    @type mode:             str
    '''
    assert mode in ['median', 'mode', 'mean'], 'Argument "mode" must be "median", "mode", or "mean"'

    if mode=='median':
        net_mode = 'Median'
        net_func = np.nanmedian
    elif mode=='mode':
        net_mode = 'Mode'
        net_func = lambda x: np.bincount(x.astype(int)).argmax()  # simple mode function for integers
    else:
        net_mode = 'Mean'
        net_func = np.nanmean

    wirescan_list = np.sort(np.unique(calibdict['wire_scan']))
    tsky_scan_list = []
    net_metric_list = []

    for wirescan in wirescan_list:
        mask = (calibdict['wire_scan']==wirescan)
        dict_scan = {key: calibdict[key][mask] for key in calibdict}
        net_metric = net_func(dict_scan['NET'])
        net_metric_list.append(net_metric)
        tsky_scan_list.append(dict_scan['T_sky'][0])  # all points in scan have same T_sky

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(tsky_scan_list, net_metric_list)
    ax.set_xlabel('T_sky (K)')
    ax.set_ylabel('%s NET (mK*sqrt(s))' % net_mode)
    ax.set_ylim(0, 20)
    ax.set_title('%s NET vs T_sky' % net_mode)
    ax.grid()

    

'''Function for skydip fit from year-long calibrations'''
def skydip(calibdict, max_NET=1e9, max_dF2K=10, weight_by='df_max', punish='bigger', outliers='discard', decrosstalk=True, diameter=10., sig=1., full_output=False, inspect=False):
    '''
    Fits the skydip: F_max(T_sky) = m T_sky + n parameters for all KIDs based on the data found in calibration
    dictionary calibdict. Allows for maximum NET and dF2K filtering for the fitting, and weighting the points by any
    attribute, punishing bigger or smaller values of it for each point.

    Returns a dictionary with fitted parameters, status of the fit, and other relevant parameters. 
    If parameter full_output==True, then it returns (0) the dictionary returns, as well as the lists of
    (1) KIDs with <30 points, (2) KIDs with high NET, (3) KIDs with high dF2K, 
    (4) KIDs with not enough points to fit, and (5) KIDs where the fit is bad (positive slope or high RMS)

    NOTE: In the output dictionary, frequency units are MHz and tone_power units are dBm

    @param calibdict:       Calibration dictionary, built or loaded from function "calDict".
    @type calibdict:        dictionary of type attribute:array
    @param max_NET:         Maximum NET to use points for fitting.
    @type max_NET:          float
    @param max_dF2K:        Maximum dF2K to use points for fitting. Extremely high dF2K values are when KIDs are
                            non-responsive or their calibration circle was way off.
    @type max_dF2K:         float
    @param weight_by:       calibdict attribute to weigh the points when fitting F_max(T_sky). Value None weighs
                            all points equally, but may lead to bad fits due to degeneracy in some attributes.
                            Default "df_max" is the f_max - f_tone; smaller values of abs(df_max) indicate that
                            we've put the tones closer to the actual responsivity peak. Weighing by this attribute
                            usually recovers the most KIDs.
    @type weight_by:        str
    @param punish:          whether to punish "smaller" or "bigger" values of attribute weight_by.
    @type punish:           str, 'smaller' or 'bigger'
    @param outliers:        keep outliers from the fit ('keep') or not ('discard')
    @type outliers:         str
    @param decrosstalk:     whether to try fitting F_max(T_sky) only with points near the modal
                            (x, y) position of the kid in the array when the fit with all points fails.
                            Default True recovers most cross-talking KIDs.
    @type decrosstalk:      bool
    @param diameter:        diameter of the circle around modal (x, y) bin to consider when decrosstalk=True.
                            Default is permissive enough to get al "good" points around modal position without
                            taking points from immediately adjacent KIDs in the physical plane of the array.
    @type diameter:         float
    @param sig:             for each KID, how many tone_power standard deviations away from mean tone_power to ignore, default=1.
                            This is a necessary filter, as the KID moves in a different F_max(T_sky) line for each tone power.
                            A big sigma (say, 1e9) eliminates this filter completely.
    @type sig:              float
    @param full_output:     if False, only returns the dictionary with fitting parameters. Otherwise, it returns (0) the dictionary,
                            and also returns the lists of (1) KIDs with <30 points, (2) KIDs with high NET, (3) KIDs with high dF2K,
                            (4) KIDs with not enough points to fit, and (5) KIDs where the fit is bad (positive slope or high RMS).
    @type full_output:      bool
    @param inspect:         whether to plot the fitting results interactively for all KIDs.
    @type inspect:          bool

    '''
        
    assert type(diameter)==float or type(diameter)==int, '"diameter" must be float or int!'
    assert punish=='bigger' or punish=='smaller', 'Argument "punish" must be "bigger" or "smaller" for values of attribute "weight_by"!'
    assert type(full_output)==bool, 'Argument "full_output" must be a bool!'
    assert outliers=='discard' or outliers=='keep', 'Argument "outliers" must be "discard" or "keep"'
    
    # Pixel size for histograms when separating by position (decrosstalk=True):
    # About 8 KIDs fit in 50 mm so one KID every ~6.3 mm.
    # Sample these scales with 2 pixels then pixel size ~ 2 mm
    # threshold = search circle radius
    # pixsize "large enough" to get mode~mean position in most populated position
    pixsize = 2  # mm

    # search radius from mode if separating by position
    threshold = diameter/2  

    if weight_by != None:
        assert weight_by in calibdict.keys(), 'Argument {} is not in the calibration dictionary. Possible values are {}'.format(weight_by, calibdict.keys())

    kidlist = np.sort(np.unique(calibdict['tone_id']))

    fitdict = dict.fromkeys(kidlist)

    info('Fitting slope and intercept: F_max(T_sky) = m T_sky + n ...')

    
    for kid in kidlist:
        # select kid
        mask = (calibdict['tone_id']==kid)
        dict_kid = {key: calibdict[key][mask] for key in calibdict}
        
        # Initial status
        t_power = np.nan
        m = np.nan
        n = np.nan
        f60K = np.nan
        wrms = np.nan
        dFdT = np.nan
        x0_pos = np.nan
        y0_pos = np.nan
        separate = 'no'
        npts = len(dict_kid['NET'])
        outliers_complex = []

        # Must have at least 30 points
        if npts<30:
            status = 'bad_less30pts'

        # If half or more or the points have NET>max_NET there's no point
        # in fitting.
        elif sum(dict_kid['NET'] > max_NET) > len(dict_kid['NET'])/2:
            npts = sum(dict_kid['NET'] > max_NET)
            status = 'bad_NET'

        # If half or more or the points have dF2K>max_dF2K there's no point
        # in fitting either.
        elif sum(dict_kid['dF2K'] > max_dF2K) > len(dict_kid['dF2K'])/2:
            npts = sum(dict_kid['dF2K'] > max_dF2K)
            status = 'bad_dF2K'

        else:
            # cut points to those with lower NET than max_NET, lower dF2K than max_dF2K, and that have positions!
            mask = (dict_kid['NET']<=max_NET) & (dict_kid['dF2K']<=max_dF2K)
            mask &= (dict_kid['x0']>-1000) & (dict_kid['y0']>-1000)
            dict_kid = {key: dict_kid[key][mask] for key in dict_kid}

            # change number of points, filtered rubbish
            npts = len(dict_kid['NET'])

            # relevant arrays
            tsky_arr = dict_kid['T_sky']
            fmax_arr = dict_kid['freq_max']
            tp_arr = dict_kid['tone_power']
            dF2K_arr = dict_kid['dF2K']

            # define direct_arr and weight_arr:
            # if no weights, all are equal
            if weight_by == None:
                weight_arr = np.ones(npts)
            # otherwise, conditions
            else:
                direct_arr = np.array(dict_kid[weight_by])

                # if used, judge only by abs(max-tone)
                if weight_by=='df_max':
                    direct_arr = np.abs(direct_arr)

                if punish=='bigger':
                    # punish high values means weight propto 1/(value)^2
                    weight_arr = 1 / direct_arr ** 2

                else:
                    # punish low values means weight propto value^2
                    weight_arr = direct_arr ** 2

            # Calculate and save positions!
            # Histogram
            xs = dict_kid['x0']
            ys = dict_kid['y0']
            xmin, xmax = np.min(xs)-threshold, np.max(xs)+threshold  # in case mode is in border
            ymin, ymax = np.min(ys)-threshold, np.max(ys)+threshold  # in case mode is in border
            binsx = int((xmax-xmin)/pixsize)
            binsy = int((ymax-ymin)/pixsize)
            if binsx == 0:
                binsx = 1
            if binsy == 0:
                binsy = 1
                    
            # 2D histogram
            H, xedges, yedges = np.histogram2d(xs, ys, bins=[binsx, binsy], range=[[xmin, xmax], [ymin, ymax]])

            # Index of maximum bin
            imax = np.unravel_index(np.argmax(H), H.shape)
            ix, iy = imax

            # Bin centers
            x_mode = 0.5 * (xedges[ix] + xedges[ix+1])
            y_mode = 0.5 * (yedges[iy] + yedges[iy+1])

            # mask near modal position
            nearmode = ((dict_kid['x0']-x_mode)**2 + (dict_kid['y0']-y_mode)**2 <= threshold**2)

            # 2D histogram
            H_2, xedges_2, yedges_2 = np.histogram2d(xs[~nearmode], ys[~nearmode], bins=[binsx, binsy], range=[[xmin, xmax], [ymin, ymax]])

            # Index of second maximum bin
            imax_2 = np.unravel_index(np.argmax(H_2), H_2.shape)
            ix_2, iy_2 = imax_2

            # Bin centers
            x_mode_2 = 0.5 * (xedges_2[ix_2] + xedges_2[ix_2+1])
            y_mode_2 = 0.5 * (yedges_2[iy_2] + yedges_2[iy_2+1])

            # mask near modal position
            nearmode_2 = ((dict_kid['x0']-x_mode_2)**2 + (dict_kid['y0']-y_mode_2)**2 <= threshold**2)

            # if the ring around the mode does not contain more points that the ring around the next modal
            # bin, then it's not the mode
            if np.sum(nearmode_2)>sum(nearmode):
                # Use second mode as central bin
                nearmode = nearmode_2
                x_mode, y_mode = x_mode_2, y_mode_2

            # SAVE MEAN OF X AND Y POSITIONS OF KID AROUND MODAL BIN!
            x0_pos = np.mean(xs[nearmode])
            y0_pos = np.mean(ys[nearmode])

            # mask near typical tone_power
            nearmeantp = (np.abs(tp_arr-np.mean(tp_arr)) <= sig * np.std(tp_arr))

            # compute the mean tone power without the outliers
            # this would be saved if the fit converges and the quality criteria are met
            t_power = np.mean(tp_arr[nearmeantp])
            
            # number of points that meet the criterium
            npts = np.sum(nearmeantp)

            # Try to fit the line to ALL points near mean TP
            # Without separating by wirescanner position firtst (AT LEAST 20 POINTS)
            if npts >= 20:
            
                m, n = np.polyfit(tsky_arr[nearmeantp], fmax_arr[nearmeantp],
                                  deg=1, w=np.sqrt(weight_arr[nearmeantp]))

                residuals = fmax_arr[nearmeantp] - (m * tsky_arr[nearmeantp] + n)
                wrms = np.sqrt(np.sum(weight_arr[nearmeantp] * residuals**2) /
                               np.sum(weight_arr[nearmeantp]))
                
                # calculate mean slope according to the selected KIDs' dF2K
                dFdT = - np.mean(1/dF2K_arr[nearmeantp]) / 1e3  # from KHz/K to MHz/K, same units as slope m

                # calculate KID freq. at 60 K Tsky
                f60K = m*60 + n  # MHz
                
                if outliers=='discard':
                    # re-fit without outliers in residuals
                    sclip = sigma_clip(residuals, sigma=3)
                    outlier = sclip.mask
                    
                    if sum(outlier)>0:
                        outlier_idx = np.where(outlier)[0]
                        outliers_complex = [tsky_arr[nearmeantp][i] + 1j*fmax_arr[nearmeantp][i]for i in outlier_idx]

                    if sum(~outlier) >= 20:
                        npts = sum(~outlier)

                        m, n = np.polyfit(tsky_arr[nearmeantp][~outlier], fmax_arr[nearmeantp][~outlier],
                                          deg=1, w=np.sqrt(weight_arr[nearmeantp][~outlier]))

                        residuals = fmax_arr[nearmeantp][~outlier] - (m * tsky_arr[nearmeantp][~outlier] + n)
                        wrms = np.sqrt(np.sum(weight_arr[nearmeantp][~outlier] * residuals**2) /
                                       np.sum(weight_arr[nearmeantp][~outlier]))

                        # calculate mean slope according to the selected KIDs' dF2K
                        dFdT = - np.mean(1/dF2K_arr[nearmeantp][~outlier]) / 1e3  # from KHz/K to MHz/K, same units as slope m

                        # calculate KID freq. at 60 K Tsky
                        f60K = m*60 + n  # MHz

                        # if resonance moves down in freq with increasing Tsky and the fit has
                        # weighted uncertainty of less than 20 KHz:
                        converged = m<0 and wrms*1e3 <= 20.

                    else:
                        npts = sum(~outlier)
                        # not enough points
                        converged = False

                else:
                    # if resonance moves down in freq with increasing Tsky and the fit has
                    # weighted uncertainty of less than 20 KHz
                    converged = m<0 and wrms*1e3 <= 20. 
                
                if converged:
                    status = 'good'

                # if it did not fit well and don't want to try separating by KID position:
                elif not converged and decrosstalk==False:
                    status = 'bad_fit'

                # if it did not fit well and want to try separating by KID position (take only near modal position)
                # NOW WE REQUIRE AT LEAST 15 POINTS, INSTEAD OF 20, AND ALLOW HIGHER WEIGHTED RMS
                else:
                    separate = 'yes'

                    # Histogram of position in array
                    #xs = dict_kid['x0']
                    #ys = dict_kid['y0']
                    #xmin, xmax = np.min(xs)-threshold, np.max(xs)+threshold  # in case mode is in border
                    #ymin, ymax = np.min(ys)-threshold, np.max(ys)+threshold  # in case mode is in border
                    #binsx = int((xmax-xmin)/pixsize)
                    #binsy = int((ymax-ymin)/pixsize)
                    #if binsx == 0:
                    #    binsx = 1
                    #if binsy == 0:
                    #    binsy = 1

                    # 2D histogram
                    #H, xedges, yedges = np.histogram2d(xs, ys, bins=[binsx, binsy], range=[[xmin, xmax], [ymin, ymax]])

                    # Index of maximum bin
                    #imax = np.unravel_index(np.argmax(H), H.shape)
                    #ix, iy = imax

                    # Bin centers
                    #x_mode = 0.5 * (xedges[ix] + xedges[ix+1])
                    #y_mode = 0.5 * (yedges[iy] + yedges[iy+1])

                    # mask near modal position
                    #nearmode = ((dict_kid['x0']-x_mode)**2 + (dict_kid['y0']-y_mode)**2 <= threshold**2)

                    # 2D histogram without first modal
                    #H_2, xedges_2, yedges_2 = np.histogram2d(xs[~nearmode], ys[~nearmode], bins=[binsx, binsy], range=[[xmin, xmax], [ymin, ymax]])

                    # Index of second maximum bin
                    #imax_2 = np.unravel_index(np.argmax(H_2), H_2.shape)
                    #ix_2, iy_2 = imax_2

                    # Bin centers
                    #x_mode_2 = 0.5 * (xedges_2[ix_2] + xedges_2[ix_2+1])
                    #y_mode_2 = 0.5 * (yedges_2[iy_2] + yedges_2[iy_2+1])

                    # mask near modal position
                    #nearmode_2 = ((dict_kid['x0']-x_mode_2)**2 + (dict_kid['y0']-y_mode_2)**2 <= threshold**2)

                    # if the ring around the mode does not contain more points that the ring around the next modal
                    # bin, then it's not the mode
                    #if np.sum(nearmode_2)>np.sum(nearmode):
                    #    # Use second mode as central bin
                    #    nearmode = nearmode_2
                    #    x_mode, y_mode = x_mode_2, y_mode_2


                    # tone power points near mode
                    tp_arr_nearmode = tp_arr[nearmode]

                    # mask near mean tone_power around modal position for unfiltered array
                    nearmeantp_global = (np.abs(tp_arr-np.mean(tp_arr_nearmode)) <= sig * np.std(tp_arr_nearmode))

                    # final mask, both near the mode and near the modal position's mean TP
                    accepted = nearmode & nearmeantp_global

                    # compute the mean tone power of points around modal position without the outliers in TP
                    # this would be saved if the fit converges and the quality criteria are met
                    t_power = np.mean(tp_arr[accepted])

                    # number of points that meet the criterium
                    npts = np.sum(accepted)

                    # Try to fit the line to all data near mean TP around modal position
                    # More permissive fit than if not separating by position!
                    if npts >= 15:
                    
                        m, n = np.polyfit(tsky_arr[accepted], fmax_arr[accepted],
                                                   deg=1, w=np.sqrt(weight_arr[accepted]))

                        residuals = fmax_arr[accepted] - (m * tsky_arr[accepted] + n)
                        wrms = np.sqrt(np.sum(weight_arr[accepted] * residuals**2) /
                                                np.sum(weight_arr[accepted]))
                        
                        # calculate mean slope according to the selected KIDs' dF2K
                        dFdT = - np.mean(1/dF2K_arr[accepted]) / 1e3  # from KHz/K to MHz/K, same units as slope m

                        # calculate KID freq. at 60 K Tsky
                        f60K = m*60 + n  # MHz
                
                        if outliers=='discard':
                            # re-fit without outliers in residuals
                            sclip = sigma_clip(residuals, sigma=3)
                            outlier = sclip.mask

                            if sum(outlier)>0:
                                outlier_idx = np.where(outlier)[0]
                                outliers_complex = [tsky_arr[accepted][i] + 1j*fmax_arr[accepted][i]for i in outlier_idx]

                            if sum(~outlier) >= 15:
                                npts = sum(~outlier)

                                m, n = np.polyfit(tsky_arr[accepted][~outlier], fmax_arr[accepted][~outlier],
                                                  deg=1, w=np.sqrt(weight_arr[accepted][~outlier]))

                                residuals = fmax_arr[accepted][~outlier] - (m * tsky_arr[accepted][~outlier] + n)
                                wrms = np.sqrt(np.sum(weight_arr[accepted][~outlier] * residuals**2) /
                                               np.sum(weight_arr[accepted][~outlier]))

                                # calculate mean slope according to the selected KIDs' dF2K
                                dFdT = - np.mean(1/dF2K_arr[accepted][~outlier]) / 1e3  # from KHz/K to MHz/K, same units as slope m

                                # calculate KID freq. at 60 K Tsky
                                f60K = m*60 + n  # MHz

                                # if resonance moves down in freq with increasing Tsky and the fit has
                                # weighted uncertainty of less than 70 KHz (more permissive):
                                converged = m<0 and wrms*1e3 <= 70

                            else:
                                npts = sum(~outlier)
                                # keep previous fit but note that after removing outliers there's not enough points
                                converged = False

                        else:
                            # if resonance moves down in freq with increasing Tsky and the fit has
                            # weighted uncertainty of less than 70 KHz (more permissive):
                            converged = m<0 and wrms*1e3 <= 70

                        if converged:
                            status = 'good'
                        else:
                            status = 'bad_fit'

                    # Not enough points after separating by position
                    else:
                        t_power = np.nan
                        m = np.nan
                        n = np.nan
                        f60K = np.nan
                        wrms = np.nan
                        dFdT = np.nan
                        x0_pos = np.nan
                        y0_pos = np.nan
                        status = 'bad_fewpts'

            # Not enough points to fit even in first stage
            # after filtering near TP and not separating
            else:
                t_power = np.nan
                m = np.nan
                n = np.nan
                f60K = np.nan
                wrms = np.nan
                dFdT = np.nan
                x0_pos = np.nan
                y0_pos = np.nan
                status = 'bad_fewpts'

        fitdict[kid] = {'tone_power':t_power, 'm':m, 'n':n, 'f60K':f60K, 'rms':wrms, 'dF/dT':dFdT, 'x0_mean':x0_pos, 'y0_mean':y0_pos, 'separate':separate, 'npts':npts, 'status':status, 'outliers':outliers_complex}

    ntotal = len(kidlist)
    nvalid = len([kid for kid in kidlist if fitdict[kid]['status']=='good'])
    ninvalid = ntotal - nvalid

    n_less30 = len([kid for kid in kidlist if fitdict[kid]['status']=='bad_less30pts'])
    n_net = len([kid for kid in kidlist if fitdict[kid]['status']=='bad_NET'])
    n_df2k = len([kid for kid in kidlist if fitdict[kid]['status']=='bad_dF2K'])
    n_fewpts = len([kid for kid in kidlist if fitdict[kid]['status']=='bad_fewpts'])
    n_badfit = len([kid for kid in kidlist if fitdict[kid]['status']=='bad_fit'])

    print('')
    info('Finished fitting: %i out of %i KIDs were succesfully characterized'%(nvalid, ntotal))
    print('')
    if max_NET > 9999:
        if max_dF2K > 9999:
            info('Of the %i bad KIDs:\n- %i have less than 30 points before filtering by anything\n- %i have half or more of their points with NET > %.1e mK sqrt(s)\n- %i have half or more of their points with dF2K > %.1e K/KHz\n- %i don\'t have enough valid points for fitting\n- %i have bad fits but should be inspected!'\
                 %(ninvalid, n_less30, n_net, max_NET, n_df2k, max_dF2K, n_fewpts, n_badfit))
        else:
            info('Of the %i bad KIDs:\n- %i have less than 30 points before filtering by anything\n- %i have half or more of their points with NET > %.1e mK sqrt(s)\n- %i have half or more of their points with dF2K > %.1f K/KHz\n- %i don\'t have enough valid points for fitting\n- %i have bad fits but should be inspected!'\
                 %(ninvalid, n_less30, n_net, max_NET, n_df2k, max_dF2K, n_fewpts, n_badfit))
    else:
        if max_dF2K > 9999:
            info('Of the %i bad KIDs:\n- %i have less than 30 points before filtering by anything\n- %i have half or more of their points with NET > %.1f mK sqrt(s)\n- %i have half or more of their points with dF2K > %.1e K/KHz\n- %i don\'t have enough valid points for fitting\n- %i have bad fits but should be inspected!'\
                 %(ninvalid, n_less30, n_net, max_NET, n_df2k, max_dF2K, n_fewpts, n_badfit))
        else:
            info('Of the %i bad KIDs:\n- %i have less than 30 points before filtering by anything\n- %i have half or more of their points with NET > %.1f mK sqrt(s)\n- %i have half or more of their points with dF2K > %.1f K/KHz\n- %i don\'t have enough valid points for fitting\n- %i have bad fits but should be inspected!'\
                 %(ninvalid, n_less30, n_net, max_NET, n_df2k, max_dF2K, n_fewpts, n_badfit))
        
        

    # Interactive plotting if inspect=True to visualize fitting:
    if inspect:
        print('')
        info('Initiating interactive plot...')
        #Message to interact with user
        msg = '\n---------------------------------------------------\n'
        msg += 'Plot next KID:                               <Enter>\n'
        msg += 'Plot previous KID:                         - <Enter>\n'
        msg += 'Go to KID Number:                   <number> <Enter>\n'
        msg += 'Cycle through "bad" KIDs:                  b <Enter>\n'
        msg += 'Exit script execution:                     q <Enter>\n'
        msg += '----------------------------------------------------\n'
        msg += 'Enter choice: '

        if weight_by==None:
            coloring=None

        else:
            if weight_by=='df_max':
                coloring = np.abs(calibdict[weight_by])
            else:
                coloring = calibdict[weight_by]
                
            vmin = np.percentile(coloring, 1)
            vmax = np.percentile(coloring, 95)

            if punish=='bigger':
                # green = less weight
                cmap = 'brg'
            else:
                # green = less weight
                cmap = 'brg_r'

        # plotting
        fig, ax = plt.subplots(1, 2)
        fig.set_size_inches(15, 8)

        pad = 10.

        currentkididx = 0
        currentkid = kidlist[0]

        badkidlist = [kid for kid in fitdict.keys() if ('bad' in fitdict[kid]['status'])]
        currentbadkididx = -1  # necessary, otherwise won't display the first bad kid when inp=='b'

        while True:
            # clear variables
            cbar = None
            cbarhist = None

            if (currentkid < kidlist[0]) or (currentkid > kidlist[-1]):
                print('\n=================================================================================')
                print('Reached maximum or minimum kid. Going back to KID 1')
                print('=================================================================================')
                currentkididx = 0
                currentkid = kidlist[0]

            # select kid
            mask = (calibdict['tone_id']==currentkid)
            dict_kid = {key: calibdict[key][mask] for key in calibdict}

            # plot previous tone placements in the back of axis 0
            mintsky, maxtsky = np.min(dict_kid['T_sky']), np.max(dict_kid['T_sky'])
            ax[0].hlines(np.unique(dict_kid['freq']), mintsky, maxtsky, linestyles=['--'], colors=['grey'], label='Previous tone placements')

            # define needed points to pass filter:
            if fitdict[currentkid]['status'] == 'bad_less30pts':
                n_needed = 30
            
            elif fitdict[currentkid]['status'] == 'bad_NET':
                n_needed = int(len(fitdict[currentkid]['NET'])/2)

            elif fitdict[currentkid]['status'] == 'bad_dF2K':
                n_needed = int(len(fitdict[currentkid]['dF2K'])/2)

            else:
                # passed filters, now only depends on whether we separated by position or not
                # either status is 'bad_fewpts', 'good' or 'bad_fit', but for all
                # npts_needed is the same.
                if fitdict[currentkid]['separate'] == 'no':
                    n_needed = 20
                else:
                    n_needed = 15


            # if KID is bad but not because of a bad fit, plot it colored by NET and see if it's worth repairing:
            if 'bad' in fitdict[currentkid]['status'] and not('fit' in fitdict[currentkid]['status']):
                coloring_aux = dict_kid['NET']
                vmin_aux = np.percentile(calibdict['NET'], 1)
                vmax_aux = np.percentile(calibdict['NET'], 95)
                sc = ax[0].scatter(dict_kid['T_sky'], dict_kid['freq_max'], marker='o', s=8, c=coloring_aux, vmin=vmin_aux, vmax=vmax_aux, cmap='brg')
                ax[0].set_ylim(np.mean(dict_kid['freq_max']) - (np.max(calibdict['freq_max'])-np.min(calibdict['freq_max']))/len(kidlist) / 2 - 0.3,
                               np.mean(dict_kid['freq_max']) + (np.max(calibdict['freq_max'])-np.min(calibdict['freq_max']))/len(kidlist) / 2)
                cbar = fig.colorbar(sc, ax=ax[0], orientation='horizontal')
                cbar.set_label('NET')
                if 'fewpts' in fitdict[currentkid]['status']:
                    fig.suptitle('F(Tsky) for KID %i: BAD! (status: %s | separate: %s | npts: %i out of %i needed to fit)'%\
                                 (currentkid, fitdict[currentkid]['status'], fitdict[currentkid]['separate'], fitdict[currentkid]['npts'], n_needed),
                                 fontsize=15)
                else:
                    fig.suptitle('F(Tsky) for KID %i: BAD! (status: %s | separate: %s | npts: %i out of %i needed to pass filter)'%\
                                 (currentkid, fitdict[currentkid]['status'], fitdict[currentkid]['separate'], fitdict[currentkid]['npts'], n_needed),
                                 fontsize=15)
                ax[1].scatter(0, 0, s=0.001)
                ax[1].text(0, 0, 'BAD KID! Is it worth saving?\nColor coding is now NET', fontdict={'fontsize':15}, horizontalalignment='center')

            # could fit, but can still be bad_fit or good
            else:
                # cut points to those with lower NET than max_NET, lower dF2K than max_dF2K, and that have positions!
                mask = (dict_kid['NET']<=max_NET) & (dict_kid['dF2K']<=max_dF2K)
                mask &= (dict_kid['x0']>-1000) & (dict_kid['y0']>-1000)
                dict_kid = {key: dict_kid[key][mask] for key in dict_kid}

                # relevant arrays
                tsky_arr = dict_kid['T_sky']
                fmax_arr = dict_kid['freq_max']
                tp_arr = dict_kid['tone_power']

                # mask near mean tone_power
                nearmeantp = (np.abs(tp_arr-np.mean(tp_arr)) <= sig * np.std(tp_arr))

                if weight_by=='df_max':
                    direct_arr = np.abs(dict_kid[weight_by])
                elif weight_by!='df_max' and weight_by!=None:
                    direct_arr = dict_kid[weight_by]

                # fit parameters
                m, n, wrms = fitdict[currentkid]['m'], fitdict[currentkid]['n'], fitdict[currentkid]['rms']
                fitlabel = r'F$_{peak}$(T$_{sky}$) = ' + '(%.3f KHz/K)'%(m*1e3) + r'T$_{sky}$ +' + '%.3f MHz'%n
                fitlabel += ' (RMS = %.1f KHz)'%(wrms*1e3)

                # Histogram
                xs = dict_kid['x0']
                ys = dict_kid['y0']
                xmin, xmax = np.min(xs)-threshold, np.max(xs)+threshold  # in case mode is in border
                ymin, ymax = np.min(ys)-threshold, np.max(ys)+threshold  # in case mode is in border
                binsx = int((xmax-xmin)/pixsize)
                binsy = int((ymax-ymin)/pixsize)
                if binsx == 0:
                    binsx = 1
                if binsy == 0:
                    binsy = 1
                    
                # 2D histogram
                H, xedges, yedges = np.histogram2d(xs, ys, bins=[binsx, binsy], range=[[xmin, xmax], [ymin, ymax]])

                # Index of maximum bin
                imax = np.unravel_index(np.argmax(H), H.shape)
                ix, iy = imax

                # Bin centers
                x_mode = 0.5 * (xedges[ix] + xedges[ix+1])
                y_mode = 0.5 * (yedges[iy] + yedges[iy+1])

                # mask near modal position
                nearmode = ((dict_kid['x0']-x_mode)**2 + (dict_kid['y0']-y_mode)**2 <= threshold**2)

                # 2D histogram
                H_2, xedges_2, yedges_2 = np.histogram2d(xs[~nearmode], ys[~nearmode], bins=[binsx, binsy], range=[[xmin, xmax], [ymin, ymax]])

                # Index of second maximum bin
                imax_2 = np.unravel_index(np.argmax(H_2), H_2.shape)
                ix_2, iy_2 = imax_2

                # Bin centers
                x_mode_2 = 0.5 * (xedges_2[ix_2] + xedges_2[ix_2+1])
                y_mode_2 = 0.5 * (yedges_2[iy_2] + yedges_2[iy_2+1])

                # mask near modal position
                nearmode_2 = ((dict_kid['x0']-x_mode_2)**2 + (dict_kid['y0']-y_mode_2)**2 <= threshold**2)

                # if the ring around the mode does not contain more points that the ring around the next modal
                # bin, then it's not the mode
                if np.sum(nearmode_2)>sum(nearmode):
                    # Use second mode as central bin
                    nearmode = nearmode_2
                    x_mode, y_mode = x_mode_2, y_mode_2

                # tone power points near mode
                tp_arr_nearmode = tp_arr[nearmode]

                # mask near mean tone_power around modal position for unfiltered array
                nearmeantp_global = (np.abs(tp_arr-np.mean(tp_arr_nearmode)) <= sig * np.std(tp_arr_nearmode))

                # final mask, both near the mode and near the modal position's mean TP
                accepted = nearmode & nearmeantp_global

                # Plotting
                
                # arrays to plot fits or circles
                tsky_toplot = np.linspace(min(tsky_arr), max(tsky_arr), 100)
                phitoplot = np.linspace(0, 2*np.pi, 100)

                # start with the histogram, which will always be there if status is good or bad_fit
                hist, _, _,img = ax[1].hist2d(xs, ys, bins=[binsx, binsy], range=[[xmin, xmax], [ymin, ymax]], cmap='autumn', norm=mcolors.LogNorm())
                cbarhist = fig.colorbar(img, ax=ax[1], orientation='horizontal')
                cbarhist.set_label('Counts')
                ax[1].axis('equal')
                ax[1].set_xlabel('Coordinate X in array')
                ax[1].set_ylabel('coordinate Y in array')
                ax[1].set_title('Position 2D Histogram for KID %i'%currentkid)

                # now the plotting of points and fits
                # first, some conditions
                if fitdict[currentkid]['separate']=='no':
                    actualmask = nearmeantp
                    linestyle = '-'
                else:
                    actualmask = accepted
                    linestyle = '--'

                if fitdict[currentkid]['status']=='good':
                    colorfit = 'green'
                    status_string = 'GOOD'
                else:
                    colorfit = 'red'
                    status_string = 'BAD FIT'

                # add outliers to mask
                if len(fitdict[currentkid]['outliers'])>0:
                    outlier_mask = ~tsky_arr==tsky_arr
                    for out in fitdict[currentkid]['outliers']:
                        outlier_mask += ((np.isclose(tsky_arr, np.real(out))) & (np.isclose(fmax_arr, np.imag(out))))
                    actualmask &= ~outlier_mask

                # plot the points ignored in fitting in white
                ax[0].scatter(tsky_arr[~actualmask], fmax_arr[~actualmask], marker='o', s=8, c='white')

                # plot the points used for fitting, green if no weights or color-coded by weight
                if weight_by==None:
                    ax[0].scatter(tsky_arr[actualmask], fmax_arr[actualmask], marker='o', s=8, c='green')
                else:
                    sc = ax[0].scatter(tsky_arr[actualmask], fmax_arr[actualmask], marker='o', s=8, c=direct_arr[actualmask], vmin=vmin, vmax=vmax, cmap=cmap)
                    cbar = fig.colorbar(sc, ax=ax[0], orientation='horizontal')
                    cbar.set_label(weight_by)

                # plot the fit, good or bad
                ax[0].plot(tsky_toplot, m*tsky_toplot+n, label=fitlabel, c=colorfit, lw=5, ls=linestyle)
                ax[0].set_ylim(np.mean(fmax_arr) - (np.max(calibdict['freq_max'])-np.min(calibdict['freq_max']))/len(kidlist) / 2 - 0.3,
                               np.mean(fmax_arr) + (np.max(calibdict['freq_max'])-np.min(calibdict['freq_max']))/len(kidlist) / 2)
                
                # Plot actual positions
                ax[1].scatter(dict_kid['x0'], dict_kid['y0'], lw=3, c='grey', alpha=0.8)

                # Plot modal bin and circle
                ax[1].scatter(x_mode, y_mode, marker='X', s=100, c='cyan', label='Mode')
                ax[1].plot(x_mode+threshold*np.cos(phitoplot), y_mode+threshold*np.sin(phitoplot), c='cyan',lw=5, label='%.3f mm-diameter ring around mode'%diameter)
                
                # Plot mean position, to be used as expected position!
                ax[1].scatter(fitdict[currentkid]['x0_mean'], fitdict[currentkid]['y0_mean'], marker='X', s=100,
                              c='green', label='Mean within ring: (%.3f, %.3f)'%(fitdict[currentkid]['x0_mean'], fitdict[currentkid]['y0_mean']))
                
                ax[1].legend();
            
                fig.suptitle('F(Tsky) for KID %i: %s (status: %s | separate: %s | npts: %i out of %i needed to fit)'%\
                             (currentkid, status_string, fitdict[currentkid]['status'], fitdict[currentkid]['separate'], fitdict[currentkid]['npts'], n_needed),
                             fontsize=15)

            ax[0].legend(loc='lower left')
            ax[0].set_xlabel('Sky temperature (K)')
            ax[0].set_xlim(mintsky, maxtsky)
            ax[0].set_ylabel('Frequency of peak (MHz)')
            ax[0].ticklabel_format(axis='y', style='plain', useOffset=False)
            plt.show()

            inp = str(raw_input(msg))
            if inp == '':
                currentkididx += 1
                currentkid = int(kidlist[currentkididx])
            elif inp =='-':
                currentkididx -=1
                currentkid = int(kidlist[currentkididx])

            # cycle through bad kids
            elif inp == 'b':
                currentbadkididx +=1 
                currentbadkid = int(badkidlist[currentbadkididx])
            
                currentkid = currentbadkid
                currentkididx = np.where(kidlist==currentkid)[0][0]

            elif inp =='q' or len(inp)>=10:
                # second condition is to avoid accidental pasting of huge text and endless loop!
                plt.close(fig)
                break
            
            else:
                try:
                    if int(inp) in kidlist:
                        currentkid = int(inp)
                        currentkididx = np.where(kidlist==currentkid)[0][0]
                        
                    else:
                        print('\n=================================================================================')
                        print('Value "{}" entered not in kid list.'.format(inp))
                        print('=================================================================================')
                except:
                    print('\n=================================================================================')
                    print('Value "{}" entered is not valid'.format(inp))
                    print('=================================================================================')

            ax[0].cla()
            ax[1].cla()
            if cbar:
                cbar.remove()
            if cbarhist:
                cbarhist.remove()

    else:
        print('')
        info('To inspect the results, re-run the function with argument inspect=True! (cycle through bad kids with input "b")')

    if full_output==False:
        return fitdict
    
    else:
        list_less30pts =       [kid for kid in kidlist if fitdict[kid]['status']=='bad_less30pts']
        list_net =             [kid for kid in kidlist if fitdict[kid]['status']=='bad_NET']
        list_df2k =            [kid for kid in kidlist if fitdict[kid]['status']=='bad_dF2K']
        list_fewpts =          [kid for kid in kidlist if fitdict[kid]['status']=='bad_fewpts']
        list_badfit =          [kid for kid in kidlist if fitdict[kid]['status']=='bad_fit']
        
        return fitdict, list_less30pts, list_net, list_df2k, list_fewpts, list_badfit