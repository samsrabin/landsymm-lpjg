#!/bin/env python
import numpy as np
import os
# import glob
#import datetime
os.chdir("/Users/Shared/GGCMI2PLUM/emulator/Sam")
from em_functions import emulate

import scipy.io as sio

emulator_dir = "../fits_yield"
co2 = 360
t = 0
w = 1
GGCM = "pDSSAT"


#%% Setup

N = sio.loadmat('nfert.remapv5e.mat')
Ncrops_PLUM = N['cropList'][0].size
cropList_PLUM = ['']*Ncrops_PLUM
PLUM_to_GGCMI_long = {}
PLUM_to_GGCMI_short = {}
print('Setting up PLUM_to_GGCMI dictionary...')
for c in np.arange(0, Ncrops_PLUM):
    thisCrop = N['cropList'][0][c][0]
    cropList_PLUM[c] = thisCrop
    if 'Oilcrops' in thisCrop \
    or 'Pulses' in thisCrop:
        if GGCM == "LPJ-GUESS":
            PLUM_to_GGCMI_long[thisCrop] = 'spring_wheat'
            PLUM_to_GGCMI_short[thisCrop] = 'swh'
        else:
            PLUM_to_GGCMI_long[thisCrop] = 'soy'
            PLUM_to_GGCMI_short[thisCrop] = 'soy'
    elif 'CerealsC3' in thisCrop:
        PLUM_to_GGCMI_long[thisCrop] = 'max_wheat'
        PLUM_to_GGCMI_short[thisCrop] = 'max'
    elif 'CerealsC4' in thisCrop:
        PLUM_to_GGCMI_long[thisCrop] = 'maize'
        PLUM_to_GGCMI_short[thisCrop] = 'mai'
    elif 'Rice' in thisCrop:
        if GGCM == "LPJ-GUESS":
            PLUM_to_GGCMI_long[thisCrop] = 'spring_wheat'
            PLUM_to_GGCMI_short[thisCrop] = 'swh'
        else:
            PLUM_to_GGCMI_long[thisCrop] = 'rice'
            PLUM_to_GGCMI_short[thisCrop] = 'ric'
    elif 'StarchyRoots' in thisCrop:
        PLUM_to_GGCMI_long[thisCrop] = 'spring_wheat'
        PLUM_to_GGCMI_short[thisCrop] = 'swh'
    elif 'ExtraCrop' in thisCrop \
    or 'Miscanthus' in thisCrop:
        PLUM_to_GGCMI_long[thisCrop] = 'NONE'
        PLUM_to_GGCMI_short[thisCrop] = 'NONE'
    else:
        raise Exception('thisCrop not recognized')
    print("    %s = %s" % (thisCrop, PLUM_to_GGCMI_long[thisCrop]))
print('Done.')

outdir = "outputs_calib"
outfile = "%s/%s.out" % (outdir, GGCM)

# Make output directory, if needed
try:
    os.makedirs(outdir)
    print("mkdir -p " + outdir)
except FileExistsError:
    # directory already exists
    pass

# Import PLUM mask and lon/lat
mask_YX = np.genfromtxt("plum/PLUM_mask.csv", delimiter=",")
lons_YX = np.genfromtxt("plum/PLUM_map_lons.csv", delimiter=",")
lats_YX = np.genfromtxt("plum/PLUM_map_lats.csv", delimiter=",")
lons = lons_YX[mask_YX == 1]
lats = lats_YX[mask_YX == 1]
lonlats = np.vstack((lons, lats))

# Begin strings for outfmt and outheader
outfmt = "%0.2f %0.2f"
outheader = "Lon Lat"
outarr_yield = lonlats


#%% Process

for c in np.arange(0,len(cropList_PLUM)):
    thisCrop_PLUM = cropList_PLUM[c]
    thisCrop_GGCMI_long = PLUM_to_GGCMI_long[thisCrop_PLUM]
    thisCrop_GGCMI_short = PLUM_to_GGCMI_long[thisCrop_PLUM]
    print(thisCrop_PLUM)
    
    if thisCrop_GGCMI_short == "NONE":
        continue
    
    # Get fertilizer application (convert kgN/m2 to kgN/ha)
    thisN = N['nfert_cYX'][c,:,:]*1e4
    
    # Emulate
    is_irr = thisCrop_PLUM[-1]=="i"
    if thisCrop_GGCMI_long == "max_wheat":
        if is_irr:
            Kw = np.load("%s/%s_%s_I.npy" % (emulator_dir, GGCM, 'winter_wheat'))
            tmpW = emulate(
                    Kw, co2, t, w, thisN
                    )
            Ks = np.load("%s/%s_%s_I.npy" % (emulator_dir, GGCM, 'spring_wheat'))
            tmpS = emulate(
                    Ks, co2, t, w, thisN
                    )
        else:
            Kw = np.load("%s/%s_%s.npy" % (emulator_dir, GGCM, 'winter_wheat'))    
            tmpW = emulate(
                Kw, co2, t, w, thisN
            )
            Ks = np.load("%s/%s_%s.npy" % (emulator_dir, GGCM, 'spring_wheat'))    
            tmpS = emulate(
                Ks, co2, t, w, thisN
            )
        tmp = np.maximum(tmpW, tmpS)
    else:
        if is_irr:
            K = np.load("%s/%s_%s_I.npy" % (emulator_dir, GGCM, thisCrop_GGCMI_long))
            tmp = emulate(
                    K, co2, t, w, thisN
                    )
        else:
            K = np.load("%s/%s_%s.npy" % (emulator_dir, GGCM, thisCrop_GGCMI_long))    
            tmp = emulate(
                K, co2, t, w, thisN
            )

    # outarr_yield comes in as (crop,cell) because vstack is the only stack that works with one-d
    # arrays like rf_10 etc., apparently. outarr_yield will be transposed for write so that each
    # row is a gridcell.
    outarr_yield = np.vstack((outarr_yield, tmp[mask_YX == 1]))

    # Update header
    outheader = (
        outheader
        + " "
        + thisCrop_PLUM
    )
    outfmt = outfmt + " %0.6f"

# Convert from tons/ha to kg/m2
outarr_yield[:][2:] = outarr_yield[:][2:] * 0.1

#%% Save in LPJ-GUESS/PLUM-readable format
np.savetxt(
    outfile,
    outarr_yield.T,
    delimiter=" ",
    fmt=outfmt,
    header=outheader,
    comments="",
)
