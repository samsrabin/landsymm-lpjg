#!/bin/env python
import numpy as np
import os
# import glob
# import datetime
from em_functions import emulate, update_out_table

emulator_dir = "../fits_yield"
co2 = 360
t = 0
w = 1
is_irrig = False
GGCM = "pDSSAT"

os.chdir("/Users/Shared/GGCMI2PLUM/emulator/Sam/")

# crop_list_long = ['maize', 'winter_wheat', 'spring_wheat']
# crop_list_short = ["mai", "wwh", "swh"]
crop_list_long = ['maize', 'winter_wheat', 'spring_wheat', 'soy', 'rice']
crop_list_short = ["mai", "wwh", "swh", "soy", "ric"]
N_list = [10, 60, 200]

outdir = "outputs_recreate_phase2"
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

for c in np.arange(0,len(crop_list_long)):
    crop_long = crop_list_long[c]
    if GGCM == "LPJ-GUESS" and (crop_long == "soy" or crop_long == "rice"):
        continue
    crop_short = crop_list_short[c]
    print(crop_long)
    # Emulate
    K = np.load("%s/%s_%s.npy" % (emulator_dir, GGCM, crop_long))
    KI = np.load("%s/%s_%s_I.npy" % (emulator_dir, GGCM, crop_long))
    rf_10 = emulate(
        K, co2, t, w, 10, "N", is_irrig
    )
    rf_60 = emulate(
        K, co2, t, w, 60, "N", is_irrig
    )
    rf_200 = emulate(
        K, co2, t, w, 200, "N", is_irrig
    )
    ir_10 = emulate(
        KI, co2, t, w, 10, "NI", is_irrig
    )
    ir_60 = emulate(
        KI, co2, t, w, 60, "NI", is_irrig
    )
    ir_200 = emulate(
        KI, co2, t, w, 200, "NI", is_irrig
    )

    # Update output table
    outarr_yield, outheader, outfmt = update_out_table(
        outarr_yield,
        outheader,
        outfmt,
        crop_short,
        rf_10,
        rf_60,
        rf_200,
        ir_10,
        ir_60,
        ir_200,
        mask_YX
    )

# Save in LPJ-GUESS/PLUM-readable format
np.savetxt(
    outfile,
    outarr_yield.T,
    delimiter=" ",
    fmt=outfmt,
    header=outheader,
    comments="",
)
