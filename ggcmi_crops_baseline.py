#!/bin/env python
from netCDF4 import Dataset
from em_functions \
    import emulate_old, emulate, update_out_table, save_out_table, import_parameter_netcdfs
import numpy as np
import os
# import glob
# import datetime
np.random.seed(1234)


def do_emulation(emulator_dir, GGCMIcrop, co2, t, w, is_irrig):

    # Get files and read parameters
    if is_irrig:
        KI = np.load('%s/%s_%s_irr.npy' % (emulator_dir, GGCM, GGCMIcrop))
    else:
        K, KI = import_parameter_netcdfs(emulator_dir, GGCM, GGCMIcrop, False)

    # Irrigated
    if is_irrig:
        ir_10 = emulate_old(KI, co2, t, 1, 10, "NI", True)
        ir_60 = emulate_old(KI, co2, t, 1, 60, "NI", True)
        ir_200 = emulate_old(KI, co2, t, 1, 200, "NI", True)
    else:
        ir_10 = emulate(KI, co2, t, 1, 10)
        ir_60 = emulate(KI, co2, t, 1, 60)
        ir_200 = emulate(KI, co2, t, 1, 200)

    # Rainfed
    if is_irrig:
        # No need to emulate irrigation for rainfed crops
        rf_10 = np.zeros(ir_10.shape)
        rf_60 = np.zeros(ir_10.shape)
        rf_200 = np.zeros(ir_10.shape)
    else:
        rf_10 = emulate(K, co2, t, w, 10)
        rf_60 = emulate(K, co2, t, w, 60)
        rf_200 = emulate(K, co2, t, w, 200)
    return(rf_10, rf_60, rf_200, ir_10, ir_60, ir_200)


def PLUMemulate(GGCM, mask_YX, outarr_yield, outarr_irrig):
    # GGCM: ['EPIC-TAMU','pDSSAT','LPJ-GUESS', 'LPJmL']
    # GGCMIcrop: ['maize' , 'rice', 'soy', 'winter_wheat', 'spring_wheat']

    # Begin strings for outfmt and outheader
    outfmt = "%0.2f %0.2f"
    outheader = "Lon Lat"

    # Define crops
    if GGCM == "LPJ-GUESS":
        GGCMIcrops = ["winter_wheat", "spring_wheat", "maize"]
    else:
        GGCMIcrops = ["winter_wheat", "spring_wheat", "maize", "soy", "rice"]

    # # Set up for tracking of which crops were missing
    # missing_climate_dict = {
    # "max_wheat": False,
    # "maize": False,
    # "rice":  False,
    # "soy":   False,
    # "spring_wheat": False,
    # }

    # prev_GGCMIcrop = "nothing"
    # prev_PLUMcrop = "nothing"
    for GGCMIcrop in GGCMIcrops:

        # Set climate
        t = 0.0
        w = 1.0

        # CO2 vars
        co2 = 360

        # Emulate the six management cases.
        if GGCM == "LPJmL" and GGCMIcrop == "soy":
            print("SKIPPING LPJML SOY (has Jim redone this yet?")
            rf_10_yield = np.empty(mask_YX.shape)
            rf_10_yield[:] = np.NaN
            rf_60_yield = rf_10_yield
            rf_200_yield = rf_10_yield
            ir_10_yield = rf_10_yield
            ir_60_yield = rf_10_yield
            ir_200_yield = rf_10_yield
        else:
            print("%s %s" % (GGCM, GGCMIcrop))
            emulator_dir_yield = "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/fits_yield"
            emulator_dir_irrig = "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/fits_irrig"
            rf_10_yield, rf_60_yield, rf_200_yield, ir_10_yield, ir_60_yield, ir_200_yield \
                = do_emulation(emulator_dir_yield, GGCMIcrop, co2, t, w, False)
            rf_10_irrig, rf_60_irrig, rf_200_irrig, ir_10_irrig, ir_60_irrig, ir_200_irrig\
                = do_emulation(emulator_dir_irrig, GGCMIcrop, co2, t, w, True)
        del t, w

        # Add values to output tables
        outarr_yield, outheader, outfmt = update_out_table(
            outarr_yield, outheader, outfmt, GGCMIcrop,
            rf_10_yield, rf_60_yield, rf_200_yield,
            ir_10_yield, ir_60_yield, ir_200_yield,
            mask_YX)
        outarr_irrig, x, y = update_out_table(
            outarr_irrig, outheader, outfmt, GGCMIcrop,
            rf_10_irrig, rf_60_irrig, rf_200_irrig,
            ir_10_irrig, ir_60_irrig, ir_200_irrig,
            mask_YX)

    return(outarr_yield, outarr_irrig, outheader, outfmt)


outdir_suffix = "20200211"

GGCMs = ["EPIC-TAMU", "LPJ-GUESS", "LPJmL", "pDSSAT"]
# GGCMs = ["LPJmL"]

# Import PLUM mask and lon/lat
plum_dir = "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/plum/"
mask_YX = np.genfromtxt(plum_dir + "PLUM_mask.csv", delimiter=",")
lons_YX = np.genfromtxt(plum_dir + "PLUM_map_lons.csv", delimiter=",")
lats_YX = np.genfromtxt(plum_dir + "PLUM_map_lats.csv", delimiter=",")
lons = lons_YX[mask_YX == 1]
lats = lats_YX[mask_YX == 1]
lonlats = np.vstack((lons,lats))

for GGCM in GGCMs:

    # Set up info about this run
    print(GGCM)
    outdir = \
        '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_GGCMIcropsBaseline_%s/%s' \
        % (outdir_suffix, GGCM)
    outfile_yield = '%s/yield.out' % (outdir)
    outfile_irrig = '%s/gsirrigation.out' % (outdir)

    # Skip if files already exist
    exist_yield = os.path.isfile(outfile_yield) or os.path.isfile(outfile_yield + ".gz")
    exist_irrig = os.path.isfile(outfile_irrig) or os.path.isfile(outfile_irrig + ".gz")
    if exist_yield and exist_irrig:
        print("Emulator outputs already exist; skipping.\n")
        continue

    # Emulate
    outarr_yield,outarr_irrig,outheader,outfmt = PLUMemulate(GGCM, mask_YX, lonlats, lonlats)

    # Convert from tons/ha to kg/m2
    outarr_yield[:][2:] = outarr_yield[:][2:] * 0.1

    # Save in PLUM-readable format; compress
    if not exist_yield:
        save_out_table(outdir, outfile_yield, outfmt, outheader, outarr_yield)
    if not exist_irrig:
        save_out_table(outdir, outfile_irrig, outfmt, outheader, outarr_irrig)

    # Compress and add "done" file for PLUM
    os.system('touch %s/done' % (outdir))

    print("\n")

















