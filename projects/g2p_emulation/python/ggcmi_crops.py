#!/bin/env python
from em_functions import try_load_climate, update_out_table, save_out_table, do_emulation
import numpy as np
import os
from netCDF4 import Dataset
# import glob
# import datetime
np.random.seed(1234)


def PLUMemulate(GCM, rcp, decade, GGCM, mask_YX, outarr_yield, outarr_irrig, do_adapt):
    """
    GCM: ['ACCESS1-0','bcc-csm1-1','BNU-ESM','CanESM2','CCSM4','CESM1-BGC','CMCC-CM',
          'CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','FGOALS-g2','GFDL-CM3','GFDL-ESM2G',
          'GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES',
          'inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM',
          'MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']
    rcp: 45 or 85
    decade: 0,1,2,3,4,5,6,7,8
        (which equates to 2011-2019 [shortened!], 2020-2029, ... 2090-2099)
    GGCM: ['EPIC-TAMU','pDSSAT','LPJ-GUESS', 'LPJmL']
    GGCMIcrop: ['maize' , 'rice', 'soy', 'winter_wheat', 'spring_wheat']
    """
    # Begin strings for outfmt and outheader
    outfmt = "%0.2f %0.2f"
    outheader = "Lon Lat"

    # Define crops
    if GGCM == "LPJ-GUESS":
        GGCMIcrops = ["winter_wheat", "spring_wheat", "maize"]
    else:
        GGCMIcrops = ["winter_wheat", "spring_wheat", "maize", "soy", "rice"]

    # Set up for tracking of which crops were missing
    missing_climate_dict = {
        "max_wheat": False,
        "maize": False,
        "rice": False,
        "soy": False,
        "spring_wheat": False,
    }

    # prev_GGCMIcrop = "nothing"
    # prev_PLUMcrop = "nothing"
    any_ok = False
    for GGCMIcrop in GGCMIcrops:

        print("%s rcp%d dec%d %s %s" % (GCM, rcp, decade, GGCM, GGCMIcrop))

        # Load climate files
        climate_dir = \
            "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/climate/emulator_gs_deltas/CMIP5"
        missing_climate = False
        t, missing_climate = try_load_climate(
            climate_dir, "tas", GCM, GGCMIcrop, rcp, missing_climate, False)
        w, missing_climate = try_load_climate(
            climate_dir, "pr", GCM, GGCMIcrop, rcp, missing_climate, False)

        # If any climate files were missing, try again in second location
        if missing_climate:
            climate_dir2 = \
                "/Volumes/WDMPP_Storage/Shared/GGCMI2PLUM_sh/emulation/" \
                + "inputs/climate/agmerra/cmip5"
            missing_climate = False
            t, missing_climate = try_load_climate(
                climate_dir2, "tas", GCM, GGCMIcrop, rcp, missing_climate, True)
            w, missing_climate = try_load_climate(
                climate_dir2, "pr", GCM, GGCMIcrop, rcp, missing_climate, True)

        # If any climate files are STILL missing, skip to the next crop.
        if missing_climate:
            missing_climate_dict[GGCMIcrop] = True
            print("    Skipping.")
            continue
        any_ok = True

        # CO2 vars
        if rcp == 45:
            co2 = [400.1285, 423.0875, 447.9455, 473.69, 497.703, 516.5865, 527.72, 532.4395,
                   536.0495]
        elif rcp == 85:
            co2 = [402.552, 432.3075, 469.135, 514.989, 572.0315, 640.299, 717.63, 801.4935,
                   890.3395]

        # Emulate the six management cases.
        emulator_dir_yield = "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/fits_yield"
        emulator_dir_irrig = "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/fits_irrig"
        rf_10_yield, rf_60_yield, rf_200_yield, ir_10_yield, ir_60_yield, ir_200_yield \
            = do_emulation(emulator_dir_yield, GGCM, GGCMIcrop, co2, t, w, False, decade,
                           do_adapt)
        rf_10_irrig, rf_60_irrig,rf_200_irrig, ir_10_irrig,ir_60_irrig, ir_200_irrig \
            = do_emulation(emulator_dir_irrig, GGCM, GGCMIcrop, co2, t, w, True, decade,
                           do_adapt)
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

    return(outarr_yield, outarr_irrig, outheader, outfmt,any_ok)


do_adapt = False
outdir_suffix = "20200310"

GCMs = ["IPSL-CM5A-MR_r1i1p1"]
# GCMs = ["IPSL-CM5A-MR","GFDL-ESM2M","MIROC5", "HadGEM2-ES"]
# GCMs = ["IPSL-CM5A-LR", "GFDL-ESM2M", "MIROC5", "HadGEM2-ES"] # ISIMIP2b selection
# GCMs = ['ACCESS1-0','bcc-csm1-1','BNU-ESM','CanESM2','CCSM4','CESM1-BGC','CMCC-CM',
#         'CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','FGOALS-g2','GFDL-CM3','GFDL-ESM2G',
#         'GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC',
#         'HadGEM2-ES','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5',
#         'MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']
# rcps = [45, 85]
rcps = [45]
decades = range(9)
# GGCMs = ["EPIC-TAMU", "LPJ-GUESS", "LPJmL", "pDSSAT"]
GGCMs = ["pDSSAT"]
# GGCMs = ["EPIC-TAMU", "LPJ-GUESS", "LPJmL"]

# Import PLUM mask and lon/lat
plum_dir = "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/plum/"
mask_YX = np.genfromtxt(plum_dir+"PLUM_mask.csv", delimiter=",")
lons_YX = np.genfromtxt(plum_dir+"PLUM_map_lons.csv", delimiter=",")
lats_YX = np.genfromtxt(plum_dir+"PLUM_map_lats.csv", delimiter=",")
lons = lons_YX[mask_YX == 1]
lats = lats_YX[mask_YX == 1]
lonlats = np.vstack((lons,lats))

for decade in decades:
    for rcp in rcps:
        for GCM in GCMs:
            for GGCM in GGCMs:

                # Set up info about this run
                y1 = 2011+10*decade
                yN = min(2099, 2020+10*decade)
                decade_str = '%d-%d' % (y1, yN)
                print("%s rcp%d %s %s..." % (GCM, rcp, decade_str, GGCM))
                # outdir = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/'\
                #     + 'outputs_GGCMIcrops_%s/%s/%s/rcp%d/%s' \
                #     % (datetime.datetime.now().strftime("%Y%m%d"), GCM, GGCM, rcp,
                #        decade_str)
                outdir = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/' \
                    + 'outputs_GGCMIcrops_%s/%s/%s/rcp%d/%s' \
                    % (outdir_suffix, GCM, GGCM, rcp, decade_str)
                outfile_yield = '%s/yield.out' % (outdir)
                outfile_irrig = '%s/gsirrigation.out' % (outdir)

                # Skip if files already exist
                exist_yield = os.path.isfile(outfile_yield) \
                    or os.path.isfile(outfile_yield + ".gz")
                exist_irrig = os.path.isfile(outfile_irrig) \
                    or os.path.isfile(outfile_irrig + ".gz")
                if exist_yield and exist_irrig:
                    print("Emulator outputs already exist; skipping.\n")
                    continue

                # Emulate
                outarr_yield, outarr_irrig, outheader, outfmt, any_ok = PLUMemulate(
                    GCM, rcp, decade, GGCM, mask_YX, lonlats, lonlats, do_adapt)

                # If nothing actually processed, skip to next GGCM
                if not any_ok:
                    print("\n")
                    continue

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

















