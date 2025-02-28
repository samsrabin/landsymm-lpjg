#!/bin/env python
import numpy as np
import os
# import glob
import datetime
from em_functions import try_load_climate, do_emulation, update_out_table

np.random.seed(1234)


def do_emulation_wheats(emulator_dir, co2, ts, ws, tw, ww, is_irrig, GGCM, decade):
    # rf_10w, rf_60w, rf_200w, ir_10w, ir_60w, ir_200w = do_emulation(
    #     emulator_dir, "winter_wheat", co2, tw, ww, is_irrig, GGCM, decade
    # )
    # rf_10s, rf_60s, rf_200s, ir_10s, ir_60s, ir_200s = do_emulation(
    #     emulator_dir, "spring_wheat", co2, ts, ws, is_irrig, GGCM, decade
    # )
    rf_10w, rf_60w, rf_200w, ir_10w, ir_60w, ir_200w \
        = do_emulation(emulator_dir, GGCM, 'winter_wheat', co2, tw, ww, is_irrig, decade,
                       False)
    rf_10s, rf_60s, rf_200s, ir_10s, ir_60s, ir_200s \
        = do_emulation(emulator_dir, GGCM, 'spring_wheat', co2, ts, ws, is_irrig, decade,
                       False)

    # Get maximum winter or spring
    if is_irrig:
        # No irrigation if rainfed
        rf_10 = np.zeros(ir_10w.shape)
        rf_60 = np.zeros(ir_10w.shape)
        rf_200 = np.zeros(ir_10w.shape)
    else:
        rf_10 = np.maximum(rf_10w, rf_10s)
        rf_60 = np.maximum(rf_60w, rf_60s)
        rf_200 = np.maximum(rf_200w, rf_200s)
    ir_10 = np.maximum(ir_10w, ir_10s)
    ir_60 = np.maximum(ir_60w, ir_60s)
    ir_200 = np.maximum(ir_200w, ir_200s)

    return (rf_10, rf_60, rf_200, ir_10, ir_60, ir_200)


def PLUMemulate(GCM, rcp, decade, GGCM, mask_YX, outarr_yield, outarr_irrig):
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

    # Define PLUM output crops
    PLUMcrops = ["CerealsC3", "CerealsC4", "Rice", "Oilcrops", "Pulses", "StarchyRoots"]

    # Define matching to GGCMI crops
    PLUM2GGCMI_dict = {
        "CerealsC3": "max_wheat",
        "CerealsC4": "maize",
        "Rice": "rice",
        "Oilcrops": "soy",
        "Pulses": "soy",
        "StarchyRoots": "spring_wheat",
    }
    if GGCM == "LPJ-GUESS":
        PLUM2GGCMI_dict["Rice"] = "spring_wheat"
        PLUM2GGCMI_dict["Oilcrops"] = "spring_wheat"
        PLUM2GGCMI_dict["Pulses"] = "spring_wheat"

    # Set up for tracking of which crops were missing
    missing_climate_dict = {
        "max_wheat": False,
        "maize": False,
        "rice": False,
        "soy": False,
        "spring_wheat": False,
    }

    prev_GGCMIcrop = "nothing"
    prev_PLUMcrop = "nothing"

    # Set up dummies so parser doesn't get upset
    rf_10_yield = np.nan
    rf_60_yield = np.nan
    rf_200_yield = np.nan
    ir_10_yield = np.nan
    ir_60_yield = np.nan
    ir_200_yield = np.nan
    rf_10_irrig = np.nan
    rf_60_irrig = np.nan
    rf_200_irrig = np.nan
    ir_10_irrig = np.nan
    ir_60_irrig = np.nan
    ir_200_irrig = np.nan

    # Loop through crops
    for PLUMcrop in PLUMcrops:

        GGCMIcrop = PLUM2GGCMI_dict[PLUMcrop]
        print("%s rcp%d dec%d %s %s %s" % (GCM, rcp, decade, GGCM, PLUMcrop, GGCMIcrop))

        # Reuse data from previous PLUMcrop if it used the same GGCMIcrop.
        # Could make this more efficient by making it so ANY previous PLUMcrop
        # would work, but this will do for now. Just arrange dict to maximize
        # consecutive GGCMIcrops.
        if (GGCMIcrop == prev_GGCMIcrop) and not (prev_PLUMcrop == "CerealsC3"):
            print(
                "    %s and %s both use emulator for %s."
                % (prev_PLUMcrop, PLUMcrop, prev_GGCMIcrop)
            )
            if missing_climate_dict[prev_GGCMIcrop]:
                print(
                    "    But %s climate is missing, so skipping %s."
                    % (GGCMIcrop, PLUMcrop)
                )
            else:
                print("    Using results from %s for %s." % (prev_PLUMcrop, PLUMcrop))
                outarr_yield, outheader, outfmt = update_out_table(
                    outarr_yield,
                    outheader,
                    outfmt,
                    PLUMcrop,
                    rf_10_yield,
                    rf_60_yield,
                    rf_200_yield,
                    ir_10_yield,
                    ir_60_yield,
                    ir_200_yield,
                    mask_YX,
                )
                outarr_irrig, x, y = update_out_table(
                    outarr_irrig,
                    outheader,
                    outfmt,
                    PLUMcrop,
                    rf_10_irrig,
                    rf_60_irrig,
                    rf_200_irrig,
                    ir_10_irrig,
                    ir_60_irrig,
                    ir_200_irrig,
                    mask_YX,
                )
            prev_GGCMIcrop = GGCMIcrop
            prev_PLUMcrop = PLUMcrop
            continue

        # Load climate files
        climate_dir = \
            "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/climate/emulator_gs_deltas/CMIP5"
        missing_climate = False
        verbose = True
        if PLUM2GGCMI_dict[PLUMcrop] == "max_wheat":
            GGCMIcrop_tmp = "winter_wheat"
            tw, missing_climate = try_load_climate(
                climate_dir, "tas", GCM, GGCMIcrop_tmp, rcp, missing_climate, verbose)
            ww, missing_climate = try_load_climate(
                climate_dir, "pr", GCM, GGCMIcrop_tmp, rcp, missing_climate, verbose)
            GGCMIcrop_tmp = "spring_wheat"
            ts, missing_climate = try_load_climate(
                climate_dir, "tas", GCM, GGCMIcrop_tmp, rcp, missing_climate, verbose)
            ws, missing_climate = try_load_climate(
                climate_dir, "pr", GCM, GGCMIcrop_tmp, rcp, missing_climate, verbose)
            del GGCMIcrop_tmp
        else:
            t, missing_climate = try_load_climate(
                climate_dir, "tas", GCM, GGCMIcrop, rcp, missing_climate, verbose)
            w, missing_climate = try_load_climate(
                climate_dir, "pr", GCM, GGCMIcrop, rcp, missing_climate, verbose)

        # If any climate files were missing, skip to the next crop.
        if missing_climate:
            missing_climate_dict[GGCMIcrop] = True
            prev_GGCMIcrop = GGCMIcrop
            prev_PLUMcrop = PLUMcrop
            print("    Skipping.")
            breakpoint()
            continue

        # CO2 vars
        if rcp == 45:
            co2 = [
                400.1285,
                423.0875,
                447.9455,
                473.69,
                497.703,
                516.5865,
                527.72,
                532.4395,
                536.0495,
            ]
        elif rcp == 85:
            co2 = [
                402.552,
                432.3075,
                469.135,
                514.989,
                572.0315,
                640.299,
                717.63,
                801.4935,
                890.3395,
            ]

        # Emulate the six management cases.
        # If CerealsC3, emulate both winter and spring wheat, then find the maximum
        # at each treatment.
        emulator_dir_yield = "inputs/fits_yield"
        emulator_dir_irrig = "inputs/fits_irrig"
        if PLUM2GGCMI_dict[PLUMcrop] == "max_wheat":
            rf_10_yield, rf_60_yield, rf_200_yield, \
                ir_10_yield, ir_60_yield, ir_200_yield = do_emulation_wheats(
                    emulator_dir_yield, co2, ts, ws, tw, ww, False, GGCM, decade
                )
            rf_10_irrig, rf_60_irrig, rf_200_irrig, \
                ir_10_irrig, ir_60_irrig, ir_200_irrig = do_emulation_wheats(
                    emulator_dir_irrig, co2, ts, ws, tw, ww, True, GGCM, decade
                )
            del ts, ws, tw, ww
        else:
            rf_10_yield, rf_60_yield, rf_200_yield, \
                ir_10_yield, ir_60_yield, ir_200_yield = do_emulation(
                    emulator_dir_yield, GGCM, GGCMIcrop, co2, t, w, False, decade, False
                )
            rf_10_irrig, rf_60_irrig, rf_200_irrig, \
                ir_10_irrig, ir_60_irrig, ir_200_irrig = do_emulation(
                    emulator_dir_irrig, GGCM, GGCMIcrop, co2, t, w, True, decade, False
                )
            del t, w

        # Add values to output tables
        outarr_yield, outheader, outfmt = update_out_table(
            outarr_yield,
            outheader,
            outfmt,
            PLUMcrop,
            rf_10_yield,
            rf_60_yield,
            rf_200_yield,
            ir_10_yield,
            ir_60_yield,
            ir_200_yield,
            mask_YX,
        )
        outarr_irrig, x, y = update_out_table(
            outarr_irrig,
            outheader,
            outfmt,
            PLUMcrop,
            rf_10_irrig,
            rf_60_irrig,
            rf_200_irrig,
            ir_10_irrig,
            ir_60_irrig,
            ir_200_irrig,
            mask_YX,
        )

        prev_GGCMIcrop = GGCMIcrop
        prev_PLUMcrop = PLUMcrop

    return (outarr_yield, outarr_irrig, outheader, outfmt)


GCMs = ["IPSL-CM5A-MR_r1i1p1"]
# GCMs = ["IPSL-CM5A-LR", "GFDL-ESM2M", "MIROC5", "HadGEM2-ES"] # ISIMIP2b selection
# GCMs = ['ACCESS1-0','bcc-csm1-1','BNU-ESM','CanESM2','CCSM4','CESM1-BGC','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','FGOALS-g2','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']
rcps = [45, 85]
decades = range(9)
GGCMs = ["EPIC-TAMU", "LPJ-GUESS", "LPJmL", "pDSSAT"]
# GGCMs = ["LPJ-GUESS", "LPJmL", "pDSSAT"]

# Which system are we on?
if os.getcwd()[0:6]=='/Users':
    os.chdir("/Users/Shared/GGCMI2PLUM_sh/emulation/")
elif os.getcwd()[0:7]=='/project':
    os.chdir("/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/")
else:
    raise NameError('System not recognized')

# Import PLUM mask and lon/lat
mask_YX = np.genfromtxt("inputs/plum/PLUM_mask.csv", delimiter=",")
lons_YX = np.genfromtxt("inputs/plum/PLUM_map_lons.csv", delimiter=",")
lats_YX = np.genfromtxt("inputs/plum/PLUM_map_lats.csv", delimiter=",")
lons = lons_YX[mask_YX == 1]
lats = lats_YX[mask_YX == 1]
lonlats = np.vstack((lons, lats))

for decade in decades:
    for rcp in rcps:
        for GCM in GCMs:
            for GGCM in GGCMs:

                # Set up info about this run
                decade_str = "%d-%d" % (2011 + 10 * decade, 2020 + 10 * decade)
                print("%s rcp%d %s %s..." % (GCM, rcp, decade_str, GGCM))
                outdir = (
                    "outputs/potYields_%s/%s/%s/rcp%d/%s"
                    % (
                        datetime.datetime.now().strftime("%Y%m%d"),
                        GCM,
                        GGCM,
                        rcp,
                        decade_str,
                    )
                )
                outfile_yield = "%s/yield.out" % (outdir)
                outfile_irrig = "%s/gsirrigation.out" % (outdir)

                # Make output directory, if needed
                try:
                    os.makedirs(outdir)
                    print("mkdir -p " + outdir)
                except FileExistsError:
                    # directory already exists
                    pass

                # Emulate
                outarr_yield, outarr_irrig, outheader, outfmt = PLUMemulate(
                    GCM, rcp, decade, GGCM, mask_YX, lonlats, lonlats
                )

                # Convert from tons/ha to kg/m2
                outarr_yield[:][2:] = outarr_yield[:][2:] * 0.1

                # Save in PLUM-readable format
                np.savetxt(
                    outfile_yield,
                    outarr_yield.T,
                    delimiter=" ",
                    fmt=outfmt,
                    header=outheader,
                    comments="",
                )
                np.savetxt(
                    outfile_irrig,
                    outarr_irrig.T,
                    delimiter=" ",
                    fmt=outfmt,
                    header=outheader,
                    comments="",
                )

                # Compress and add "done" file for PLUM
                os.system("gzip %s" % (outfile_yield))
                os.system("gzip %s" % (outfile_irrig))
                os.system("touch %s/done" % (outdir))
