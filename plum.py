#!/bin/env python
import numpy as np;
import os
import glob
np.random.seed(1234)
def emulate(K, c, t, w, n, case):
    #Case: NN- no nitrogen, NNI - No nitrogen irr, N- with nitrogen, NI- with nitrogen irr
    if (case == 'NN'):
        Y = (K[0,:,:]              +
        K[1,:,:]  * c          +
        K[2,:,:]  * t          +
        K[3,:,:]  * w          +
        K[4,:,:]  * c ** 2     +
        K[5,:,:]  * t ** 2     +
        K[6,:,:]  * w ** 2     +
        K[7,:,:]  * c *  w     +
        K[8,:,:]  * t *  w     +
        K[9,:,:]  * t ** 3     +
        K[10,:,:] * w ** 3     +
        K[11,:,:] * t ** 2 * w +
        K[12,:,:] * w ** 2 * t  )

    elif (case == 'NNI'):
        Y = (K[0,:,:]              +
        K[1,:,:]  * c          +
        K[2,:,:]  * t          +
        K[3,:,:]  * c ** 2     +
        K[4,:,:]  * t ** 2     +
        K[5,:,:]  * t ** 3      )

    elif (case == 'N'):
        Y = (K[0,:,:]              +
        K[1,:,:]  * c          +
        K[2,:,:]  * t          +
        K[3,:,:]  * w          +
        K[4,:,:]  * n          +
        K[5,:,:]  * c ** 2     +
        K[6,:,:]  * t ** 2     +
        K[7,:,:]  * w ** 2     +
        K[8,:,:]  * n ** 2     +
        K[9,:,:]  * c *  w     +
        K[10,:,:] * c *  n     +
        K[11,:,:] * t *  w     +
        K[12,:,:] * t *  n     +
        K[13,:,:] * w *  n     +
        K[14,:,:] * t ** 3     +
        K[15,:,:] * w ** 3     +
        K[16,:,:] * t *  w * n +
        K[17,:,:] * t ** 2 * w +
        K[18,:,:] * w ** 2 * t +
        K[19,:,:] * w ** 2 * n +
        K[20,:,:] * n ** 2 * c +
        K[21,:,:] * n ** 2 * t +
        K[22,:,:] * n ** 2 * w  )

    elif (case == 'NI'):
        Y = (K[0,:,:]          +
        K[1,:,:]  * c          +
        K[2,:,:]  * t          +
        K[3,:,:]  * n          +
        K[4,:,:]  * c ** 2     +
        K[5,:,:]  * t ** 2     +
        K[6,:,:]  * n ** 2     +
        K[7,:,:]  * c *  n     +
        K[8,:,:]  * t *  n     +
        K[9,:,:]  * t ** 3     +
        K[10,:,:] * n ** 2 * c +
        K[11,:,:] * n ** 2 * t  )

    elif (case == 'NL'):
        Y = (K[0,:,:]          +
        K[1,:,:]  * c          +
        K[2,:,:]  * t          +
        K[3,:,:]  * w          +
        K[4,:,:]  * n          +
        K[5,:,:]  * c ** 2     +
        K[6,:,:]  * t ** 2     +
        K[7,:,:]  * w ** 2     +
        K[8,:,:]  * n ** 2     +
        K[9,:,:]  * c *  w     +
        K[10,:,:] * c *  n     +
        K[11,:,:] * t *  w     +
        K[12,:,:] * t *  n     +
        K[13,:,:] * w *  n     +
        K[14,:,:] * t ** 3     +
        K[15,:,:] * w ** 3     +
        K[16,:,:] * t *  w * n +
        K[17,:,:] * t ** 2 * w +
        K[18,:,:] * w ** 2 * t +
        K[19,:,:] * w ** 2 * n +
        K[20,:,:] * n ** 2 * c +
        K[21,:,:] * n ** 2 * t +
        K[22,:,:] * n ** 2 * w +
        K[23,:,:] * n ** 3      )

    elif (case == 'NIL'):
        Y = (K[0,:,:]          +
        K[1,:,:]  * c          +
        K[2,:,:]  * t          +
        K[3,:,:]  * n          +
        K[4,:,:]  * c ** 2     +
        K[5,:,:]  * t ** 2     +
        K[6,:,:]  * n ** 2     +
        K[7,:,:]  * c *  n     +
        K[8,:,:]  * t *  n     +
        K[9,:,:]  * t ** 3     +
        K[10,:,:] * n ** 2 * c +
        K[11,:,:] * n ** 2 * t +
        K[12,:,:] * n ** 3      )

    Y       = np.nan_to_num(Y)
    Y[Y<0]  = 0
    Y[Y>30] = 30
    return(Y)


def try_load_climate(climate_dir, var, GGCM, GGCMIcrop, rcp, missing_climate):
    # If no previous climate file was missing, load climate file, if it exists.
    # Otherwise, change missing_climate to True.
    in_array = -1
    if not missing_climate:
        climate_file = '%s/%s_%s_%s_rf_%d.npy'%(climate_dir, var, GCM, GGCMIcrop, rcp)
        try:
            in_array = np.load(climate_file)
        except FileNotFoundError:
            print('    File not found: %s'%(climate_file))
            missing_climate = True
    return(in_array, missing_climate)


def update_out_table(outarray, outheader, outfmt, PLUMcrop, rf_10, rf_60, rf_200, ir_10, ir_60, ir_200, mask_YX):
    # outarray comes in as (crop,cell) because vstack is the only stack that works with one-d
    # arrays like rf_10 etc., apparently. outarray will be transposed for write so that each
    # row is a gridcell.
    outarray = np.vstack((outarray, rf_10[mask_YX==1]))
    outarray = np.vstack((outarray, rf_60[mask_YX==1]))
    outarray = np.vstack((outarray, rf_200[mask_YX==1]))
    outarray = np.vstack((outarray, ir_10[mask_YX==1]))
    outarray = np.vstack((outarray, ir_60[mask_YX==1]))
    outarray = np.vstack((outarray, ir_200[mask_YX==1]))
    
    # Update header
    outheader = outheader + " " + PLUMcrop + "010 " + PLUMcrop + "060 " + PLUMcrop + "200 " + PLUMcrop + "i010 " + PLUMcrop + "i060 " + PLUMcrop + "i200"
    outfmt = outfmt + " %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f"

    return(outarray, outheader, outfmt)



def PLUMemulate(GCM, rcp, decade, GGCM, mask_YX, outarray):
    # GCM: ['ACCESS1-0','bcc-csm1-1','BNU-ESM','CanESM2','CCSM4','CESM1-BGC','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','FGOALS-g2',
    #      'GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A-LR',
    #      'IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']
    # rcp: 45 or 85
    # decade: 0,1,2,3,4,5,6,7,8  (which equates to 2010-2020, 2030-2040, ... 2090-2100)
    # GGCM: ['EPIC-TAMU','pDSSAT','LPJ-GUESS', 'LPJmL']
    # GGCMIcrop: ['maize' , 'rice', 'soy', 'winter_wheat', 'spring_wheat']

    # Begin strings for outfmt and outheader
    outfmt="%0.2f %0.2f"
    outheader = "Lon Lat"

    # Define PLUM output crops
    PLUMcrops = ["CerealsC3", "CerealsC4", "Rice", "Oilcrops", "Pulses", "StarchyRoots"];

    # Define matching to GGCMI crops
    PLUM2GGCMI_dict={
    "CerealsC3":    "max_wheat", # Should not ever actually be used
    "CerealsC4":    "maize",
    "Rice":         "rice",
    "Oilcrops":     "soy",
    "Pulses":       "soy",
    "StarchyRoots": "spring_wheat"
    }
    if GGCM == "LPJ-GUESS":
        PLUM2GGCMI_dict["Rice"] = "spring_wheat"
        PLUM2GGCMI_dict["Oilcrops"] = "spring_wheat"
        PLUM2GGCMI_dict["Pulses"] = "spring_wheat"
        PLUM2GGCMI_dict["StarchyRoots"] = "spring_wheat"

    # Set up for tracking of which crops were missing
    missing_climate_dict = {
    "max_wheat": False,
    "maize": False,
    "rice":  False,
    "soy":   False,
    "spring_wheat": False,
    }

    prev_GGCMIcrop = "nothing"
    prev_PLUMcrop = "nothing"
    for PLUMcrop in PLUMcrops:

        GGCMIcrop = PLUM2GGCMI_dict[PLUMcrop]
        print("%s rcp%d dec%d %s %s %s"%(GCM, rcp, decade, GGCM, PLUMcrop, GGCMIcrop))

        # Reuse data from previous PLUMcrop if it used the same GGCMIcrop.
        # Could make this more efficient by making it so ANY previous PLUMcrop
        # would work, but this will do for now. Just arrange dict to maximize
        # consecutive GGCMIcrops.
        if ((GGCMIcrop == prev_GGCMIcrop) and not (prev_PLUMcrop == "CerealsC3")):
            print('    %s and %s both use emulator for %s.'%(prev_PLUMcrop, PLUMcrop, prev_GGCMIcrop))
            if missing_climate_dict[prev_GGCMIcrop]:
                print('    But %s climate is missing, so skipping %s.'%(GGCMIcrop, PLUMcrop))
            else:
                print('    Using results from %s for %s.'%(prev_PLUMcrop, PLUMcrop))
                outarray,outheader,outfmt = update_out_table(outarray, outheader, outfmt, PLUMcrop, rf_10, rf_60, rf_200, ir_10, ir_60, ir_200, mask_YX)
            prev_GGCMIcrop = GGCMIcrop
            prev_PLUMcrop = PLUMcrop
            continue

        # Load Emulator params
        emulator_dir = "/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/crop"
        if PLUMcrop == "CerealsC3":
            GGCMIcrop_tmp = "winter_wheat"
            Kw  = np.load('%s/%s_%s.npy'%(emulator_dir, GGCM, GGCMIcrop_tmp))
            KIw = np.load('%s/%s_%s_I.npy'%(emulator_dir, GGCM, GGCMIcrop_tmp))
            GGCMIcrop_tmp = "spring_wheat"
            Ks  = np.load('%s/%s_%s.npy'%(emulator_dir, GGCM, GGCMIcrop_tmp))
            KIs = np.load('%s/%s_%s_I.npy'%(emulator_dir, GGCM, GGCMIcrop_tmp))
            del GGCMIcrop_tmp
        else:
            K  = np.load('%s/%s_%s.npy'%(emulator_dir, GGCM, GGCMIcrop))
            KI = np.load('%s/%s_%s_I.npy'%(emulator_dir, GGCM, GGCMIcrop))

        # Load climate files
        #climate_dir = "/project/ggcmi/AgMIP.output/Jim_Emulator/agmerra/cmip5"
        climate_dir = "/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/climate"
        missing_climate = False
        if PLUMcrop == "CerealsC3":
            GGCMIcrop_tmp = "winter_wheat"
            tw,missing_climate = try_load_climate(climate_dir, "tas", GGCM, GGCMIcrop_tmp, rcp, missing_climate)
            ww,missing_climate = try_load_climate(climate_dir, "pr", GGCM, GGCMIcrop_tmp, rcp, missing_climate)
            GGCMIcrop_tmp = "spring_wheat"
            ts,missing_climate = try_load_climate(climate_dir, "tas", GGCM, GGCMIcrop_tmp, rcp, missing_climate)
            ws,missing_climate = try_load_climate(climate_dir, "pr", GGCM, GGCMIcrop_tmp, rcp, missing_climate)
            del GGCMIcrop_tmp
        else:
            t,missing_climate = try_load_climate(climate_dir, "tas", GGCM, GGCMIcrop, rcp, missing_climate)
            w,missing_climate = try_load_climate(climate_dir, "pr", GGCM, GGCMIcrop, rcp, missing_climate)

        # If any climate files were missing, skip to the next crop.
        if missing_climate:
            missing_climate_dict[GGCMIcrop] = True
            prev_GGCMIcrop = GGCMIcrop
            prev_PLUMcrop = PLUMcrop
            print("    Skipping.")
            continue

        # CO2 vars
        if rcp == 45: co2 =[400.1285, 423.0875, 447.9455, 473.69, 497.703, 516.5865, 527.72, 532.4395, 536.0495]
        elif rcp == 85: co2 = [402.552, 432.3075, 469.135, 514.989, 572.0315, 640.299, 717.63, 801.4935, 890.3395]

        # Emulate the six management cases. If CerealsC3, emulate both winter and spring wheat,
        # then find the maximum at each treatment.
        if PLUMcrop == "CerealsC3":
            rf_10w  = emulate(Kw,  co2[decade], tw[decade,:,:],  ww[decade,:,:], 10,  'N')
            rf_10s  = emulate(Ks,  co2[decade], ts[decade,:,:],  ws[decade,:,:], 10,  'N')
            rf_10   = np.maximum(rf_10w, rf_10s)
            rf_60w  = emulate(Kw,  co2[decade], tw[decade,:,:],  ww[decade,:,:], 60,  'N')
            rf_60s  = emulate(Ks,  co2[decade], ts[decade,:,:],  ws[decade,:,:], 60,  'N')
            rf_60   = np.maximum(rf_60w, rf_60s)
            rf_200w  = emulate(Kw,  co2[decade], tw[decade,:,:],  ww[decade,:,:], 200,  'N')
            rf_200s  = emulate(Ks,  co2[decade], ts[decade,:,:],  ws[decade,:,:], 200,  'N')
            rf_200   = np.maximum(rf_200w, rf_200s)
            ir_10w  = emulate(Kw,  co2[decade], tw[decade,:,:],  1, 10,  'N')
            ir_10s  = emulate(Ks,  co2[decade], ts[decade,:,:],  1, 10,  'N')
            ir_10   = np.maximum(ir_10w, ir_10s)
            ir_60w  = emulate(Kw,  co2[decade], tw[decade,:,:],  1, 60,  'N')
            ir_60s  = emulate(Ks,  co2[decade], ts[decade,:,:],  1, 60,  'N')
            ir_60   = np.maximum(ir_60w, ir_60s)
            ir_200w  = emulate(Kw,  co2[decade], tw[decade,:,:],  1, 200,  'N')
            ir_200s  = emulate(Ks,  co2[decade], ts[decade,:,:],  1, 200,  'N')
            ir_200   = np.maximum(ir_200w, ir_200s)
        else:
            rf_10  = emulate(K,  co2[decade], t[decade,:,:],  w[decade,:,:], 10,  'N')
            rf_60  = emulate(K,  co2[decade], t[decade,:,:],  w[decade,:,:], 60,  'N')
            rf_200 = emulate(K,  co2[decade], t[decade,:,:],  w[decade,:,:], 200, 'N')
            ir_10  = emulate(KI, co2[decade], t[decade,:,:], 1, 10,  'NI')
            ir_60  = emulate(KI, co2[decade], t[decade,:,:], 1, 60,  'NI')
            ir_200 = emulate(KI, co2[decade], t[decade,:,:], 1, 200, 'NI')

        # outarray comes in as (crop,cell) because vstack is the only stack that works with one-d
        # arrays like rf_10 etc., apparently. outarray will be transposed for write so that each
        # row is a gridcell.
        outarray = np.vstack((outarray, rf_10[mask_YX==1]))
        outarray = np.vstack((outarray, rf_60[mask_YX==1]))
        outarray = np.vstack((outarray, rf_200[mask_YX==1]))
        outarray = np.vstack((outarray, ir_10[mask_YX==1]))
        outarray = np.vstack((outarray, ir_60[mask_YX==1]))
        outarray = np.vstack((outarray, ir_200[mask_YX==1]))

        # Update header
        outheader = outheader + " " + PLUMcrop + "010 " + PLUMcrop + "060 " + PLUMcrop + "200 " + PLUMcrop + "i010 " + PLUMcrop + "i060 " + PLUMcrop + "i200"
        outfmt = outfmt + " %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f"

        prev_GGCMIcrop = GGCMIcrop
        prev_PLUMcrop = PLUMcrop

    return(outarray, outheader, outfmt)

GCMs = ["IPSL-CM5A-MR"]
rcps = [45, 85]
decades = range(9)
GGCMs = ["EPIC-TAMU", "LPJ-GUESS", "LPJmL", "pDSSAT"]

# Import PLUM mask and lon/lat
plum_dir = "/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/plum/"
mask_YX = np.genfromtxt(plum_dir+"PLUM_mask.csv", delimiter=",")
lons_YX = np.genfromtxt(plum_dir+"PLUM_map_lons.csv", delimiter=",")
lats_YX = np.genfromtxt(plum_dir+"PLUM_map_lats.csv", delimiter=",")
lons = lons_YX[mask_YX==1]
lats = lats_YX[mask_YX==1]
lonlats = np.vstack((lons,lats))

for decade in decades:
    for rcp in rcps:
        for GCM in GCMs:
            for GGCM in GGCMs:

                # Set up info about this run
                decade_str = str(2011+10*decade) + "-" + str(2020+10*decade)
                print("%s rcp%d %s %s..."%(GCM, rcp, decade_str, GGCM))
                outdir = "/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/yields.nowheatIR/" + GCM + "/rcp" + str(rcp) + "/" + GGCM + "/" + decade_str + "/"
                outfile = outdir + "yield.nowheatIR." + decade_str + ".out"

                # Make output directory, if needed
                try:
                    os.makedirs(outdir)
                    print("mkdir -p " + outdir)
                except FileExistsError:
                    # directory already exists
                    pass

                # Emulate
                outarray,outheader,outfmt = PLUMemulate(GCM, rcp, decade, GGCM, mask_YX, lonlats)

                # Save in PLUM-readable format
                np.savetxt(outfile, outarray.T,
                delimiter=" ",
                fmt=outfmt,
                header=outheader,
                comments="")

















