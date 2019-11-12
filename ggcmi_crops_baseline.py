#!/bin/env python
import numpy as np;
import os
import glob
import datetime
np.random.seed(1234)
def emulate(K, c, t, w, n, case, is_irrig):
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

    # Limit output values
    Y[Y<0]  = 0
    if not is_irrig: Y[Y>30] = 30
    return(Y)


def do_emulation(emulator_dir, GGCMIcrop, co2, t, w, is_irrig):

    # Irrigated
    if is_irrig:
        KI = np.load('%s/%s_%s_irr.npy'%(emulator_dir, GGCM, GGCMIcrop))
    else:
        KI = np.load('%s/%s_%s_I.npy'%(emulator_dir, GGCM, GGCMIcrop))
    ir_10  = emulate(KI, co2, t, 1, 10,  'NI', is_irrig)
    ir_60  = emulate(KI, co2, t, 1, 60,  'NI', is_irrig)
    ir_200 = emulate(KI, co2, t, 1, 200, 'NI', is_irrig)

    # Rainfed
    if is_irrig:
        # No need to emulate irrigation for rainfed crops
        rf_10 = np.zeros(ir_10.shape)
        rf_60 = np.zeros(ir_10.shape)
        rf_200 = np.zeros(ir_10.shape)
    else:
        K  = np.load('%s/%s_%s.npy'%(emulator_dir, GGCM, GGCMIcrop))
        rf_10  = emulate(K,  co2, t,  w, 10,  'N', is_irrig)
        rf_60  = emulate(K,  co2, t,  w, 60,  'N', is_irrig)

        rf_200 = emulate(K,  co2, t,  w, 200, 'N', is_irrig)
    return(rf_10, rf_60, rf_200, ir_10, ir_60, ir_200)


def do_emulation_wheats(emulator_dir, co2, ts, ws, tw, ww, is_irrig):
    rf_10w,rf_60w,rf_200w,ir_10w,ir_60w,ir_200w = do_emulation(emulator_dir, "winter_wheat", co2, tw, ww, is_irrig)
    rf_10s,rf_60s,rf_200s,ir_10s,ir_60s,ir_200s = do_emulation(emulator_dir, "spring_wheat", co2, ts, ws, is_irrig)

    # Get maximum winter or spring
    if is_irrig:
        # No irrigation if rainfed
        rf_10 = np.zeros(ir_10w.shape)
        rf_60 = np.zeros(ir_10w.shape)
        rf_200 = np.zeros(ir_10w.shape)
    else:
        rf_10   = np.maximum(rf_10w, rf_10s)
        rf_60   = np.maximum(rf_60w, rf_60s)
        rf_200   = np.maximum(rf_200w, rf_200s)
    ir_10   = np.maximum(ir_10w, ir_10s)
    ir_60   = np.maximum(ir_60w, ir_60s)
    ir_200   = np.maximum(ir_200w, ir_200s)

    return(rf_10, rf_60, rf_200, ir_10, ir_60, ir_200)


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



def PLUMemulate(GGCM, mask_YX, outarr_yield, outarr_irrig):
    # GGCM: ['EPIC-TAMU','pDSSAT','LPJ-GUESS', 'LPJmL']
    # GGCMIcrop: ['maize' , 'rice', 'soy', 'winter_wheat', 'spring_wheat']

    # Begin strings for outfmt and outheader
    outfmt="%0.2f %0.2f"
    outheader = "Lon Lat"

    # Define crops
    if GGCM == "LPJ-GUESS":
        GGCMIcrops = ["winter_wheat", "spring_wheat", "maize"];
    else:
        GGCMIcrops = ["winter_wheat", "spring_wheat", "maize", "soy", "rice"];

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
    for GGCMIcrop in GGCMIcrops:

        print("%s %s"%(GGCM, GGCMIcrop))

        # Set climate
        t = 0.0
        w = 1.0

        # CO2 vars
        co2 = 360

        # Emulate the six management cases.
        emulator_dir_yield = "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/fits_yield"
        emulator_dir_irrig = "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/fits_irrig"
        rf_10_yield,rf_60_yield,rf_200_yield,ir_10_yield,ir_60_yield,ir_200_yield = do_emulation(emulator_dir_yield,GGCMIcrop, co2, t, w, False)
        rf_10_irrig,rf_60_irrig,rf_200_irrig,ir_10_irrig,ir_60_irrig,ir_200_irrig = do_emulation(emulator_dir_irrig,GGCMIcrop, co2, t, w, True)
        del t, w

        # Add values to output tables
        outarr_yield,outheader,outfmt = update_out_table(outarr_yield, outheader, outfmt, GGCMIcrop, rf_10_yield, rf_60_yield, rf_200_yield, ir_10_yield, ir_60_yield, ir_200_yield, mask_YX)
        outarr_irrig,x,y = update_out_table(outarr_irrig, outheader, outfmt, GGCMIcrop, rf_10_irrig, rf_60_irrig, rf_200_irrig, ir_10_irrig, ir_60_irrig, ir_200_irrig, mask_YX)


    return(outarr_yield, outarr_irrig, outheader, outfmt)
    
    
outdir_suffix = "20191008"

GGCMs = ["EPIC-TAMU", "LPJ-GUESS", "LPJmL", "pDSSAT"]
#GGCMs = ["LPJmL"]

# Import PLUM mask and lon/lat
plum_dir = "/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/plum/"
mask_YX = np.genfromtxt(plum_dir+"PLUM_mask.csv", delimiter=",")
lons_YX = np.genfromtxt(plum_dir+"PLUM_map_lons.csv", delimiter=",")
lats_YX = np.genfromtxt(plum_dir+"PLUM_map_lats.csv", delimiter=",")
lons = lons_YX[mask_YX==1]
lats = lats_YX[mask_YX==1]
lonlats = np.vstack((lons,lats))

for GGCM in GGCMs:

    # Set up info about this run
    print(GGCM)
    outdir = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_GGCMIcropsBaseline_%s/%s'%(outdir_suffix, GGCM)
    outfile_yield = '%s/yield.out'%(outdir)
    outfile_irrig = '%s/gsirrigation.out'%(outdir)
    
    # Skip if files already exist
    exist_yield = os.path.isfile(outfile_yield) or os.path.isfile(outfile_yield + ".gz")
    exist_irrig = os.path.isfile(outfile_irrig) or os.path.isfile(outfile_irrig + ".gz")
    if exist_yield and exist_irrig:
        print("Emulator outputs already exist; skipping.\n")
        continue

    # Emulate
    outarr_yield,outarr_irrig,outheader,outfmt = PLUMemulate(GGCM, mask_YX, lonlats, lonlats)
    
    # Make output directory, if needed
    try:
        os.makedirs(outdir)
        print("mkdir -p " + outdir)
    except FileExistsError:
        # directory already exists
        pass

    # Convert from tons/ha to kg/m2
    outarr_yield[:][2:] = outarr_yield[:][2:] * 0.1

    # Save in PLUM-readable format; compress
    if not exist_yield:
        np.savetxt(outfile_yield, outarr_yield.T, delimiter=" ", fmt=outfmt, header=outheader, comments="")
        os.system('gzip %s'%(outfile_yield))
    if not exist_irrig:
        np.savetxt(outfile_irrig, outarr_irrig.T, delimiter=" ", fmt=outfmt, header=outheader, comments="")
        os.system('gzip %s'%(outfile_irrig))
        
    # Compress and add "done" file for PLUM
    os.system('touch %s/done'%(outdir))
    
    print("\n")

















