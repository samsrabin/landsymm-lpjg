#!/bin/env python
import numpy as np;
import os
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

def PLUMemulate(GCM, rcp, decade, GGCM, mask_YX, outarray):
    # GCM: ['ACCESS1-0','bcc-csm1-1','BNU-ESM','CanESM2','CCSM4','CESM1-BGC','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','FGOALS-g2',
    #      'GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A-LR',
    #      'IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']
    # rcp: 45 or 85
    # decade: 0,1,2,3,4,5,6,7,8  (which equates to 2010-2020, 2030-2040, ... 2090-2100)
    # GGCM: ['EPIC-TAMU','pDSSAT','LPJ-GUESS', 'LPJmL']
    # crop: ['maize' , 'rice', 'soy', 'winter_wheat', 'spring_wheat']

    #if ((GGCM == 'LPJ-GUESS') & (crop == 'rice')): print('error: LPJ-GUESS rice is not availible')
    #elif ((GGCM == 'LPJ-GUESS') & (crop == 'soy')): print('error: LPJ-GUESS soy is not availible')
    #else: pass

    # Begin strings for outfmt and outheader
    outfmt="%0.2f %0.2f"
    outheader = "Lon Lat"

    # Define PLUM output crops
    PLUMcrops = ["CerealsC3", "CerealsC4", "Rice", "Oilcrops", "Pulses", "StarchyRoots"];

    # Define matching to GGCMI crops
    PLUM2GGCMI_dict={
    #"CerealsC3":    "max_wheat",
    "CerealsC3":    "winter_wheat",
    "CerealsC4":    "maize",
    "Rice":        "rice",
    "Oilcrops":    "soy",
    "Pulses":    "soy",
    "StarchyRoots": "spring_wheat"
    }
    if GGCM == "LPJ-GUESS":
        PLUM2GGCMI_dict["Rice"] = "spring_wheat"
        PLUM2GGCMI_dict["Oilcrops"] = "spring_wheat"
        PLUM2GGCMI_dict["StarchyRoots"] = "spring_wheat"

#    # Define crop name list and dictionary for output crop names
#    if GGCM == "LPJ-GUESS":
#        crops = ["maize", "winter_wheat", "spring_wheat"]
#    else:
#        crops = ["maize" , "rice", "soy", "winter_wheat", "spring_wheat"]
#    
#    PLUMcrop_dict={
#        "maize":    "maiz",
#        "rice":        "rice",
#        "soy":        "soya",
#        "winter_wheat": "wwhe",
#        "spring_wheat": "swhe",
#    }

    for PLUMcrop in PLUMcrops:

        print("%s rcp%d dec%d %s %s"%(GCM, rcp, decade, GGCM, PLUMcrop))
        crop = PLUM2GGCMI_dict[PLUMcrop]

        # Could make this more efficient by not redoing emulation for two PLUMcrops that use the same GGCMI crop

        # Load Emulator params
        if PLUMcrop == "CerealsC3":
            tmp = "winter_wheat"
            Kw  = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/crop/%s_%s.npy'%(GGCM, tmp))
            KIw = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/crop/%s_%s_I.npy'%(GGCM, tmp))
            tmp = "spring_wheat"
            Ks  = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/crop/%s_%s.npy'%(GGCM, tmp))
            KIs = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/crop/%s_%s_I.npy'%(GGCM, tmp))
        else:
            K  = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/crop/%s_%s.npy'%(GGCM, crop))
            KI = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/crop/%s_%s_I.npy'%(GGCM, crop))

        # Load climate files
        t  = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/climate/rcp%s/tas_%s_%s_rf.npy'%(rcp, GCM, crop))
        if ((crop == 'winter_wheat') | (crop == 'spring_wheat')):
            tI = t
        else:
            tI = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/climate/rcp%s/tas_%s_%s_ir.npy'%(rcp, GCM, crop))
        w  = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/climate/rcp%s/pr_%s_%s_rf.npy'%(rcp, GCM, crop))

        # CO2 vars
        if rcp == 85: co2 = [402.552, 432.3075, 469.135, 514.989, 572.0315, 640.299, 717.63, 801.4935, 890.3395]
        elif rcp == 45: co2 =[400.1285, 423.0875, 447.9455, 473.69, 497.703, 516.5865, 527.72, 532.4395, 536.0495]

        # Emulate the six management cases. If CerealsC3, emulate both winter and spring wheat,
        # then find the maximum at each treatment.
        if PLUMcrop == "CerealsC3":
            rf_10w  = emulate(Kw,  co2[decade], t[decade,:,:],  w[decade,:,:], 10,  'N')
            rf_10s  = emulate(Ks,  co2[decade], t[decade,:,:],  w[decade,:,:], 10,  'N')
            rf_10   = np.maximum(rf_10w, rf_10s)
            rf_60w  = emulate(Kw,  co2[decade], t[decade,:,:],  w[decade,:,:], 60,  'N')
            rf_60s  = emulate(Ks,  co2[decade], t[decade,:,:],  w[decade,:,:], 60,  'N')
            rf_60   = np.maximum(rf_60w, rf_60s)
            rf_200w  = emulate(Kw,  co2[decade], t[decade,:,:],  w[decade,:,:], 200,  'N')
            rf_200s  = emulate(Ks,  co2[decade], t[decade,:,:],  w[decade,:,:], 200,  'N')
            rf_200   = np.maximum(rf_200w, rf_200s)
            ir_10w  = emulate(Kw,  co2[decade], t[decade,:,:],  w[decade,:,:], 10,  'N')
            ir_10s  = emulate(Ks,  co2[decade], t[decade,:,:],  w[decade,:,:], 10,  'N')
            ir_10   = np.maximum(ir_10w, ir_10s)
            ir_60w  = emulate(Kw,  co2[decade], t[decade,:,:],  w[decade,:,:], 60,  'N')
            ir_60s  = emulate(Ks,  co2[decade], t[decade,:,:],  w[decade,:,:], 60,  'N')
            ir_60   = np.maximum(ir_60w, ir_60s)
            ir_200w  = emulate(Kw,  co2[decade], t[decade,:,:],  w[decade,:,:], 200,  'N')
            ir_200s  = emulate(Ks,  co2[decade], t[decade,:,:],  w[decade,:,:], 200,  'N')
            ir_200   = np.maximum(ir_200w, ir_200s)
        else:
            rf_10  = emulate(K,  co2[decade], t[decade,:,:],  w[decade,:,:], 10,  'N')
            rf_60  = emulate(K,  co2[decade], t[decade,:,:],  w[decade,:,:], 60,  'N')
            rf_200 = emulate(K,  co2[decade], t[decade,:,:],  w[decade,:,:], 200, 'N')
            ir_10  = emulate(KI, co2[decade], tI[decade,:,:], 1, 10,  'NI')
            ir_60  = emulate(KI, co2[decade], tI[decade,:,:], 1, 60,  'NI')
            ir_200 = emulate(KI, co2[decade], tI[decade,:,:], 1, 200, 'NI')
        outarray = np.vstack((outarray, rf_10[mask_YX==1]))
        outarray = np.vstack((outarray, rf_60[mask_YX==1]))
        outarray = np.vstack((outarray, rf_200[mask_YX==1]))
        outarray = np.vstack((outarray, ir_10[mask_YX==1]))
        outarray = np.vstack((outarray, ir_60[mask_YX==1]))
        outarray = np.vstack((outarray, ir_200[mask_YX==1]))

        # Update header
#        PLUMcrop = PLUMcrop_dict.get(crop)
        outheader = outheader + " " + PLUMcrop + "010 " + PLUMcrop + "060 " + PLUMcrop + "200 " + PLUMcrop + "i010 " + PLUMcrop + "i060 " + PLUMcrop + "i200"
        outfmt = outfmt + " %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f"

    #return(ir_10, ir_60, ir_200, rf_10, rf_60, rf_200)
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
                #outfile = outdir + GCM + "_rcp" + str(rcp) + "_" + GGCM + "_" + crop + "_" + decade_str + ".csv"
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

















