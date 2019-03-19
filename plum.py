#!/bin/env python
import numpy as np; from netCDF4 import Dataset
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
		Y = (K[0,:,:]              +
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
			K[22,:,:] * n ** 2 * w +
			K[23,:,:] * n ** 3      )

	elif (case == 'NIL'):
		Y = (K[0,:,:]              +
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

def PLUMemulate(GCM, rcp, decade, GGCM, crop):
	# GCM: ['ACCESS1-0','bcc-csm1-1','BNU-ESM','CanESM2','CCSM4','CESM1-BGC','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','FGOALS-g2',
			#'GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A-LR',
			#'IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']
	# rcp: 45 or 85
	# decade: 0,1,2,3,4,5,6,7,8  (which equates to 2010-2020, 2030-2040, ... 2090-2100)
	# GGCM: ['EPIC-TAMU','pDSSAT','LPJ-GUESS', 'LPJmL']
	# crop: ['maize' , 'rice', 'soy', 'winter_wheat', 'spring_wheat']

	if ((GGCM == 'LPJ-GUESS') & (crop == 'rice')): print('error: LPJ-GUESS rice is not availible')
	elif ((GGCM == 'LPJ-GUESS') & (crop == 'soy')): print('error: LPJ-GUESS soy is not availible')
	else: pass

	print(GCM + str(rcp) + str(decade) + GGCM + crop)

	# Load Emulator params
	print('Loading emulator parameters...')
	K  = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/crop/%s_%s.npy'%(GGCM, crop))
	KI = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/crop/%s_%s_I.npy'%(GGCM, crop))

	# Load climate files
	print('Loading climate files...')
	t  = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/climate/rcp%s/tas_%s_%s_rf.npy'%(rcp, GCM, crop))
	if ((crop == 'winter_wheat') | (crop == 'spring_wheat')): 
		tI = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/climate/rcp%s/tas_%s_%s_ir.npy'%(rcp, GCM, crop))
	else: tI = t
	w  = np.load('/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/climate/rcp%s/tas_%s_%s_rf.npy'%(rcp, GCM, crop))

	# CO2 vars
	if rcp == 85: co2 = [402.552, 432.3075, 469.135, 514.989, 572.0315, 640.299, 717.63, 801.4935, 890.3395]
	elif rcp == 45: co2 =[400.1285, 423.0875, 447.9455, 473.69, 497.703, 516.5865, 527.72, 532.4395, 536.0495]

	#for decade in range(9): might want to loop over decades at this point
	# Emulate the six management cases
	print('Emulating RF N1/3...')
	rf_10  = emulate(K,  co2[decade], t[decade,:,:],  w[decade,:,:], 10,  'N')
	print('Emulating RF N2/3...')
	rf_60  = emulate(K,  co2[decade], t[decade,:,:],  w[decade,:,:], 60,  'N')
	print('Emulating RF N3/3...')
	rf_200 = emulate(K,  co2[decade], t[decade,:,:],  w[decade,:,:], 200, 'N')
	print('Emulating IR N1/3...')
	ir_10  = emulate(KI, co2[decade], tI[decade,:,:], 1,             10,  'NI')
	print('Emulating IR N2/3...')
	ir_60  = emulate(KI, co2[decade], tI[decade,:,:], 1,             60,  'NI')
	print('Emulating IR N3/3...')
	ir_200 = emulate(KI, co2[decade], tI[decade,:,:], 1,             200, 'NI')

	return(ir_10, ir_60, ir_200, rf_10, rf_60, rf_200)

### Random test ###
GCM = 'ACCESS1-0'
rcp = 85
decade = 8
GGCM = 'pDSSAT'
crop = 'maize'

# Make output directory, if needed
outdir = "/project/ggcmi/AgMIP.output/Jim_Emulator/Sam/yields/" + GCM + "/rcp" + str(rcp) + "/" + GGCM + "/" + crop + "/"
try:
    os.makedirs(outdir)
    print("mkdir -p " + outdir)
except FileExistsError:
    # directory already exists
    pass

# Emulate
ir_10,ir_60,ir_200,rf_10,rf_60,rf_200 = PLUMemulate(GCM, rcp, decade, GGCM, crop)

# Save
outfile_base = outdir + GCM + "_rcp" + str(rcp) + "_" + GGCM + "_"
print(outfile_base)
#np.savetxt(outfile_base+crop+str(decade)+".csv", ir_10, delimiter=',')















