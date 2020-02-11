#!/bin/env python
import numpy as np
import os


def emulate_old(K, c, t, w, n, case, is_irrig):
    '''
    Cases:
        NN: no nitrogen
        NNI: no nitrogen, irrigated
        N: with nitrogen
        NI: with nitrogen, irrigated
    '''
    if case == "NN":
        Y = (
            K[0, :, :]
            + K[1, :, :] * c
            + K[2, :, :] * t
            + K[3, :, :] * w
            + K[4, :, :] * c ** 2
            + K[5, :, :] * t ** 2
            + K[6, :, :] * w ** 2
            + K[7, :, :] * c * w
            + K[8, :, :] * t * w
            + K[9, :, :] * t ** 3
            + K[10, :, :] * w ** 3
            + K[11, :, :] * t ** 2 * w
            + K[12, :, :] * w ** 2 * t
        )

    elif case == "NNI":
        Y = (
            K[0, :, :]
            + K[1, :, :] * c
            + K[2, :, :] * t
            + K[3, :, :] * c ** 2
            + K[4, :, :] * t ** 2
            + K[5, :, :] * t ** 3
        )

    elif case == "N":
        Y = (
            K[0, :, :]
            + K[1, :, :] * c
            + K[2, :, :] * t
            + K[3, :, :] * w
            + K[4, :, :] * n
            + K[5, :, :] * c ** 2
            + K[6, :, :] * t ** 2
            + K[7, :, :] * w ** 2
            + K[8, :, :] * n ** 2
            + K[9, :, :] * c * w
            + K[10, :, :] * c * n
            + K[11, :, :] * t * w
            + K[12, :, :] * t * n
            + K[13, :, :] * w * n
            + K[14, :, :] * t ** 3
            + K[15, :, :] * w ** 3
            + K[16, :, :] * t * w * n
            + K[17, :, :] * t ** 2 * w
            + K[18, :, :] * w ** 2 * t
            + K[19, :, :] * w ** 2 * n
            + K[20, :, :] * n ** 2 * c
            + K[21, :, :] * n ** 2 * t
            + K[22, :, :] * n ** 2 * w
        )

    elif case == "NI":
        Y = (
            K[0, :, :]
            + K[1, :, :] * c
            + K[2, :, :] * t
            + K[3, :, :] * n
            + K[4, :, :] * c ** 2
            + K[5, :, :] * t ** 2
            + K[6, :, :] * n ** 2
            + K[7, :, :] * c * n
            + K[8, :, :] * t * n
            + K[9, :, :] * t ** 3
            + K[10, :, :] * n ** 2 * c
            + K[11, :, :] * n ** 2 * t
        )

    elif case == "NL":
        Y = (
            K[0, :, :]
            + K[1, :, :] * c
            + K[2, :, :] * t
            + K[3, :, :] * w
            + K[4, :, :] * n
            + K[5, :, :] * c ** 2
            + K[6, :, :] * t ** 2
            + K[7, :, :] * w ** 2
            + K[8, :, :] * n ** 2
            + K[9, :, :] * c * w
            + K[10, :, :] * c * n
            + K[11, :, :] * t * w
            + K[12, :, :] * t * n
            + K[13, :, :] * w * n
            + K[14, :, :] * t ** 3
            + K[15, :, :] * w ** 3
            + K[16, :, :] * t * w * n
            + K[17, :, :] * t ** 2 * w
            + K[18, :, :] * w ** 2 * t
            + K[19, :, :] * w ** 2 * n
            + K[20, :, :] * n ** 2 * c
            + K[21, :, :] * n ** 2 * t
            + K[22, :, :] * n ** 2 * w
            + K[23, :, :] * n ** 3
        )

    elif case == "NIL":
        Y = (
            K[0, :, :]
            + K[1, :, :] * c
            + K[2, :, :] * t
            + K[3, :, :] * n
            + K[4, :, :] * c ** 2
            + K[5, :, :] * t ** 2
            + K[6, :, :] * n ** 2
            + K[7, :, :] * c * n
            + K[8, :, :] * t * n
            + K[9, :, :] * t ** 3
            + K[10, :, :] * n ** 2 * c
            + K[11, :, :] * n ** 2 * t
            + K[12, :, :] * n ** 3
        )

    Y = np.nan_to_num(Y)

    # Limit output values
    Y[Y < 0] = 0
    if not is_irrig:
        Y[Y > 30] = 30
    return Y


def emulate(K, C, T, W, N):
    if K.shape[0] == 34:
        Y = (K[0,:,:] + K[1,:,:]*C + K[2,:,:]*T + K[3,:,:]*W + K[4,:,:]*N + K[5,:,:]*C**2
             + K[6,:,:]*C*T + K[7,:,:]*C*W + K[8,:,:]*C*N + K[9,:,:]*T**2 + K[10,:,:]*T*W
             + K[11,:,:]*T*N + K[12,:,:]*W**2 + K[13,:,:]*W*N + K[14,:,:]*N**2
             + K[15,:,:]*C**3 + K[16,:,:]*C**2*T + K[17,:,:]*C**2*W + K[18,:,:]*C**2*N
             + K[19,:,:]*C*T**2 + K[20,:,:]*C*T*W + K[21,:,:]*C*T*N + K[22,:,:]*C*W**2
             + K[23,:,:]*C*W*N + K[24,:,:]*C*N**2 + K[25,:,:]*T**3 + K[26,:,:]*T**2*W
             + K[27,:,:]*T**2*N + K[28,:,:]*T*W**2 + K[29,:,:]*T*W*N + K[30,:,:]*T*N**2
             + K[31,:,:]*W**3 + K[32,:,:]*W**2*N + K[33,:,:]*W*N**2)

    elif K.shape[0] == 20:
        Y = (K[0,:,:] + K[1,:,:]*C + K[2,:,:]*T + K[3,:,:]*W + K[4,:,:]*C**2 + K[5,:,:]*C*T
             + K[6,:,:]*C*W + K[7,:,:]*T**2 + K[8,:,:]*T*W + K[9,:,:]*W**2 + K[10,:,:]*C**3
             + K[11,:,:]*C**2*T + K[12,:,:]*C**2*W + K[13,:,:]*C*T**2 + K[14,:,:]*C*T*W
             + K[15,:,:]*C*W**2 + K[16,:,:]*T**3 + K[17,:,:]*T**2*W + K[18,:,:]*T*W**2
             + K[19,:,:]*W**3)

    elif K.shape[0] == 19:
        Y = (K[0,:,:] + K[1,:,:]*C + K[2,:,:]*T + K[3,:,:]*N + K[4,:,:]*C**2 + K[5,:,:]*C*T
             + K[6,:,:]*C*N + K[7,:,:]*T**2 + K[8,:,:]*T*N + K[9,:,:]*N**2 + K[10,:,:]*C**3
             + K[11,:,:]*C**2*T + K[12,:,:]*C**2*N + K[13,:,:]*C*T**2 + K[14,:,:]*C*T*N
             + K[15,:,:]*C*N**2 + K[16,:,:]*T**3 + K[17,:,:]*T**2*N + K[18,:,:]*T*N**2)

    elif K.shape[0] == 10:
        Y = (K[0,:,:] + K[1,:,:]*C + K[2,:,:]*T + K[3,:,:]*C**2 + K[4,:,:]*C*T
             + K[5,:,:]*T**2 + K[6,:,:]*C**3 + K[7,:,:]*C**2*T + K[8,:,:]*C*T**2
             + K[9,:,:]*T**3)
    else:
        raise Exception("K.shape[0] not recognized: %d" % K.shape[0])

    Y = np.nan_to_num(Y)
    Y[Y < 0.01] = 0
    return Y


def do_emulation(emulator_dir, GGCMIcrop, co2, t, w, is_irrig, GGCM, decade):
    if is_irrig:
        KI = np.load("%s/%s_%s_irr.npy" % (emulator_dir, GGCM, GGCMIcrop))
    else:
        KI = np.load("%s/%s_%s_I.npy" % (emulator_dir, GGCM, GGCMIcrop))
    ir_10 = emulate(KI, co2[decade], t[decade, :, :], 1, 10)
    ir_60 = emulate(KI, co2[decade], t[decade, :, :], 1, 60)
    ir_200 = emulate(KI, co2[decade], t[decade, :, :], 1, 200)

    # Rainfed
    if is_irrig:
        # No need to emulate irrigation for rainfed crops
        rf_10 = np.zeros(ir_10.shape)
        rf_60 = np.zeros(ir_10.shape)
        rf_200 = np.zeros(ir_10.shape)
    else:
        K = np.load("%s/%s_%s.npy" % (emulator_dir, GGCM, GGCMIcrop))
        rf_10 = emulate(
            K, co2[decade], t[decade, :, :], w[decade, :, :], 10
        )
        rf_60 = emulate(
            K, co2[decade], t[decade, :, :], w[decade, :, :], 60
        )
        rf_200 = emulate(
            K, co2[decade], t[decade, :, :], w[decade, :, :], 200
        )
    return (rf_10, rf_60, rf_200, ir_10, ir_60, ir_200)


def update_out_table(
    outarray,
    outheader,
    outfmt,
    PLUMcrop,
    rf_10,
    rf_60,
    rf_200,
    ir_10,
    ir_60,
    ir_200,
    mask_YX
):
    # outarray comes in as (crop,cell) because vstack is the only stack that works with one-d
    # arrays like rf_10 etc., apparently. outarray will be transposed for write so that each
    # row is a gridcell.
    outarray = np.vstack((outarray, rf_10[mask_YX == 1]))
    outarray = np.vstack((outarray, rf_60[mask_YX == 1]))
    outarray = np.vstack((outarray, rf_200[mask_YX == 1]))
    outarray = np.vstack((outarray, ir_10[mask_YX == 1]))
    outarray = np.vstack((outarray, ir_60[mask_YX == 1]))
    outarray = np.vstack((outarray, ir_200[mask_YX == 1]))

    # Update header
    outheader = (
        outheader + " "
        + PLUMcrop + "010 "
        + PLUMcrop + "060 "
        + PLUMcrop + "200 "
        + PLUMcrop + "i010 "
        + PLUMcrop + "i060 "
        + PLUMcrop + "i200"
    )
    outfmt = outfmt + " %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f"

    return (outarray, outheader, outfmt)


def save_out_table(
        outdir,
        outfile,
        outfmt,
        outheader,
        outarr):

    # Make output directory, if needed
    try:
        os.makedirs(outdir)
        print("mkdir -p " + outdir)
    except FileExistsError:
        # directory already exists
        pass

    # Save in PLUM-readable format; compress
    np.savetxt(outfile, outarr.T, delimiter=" ", fmt=outfmt, header=outheader, comments="")
    os.system('gzip %s' % (outfile))
