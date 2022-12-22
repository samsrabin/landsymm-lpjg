#!/bin/bash

MODELLIST=$(ssh rabin@midway.rcc.uchicago.edu ls -d /project/ggcmi/AgMIP.output/*/phase2 | sed "s@/project/ggcmi/AgMIP.output/@@" | sed "s@/phase2@@")

for M in ${MODELLIST[@]}; do
    CROPLIST=$(ssh rabin@midway.rcc.uchicago.edu ls /project/ggcmi/AgMIP.output/${M}/phase2/)
    for C in ${CROPLIST[@]}; do
        FILELIST=$(rsync -vtn -e ssh rabin@midway.rcc.uchicago.edu:/project/ggcmi/AgMIP.output/${M}/phase2/${C}/A0/yield/*_global_annual_1980_2010_C360_T0_*_N*_A0.nc4 ${M}/phase2/${C}/A0/yield/ | grep "C360" | grep "W[0i]")
        if [[ $FILELIST == "" ]]; then
            echo "${M}/${C} already up-to-date; skipping."
        else
            ./get_yields_phase2.sh $M $C
        fi
    done
done

exit 0