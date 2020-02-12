#!/bin/bash

MODEL=$1
CROP=$2

DIR_LOCAL=${MODEL}/phase2/${CROP}/A0/yield
mkdir -p $DIR_LOCAL
DIR_REMOTE=/project/ggcmi/AgMIP.output/${DIR_LOCAL}

USERNAME=rabin
REMOTE=midway.rcc.uchicago.edu
TRANSFERTXT0="$USERNAME@$REMOTE:${DIR_REMOTE}/*_global_annual_1980_2010_C360_T0_*_N*_A0.nc4 ${DIR_LOCAL}/"
#echo $TRANSFERTXT0
#exit
FILELIST=$(rsync -vtn -e ssh $TRANSFERTXT0 | grep "C360" | grep "W[0i]")

###rsync -vtn -include='*_global_annual_1980_2010_C360_T0_W0_N*_A0.nc4' --include='*_global_annual_1980_2010_C360_T0_Winf_N*_A0.nc4' -e ssh rabin@midway.rcc.uchicago.edu:/project/ggcmi/AgMIP.output/LPJ-GUESS/phase2/maize/A0/yield/ LPJ-GUESS/phase2/maize/A0/yield/

for F in ${FILELIST[@]}; do

    # Do we need to rsync this file? This does a dry run and counts the number of files that would be transferred. (Should be either 0 or 1.)
    TRANSFERTXT="$USERNAME@$REMOTE:${DIR_REMOTE}/${F} ${DIR_LOCAL}/"
#    echo $TRANSFERTXT
    NOTDONEYET=$(rsync -vtn -e ssh $TRANSFERTXT | grep "C360" | wc -l)
#    echo $F
#    echo $NOTDONEYET
    NTRIES=0

    if [[ ${NOTDONEYET} -eq 1 ]]; then
       echo "Starting $F"
       while [[ ${NOTDONEYET} -eq 1 ]]; do
          NTRIES=$((NTRIES+1))
          if [[ ${NTRIES} -gt 1 ]]; then
             echo "Starting try ${NTRIES} (${THEPROBLEM}) of $F"
          fi
          set +e
          rsync -t --progress --timeout=$TIMEOUT --bwlimit=$BWLIMIT --partial --append -e ssh $TRANSFERTXT
    ###      rsync -vt --progress --timeout=$TIMEOUT --bwlimit=$BWLIMIT --partial --append -e ssh $TRANSFERTXT
          if [[ $? -eq 30 ]]; then
            THEPROBLEM=timeout
          elif [[ $? -gt 0 ]]; then
            THEPROBLEM="error $?"
          else
            THEPROBLEM=other
          fi
          set -e
          NOTDONEYET=$(rsync -vtn -e ssh $TRANSFERTXT | grep "C360" | wc -l)
       done
       echo "Done with $F"
    else
       echo "Skipping $F"
    fi

done

exit 0