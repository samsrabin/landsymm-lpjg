#!/bin/bash
set -e

d1=/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/MATLAB_work
d2=/Users/sam/Documents/git_repos/g2p_emulation

for f in make_justCPB.m make_himelo_filled.m make_himelo.m removeCols_addZeros.m get_remapv2_keys.m; do

	if [[ $(find . -name "${f}" | wc -l ) -lt 2 ]]; then
		continue
	fi
	
	f1=$(find $d1 -name $f)
	f2=$(find $d2 -name $f)
	echo $f
	ls -lh $f1
	ls -lh $f2
	echo " "
done

exit 0