#!/bin/bash
set -e

thisdir=$1
if [[ "${thisdir}" == "" ]]; then
	echo "You must provide a directory."
	exit 1
fi

cd ${thisdir}
filelist=($(find . -name "*.out" | sort))
Nfiles=${#filelist[@]}

counter=0
for f in "${filelist[@]}"; do
	counter=$(( counter + 1 ))
	echo "${f} (${counter} of ${Nfiles})..."
	gzip ${f}
done

echo "Done!"

exit 0
