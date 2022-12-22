#!/bin/bash

thisDir=$1
LUdir=$(grep 'param "file_lu"' $thisDir/landcover.ins \
    | grep -v -e "^[[:blank:]]!" \
    | sed 's@param "file_lu" (str "/project/fh1-project-lpjgpi/lr8247/PLUM/input/PLUMouts_2011-2100/@@' \
    | sed 's@param "file_lu" (str "/project/fh1-project-lpjgpi/lr8247/input/PLUMouts_2011-2100/@@' \
    | sed 's@/landcover.txt")@@')
echo $LUdir