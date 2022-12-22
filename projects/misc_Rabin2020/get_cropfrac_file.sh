#!/bin/bash

thisDir=$1
CFfile=$(grep 'param "file_lucrop"' "$thisDir/landcover.ins" \
    | grep -v -e "^[[:blank:]]!" \
    | sed 's@  @ @' \
    | sed 's@param "file_lucrop" (str "/project/fh1-project-lpjgpi/lr8247/PLUM/input/PLUMouts_2011-2100/@@' \
    | sed 's@param "file_lucrop" (str "/project/fh1-project-lpjgpi/lr8247/input/PLUMouts_2011-2100/@@' \
    | sed 's@")@@')
echo $CFfile