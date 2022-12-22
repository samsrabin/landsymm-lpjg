#!/bin/bash

thisDir=$1
LUfile=$(grep 'param "file_nfert"' "$thisDir/landcover.ins" \
    | grep -v -e "^[[:blank:]]!" \
    | sed 's@param "file_nfert" (str "/project/fh1-project-lpjgpi/lr8247/PLUM/input/PLUMouts_2011-2100/@@' \
    | sed 's@param "file_nfert" (str "/project/fh1-project-lpjgpi/lr8247/input/PLUMouts_2011-2100/@@' \
    | sed 's@param "file_nfert" (str "/project/fh1-project-lpjgpi/lr8247/input/Nfert/@@' \
    | sed 's@param "file_nfert" (str "/project/fh1-project-lpjgpi/lr8247/PLUM/input/Nfert/@@' \
    | sed 's@")@@')
echo $LUfile