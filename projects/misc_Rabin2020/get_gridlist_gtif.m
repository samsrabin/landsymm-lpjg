[~,R] = geotiffread('/Users/sam/Geodata/General/Countries RASTER/countries_raster_halfdeg.tif') ;
outDir = addslashifneeded('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work/') ;



load([outDir 'gridlist_mask.mat']) ;
gridlist_mask_dbl = double(gridlist_mask) ;
geotiffwrite([outDir 'gridlist_mask.tif'],flipud(gridlist_mask_dbl),R) ;