function in_x = remap_combine_Nfert(thisCrop_in, croparea_in, nfert_dir, gridlist)

if ~strcmp(thisCrop_in, 'combined_sugars')
    error('Not recognized: %s', thisCrop_in)
end

cropfrac_is_beet = contains(croparea_in.varNames, 'Sugarbeet') ;
if sum(cropfrac_is_beet) ~= 2
    error('Expected 2 croparea_in.varNames containing Sugarbeet; found %d', ...
        sum(cropfrac_is_beet))
end
area_beet_x = sum(croparea_in.garr_xv(:,cropfrac_is_beet), 2) ;
cropfrac_is_cane = contains(croparea_in.varNames, 'Sugarcane') ;
if sum(cropfrac_is_cane) ~= 2
    error('Expected 2 croparea_in.varNames containing Sugarcane; found %d', ...
        sum(cropfrac_is_cane))
end
area_cane_x = sum(croparea_in.garr_xv(:,cropfrac_is_cane), 2) ;

area_sugar_x = area_beet_x + area_cane_x ;
area_sugar_x(area_sugar_x==0) = 1 ; % avoid /0
wt_beet_x = area_beet_x ./ area_sugar_x ;
wt_cane_x = area_cane_x ./ area_sugar_x ;

file_N_beet = fullfile(nfert_dir, ...
    sprintf('agmip_%s_apprate_fill_NPK_0.5.nc4', ...
    'sugarbeet')) ;
beet_YX = flipud(transpose(ncread(file_N_beet, 'Napprate'))) ;
file_N_cane = fullfile(nfert_dir, ...
    sprintf('agmip_%s_apprate_fill_NPK_0.5.nc4', ...
    'sugarcane')) ;
cane_YX = flipud(transpose(ncread(file_N_cane, 'Napprate'))) ;

beet_x = beet_YX(gridlist.list2map) ;
cane_x = cane_YX(gridlist.list2map) ;

in_x = beet_x.*wt_beet_x + cane_x.*wt_cane_x ;

end