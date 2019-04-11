list1 = {'NATURAL','CROPLAND','PASTURE','BARREN'} ;
list2 = {'ntrl','crop','past','bare'} ;

tmp = nansum(nansum(land_area_YX)) / nansum(nansum(gcel_area_YX)) ;

for L = 1:length(list1)
    thisLU = list1{L} ;
    thisLU_short = list2{L} ;
    for r = 2%1:Nruns
        disp(runList{r})
        this2010_YX = maps_LU_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,thisLU),end) ;
        this2011_YX = maps_LU_d1.maps_YXvyr(:,:,strcmp(maps_LU_d1.varNames,thisLU),1,r) ;
        disp('From maps:')
        fprintf('\t2010 %s: %0.3e\n', thisLU, nansum(nansum(this2010_YX.*gcel_area_YX))) ;
        fprintf('\t2011 %s: %0.3e\n', thisLU, nansum(nansum(this2011_YX.*gcel_area_YX))) ;
        disp('From TS:')
        fprintf('\t2010 %s: %0.3e\n', thisLU, eval(sprintf('ts_LUarea_%s_bl(end)',thisLU_short))) ;
        fprintf('\t2011 %s: %0.3e\n', thisLU, eval(sprintf('ts_LUarea_%s_yr(1,r)',thisLU_short))) ;
        disp(' ')
    end
end
    

% shademap(maps_LU_d1.maps_YXvyr(:,:,1,1,3)>0);
% title('This should NOT have cropland everywhere!')