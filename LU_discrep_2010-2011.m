thisLU = 'NATURAL' ;

for r = 1:Nruns
    disp(runList{r})
    this2010_YX = maps_LU_d9.maps_YXvyr(:,:,strcmp(maps_LU_d9.varNames,thisLU),end,r) ;
    this2011_YX = maps_LU_d1.maps_YXvyr(:,:,strcmp(maps_LU_d1.varNames,thisLU),1,r) ;
    disp('From maps:')
    fprintf('\t2010 %s: %0.3e\n', thisLU, nansum(nansum(this2010_YX.*land_area_YX))) ;
    fprintf('\t2011 %s: %0.3e\n', thisLU, nansum(nansum(this2011_YX.*land_area_YX))) ;
    disp('From TS:')
    fprintf('\t2010 %s: %0.3e\n', thisLU, ts_LUarea_ntrl_bl(end)) ;
    fprintf('\t2011 %s: %0.3e\n', thisLU, ts_LUarea_ntrl_yr(1,r)) ;
    disp(' ')
end
    

shademap(maps_LU_d1.maps_YXvyr(:,:,1,1,3)>0);
title('This should NOT have cropland everywhere!')