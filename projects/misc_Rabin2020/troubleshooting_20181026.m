%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Troubleshooting 2018-10-26 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisCell = [74.75 62.75] ;


%% Setup

lons_map = 0.25 + repmat(-180:0.5:179.5,[360 1]) ;
lats_map = 0.25 + repmat((-90:0.5:89.5)',[1 720]) ;

thisLon = thisCell(1) ;
thisLat = thisCell(2) ;

[i j] = find(lons_map==thisLon & lats_map==thisLat) ;

file_lu = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/LU.remapv6.20180214.ecFertIrr0.setaside0103.m4.someOfEachCrop.txt' ;
file_cf = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/cropfracs.remapv6.20180214.ecFertIrr0.setaside0103.m4.someOfEachCrop.txt' ;


%% Import

lu = lpjgu_matlab_readTable(file_lu) ;
cf = lpjgu_matlab_readTable(file_cf) ;

yearList = unique(lu.Year) ;
if ~isequal(yearList,unique(cf.Year))
    error('Intersect years')
end
Nyears = length(yearList) ;


%% List values

thisLU = 'CROPLAND' ;
thisCrop = 'Miscanthus' ;

thisLU_i = find(strcmp(lu.Properties.VariableNames,thisLU)) ;
thisCrop_i = find(strcmp(cf.Properties.VariableNames,thisCrop)) ;

for y = 1:Nyears
    thisYear = yearList(y) ;
    thisArea_thisLU = table2array(lu(lu.Lon==thisLon & lu.Lat==thisLat & lu.Year==thisYear,thisLU_i)) ;
    thisArea_thisCrop = table2array(cf(cf.Lon==thisLon & cf.Lat==thisLat & cf.Year==thisYear,thisCrop_i)) ;
    fprintf('%d: CROP=%0.4e, %s=%0.4e, combined=%0.4e\n',...
        thisYear, thisArea_thisLU, thisCrop, thisArea_thisCrop, thisArea_thisLU*thisArea_thisCrop)
end