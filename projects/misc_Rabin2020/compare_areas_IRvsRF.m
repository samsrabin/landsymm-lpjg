%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare irrigated vs. rainfed areas %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lu_file = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
cf_file = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.Misc0s.txt' ;


%% Import

% LPJ-GUESS input files
lu = lpjgu_matlab_readTable_then2map(lu_file) ;
cf = lpjgu_matlab_readTable_then2map(cf_file) ;

% Import land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd


%% Process

ca_YXcy = (lu.maps_YXvy(:,:,strcmp(lu.varNames,'CROPLAND'),:) .* repmat(land_area_YX,[1 1 1 length(lu.yearList)])) ...
    .* repmat(cf.maps_YXv,[1 1 1 length(lu.yearList)]) ;
isirr = ~cellfun(@isempty,regexp(cf.varNames,'i$')) ;

ca_ir = squeeze(nansum2(ca_YXcy(:,:,isirr,:),1:3)) ;
ca_rf = squeeze(nansum2(ca_YXcy(:,:,~isirr,:),1:3)) ;


%% Plot

figure ;
plot(lu.yearList,[ca_rf ca_ir]) ;
legend({'Rainfed','Irrigated'})