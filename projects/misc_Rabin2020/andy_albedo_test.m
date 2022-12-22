
baresoil_albedo_file = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work/soilmap.txt' ;

fpc_file = '/Users/sam/Geodata/Andy_albedo_test/MAgPIE_BIOADAFF/fpc.out.gz' ;
snowdepth_file = '/Users/sam/Geodata/Andy_albedo_test/MAgPIE_BIOADAFF/msnow.out.gz' ;
% LU_file = '/Users/sam/Geodata/Andy_albedo_test/MAgPIE_BIOADAFF/MAgPIE_MIT-BIO-AD-AFF-RCP26-co2_LPJG_landuse_1901_2100_midpoint.txt.gz' ;

% fpc_file = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_remap6p3/output-2018-12-01-082125/fpc.out' ;
% snowdepth_file = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_remap6p3/output-2018-12-01-082125/msnowdepth.out' ;
LU_file = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6p3/LU.remapv6p3.txt' ;


%% Read baresoil albedo
disp('Reading baresoil_albedo...')

baresoil_albedo = lpjgu_matlab_readTable_then2map(baresoil_albedo_file,'force_mat_save',true,'verbose',false) ;
baresoil_albedo_YX = baresoil_albedo.maps_YXv ;
clear baresoil_albedo

disp('Done.')


%% Land area
disp('Reading land area...')

% Import land area (km2 to m2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;

%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
tmp = gcel_area_YXqd(:,1:2:1440) + gcel_area_YXqd(:,2:2:1440) ;
gcel_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;

% Convert to m2
land_area_YX = land_area_YX*1e6 ;
gcel_area_YX = gcel_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd

disp('Done.')


%% Read LU
disp('Reading LU...')

LU = lpjgu_matlab_readTable_then2map(LU_file,'force_mat_save',true) ;

LU.maps_YXvy(:,:,:,LU.yearList<2000 | LU.yearList>2009) = [] ;
LU.yearList(LU.yearList<2000 | LU.yearList>2009) = [] ;

disp('Done.')


%% Read FPC
disp('Reading and processing FPC...')

fpc_orig = lpjgu_matlab_readTable_then2map(fpc_file,'force_mat_save',true) ;
fpc_orig.maps_YXvy(:,:,:,fpc_orig.yearList<2000 | fpc_orig.yearList>2009) = [] ;
fpc_orig.yearList(fpc_orig.yearList<2000 | fpc_orig.yearList>2009) = [] ;

% [fpc, ~, ~] = CrOp_and_CrOpi(fpc, 'fpc', cropTypes, merge_or_replace) ;

% Normalize to 1
fpc = fpc_orig ;
fpc_total_YX1y = fpc.maps_YXvy(:,:,strcmp(fpc.varNames,'Total'),:) ;
isbad_YX1y = fpc_total_YX1y > 1 ;
i = 0 ;
while any(isbad_YX1y(:))
    i = i + 1 ;
    if i > 5
        error('Too many iterations!')
    end
    isbad_YXvy = repmat(isbad_YX1y,[1 1 size(fpc.maps_YXvy,3) 1]) ;
    fpc_total_YXvy = repmat(fpc_total_YX1y,[1 1 size(fpc.maps_YXvy,3) 1]) ;
    fpc.maps_YXvy(isbad_YXvy) = fpc.maps_YXvy(isbad_YXvy) ./ fpc_total_YXvy(isbad_YXvy) ;
    clear fpc_total_YXvy
    fpc_total_YX1y = fpc.maps_YXvy(:,:,strcmp(fpc.varNames,'Total'),:) ;
    isbad_YX1y = fpc_total_YX1y > 1 ;
end
clear fpc_total_YX1y isbad* i


disp('Done.')


%% Read snow depth
disp('Reading snowdepth...')

snowdepth = lpjgu_matlab_readTable_then2map(snowdepth_file,'force_mat_save',true) ;
snowdepth.maps_YXvy(:,:,:,snowdepth.yearList<2000 | snowdepth.yearList>2009) = [] ;
snowdepth.yearList(snowdepth.yearList<2000 | snowdepth.yearList>2009) = [] ;

disp('Done.')


%% Get albedo: original method

[albedo1_jan_YXy,albedo1_jul_YXy] = ...
    get_albedo(fpc, snowdepth, LU, baresoil_albedo_YX, land_area_YX, fpc.varNames) ;


%% Get albedo: method 2

% [albedo2_jan_YXy,albedo2_jul_YXy] = ...
%     get_albedo_nobare(fpc, snowdepth, LU, baresoil_albedo_YX, land_area_YX, fpc.varNames) ;


%% Get albedo: method 3

[albedo3_jan_YXy,albedo3_jul_YXy] = ...
    get_albedo3(fpc, snowdepth, LU, baresoil_albedo_YX, land_area_YX, fpc.varNames) ;


%%

LU_tmp = LU ;
LU_tmp.maps_YXvy = ...
    LU.maps_YXvy(:,:,contains(LU.varNames,{'CROPLAND','PASTURE','NATURAL'}),:) ...
    ./ sum(LU.maps_YXvy(:,:,contains(LU.varNames,{'CROPLAND','PASTURE','NATURAL'}),:),3) ;
LU_tmp.varNames(strcmp(LU_tmp.varNames,'BARREN')) = [] ;
    

[albedo3_jan_YXy,albedo3_jul_YXy] = ...
    get_albedo3(fpc, snowdepth, LU_tmp, baresoil_albedo_YX, land_area_YX, fpc.varNames) ;




