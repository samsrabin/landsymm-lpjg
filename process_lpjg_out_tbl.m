%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Work with LPJ-GUESS outputs for impacts paper %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inDir = '/Users/Shared/PLUM/trunk_runs/PLUM2LPJG_26_s1_2011-2100/output-2017-11-21-122520' ;
% inDir = '/Users/Shared/PLUM/trunk_runs/PLUM2LPJG_45_s1_2011-2100/output-2017-11-22-052523' ;
% inDir = '/Users/Shared/PLUM/trunk_runs/PLUM2LPJG_60_s1_2011-2100/output-2017-11-22-120456' ;

%% Setup

inDir = addslashifneeded(inDir) ;

if strcmp(inDir,'/Users/Shared/PLUM/trunk_runs/PLUM2LPJG_26_s1_2011-2100/output-2017-11-21-122520/')
    LUfile = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/26_s1.forLPJG.MATLAB/landcover.txt' ;
    cropfile = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/26_s1.forLPJG.MATLAB/cropfractions.txt' ;
else
    error('Specify LUfile and cropfile for this inDir!')
end

cd '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work'
addpath(genpath(pwd))


%% Import

% Carbon: Vegetation, soil/litter, product, total (kgC/m2)
cpool = lpjgu_matlab_readTable([inDir 'cpool.out'],'do_save_mat',true) ;
cpool_map = lpjgu_matlab_readTable_then2map([inDir 'cpool.out'],'force_mat_save',true) ;

% Get mapping info
list2map = cpool_map.list_to_map ;
lons = cpool.Lon(cpool.Year==min(cpool.Year)) ;
lats = cpool.Lat(cpool.Year==min(cpool.Year)) ;

% Evapotranspiration (mm/yr)
aaet = lpjgu_matlab_readTable_then2map([inDir 'aaet.out'],'do_save_mat',true) ;

% Runoff (mm/yr)
tot_runoff = lpjgu_matlab_readTable_then2map([inDir 'tot_runoff.out'],'do_save_mat',true) ;
% Runoff (mm/month)
mon_runoff = lpjgu_matlab_readTable_then2map([inDir 'mrunoff.out'],'do_save_mat',true) ;

% BVOCs (isoprene, monoterpenes) (mgC/m2)
aiso = lpjgu_matlab_readTable_then2map([inDir 'aiso.out'],'do_save_mat',true) ;
amon = lpjgu_matlab_readTable_then2map([inDir 'amon.out'],'do_save_mat',true) ;

% N flux (kgN/ha)
nflux = lpjgu_matlab_readTable_then2map([inDir 'nflux.out'],'do_save_mat',true) ;

% Foliar projective cover
fpc = lpjgu_matlab_readTable_then2map([inDir 'fpc.out'],'do_save_mat',true) ;
is_ntrl_tree = strcmp(fpc.varNames,'BNE') | strcmp(fpc.varNames,'BINE') ...
             | strcmp(fpc.varNames,'BNS') | strcmp(fpc.varNames,'TeNE') ...
             | strcmp(fpc.varNames,'TeBS') | strcmp(fpc.varNames,'IBS') ...
             | strcmp(fpc.varNames,'TeBE') | strcmp(fpc.varNames,'TrBe') ...
             | strcmp(fpc.varNames,'TrIBE') | strcmp(fpc.varNames,'TrBR') ;
is_ntrl_grass = strcmp(fpc.varNames,'C3G') | strcmp(fpc.varNames,'C4G') ;
is_ntrl = is_ntrl_tree | is_ntrl_grass ;
% Normalize to 1
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

% Import land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calibration_for_PLUM/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
land_area_YX_m2 = land_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd

% Import land use and crop fractions
LU = lpjgu_matlab_readTable_then2map(LUfile,'do_save_mat',true) ;
if ~isequal(LU.yearList,fpc.yearList)
    error('Rework this so yearLists match!')
end
cropfracs = lpjgu_matlab_readTable_then2map(cropfile,'do_save_mat',true) ;
if ~isequal(cropfracs.yearList,fpc.yearList)
    error('Rework this so yearLists match!')
end

% Import bare-soil albedo
baresoil_albedo_file = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work/soilmap.txt' ;
baresoil_albedo = lpjgu_matlab_readTable_then2map(baresoil_albedo_file,'force_mat_save',true) ;
baresoil_albedo_YX = baresoil_albedo.maps_YXv ;
clear baresoil_albedo

disp('Done.')



%% Get 2091-2100 values

y1 = 2091 ;
yN = 2100 ;
theseYears = (cpool.yearList>=y1 & cpool.yearList<=yN) ;

disp('Carbon:')

tmp = mean(cpool.maps_YXvy(:,:,strcmp(cpool.varNames,'VegC'),theseYears),4) .* land_area_YX_m2 ;
tmp = sum(tmp(~isnan(tmp))) * 1e-3*1e-9 ;
fprintf('%.1f\t%s\n',tmp,'Vegetation C (Gt)')
clear tmp

tmp1 = mean(cpool.maps_YXvy(:,:,strcmp(cpool.varNames,'LitterC'),theseYears),4) .* land_area_YX_m2 ;
tmp2 = mean(cpool.maps_YXvy(:,:,strcmp(cpool.varNames,'SoilC'),theseYears),4) .* land_area_YX_m2 ;
tmp = (sum(tmp1(~isnan(tmp1))) + sum(tmp2(~isnan(tmp2)))) * 1e-3*1e-9 ;
fprintf('%.1f\t%s\n',tmp,'Soil+Litter C (Gt)')
clear tmp

tmp = mean(cpool.maps_YXvy(:,:,strcmp(cpool.varNames,'HarvSlowC'),theseYears),4) .* land_area_YX_m2 ;
tmp = sum(tmp(~isnan(tmp))) * 1e-3*1e-9 ;
fprintf('%.1f\t%s\n',tmp,'Product C (Gt)')
clear tmp

tmp = mean(cpool.maps_YXvy(:,:,strcmp(cpool.varNames,'Total'),theseYears),4) .* land_area_YX_m2 ;
tmp = sum(tmp(~isnan(tmp))) * 1e-3*1e-9 ;
fprintf('%.1f\t%s\n',tmp,'Total C (Gt)')
clear tmp

disp(' ')
disp('Water:')

tmp = mean(aaet.maps_YXvy(:,:,strcmp(aaet.varNames,'Total'),theseYears),4) .* land_area_YX ;
tmp = sum(tmp(~isnan(tmp))) * 1e-3*1e-3*1e-3 ;
fprintf('%.1f\t%s\n',tmp,'Annual evapotranspiration (1000 km3/yr)')
clear tmp

tmp = mean(tot_runoff.maps_YXvy(:,:,strcmp(tot_runoff.varNames,'Total'),theseYears),4) .* land_area_YX ;
tmp = sum(tmp(~isnan(tmp))) * 1e-3*1e-3*1e-3 ;
fprintf('%.1f\t%s\n',tmp,'Annual runoff (1000 km3/yr)')
clear tmp

tmp = max(mean(mon_runoff.maps_YXvy(:,:,:,theseYears),4),[],3) .* land_area_YX ;
tmp = sum(tmp(~isnan(tmp))) * 1e-3*1e-3*1e-3 ;
fprintf('%.1f\t%s\n',tmp,'Peak monthly runoff (1000 km3/month)')
clear tmp

disp(' ')
disp('N loss:')

% kgN/ha/yr --> TgN/yr
tmp = mean(nflux.maps_YXvy(:,:,strcmp(nflux.varNames,'flux'),theseYears),4) .* (land_area_YX_m2*1e-4) ;
tmp = sum(tmp(~isnan(tmp))) * 1e-9 ;
fprintf('%.1f\t%s\n',tmp,'Gaseous (TgN/yr)')
tmp = mean(nflux.maps_YXvy(:,:,strcmp(nflux.varNames,'leach'),theseYears),4) .* (land_area_YX_m2*1e-4) ;
tmp = sum(tmp(~isnan(tmp))) * 1e-9 ;
fprintf('%.1f\t%s\n',tmp,'Dissolved (TgN/yr)')
tmp = sum(mean(nflux.maps_YXvy(:,:,strcmp(nflux.varNames,'harvest')|strcmp(nflux.varNames,'Slow_h'),theseYears),4),3) .* (land_area_YX_m2*1e-4) ;
tmp = sum(tmp(~isnan(tmp))) * 1e-9 ;
fprintf('%.1f\t%s\n',tmp,'Harvest (TgN/yr)')
tmp = mean(nflux.maps_YXvy(:,:,strcmp(nflux.varNames,'LU_ch'),theseYears),4) .* (land_area_YX_m2*1e-4) ;
tmp = sum(tmp(~isnan(tmp))) * 1e-9 ;
fprintf('%.1f\t%s\n',tmp,'Land-use change (TgN/yr)')
clear tmp

disp(' ')
disp('BVOCs:')

tmp = mean(aiso.maps_YXvy(:,:,strcmp(aiso.varNames,'Total'),theseYears),4) .* land_area_YX_m2 ;
tmp = sum(tmp(~isnan(tmp))) * 1e-15 ;
fprintf('%.1f\t%s\n',tmp,'Isoprene (TgC/yr)')
clear tmp

tmp = mean(amon.maps_YXvy(:,:,strcmp(amon.varNames,'Total'),theseYears),4) .* land_area_YX_m2 ;
tmp = sum(tmp(~isnan(tmp))) * 1e-15 ;
fprintf('%.1f\t%s\n',tmp,'Monoterpene (TgC/yr)')
clear tmp

disp(' ')
get_albedo(fpc,LU,cropfracs,baresoil_albedo_YX,land_area_YX);


