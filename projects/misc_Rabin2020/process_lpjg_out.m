%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Work with LPJ-GUESS outputs for impacts paper %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% topDir = '/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/' ; topDir = addslashifneeded(topDir) ;
% inDir = [topDir 'PLUM2LPJG_26_s1_2011-2100/output-2017-11-21-122520'] ;
% inDir = [topDir 'PLUM2LPJG_45_s1_2011-2100/output-2017-11-22-052523'] ;
% inDir = [topDir 'PLUM2LPJG_60_s1_2011-2100/output-2017-11-22-120456'] ;
% inDir = [topDir 'PLUM2LPJG_85_s1_2011-2100/output-2017-11-22-114918_few'] ;

topDir = '/Users/Shared/PLUM/trunk_runs' ; topDir = addslashifneeded(topDir) ;
% inDir = [topDir 'PLUM2LPJG_26_s1_1850-2010/output-2017-11-24-175840'] ;
inDir = [topDir 'PLUM2LPJG_26_s1_2011-2100/output-2017-11-25-102553'] ;
% inDir = [topDir 'PLUM2LPJG_45_s1_2011-2100/output-2017-11-25-110709'] ;
% inDir = [topDir 'PLUM2LPJG_60_s1_2011-2100/output-2017-11-25-104811'] ;
% inDir = [topDir 'PLUM2LPJG_85_s1_2011-2100/output-2017-11-25-112417'] ;


%% Setup

inDir = addslashifneeded(inDir) ;

if strcmp(inDir,[topDir 'PLUM2LPJG_26_s1_1850-2010/output-2017-11-24-175840/'])
    LUfile = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
    cropfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/LU/crop_fractions_fromMIRCA_PLUM_1901_2014_20160906.MP.noIrr.plus_Nfert0-200-1000xIrr_factorial_empties.inclZeroFert.out' ;
elseif strcmp(inDir,'/Users/Shared/PLUM/trunk_runs/PLUM2LPJG_26_s1_2011-2100/output-2017-11-21-122520/') ...
    || strcmp(inDir,[topDir 'PLUM2LPJG_26_s1_2011-2100/output-2017-11-25-102553/'])
    LUfile = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/26_s1.forLPJG.MATLAB/landcover.txt' ;
    cropfile = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/26_s1.forLPJG.MATLAB/cropfractions.txt' ;
elseif strcmp(inDir,'/Users/Shared/PLUM/trunk_runs/PLUM2LPJG_45_s1_2011-2100/output-2017-11-22-052523/') ...
    || strcmp(inDir,[topDir 'PLUM2LPJG_45_s1_2011-2100/output-2017-11-25-110709'])
    LUfile = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/45.forLPJG.MATLAB/landcover.txt' ;
    cropfile = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/45.forLPJG.MATLAB/cropfractions.txt' ;
elseif strcmp(inDir,'/Users/Shared/PLUM/trunk_runs/PLUM2LPJG_60_s1_2011-2100/output-2017-11-22-120456/') ...
    || strcmp(inDir,[topDir 'PLUM2LPJG_60_s1_2011-2100/output-2017-11-25-104811'])
    LUfile = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/60.forLPJG.MATLAB/landcover.txt' ;
    cropfile = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/60.forLPJG.MATLAB/cropfractions.txt' ;
elseif strcmp(inDir,'/Users/Shared/PLUM/trunk_runs/PLUM2LPJG_85_s1_2011-2100/output-2017-11-22-114918_few/') ...
    || strcmp(inDir,[topDir 'PLUM2LPJG_85_s1_2011-2100/output-2017-11-25-112417/'])
    LUfile = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/85.forLPJG.MATLAB/landcover.txt' ;
    cropfile = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/85.forLPJG.MATLAB/cropfractions.txt' ;
else
    error('Specify LUfile and cropfile for this inDir!')
end

addpath(genpath(landsymm_lpjg_path()))


%% Import

% Carbon: Vegetation, soil/litter, product, total (kgC/m2)
cpool = lpjgu_matlab_readTable_then2map([inDir 'cpool.out'],'force_mat_save',true) ;
% cpool_table = lpjgu_matlab_readTable([inDir 'cpool.out'],'do_save_mat',true) ;

% Evapotranspiration (mm/yr)
aaet = lpjgu_matlab_readTable_then2map([inDir 'aaet.out'],'force_mat_save',true) ;

% Runoff (mm/yr)
tot_runoff = lpjgu_matlab_readTable_then2map([inDir 'tot_runoff.out'],'force_mat_save',true) ;
% Runoff (mm/month)
mon_runoff = lpjgu_matlab_readTable_then2map([inDir 'mrunoff.out'],'force_mat_save',true) ;

% Monthly snow depth (m)
snowdepth = lpjgu_matlab_readTable_then2map([inDir 'msnowdepth.out'],'force_mat_save',true) ;

% BVOCs (isoprene, monoterpenes) (mgC/m2)
aiso = lpjgu_matlab_readTable_then2map([inDir 'aiso.out'],'force_mat_save',true) ;
amon = lpjgu_matlab_readTable_then2map([inDir 'amon.out'],'force_mat_save',true) ;

% N flux (kgN/ha)
nflux = lpjgu_matlab_readTable_then2map([inDir 'nflux.out'],'force_mat_save',true) ;

% Foliar projective cover
fpc = lpjgu_matlab_readTable_then2map([inDir 'fpc.out'],'force_mat_save',true) ;
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
LU = lpjgu_matlab_readTable_then2map(LUfile,'force_mat_save',true) ;
if ~isequal(LU.yearList,fpc.yearList)
    if min(LU.yearList) <= min(fpc.yearList) && max(LU.yearList) >= max(fpc.yearList)
        LU.maps_YXvy = LU.maps_YXvy(:,:,:,LU.yearList>= min(fpc.yearList) & LU.yearList>= max(fpc.yearList)) ;
        LU.yearList = transpose(LU.yearList(1):max(fpc.yearList)) ;
    end
    if ~isequal(squeeze(LU.yearList),squeeze(fpc.yearList))
        error('Rework this so yearLists match!')
    end
end
cropfracs = lpjgu_matlab_readTable_then2map(cropfile,'force_mat_save',true) ;
if ~isequal(cropfracs.yearList,fpc.yearList)
    if min(cropfracs.yearList) <= min(fpc.yearList) && max(cropfracs.yearList) >= max(fpc.yearList)
        cropfracs.maps_YXvy = cropfracs.maps_YXvy(:,:,:,cropfracs.yearList>=min(fpc.yearList) & cropfracs.yearList>= max(fpc.yearList)) ;
        cropfracs.yearList = transpose(cropfracs.yearList(1):max(fpc.yearList)) ;
    elseif min(cropfracs.yearList) > min(fpc.yearList) && max(cropfracs.yearList) == max(fpc.yearList) && isequaln(cropfracs.maps_YXvy(:,:,:,1),cropfracs.maps_YXvy(:,:,:,end))
        Nmissing = length(fpc.yearList(1):cropfracs.yearList(1)-1) ;
        cropfracs.maps_YXvy = cat(4,repmat(cropfracs.maps_YXvy(:,:,:,1),[1 1 1 Nmissing]),cropfracs.maps_YXvy(:,:,:,cropfracs.yearList>=min(fpc.yearList) & cropfracs.yearList>= max(fpc.yearList))) ;
        cropfracs.yearList = fpc.yearList ;
    end
    if ~isequal(cropfracs.yearList,fpc.yearList)
        error('Rework this so yearLists match!')
    end
end

% Import bare-soil albedo
baresoil_albedo_file = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work/soilmap.txt' ;
baresoil_albedo = lpjgu_matlab_readTable_then2map(baresoil_albedo_file,'force_mat_save',true) ;
baresoil_albedo_YX = baresoil_albedo.maps_YXv ;
clear baresoil_albedo

% Misc
yearList = cpool.yearList ;
Nyears = length(yearList) ;

disp('Done.')


%% Get time series

disp('Getting time series...')

cf_cpool = 1e-3*1e-9 ;
cpool_ts.VegC = getTS(cpool,'VegC',land_area_YX_m2) * cf_cpool ;
cpool_ts.LitterC_SoilC = getTS(cpool,{'LitterC','SoilC'},land_area_YX_m2) * cf_cpool ;
cpool_ts.HarvSlowC = getTS(cpool,'HarvSlowC',land_area_YX_m2) * cf_cpool ;
cpool_ts.Total = getTS(cpool,'Total',land_area_YX_m2) * cf_cpool ;

cf_water = 1e-3*1e-3*1e-3 ;
aaet_ts = getTS(aaet,'Total',land_area_YX) * cf_water ;
tot_runoff_ts = getTS(tot_runoff,'Total',land_area_YX) * cf_water ;

cf_nflux = 1e-9 ;
nflux_ts.flux = getTS(nflux,'flux',land_area_YX_m2*1e-4) * cf_nflux ;
nflux_ts.leach = getTS(nflux,'leach',land_area_YX_m2*1e-4) * cf_nflux ;
nflux_ts.harvest = getTS(nflux,'harvest',land_area_YX_m2*1e-4) * cf_nflux ;
nflux_ts.LU_ch = getTS(nflux,'LU_ch',land_area_YX_m2*1e-4) * cf_nflux ;

cf_bvoc = 1e-15 ;
aiso_ts = getTS(aiso,'Total',land_area_YX_m2) * cf_bvoc ;
amon_ts = getTS(amon,'Total',land_area_YX_m2) * cf_bvoc ;

save([inDir 'timeseries.mat'],'cpool_ts','aaet_ts','tot_runoff_ts',...
                              'nflux_ts','aiso_ts','amon_ts') ;

disp('Done.')


%% Get last-decade maps

disp('Getting last-decade maps...')

cpool_ld = cpool ;
cpool_ld.maps_YXvy = cpool_ld.maps_YXvy(end-9:end) ;
cpool_ld.yearList = cpool_ld.yearList(end-9:end) ;
save([inDir 'last_decade.mat'],'cpool_ld') ;
clear *_ld

aaet_ld = aaet ;
aaet_ld.maps_YXvy = aaet_ld.maps_YXvy(end-9:end) ;
aaet_ld.yearList = aaet_ld.yearList(end-9:end) ;
save([inDir 'last_decade.mat'],'aaet_ld','-append') ;
clear *_ld

tot_runoff_ld = tot_runoff ;
tot_runoff_ld.maps_YXvy = tot_runoff_ld.maps_YXvy(end-9:end) ;
tot_runoff_ld.yearList = tot_runoff_ld.yearList(end-9:end) ;
save([inDir 'last_decade.mat'],'tot_runoff_ld','-append') ;
clear *_ld

mon_runoff_ld = mon_runoff ;
mon_runoff_ld.maps_YXvy = mon_runoff_ld.maps_YXvy(end-9:end) ;
mon_runoff_ld.yearList = mon_runoff_ld.yearList(end-9:end) ;
save([inDir 'last_decade.mat'],'mon_runoff_ld','-append') ;
clear *_ld

snowdepth_ld = snowdepth ;
snowdepth_ld.maps_YXvy = snowdepth_ld.maps_YXvy(end-9:end) ;
snowdepth_ld.yearList = snowdepth_ld.yearList(end-9:end) ;
save([inDir 'last_decade.mat'],'snowdepth_ld','-append') ;
clear *_ld

aiso_ld = aiso ;
aiso_ld.maps_YXvy = aiso_ld.maps_YXvy(end-9:end) ;
aiso_ld.yearList = aiso_ld.yearList(end-9:end) ;
save([inDir 'last_decade.mat'],'aiso_ld','-append') ;
clear *_ld

amon_ld = amon ;
amon_ld.maps_YXvy = amon_ld.maps_YXvy(end-9:end) ;
amon_ld.yearList = amon_ld.yearList(end-9:end) ;
save([inDir 'last_decade.mat'],'amon_ld','-append') ;
clear *_ld

nflux_ld = nflux ;
nflux_ld.maps_YXvy = nflux_ld.maps_YXvy(end-9:end) ;
nflux_ld.yearList = nflux_ld.yearList(end-9:end) ;
save([inDir 'last_decade.mat'],'nflux_ld','-append') ;
clear *_ld

fpc_ld = fpc ;
fpc_ld.maps_YXvy = fpc_ld.maps_YXvy(end-9:end) ;
fpc_ld.yearList = fpc_ld.yearList(end-9:end) ;
save([inDir 'last_decade.mat'],'fpc_ld','-append') ;
clear *_ld

LU_ld = LU ;
LU_ld.maps_YXvy = LU_ld.maps_YXvy(end-9:end) ;
LU_ld.yearList = LU_ld.yearList(end-9:end) ;
save([inDir 'last_decade.mat'],'LU_ld','-append') ;
clear *_ld

disp('Done.')


%% Get one decade's values

y1 = max(cpool.yearList) - 9 ;
yN = max(cpool.yearList) ;
disp([num2str(y1) '-' num2str(yN)])
disp(' ')
theseYears = (cpool.yearList>=y1 & cpool.yearList<=yN) ;

disp('Carbon:')
fprintf('%.1f\t%s\n',mean(cpool_ts.VegC(theseYears)),'Vegetation C (Gt)')
fprintf('%.1f\t%s\n',mean(cpool_ts.LitterC_SoilC(theseYears)),'Soil+Litter C (Gt)')
fprintf('%.1f\t%s\n',mean(cpool_ts.HarvSlowC(theseYears)),'Product C (Gt)')
fprintf('%.1f\t%s\n',mean(cpool_ts.Total(theseYears)),'Total C (Gt)')

disp(' ')
disp('Water:')
fprintf('%.1f\t%s\n',mean(aaet_ts(theseYears)),'Annual evapotranspiration (1000 km3/yr)')
fprintf('%.1f\t%s\n',mean(tot_runoff_ts(theseYears)),'Annual runoff (1000 km3/yr)')
tmp = max(mean(mon_runoff.maps_YXvy(:,:,:,theseYears),4),[],3) .* land_area_YX ;
tmp = sum(tmp(~isnan(tmp))) * cf_water ;
fprintf('%.1f\t%s\n',tmp,'Peak monthly runoff (1000 km3/month)')
clear tmp

disp(' ')
disp('N loss:')

% kgN/ha/yr --> TgN/yr
fprintf('%.1f\t%s\n',mean(nflux_ts.flux(theseYears)),'Gaseous (TgN/yr)')
fprintf('%.1f\t%s\n',mean(nflux_ts.leach(theseYears)),'Dissolved (TgN/yr)')
fprintf('%.1f\t%s\n',mean(nflux_ts.harvest(theseYears)),'Harvest (TgN/yr)')
fprintf('%.1f\t%s\n',mean(nflux_ts.LU_ch(theseYears)),'Land-use change (TgN/yr)')

disp(' ')
disp('BVOCs:')
fprintf('%.1f\t%s\n',mean(aiso_ts(theseYears)),'Isoprene (TgC/yr)')
fprintf('%.1f\t%s\n',mean(amon_ts(theseYears)),'Monoterpene (TgC/yr)')

disp(' ')
get_albedo(fpc,LU,cropfracs,baresoil_albedo_YX,land_area_YX);


%% Plot time series

%%% Options %%%%%%%%%%%
nsubplot = [4 3] ;
spacing = [0.075 0.05] ;   % [v h]
fontsize = 14 ;
%%%%%%%%%%%%%%%%%%%%%%%

% Setup
ny = nsubplot(1) ;
nx = nsubplot(2) ;

figure('Color','w','Position',figurePos) ;

% Total C
subplot_tight(ny,nx,1,spacing) ;
plot(yearList, cpool_ts.Total)
title('Total C')
ylabel('GtC')
set(gca,'FontSize',fontsize)

% Vegetation C
subplot_tight(ny,nx,2,spacing) ;
plot(yearList, cpool_ts.VegC)
title('Vegetation C')
ylabel('GtC')
set(gca,'FontSize',fontsize)

% Soil+Litter C
subplot_tight(ny,nx,3,spacing) ;
plot(yearList, cpool_ts.LitterC_SoilC)
title('Soil+Litter C')
ylabel('GtC')
set(gca,'FontSize',fontsize)

% January albedo

% July albedo

% Evapotranspiration
subplot_tight(ny,nx,6,spacing) ;
plot(yearList, aaet_ts)
% tmp = movmean(aaet_ts,5) ;
% plot(yearList(3:end-2),tmp(3:end-2))
title('Evapotranspiration')
ylabel('1000 km^3')
set(gca,'FontSize',fontsize)

% Runoff
subplot_tight(ny,nx,7,spacing) ;
plot(yearList, tot_runoff_ts)
title('Runoff')
ylabel('1000 km^3')
set(gca,'FontSize',fontsize)

% Peak monthly runoff

% Crop production (Ecal/yr)

% N loss
subplot_tight(ny,nx,10,spacing) ;
plot(yearList, nflux_ts.flux+nflux_ts.leach+nflux_ts.harvest+nflux_ts.LU_ch)
title('N loss')
ylabel('TgN')
set(gca,'FontSize',fontsize)

% Isoprene emissions
subplot_tight(ny,nx,11,spacing) ;
plot(yearList, aiso_ts)
title('Isoprene emissions')
ylabel('TgC')
set(gca,'FontSize',fontsize)

% Monoterpene emissions
subplot_tight(ny,nx,12,spacing) ;
plot(yearList, amon_ts)
title('Monoterpene emissions')
ylabel('TgC')
set(gca,'FontSize',fontsize)




