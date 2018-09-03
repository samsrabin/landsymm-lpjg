%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing harmonized PLUM land use trajectory %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisVer = '' ;
% thisVer = 'v2.' ;

runList = {'SSP1.v10.s1' ;
            'SSP3.v10.s1' ;
            'SSP4.v10.s1' ;
            'SSP5.v10.s1'} ;

base_year = 2010 ;

yearList = 2011:2100 ;

norm2extra = 0.177 ;

out_dir = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/harmonization_figs_20180831/' ;


%% Setup

topDir = addslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG') ;

PLUM_base_in = [topDir addslashifneeded(runList{1}) '2010/'] ;

Nyears = length(yearList) ;
Nruns = length(runList) ;

addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;

% % Make lower-left lat/lon map (for compat. with PLUM style)
% lons_map_2deg = repmat(-180:2:178,[90 1]) ;
% lats_map_2deg = repmat((-90:2:88)',[1 180]) ;
% lons_map = repmat(-180:0.5:179.5,[360 1]) ;
% lats_map = repmat((-90:0.5:89.5)',[1 720]) ;


%% Import and aggregate reference LU map

disp('Importing reference LU map...')

% Get LUH2 land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
landArea_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear *_YXqd

% Import LUH2
luh2_file = '/Users/Shared/PLUM/input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
luh2 = lpjgu_matlab_readTable_then2map(luh2_file) ;
if ~isempty(find(luh2.maps_YXvy(:,:,contains(luh2.varNames,{'URBAN','PEATLAND'}),:)>0,1))
    error('This code is not designed to handle LUH2 inputs with any URBAN or PEATLAND area!')
end
luh2.maps_YXvy = luh2.maps_YXvy(:,:,~contains(luh2.varNames,{'URBAN','PEATLAND'}),:) ;
luh2.varNames = luh2.varNames(~contains(luh2.varNames,{'URBAN','PEATLAND'})) ;
yearList_luh2 = luh2.yearList ;
Nyears_luh2 = length(yearList_luh2) ;

% Harmonize LUH2 mask and PLUM mask
file_in = [PLUM_base_in 'LandCoverFract.txt'] ;
S = lpjgu_matlab_readTable_then2map(file_in,'verboseIfNoMat',false,'force_mat_nosave',true) ;
mask_YX = isnan(S.maps_YXv(:,:,1)) ...
    | sum(S.maps_YXv(:,:,contains(S.varNames,{'CROPLAND','PASTURE','NATURAL'})),3)==0 ...
    | landArea_YX==0 ;
landArea_YX(mask_YX) = 0 ;
clear S

% Get repmat 0.5º land area
LUnames = luh2.varNames ;
Nlu = length(LUnames) ;
landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;

% Convert from "fraction of land" to "land area (km2)"
% (This also masks where needed due to harmonization of LUH2+PLUM masks)
luh2.maps_YXvy = luh2.maps_YXvy .* repmat(landArea_YXv,[1 1 1 Nyears_luh2]) ;

% Import LUH2 base_year crop fractions
luh2_cropf = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/remaps_v4/cropfracs.remapv4.20180214.cgFertIrr0.setaside0103.m0.txt',...
    'verboseIfNoMat',false,'force_mat_save',true) ;
% Combine CC3G and CC4G into ExtraCrop
if any(strcmp(luh2_cropf.varNames,'CC3G')) && any(strcmp(luh2_cropf.varNames,'CC4G'))
    luh2_cropf.maps_YXv(:,:,strcmp(luh2_cropf.varNames,'CC3G')) = ...
        luh2_cropf.maps_YXv(:,:,strcmp(luh2_cropf.varNames,'CC3G')) ...
      + luh2_cropf.maps_YXv(:,:,strcmp(luh2_cropf.varNames,'CC4G')) ;
    luh2_cropf.maps_YXv(:,:,strcmp(luh2_cropf.varNames,'CC4G')) = [] ;
    luh2_cropf.varNames(strcmp(luh2_cropf.varNames,'CC4G')) = [] ;
    luh2_cropf.varNames{strcmp(luh2_cropf.varNames,'CC3G')} = 'ExtraCrop' ;
end
% Get previous crop types
tmp = luh2_cropf.varNames ;
for c = 1:length(luh2_cropf.varNames)
    thisCrop = luh2_cropf.varNames{c} ;
    thisCropI = [thisCrop 'i'] ;
    tmp(strcmp(tmp,thisCropI)) = [] ;
end
LPJGcrops = tmp ;
Ncrops_lpjg = length(LPJGcrops) ;
Ncrops = Ncrops_lpjg ;
Nagri = Ncrops_lpjg + 1 ;
clear tmp
% Combine irrigated and rainfed
any_irrigated = ~isequal(luh2_cropf.varNames,LPJGcrops) ;
if any_irrigated
    tmp.varNames = LPJGcrops ;
    tmp.maps_YXv = nan(size(luh2_cropf.maps_YXv,1), ...
                        size(luh2_cropf.maps_YXv,2), ...
                        Ncrops_lpjg) ;
    for c = 1:Ncrops_lpjg
        thisCrop = LPJGcrops{c} ;
        thisCropI = [thisCrop 'i'] ;
        is_thisCrop  = strcmp(luh2_cropf.varNames,thisCrop) ;
        is_thisCropI = strcmp(luh2_cropf.varNames,thisCropI) ;
        if any(is_thisCropI)
            tmp.maps_YXv(:,:,c) = luh2_cropf.maps_YXv(:,:,is_thisCrop) ...
                                + luh2_cropf.maps_YXv(:,:,is_thisCropI) ;
        else
            tmp.maps_YXv(:,:,c) = luh2_cropf.maps_YXv(:,:,is_thisCrop) ;
        end
    end
    luh2_cropf = tmp ;
    clear tmp
end

% Convert luh2 to have individual crops instead of CROPLAND
if ~exist('luh2_orig','var')
    luh2_orig = luh2 ;
end
luh2_cropa = luh2_cropf ;
luh2_cropa.maps_YXvy = repmat(luh2_cropf.maps_YXv,[1 1 1 Nyears_luh2]) .* luh2_orig.maps_YXvy(:,:,strcmp(luh2_orig.varNames,'CROPLAND'),:) ;
clear luh2
luh2.varNames = [luh2_cropa.varNames luh2_orig.varNames(~strcmp(luh2_orig.varNames,'CROPLAND'))] ;
luh2.maps_YXvy = cat(3,luh2_cropa.maps_YXvy,luh2_orig.maps_YXvy(:,:,~strcmp(luh2_orig.varNames,'CROPLAND'),:)) ;

% Get info
LUnames = luh2.varNames ;
Nlu = length(LUnames) ;
notBare = ~strcmp(LUnames,'BARREN') ;
isAgri = ~strcmp(LUnames,'NATURAL') & notBare ;
LUnames_agri = LUnames(isAgri) ;
isAgri_isPast = strcmp(LUnames_agri,'PASTURE') ;
isCrop = ~strcmp(LUnames,'NATURAL') & ~strcmp(LUnames,'PASTURE') & notBare ;

% Get time series
luh2.maps_YXvy(isnan(luh2.maps_YXvy)) = 0 ;
ts_luh2_cy = squeeze(sum(sum(luh2.maps_YXvy(:,:,isCrop,:),2),1)) ;

% Get LUH2 fraction that is vegetated, barren
luh2_vegd_YX = sum(luh2.maps_YXvy(:,:,notBare,end),3) ;
luh2_vegdFrac_YX = luh2_vegd_YX ./ landArea_YX ;
luh2_bare_YX = luh2.maps_YXvy(:,:,strcmp(LUnames,'BARREN'),end) ;
luh2_bareFrac_YX = luh2_bare_YX ./ landArea_YX ;
             
% % Get lats/lons
% list2map_2deg = find(landArea_2deg_YX>0) ;
% lons_2deg = lons_map_2deg(list2map_2deg) ;
% lats_2deg = lats_map_2deg(list2map_2deg) ;
% list2map = find(landArea_YX>0) ;
% lons = lons_map(list2map) ;
% lats = lats_map(list2map) ;

% Get PLUM crop types
file_in = [PLUM_base_in 'LandUse.txt'] ;
T = lpjgu_matlab_readTable(file_in,'verboseIfNoMat',false,'dont_save_MAT',true) ;
PLUMcrops = T.Properties.VariableNames ;
PLUMcrops = PLUMcrops(...
    contains(PLUMcrops,'_A') ...
    & ~contains(PLUMcrops,'ruminants') ...
    & ~contains(PLUMcrops,'monogastrics') ...
    & ~contains(PLUMcrops,'pasture') ...
    ) ;
PLUMcrops = strrep(PLUMcrops,'_A','') ;
Ncfts_plum = length(PLUMcrops) ;
clear T

% Get PLUM-to-LPJG mapping
PLUMtoLPJG = {} ;
PLUMtoLPJG{strcmp(LPJGcrops,'CerealsC3')} = 'wheat' ;
PLUMtoLPJG{strcmp(LPJGcrops,'CerealsC4')} = 'maize' ;
PLUMtoLPJG{strcmp(LPJGcrops,'Rice')} = 'rice' ;
PLUMtoLPJG{strcmp(LPJGcrops,'Oilcrops')} = 'oilcrops' ;
PLUMtoLPJG{strcmp(LPJGcrops,'Pulses')} = 'pulses' ;
PLUMtoLPJG{strcmp(LPJGcrops,'StarchyRoots')} = 'starchyRoots' ;
PLUMtoLPJG{strcmp(LPJGcrops,'Miscanthus')} = 'energycrops' ;
PLUMtoLPJG{strcmp(LPJGcrops,'ExtraCrop')} = 'setaside' ;

% Check that every PLUM crop is mapped
checkPLUMmap(PLUMtoLPJG,PLUMcrops) ;

disp('Done importing reference LU map and doing other setup.')



%% Import PLUM (original + harmonized)

is2deg = false ;
bareFrac_y0_YX = luh2_bareFrac_YX ;

disp('Setting up PLUM*_YXvyr arrays...')
PLUMorig_YXvyr = nan(360,720,Nlu,Nyears,Nruns,'single') ;
PLUMharm_YXvyr = nan(360,720,Nlu,Nyears,Nruns,'single') ;

for r = 1:Nruns
    thisRun = removeslashifneeded(runList{r}) ;
    
    % Original
    fprintf('Importing %s...\n', thisRun) ;
    S_out = PLUMharm_pp_readPLUM([topDir thisRun], base_year, yearList, ...
        landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
        is2deg, bareFrac_y0_YX, norm2extra) ;
    PLUMorig_YXvyr(:,:,:,:,r) = S_out.maps_YXvy ;
    clear S_out
    
    % Harmonized
    fprintf('Importing %s.harm...\n', thisRun) ;
    S_out = PLUMharm_pp_readPLUM([topDir thisRun '.harm'],base_year,yearList, ...
        landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
        is2deg, [], 0) ;
    PLUMharm_YXvyr(:,:,:,:,r) = S_out.maps_YXvy ;
    clear S_out
    
end

disp('Done reading PLUM.')


%% Time series of LUs

ts_luh2_cy = squeeze(nansum(nansum(luh2.maps_YXvy,1),2)) ;
ts_luh2_cy = cat(1,sum(ts_luh2_cy(isCrop,:),1),ts_luh2_cy(~isCrop,:)) ;
ts_orig_cyr = squeeze(nansum(nansum(PLUMorig_YXvyr,1),2)) ;
ts_orig_cyr = cat(1,sum(ts_orig_cyr(isCrop,:,:),1),ts_orig_cyr(~isCrop,:,:)) ;
ts_harm_cyr = squeeze(nansum(nansum(PLUMharm_YXvyr,1),2)) ;
ts_harm_cyr = cat(1,sum(ts_harm_cyr(isCrop,:,:),1),ts_harm_cyr(~isCrop,:,:)) ;
combinedLUs = [{'CROPLAND'} LUnames(~isCrop)] ;

spacing = [0.05 0.1] ;

figure('Color','w','Position',figurePos)

incl_luh2_years = yearList_luh2>=1950 & yearList_luh2<=2010 ;

for v = 1:length(combinedLUs)
    subplot_tight(2,2,v,spacing) ;
    plot(yearList_luh2(incl_luh2_years),ts_luh2_cy(v,incl_luh2_years)*1e-6,'-k','LineWidth',2) ;
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList,squeeze(ts_harm_cyr(v,:,:))*1e-6,'-','LineWidth',1)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList,squeeze(ts_orig_cyr(v,:,:))*1e-6,'--','LineWidth',1)
    hold off
    title(combinedLUs{v})
    set(gca,'FontSize',14)
    ylabel('Million km2')
    legend({'LUH2','SSP1','SSP3','SSP4','SSP5'},'Location','NorthWest')
end

% Save
export_fig([out_dir 'timeSeries_landUse.pdf']) ;


%% Time series of crops

ts_luh2_cy = squeeze(nansum(nansum(luh2.maps_YXvy,1),2)) ;
ts_orig_cyr = squeeze(nansum(nansum(PLUMorig_YXvyr,1),2)) ;
ts_harm_cyr = squeeze(nansum(nansum(PLUMharm_YXvyr,1),2)) ;

spacing = [0.05 0.1] ;

figure('Color','w','Position',figurePos)

incl_luh2_years = yearList_luh2>=1950 & yearList_luh2<=2010 ;

for v = 1:Ncrops
    subplot_tight(4,2,v,spacing) ;
    plot(yearList_luh2(incl_luh2_years),ts_luh2_cy(v,incl_luh2_years)*1e-6,'-k','LineWidth',2) ;
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList,squeeze(ts_harm_cyr(v,:,:))*1e-6,'-','LineWidth',1)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList,squeeze(ts_orig_cyr(v,:,:))*1e-6,'--','LineWidth',1)
    hold off
    title(LPJGcrops{v})
    set(gca,'FontSize',14)
    ylabel('Million km2')
    legend({'LUH2','SSP1','SSP3','SSP4','SSP5'},'Location','NorthWest')
end

% Save
export_fig([out_dir 'timeSeries_crops.pdf']) ;


%% Maps: At three years

theseYears = [2011 2050 2100] ;
spacing = [0.01 0.025] ;
cbar_loc = 'SouthOutside' ;
y1 = 60 ;
fontSize = 14 ;
png_res = 150 ;

for r = 1:Nruns
    thisRun = runList{r} ;
    for v = 1:Nlu
        figure('Color','w','Position',figurePos) ;
        thisLU = LUnames{v} ;
        for y = 1:3
            thisYear = theseYears(y) ;
            h1 = subplot_tight(2,3,y,spacing) ;
            tmp = PLUMorig_YXvyr(:,:,v,yearList==thisYear,r) ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s orig: %s, %d',thisRun,thisLU,thisYear)) ;
            set(gca,'FontSize',fontSize)
            h2 = subplot_tight(2,3,y+3,spacing) ;
            tmp = PLUMharm_YXvyr(:,:,v,yearList==thisYear,r) ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s harm: %s, %d (km^2)',thisRun,thisLU,thisYear)) ;
            new_caxis = [0 max([caxis(h1) caxis(h2)])] ;
            caxis(h1,new_caxis) ;
            caxis(h2,new_caxis) ;
            set(gca,'FontSize',fontSize)
        end
        export_fig([out_dir 'maps_' thisLU '_' strrep(num2str(theseYears),'  ','-') '_' thisRun '.png'],['-r' num2str(png_res)]) ;
        close
    end
end


%% Maps: Diffs between orig and harm at 3 years

theseYears = [2011 2050 2100] ;
spacing = [0.01 0.025] ;
cbar_loc = 'SouthOutside' ;
y1 = 60 ;
fontSize = 14 ;
png_res = 150 ;
thisPos = [1         500        1440         305] ;

for r = 1:Nruns
    thisRun = runList{r} ;
    for v = 1:Nlu
        figure('Color','w','Position',thisPos) ;
        thisLU = LUnames{v} ;
        for y = 1:3
            thisYear = theseYears(y) ;
            h1 = subplot_tight(1,3,y,spacing) ;
            tmp1 = PLUMorig_YXvyr(:,:,v,yearList==thisYear,r) ;
            tmp2 = PLUMharm_YXvyr(:,:,v,yearList==thisYear,r) ;
            tmp = tmp2 - tmp1 ;
%             tmp = tmp2/sum(tmp2(:)) - tmp1/sum(tmp1(:)) ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colormap(flipud(brewermap(64,'rdbu_ssr'))) ;
            caxis([-1 1]*max(abs(caxis))) ;
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s orig: %s, %d',thisRun,thisLU,thisYear)) ;
            set(gca,'FontSize',fontSize)
        end
        export_fig([out_dir 'mapsOHdiffs_' thisLU '_' strrep(num2str(theseYears),'  ','-') '_' thisRun '.png'],['-r' num2str(png_res)]) ;
        close
    end
end




%% Maps: Differences between two pairs of years

theseYears = [2011 2050 2100] ;
spacing = [0.01 0.025] ;
cbar_loc = 'SouthOutside' ;
y1 = 60 ;
fontSize = 14 ;
png_res = 150 ;

for r = 1:Nruns
    thisRun = runList{r} ;
    for v = 1:Nlu
        figure('Color','w','Position',figurePos) ;
        thisLU = LUnames{v} ;
        for y = 1:2
            thisYear1 = theseYears(y) ;
            thisYear2 = theseYears(y+1) ;
            h1 = subplot_tight(2,2,y,spacing) ;
            tmp1 = PLUMorig_YXvyr(:,:,v,yearList==thisYear1,r) ;
            tmp2 = PLUMorig_YXvyr(:,:,v,yearList==thisYear2,r) ;
            tmp = tmp2 - tmp1 ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colormap(brewermap(64,'rdbu_ssr')) ;
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s orig: %s, %d-%d',thisRun,thisLU,thisYear1,thisYear2)) ;
            set(gca,'FontSize',fontSize)
            h2 = subplot_tight(2,2,y+2,spacing) ;
            tmp1 = PLUMharm_YXvyr(:,:,v,yearList==thisYear1,r) ;
            tmp2 = PLUMharm_YXvyr(:,:,v,yearList==thisYear2,r) ;
            tmp = tmp2 - tmp1 ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colormap(brewermap(64,'rdbu_ssr')) ;
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s harm: %s, %d-%d',thisRun,thisLU,thisYear1,thisYear2)) ;
            new_caxis = [-1 1] * max([caxis(h1) caxis(h2)]) ;
            caxis(h1,new_caxis) ;
            caxis(h2,new_caxis) ;
            set(gca,'FontSize',fontSize)
        end
        export_fig([out_dir 'mapsChgs_' thisLU '_' strrep(num2str(theseYears),'  ','-') '_' thisRun '.png'],['-r' num2str(png_res)]) ;
        close
    end
end










