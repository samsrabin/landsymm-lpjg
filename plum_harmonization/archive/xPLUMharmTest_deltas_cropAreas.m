%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Are we getting the deltas we expect? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisVer = '' ;
% thisVer = '2deg.' ;
% thisVer = 'orig.' ;

base_year = 2010 ;
year0 = 2010 ;
% lastYear = 2015 ;

do_debug = false ;

PLUM_in_top = removeslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1') ;

outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;

conserv_tol_pct = 0.1 ;
norm2extra = 0.177 ;


%% Setup

PLUM_out_top = addslashifneeded([PLUM_in_top '.harm/']) ;
unix(['mkdir -p ' PLUM_out_top]) ;
PLUM_base_in = addslashifneeded([PLUM_in_top '/' num2str(base_year)]) ;
PLUM_in_top = addslashifneeded(PLUM_in_top) ;

addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;

% Make lower-left lat/lon map (for compat. with PLUM style)
lons_map_2deg = repmat(-180:2:178,[90 1]) ;
lats_map_2deg = repmat((-90:2:88)',[1 180]) ;
lons_map = repmat(-180:0.5:179.5,[360 1]) ;
lats_map = repmat((-90:0.5:89.5)',[1 720]) ;

year1 = year0 + 1 ;
year2 = year1 + 1 ;
year0s = num2str(year0) ;
year1s = num2str(year1) ;
year2s = num2str(year2) ;


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

% Import LUH2 base_year
luh2_file = '/Users/Shared/PLUM/input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
luh2 = lpjgu_matlab_readTable_then2map(luh2_file) ;
if ~isempty(find(luh2.maps_YXvy(:,:,contains(luh2.varNames,{'URBAN','PEATLAND'}),:)>0,1))
    error('This code is not designed to handle LUH2 inputs with any URBAN or PEATLAND area!')
end
luh2_base.maps_YXv = luh2.maps_YXvy(:,:,~contains(luh2.varNames,{'URBAN','PEATLAND'}),luh2.yearList==base_year) ;
luh2_base.varNames = luh2.varNames(~contains(luh2.varNames,{'URBAN','PEATLAND'})) ;
luh2_base.list_to_map = luh2.list_to_map ;
clear luh2

% Harmonize LUH2 mask and PLUM mask
file_in = [PLUM_base_in 'LandCoverFract.txt'] ;
S = lpjgu_matlab_readTable_then2map(file_in,'verboseIfNoMat',false,'force_mat_nosave',true) ;
mask_YX = isnan(S.maps_YXv(:,:,1)) ...
    | sum(S.maps_YXv(:,:,contains(S.varNames,{'CROPLAND','PASTURE','NATURAL'})),3)==0 ...
    | landArea_YX==0 ;
landArea_YX(mask_YX) = 0 ;
clear S

% Get repmat 0.5º land area
LUnames = luh2_base.varNames ;
Nlu = length(LUnames) ;
landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;

% Convert from "fraction of land" to "land area (km2)"
% (This also masks where needed due to harmonization of LUH2+PLUM masks)
luh2_base.maps_YXv = luh2_base.maps_YXv .* landArea_YXv ;

% Import LUH2 base_year crop fractions
% remaps_v4 has setAside and unhandledCrops in CC3G and CC4G. In the next
% step, we will combine CC3G and CC4G into ExtraCrop. (Was called SetAside,
% but renamed to avoid confusion, considering that this is not just
% setAside but also unhandledCrops.)
luh2_base_cropf = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/remaps_v4/cropfracs.remapv4.20180214.cgFertIrr0.setaside0103.m0.txt',...
    'verboseIfNoMat',false,'force_mat_save',true) ;
% Get just base year, if needed
if isfield(luh2_base_cropf,'yearList')
    tmp = luh2_base_cropf.maps_YXvy(:,:,luh2_base_cropf.yearList==base_year) ;
    luh2_base_cropf = rmfield(luh2_base_cropf,{'maps_YXvy','yearList'}) ;
    luh2_base_cropf.maps_YXv = tmp ;
    clear tmp
end
% Combine CC3G and CC4G into ExtraCrop
if any(strcmp(luh2_base_cropf.varNames,'CC3G')) && any(strcmp(luh2_base_cropf.varNames,'CC4G'))
    luh2_base_cropf.maps_YXv(:,:,strcmp(luh2_base_cropf.varNames,'CC3G')) = ...
        luh2_base_cropf.maps_YXv(:,:,strcmp(luh2_base_cropf.varNames,'CC3G')) ...
      + luh2_base_cropf.maps_YXv(:,:,strcmp(luh2_base_cropf.varNames,'CC4G')) ;
    luh2_base_cropf.maps_YXv(:,:,strcmp(luh2_base_cropf.varNames,'CC4G')) = [] ;
    luh2_base_cropf.varNames(strcmp(luh2_base_cropf.varNames,'CC4G')) = [] ;
    luh2_base_cropf.varNames{strcmp(luh2_base_cropf.varNames,'CC3G')} = 'ExtraCrop' ;
end
% Get previous crop types
tmp = luh2_base_cropf.varNames ;
for c = 1:length(luh2_base_cropf.varNames)
    thisCrop = luh2_base_cropf.varNames{c} ;
    thisCropI = [thisCrop 'i'] ;
    tmp(strcmp(tmp,thisCropI)) = [] ;
end
LPJGcrops = tmp ;
Ncrops_lpjg = length(LPJGcrops) ;
Nagri = Ncrops_lpjg + 1 ;
clear tmp
% Combine irrigated and rainfed
any_irrigated = ~isequal(luh2_base_cropf.varNames,LPJGcrops) ;
if any_irrigated
    tmp.varNames = LPJGcrops ;
    tmp.maps_YXv = nan(size(luh2_base_cropf.maps_YXv,1), ...
                       size(luh2_base_cropf.maps_YXv,2), ...
                       Ncrops_lpjg) ;
    for c = 1:Ncrops_lpjg
        thisCrop = LPJGcrops{c} ;
        thisCropI = [thisCrop 'i'] ;
        is_thisCrop  = strcmp(luh2_base_cropf.varNames,thisCrop) ;
        is_thisCropI = strcmp(luh2_base_cropf.varNames,thisCropI) ;
        if any(is_thisCropI)
            tmp.maps_YXv(:,:,c) = luh2_base_cropf.maps_YXv(:,:,is_thisCrop) ...
                                + luh2_base_cropf.maps_YXv(:,:,is_thisCropI) ;
        else
            tmp.maps_YXv(:,:,c) = luh2_base_cropf.maps_YXv(:,:,is_thisCrop) ;
        end
    end
    luh2_base_cropf = tmp ;
    clear tmp
end

% Convert luh2_base to have individual crops instead of CROPLAND
if ~exist('luh2_base_orig','var')
    luh2_base_orig = luh2_base ;
end
luh2_base_cropa = luh2_base_cropf ;
luh2_base_cropa.maps_YXv = luh2_base_cropf.maps_YXv .* luh2_base_orig.maps_YXv(:,:,strcmp(luh2_base_orig.varNames,'CROPLAND')) ;
clear luh2_base
luh2_base.varNames = [luh2_base_cropa.varNames luh2_base_orig.varNames(~strcmp(luh2_base_orig.varNames,'CROPLAND'))] ;
luh2_base.maps_YXv = cat(3,luh2_base_cropa.maps_YXv,luh2_base_orig.maps_YXv(:,:,~strcmp(luh2_base_orig.varNames,'CROPLAND'))) ;

% Get info
LUnames = luh2_base.varNames ;
Nlu = length(LUnames) ;
notBare = ~strcmp(LUnames,'BARREN') ;
isAgri = ~strcmp(LUnames,'NATURAL') & notBare ;
LUnames_agri = LUnames(isAgri) ;
isAgri_isPast = strcmp(LUnames_agri,'PASTURE') ;

% Convert NaNs to zeros for addition
luh2_base.maps_YXv(isnan(luh2_base.maps_YXv)) = 0 ;

% Get LUH2 fraction that is vegetated, barren
luh2_vegd_YX = sum(luh2_base.maps_YXv(:,:,notBare),3) ;
luh2_vegdFrac_YX = luh2_vegd_YX ./ landArea_YX ;
luh2_bare_YX = luh2_base.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
luh2_bareFrac_YX = luh2_bare_YX ./ landArea_YX ;

% Aggregate from 0.5º to 2º
tmp = luh2_base.maps_YXv(:,1:4:end,:) ...
    + luh2_base.maps_YXv(:,2:4:end,:) ...
    + luh2_base.maps_YXv(:,3:4:end,:) ...
    + luh2_base.maps_YXv(:,4:4:end,:) ;
luh2_base_2deg.maps_YXv = tmp(1:4:end,:,:) ...
                        + tmp(2:4:end,:,:) ...
                        + tmp(3:4:end,:,:) ...
                        + tmp(4:4:end,:,:) ;
clear tmp
luh2_base_2deg.varNames = luh2_base.varNames ;
tmp = landArea_YX(:,1:4:end) ...
    + landArea_YX(:,2:4:end) ...
    + landArea_YX(:,3:4:end) ...
    + landArea_YX(:,4:4:end) ;
landArea_2deg_YX = tmp(1:4:end,:) ...
                 + tmp(2:4:end,:) ...
                 + tmp(3:4:end,:) ...
                 + tmp(4:4:end,:) ;
clear tmp

% Get LUH2 2-deg fraction that is vegetated, barren
luh2_vegd_2deg_YX = sum(luh2_base_2deg.maps_YXv(:,:,notBare),3) ;
luh2_vegdFrac_2deg_YX = luh2_vegd_2deg_YX ./ landArea_2deg_YX ;
luh2_bare_2deg_YX = luh2_base_2deg.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
luh2_bareFrac_2deg_YX = luh2_bare_2deg_YX ./ landArea_2deg_YX ;
             
% Get lats/lons
list2map_2deg = find(landArea_2deg_YX>0) ;
lons_2deg = lons_map_2deg(list2map_2deg) ;
lats_2deg = lats_map_2deg(list2map_2deg) ;
list2map = find(landArea_YX>0) ;
lons = lons_map(list2map) ;
lats = lats_map(list2map) ;

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



%% Import PLUM

disp('Importing PLUM...')

is2deg = strcmp(thisVer,'2deg.') ;

if is2deg
    bareFrac_y0_YX = luh2_bareFrac_2deg_YX ;
    thisLandArea_YX = landArea_2deg_YX ;
    Yincl = 15:90 ;
    thisluh2_2010_YXv = luh2_base_2deg.maps_YXv ;
else
    bareFrac_y0_YX = luh2_bareFrac_YX ;
    thisLandArea_YX = landArea_YX ;
    thisluh2_2010_YXv = luh2_base.maps_YXv ;
    Yincl = 60:360 ;
end

% Originals (only 0.5-degree possible)
if is2deg
    [~,PLUM_orig_year0] = PLUMharm_processPLUMin_areaCrops(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year0s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, bareFrac_y0_YX, PLUMtoLPJG, LPJGcrops, norm2extra) ;
    [~,PLUM_orig_year1] = PLUMharm_processPLUMin_areaCrops(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year1s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, bareFrac_y0_YX, PLUMtoLPJG, LPJGcrops, norm2extra) ;
    [~,PLUM_orig_year2] = PLUMharm_processPLUMin_areaCrops(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year2s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, bareFrac_y0_YX, PLUMtoLPJG, LPJGcrops, norm2extra) ;
else
    [PLUM_orig_year0,~] = PLUMharm_processPLUMin_areaCrops(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year0s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, bareFrac_y0_YX, PLUMtoLPJG, LPJGcrops, norm2extra) ;
    [PLUM_orig_year1,~] = PLUMharm_processPLUMin_areaCrops(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year1s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, bareFrac_y0_YX, PLUMtoLPJG, LPJGcrops, norm2extra) ;
    [PLUM_orig_year2,~] = PLUMharm_processPLUMin_areaCrops(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year2s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, bareFrac_y0_YX, PLUMtoLPJG, LPJGcrops, norm2extra) ;
end

% Harmonized
if is2deg
    [~,PLUM_harm_year1] = PLUMharm_processPLUMin_areaCrops(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm/' year1s '/LandCoverFract.base2010.' thisVer 'txt'],...
        thisLandArea_YX, LUnames, [], PLUMtoLPJG, LPJGcrops, 0) ;
    [~,PLUM_harm_year2] = PLUMharm_processPLUMin_areaCrops(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm/' year2s '/LandCoverFract.base2010.' thisVer 'txt'],...
        thisLandArea_YX, LUnames, [], PLUMtoLPJG, LPJGcrops, 0) ;
else
    [PLUM_harm_year1,~] = PLUMharm_processPLUMin_areaCrops(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm/' year1s '/LandCoverFract.base2010.' thisVer 'txt'],...
        thisLandArea_YX, LUnames, [], PLUMtoLPJG, LPJGcrops, 0) ;
    [PLUM_harm_year2,~] = PLUMharm_processPLUMin_areaCrops(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm/' year2s '/LandCoverFract.base2010.' thisVer 'txt'],...
        thisLandArea_YX, LUnames, [], PLUMtoLPJG, LPJGcrops, 0) ;
end

disp('Done importing PLUM.')


%% Map transitions for 2010 to 2011, including original L2010-->P2011
spacing = 0.025 ;
thisPos = figurePos ;

if year0~=2010
    error('This doesn''t make sense to do without year0==2010.')
end

figure('Color','w','Position',thisPos) ;

h = cell(12,1) ;
for v = 1:4
    h{v} = subplot_tight(3,4,v,spacing) ;
    if v==1
        vv = 1:8 ;
    else
        vv = v + 7 ;
    end
    thisDiff_YX = sum(PLUM_orig_year1.maps_YXv(:,:,vv),3) ...
                - sum(PLUM_orig_year0.maps_YXv(:,:,vv),3) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    if v==1
        title('Original year1-year0: CROPLAND')
    else
        title(['Original year1-year0: ' LUnames{vv}])
    end
    set(gca,'FontSize',14)
end
for v = 1:4
    h{4+v} = subplot_tight(3,4,4+v,spacing) ;
    if v==1
        vv = 1:8 ;
    else
        vv = v + 7 ;
    end
    thisDiff_YX = sum(PLUM_harm_year1.maps_YXv(:,:,vv),3) ...
                - sum(thisluh2_2010_YXv(:,:,vv),3) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    if v==1
        title('Harmonized year1-year0: CROPLAND')
    else
        title(['Harmonized year1-year0: ' LUnames{vv}])
    end
    set(gca,'FontSize',14)
end
for v = 1:4
    h{8+v} = subplot_tight(3,4,8+v,spacing) ;
    if v==1
        vv = 1:8 ;
    else
        vv = v + 7 ;
    end
    thisDiff_YX = sum(PLUM_orig_year1.maps_YXv(:,:,vv),3) ...
                - sum(thisluh2_2010_YXv(:,:,vv),3) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    if v==1
        title('Oyear1-Lyear0: CROPLAND')
    else
        title(['Oyear1-Lyear0: ' LUnames{vv}])
    end
    set(gca,'FontSize',14)
end
for v = 1:4
    new_caxis = get_combined_caxis(h{v},h{4+v}) ;
    caxis(h{v},new_caxis) ;
    caxis(h{4+v},new_caxis) ;
end


%% Crop maps: LUH2_2010 vs. Harmonized_2011
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

if year0~=2010
    error('This doesn''t make sense to do unless year0==2010.')
end

figure('Color','w','Position',thisPos) ;

h = cell(2,1) ;
for c = 1:Nagri
    plot1 = (c-1)*2 + 1 ;
    h1 = subplot_tight(4,4,plot1,spacing) ;
    thisMap_YX = thisluh2_2010_YXv(:,:,c) ;
    thisMap_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisMap_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','EastOutside')
    title(['LUH2 ' year1s ': ' LUnames{c}])
    set(gca,'FontSize',14)
    clear h
    
    plot2 = plot1 + 1 ;
    h2 = subplot_tight(4,4,plot2,spacing) ;
    thisMap_YX = PLUM_harm_year1.maps_YXv(:,:,c) ;
    thisMap_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisMap_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','EastOutside')
    title(['Harmonized ' year1s ': ' LUnames{c}])
    set(gca,'FontSize',14)
    clear h
    
    new_caxis = get_combined_caxis(h1,h2) ;
    caxis(h1,new_caxis) ;
    caxis(h2,new_caxis) ;
end


%% Map difference in transitions for year0-year1
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

for v = 1:(Nagri+1)
    subplot_tight(3,4,v,spacing)
    thisDiffO_YX = PLUM_orig_year1.maps_YXv(:,:,v) ...
                 - PLUM_orig_year0.maps_YXv(:,:,v) ;
    thisDiffH_YX = PLUM_harm_year1.maps_YXv(:,:,v) ...
                 - thisluh2_2010_YXv(:,:,v) ;
    thisDiff_YX = thisDiffH_YX - thisDiffO_YX ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['(H' year1s '-Lyear0)-(O' year1s '-Oyear0): ' LUnames{v}])
    set(gca,'FontSize',14)
end






%% Map transitions for year0 to year1
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

h = cell(8,1) ;
for v = 1:4
    h{v} = subplot_tight(2,4,v,spacing) ;
    thisDiff_YX = PLUM_orig_year1.maps_YXv(:,:,v) ...
                - PLUM_orig_year0.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Original ' year1s '-' year0s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    h{4+v} = subplot_tight(2,4,4+v,spacing) ;
    thisDiff_YX = PLUM_harm_year1.maps_YXv(:,:,v) ...
                - thisluh2_2010_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized ' year1s '-' year0s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
% for v = 1:4
%     new_caxis = get_combined_caxis(h{v},h{4+v}) ;
%     caxis(h{v},new_caxis) ;
%     caxis(h{4+v},new_caxis) ;
% end




%% Map transitions for year1 to year2
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

h = cell(8,1) ;
for v = 1:4
    h{v} = subplot_tight(2,4,v,spacing) ;
    thisDiff_YX = PLUM_orig_year2.maps_YXv(:,:,v) ...
                - PLUM_orig_year1.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Original ' year2s '-' year1s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    h{4+v} = subplot_tight(2,4,4+v,spacing) ;
    thisDiff_YX = PLUM_harm_year2.maps_YXv(:,:,v) ...
                - PLUM_harm_year1.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized ' year2s '-' year1s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    new_caxis = get_combined_caxis(h{v},h{4+v}) ;
    caxis(h{v},new_caxis) ;
    caxis(h{4+v},new_caxis) ;
end


%% Map difference in transitions for year1-year2
spacing = 0.025 ;

thisPos = figurePos ;

figure('Color','w','Position',thisPos) ;

for v = 1:4
    subplot_tight(2,2,v,spacing)
    thisDiffO_YX = PLUM_orig_year2.maps_YXv(:,:,v) ...
                 - PLUM_orig_year1.maps_YXv(:,:,v) ;
    thisDiffH_YX = PLUM_harm_year2.maps_YXv(:,:,v) ...
                 - PLUM_harm_year1.maps_YXv(:,:,v) ;
    thisDiff_YX = thisDiffH_YX - thisDiffO_YX ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['(H' year2s '-H' year1s ')-(O' year2s '-O' year1s '): ' LUnames{v}])
    set(gca,'FontSize',14)
end


%% Maps in year1
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

h = cell(8,1) ;
for v = 1:4
    h{v} = subplot_tight(2,4,v,spacing) ;
    thisDiff_YX = PLUM_orig_year1.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Original ' year1s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    h{4+v} = subplot_tight(2,4,4+v,spacing) ;
    thisDiff_YX = PLUM_harm_year1.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized ' year1s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    new_caxis = get_combined_caxis(h{v},h{4+v}) ;
    caxis(h{v},new_caxis) ;
    caxis(h{4+v},new_caxis) ;
end


%% Maps in year2
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

h = cell(8,1) ;
for v = 1:4
    h{v} = subplot_tight(2,4,v,spacing) ;
    thisDiff_YX = PLUM_orig_year2.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Original year2: ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    h{4+v} = subplot_tight(2,4,4+v,spacing) ;
    thisDiff_YX = PLUM_harm_year2.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized year2: ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    new_caxis = get_combined_caxis(h{v},h{4+v}) ;
    caxis(h{v},new_caxis) ;
    caxis(h{4+v},new_caxis) ;
end

%% Map differences in year1
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

for v = 1:4
    subplot_tight(1,4,v,spacing)
    thisDiff_YX = PLUM_harm_year1.maps_YXv(:,:,v) ...
                - PLUM_orig_year1.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized-Original year1: ' LUnames{v}])
    set(gca,'FontSize',14)
end


%% Map differences in year2
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

for v = 1:4
    subplot_tight(1,4,v,spacing)
    thisDiff_YX = PLUM_harm_year2.maps_YXv(:,:,v) ...
                - PLUM_orig_year2.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized-Original year2: ' LUnames{v}])
    set(gca,'FontSize',14)
end


%% Map "difference in differences" of year1 vs year2
% How persistent are the differences induced by the harmonization to year0?
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

for v = 1:4
    subplot_tight(1,4,v,spacing)
    thisDiff1_YX = PLUM_harm_year1.maps_YXv(:,:,v) ...
                 - PLUM_orig_year1.maps_YXv(:,:,v) ;
    thisDiff2_YX = PLUM_harm_year2.maps_YXv(:,:,v) ...
                 - PLUM_orig_year2.maps_YXv(:,:,v) ;
    thisDiff_YX = thisDiff2_YX - thisDiff1_YX ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['(H-O year2) - (H-O year1): ' LUnames{v}])
    set(gca,'FontSize',14)
end

