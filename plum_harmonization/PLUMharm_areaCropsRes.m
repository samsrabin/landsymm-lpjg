%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LUH1-style harmonization for PLUM outputs, at cropType level %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_year = 2010 ;
year1 = 2011 ;
yearN = 2020 ;

% Directory for PLUM outputs
PLUM_in_toptop = {...
                  '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1' ;
%                   '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v10.s1' ;
%                   '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP4.v10.s1' ;
%                   '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v10.s1' ;
                  } ;
              
% Save?
save_halfDeg = true ;
save_2deg = false ;

% Coordinates of 2-degree cell to debug (leave empty for no debug)
debugIJ_2deg = [] ;
% debugIJ_2deg = [67 42] ;
% debugIJ_2deg = [73 122] ;
% debugIJ_2deg = [23 176] ;
% debugIJ_2deg = [28 149] ;

% Method for inpaint_nans()
inpaint_method = 0 ;

% Use "latest PLUM management" in cells that don't have thisCrop thisYear
% but did in a previous year? FALSE = rely solely on interpolation of
% thisYear's thisMgmt.
useLatestPLUMmgmt = false ;




%% Setup

addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;

% Output file details
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;

% PLUM protected area info
PLUM_file_res_prot = '/Users/Shared/PLUM/input/protected_areas/protected_areas_with_points.txt' ;
PLUM_file_res_terr = '/Users/Shared/PLUM/input/protected_areas/maxcropfrac2.txt' ;

% How much divergence from PLUM-original transition is acceptable? (%)
conserv_tol_pct = 0.2 ;

% The fraction of real crops (by area) not included in PLUM's 6
% non-energyCrops commodities. Multiply PLUM area by (1-norm2extra) to
% get area of the crop itself without the norm2extra fraction included.
norm2extra = 0.177 ;

% Make lower-left lat/lon map (for compat. with PLUM style)
lons_map_2deg = repmat(-180:2:178,[90 1]) ;
lats_map_2deg = repmat((-90:2:88)',[1 180]) ;
lons_map = repmat(-180:0.5:179.5,[360 1]) ;
lats_map = repmat((-90:0.5:89.5)',[1 720]) ;

% Debug?
do_debug = ~isempty(debugIJ_2deg) ;
clear db2i db2j
if do_debug
    db2i = debugIJ_2deg(1) ;
    db2j = debugIJ_2deg(2) ;
end
debug_header = ':\t\tArea\t\tFrac\n' ;


%% Import reference data

disp('Importing reference data...')

% Read fraction of VEGETATED protected by...
%%% Terrain
resFrac_terr = lpjgu_matlab_readTable_then2map(PLUM_file_res_terr,'verboseIfNoMat',false,'force_mat_nosave',true) ;
resFrac_terr_YX = 1 - resFrac_terr.maps_YXv ;
resFrac_terr_YX(isnan(resFrac_terr_YX)) = 0 ;
clear resFrac_terr
%%% Slope
fid = fopen(PLUM_file_res_prot, 'rt') ;
resFrac_prot = textscan(fid, ['%f' repmat(' %f',[1 720-1])], 'headerLines', 6) ;
fclose(fid) ;
resFrac_prot_YX = flipud(cell2mat(resFrac_prot)) ;
clear resFrac_prot
resFrac_prot_YX(resFrac_prot_YX==-9999) = 0 ;

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

% Harmonize masks
PLUM_base_in = addslashifneeded([addslashifneeded(PLUM_in_toptop{1}) num2str(base_year)]) ;
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

% Import base_year crop fractions
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
% Get "irrigation intensity" as fraction of thisCrop that is irrigated
% (later, will interpolate to cells without any of thisCrop)
luh2_base_irrig.varNames = LPJGcrops ;
luh2_base_irrig.maps_YXv = nan([size(landArea_YX) Ncrops_lpjg],'single') ;
for c = 1:Ncrops_lpjg
    thisCrop = LPJGcrops{c} ;
    if strcmp(thisCrop,'ExtraCrop')
        luh2_base_irrig.maps_YXv(:,:,c) = zeros(size(landArea_YX)) ;
        continue
    end
    thisCropI = [thisCrop 'i'] ;
    thisCropR_YX = luh2_base_cropf.maps_YXv(:,:,strcmp(luh2_base_cropf.varNames,thisCrop)) ;
    thisCropI_YX = luh2_base_cropf.maps_YXv(:,:,strcmp(luh2_base_cropf.varNames,thisCropI)) ;
    thisCrop_YX = thisCropR_YX + thisCropI_YX ;
    tmp_YX = thisCropI_YX ./ thisCrop_YX ;
    tmp_YX(thisCrop_YX==0) = NaN ;
    luh2_base_irrig.maps_YXv(:,:,c) = tmp_YX ;
end
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
    clear luh2_base_cropf
    luh2_base_cropf = tmp ;
    clear tmp
end

% Read baseline fertilization
% (later, will interpolate to cells without any of thisCrop)
luh2_base_nfert = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/remaps_v4/nfert.remapv4.20180214.cgFertIrr0.setaside0103.m0.txt',...
    'verboseIfNoMat',false,'force_mat_save',true) ;
% Get just base year, if needed
if isfield(luh2_base_nfert,'yearList')
    tmp = luh2_base_nfert.maps_YXvy(:,:,luh2_base_nfert.yearList==base_year) ;
    luh2_base_nfert = rmfield(luh2_base_nfert,{'maps_YXvy','yearList'}) ;
    luh2_base_nfert.maps_YXv = tmp ;
    clear tmp
end
% Add ExtraCrop, if needed
if any(strcmp(LPJGcrops,'ExtraCrop')) && ~any(strcmp(luh2_base_nfert.varNames,'ExtraCrop'))
    luh2_base_nfert.varNames = [luh2_base_nfert.varNames {'ExtraCrop'}] ;
    luh2_base_nfert.maps_YXv = cat(3,luh2_base_nfert.maps_YXv,zeros(size(landArea_YX))) ;
end
[~,IA,IB] = intersect(LPJGcrops,luh2_base_nfert.varNames,'stable') ;
if ~isequal(IA',1:Ncrops_lpjg)
    error('Issue with crop compatibility between LPJGcrops and baseline Nfert.')
end
luh2_base_nfert.maps_YXv = luh2_base_nfert.maps_YXv(:,:,IB) ;
luh2_base_nfert.varNames = LPJGcrops ;

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
maxLength_LUnames = max(cellfun(@length,LUnames)) ;
Nlu = length(LUnames) ;
notBare = ~strcmp(LUnames,'BARREN') ;
isAgri = ~strcmp(LUnames,'NATURAL') & notBare ;
LUnames_agri = LUnames(isAgri) ;
isAgri_isPast = strcmp(LUnames_agri,'PASTURE') ;
isCrop = ~strcmp(LUnames,'NATURAL') & ~strcmp(LUnames,'PASTURE') & notBare ;

% Convert NaNs to zeros for addition
luh2_base.maps_YXv(isnan(luh2_base.maps_YXv)) = 0 ;

% Get LUH2 fraction that is vegetated, barren
luh2_vegd_YX = sum(luh2_base.maps_YXv(:,:,notBare),3) ;
luh2_vegdFrac_YX = luh2_vegd_YX ./ landArea_YX ;
luh2_bare_YX = luh2_base.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
luh2_bareFrac_YX = luh2_bare_YX ./ landArea_YX ;

% Aggregate from 0.5º to 2º
luh2_base_2deg.maps_YXv = PLUMharm_aggregate(luh2_base.maps_YXv,0.5,2) ;
luh2_base_2deg.varNames = luh2_base.varNames ;
luh2_base_nfert_2deg.maps_YXv = PLUMharm_aggregate_mgmt(luh2_base_nfert.maps_YXv,luh2_base.maps_YXv(:,:,isCrop),0.5,2) ;
luh2_base_nfert_2deg.varNames = LPJGcrops ;
luh2_base_irrig_2deg.maps_YXv = PLUMharm_aggregate_mgmt(luh2_base_irrig.maps_YXv,luh2_base.maps_YXv(:,:,isCrop),0.5,2) ;
luh2_base_irrig_2deg.varNames = LPJGcrops ;
landArea_2deg_YX = PLUMharm_aggregate(landArea_YX,0.5,2) ;

% Interpolate management inputs
luh2_base_nfert_2deg.maps_YXv = PLUMharm_interpolateMgmt(...
    luh2_base_nfert_2deg.maps_YXv, luh2_base.maps_YXv(:,:,isCrop), landArea_2deg_YX,...
    LPJGcrops, inpaint_method) ;
luh2_base_irrig_2deg.maps_YXv = PLUMharm_interpolateMgmt(...
    luh2_base_irrig_2deg.maps_YXv, luh2_base.maps_YXv(:,:,isCrop), landArea_2deg_YX,...
    LPJGcrops, inpaint_method) ;

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

disp('Done importing reference data.')


%% Do it

% The years we want to produce PLUM outputs for (will begin with
% transitions from years(1)-1 to years(1)
years = year1:yearN ;
if year1 <= base_year
    error('year1 (%d) must be > base_year (%d)!\n', year1, base_year)
elseif year1 > yearN
    error('year1 (%d) must be <= yearN (%d)!\n', year1, yearN)
end

for d = 1:length(PLUM_in_toptop)

    %%%%%%%%%%%%%
    %%% Setup %%%
    %%%%%%%%%%%%%

    % Get directories
    PLUM_in_top = removeslashifneeded(PLUM_in_toptop{d}) ;
    disp(PLUM_in_top)
    PLUM_out_top = addslashifneeded([PLUM_in_top '.harm/']) ;
    unix(['mkdir -p ' PLUM_out_top]) ;
    PLUM_in_top = addslashifneeded(PLUM_in_top) ;

    % Read fraction of VEGETATED protected by...
    %%% PLUM's minimum natural fraction ("rate")
    [~,r] = unix(['grep "MIN_NATURAL_RATE" ' PLUM_in_top 'config.properties | sed "s/MIN_NATURAL_RATE=//"  | tr -d ''\n''']) ;
    min_natural_rate = str2double(r) ;
    resFrac_minN_YX = min_natural_rate * ones(size(resFrac_prot_YX)) ;
    %%% Total
    resFrac_YX = 1 - (1-resFrac_terr_YX).*(1-resFrac_prot_YX)*(1-min_natural_rate) ;
    % Convert to area
    resArea_YX = resFrac_YX .* luh2_vegd_YX ;
    resArea_2deg_YX = PLUMharm_aggregate(resArea_YX,0.5,2) ;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Process PLUM outputs %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Debugging output
    if do_debug
        fprintf('Center Lon %.2f, Lat %.2f, land area %0.4g\n\n',1+lons_map_2deg(db2i,db2j),1+lats_map_2deg(db2i,db2j),round(landArea_2deg_YX(db2i,db2j),4)) ;
        PLUMharm_debugOut('LUH2',luh2_base_2deg.maps_YXv,landArea_2deg_YX,debugIJ_2deg,LUnames)
    end

    bareFrac_y0_YX = luh2_bareFrac_YX ;
    if useLatestPLUMmgmt
        latestPLUMin_nfert_2deg_YXv = -1*ones([size(landArea_2deg_YX) Ncrops_lpjg]) ;
        latestPLUMin_irrig_2deg_YXv = -1*ones([size(landArea_2deg_YX) Ncrops_lpjg]) ;
    else
        latestPLUMin_nfert_2deg_YXv = [] ;
        latestPLUMin_irrig_2deg_YXv = [] ;
    end
    clear out_* in_*
    for y = 1:length(years)
        
        thisYear = years(y) ;
        disp(num2str(thisYear)) ;
        tic ;
        
        % Import previous harmonized year, if needed
        if thisYear-1 == base_year
            out_y0 = luh2_base ;
            disp('out_y0 from luh2')
            out_y0_2deg = luh2_base_2deg ;
            out_y0_agri_YXv = luh2_base.maps_YXv(:,:,isAgri) ;
            out_y0_2deg_agri_YXv = luh2_base_2deg.maps_YXv(:,:,isAgri) ;
            out_y0_nfert_2deg = luh2_base_nfert_2deg ;
            out_y0_irrig_2deg = luh2_base_irrig_2deg ;
        elseif y==1
            file_in = [removeslashifneeded(PLUM_in_top) '.harm/' num2str(thisYear-1) 'post.base' num2str(base_year) '.mat'] ;
            disp(['*y0* from ' file_in])
            load(file_in) ;
        end
        
        % Import for debugging redistribution
        if do_debug
            if y==1
                file_in = [PLUM_base_in '/LandCoverFract.txt'] ;
                [~, ~, ~, in_y0orig_2deg, in_y0orig_nfert_2deg, in_y0orig_irrig_2deg, ...
                    ~, ~] = ...
                    PLUMharm_processPLUMin_areaCrops(file_in, landArea_YX, landArea_2deg_YX, ...
                    LUnames, bareFrac_y0_YX, ...
                    latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
                    PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method) ;
            end
        else
            diffO_YXv = [] ;
            in_y0orig_2deg = [] ;
            in_y1orig_2deg = [] ;
        end

        % Import and process previous year, if needed
        if ~exist('in_y0','var')
            file_in = [PLUM_in_top num2str(thisYear-1) '/LandCoverFract.txt'] ;
            disp(['in_y0 from ' file_in])
            [in_y0, ~, ~, in_y0_2deg, in_y0_nfert_2deg, in_y0_irrig_2deg, ...
                latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv] = ...
                PLUMharm_processPLUMin_areaCrops(file_in, landArea_YX, landArea_2deg_YX, ...
                LUnames, bareFrac_y0_YX, latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
                PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method) ;
            bareFrac_y0_YX = in_y0.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ./ landArea_YX ;
            in_y0_agri_YXv = in_y0.maps_YXv(:,:,isAgri) ;
            in_y0_2deg_agri_YXv = in_y0_2deg.maps_YXv(:,:,isAgri) ;
            in_y0_2deg_vegd_YX = sum(in_y0_2deg.maps_YXv(:,:,notBare),3) ;
        end

        % Import this year and convert to area
        file_in = [PLUM_in_top num2str(thisYear) '/LandCoverFract.txt'] ;
%         disp(['in_y1 from ' file_in])
        [in_y1, ~, ~, in_y1_2deg, in_y1_nfert_2deg, in_y1_irrig_2deg, ...
            latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv] = ...
            PLUMharm_processPLUMin_areaCrops(file_in, landArea_YX, landArea_2deg_YX, ...
            LUnames, bareFrac_y0_YX, latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
            PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method) ;
        in_y1_agri_YXv = in_y1.maps_YXv(:,:,isAgri) ;
        in_y1_2deg_agri_YXv = in_y1_2deg.maps_YXv(:,:,isAgri) ;
        in_y1_ntrl_YX = in_y1.maps_YXv(:,:,strcmp(in_y1.varNames,'NATURAL')) ;
        in_y1_bare_YX = in_y1.maps_YXv(:,:,strcmp(in_y1.varNames,'BARREN')) ;
        in_y1_vegd_YX = sum(in_y1.maps_YXv(:,:,notBare),3) ;
        in_y1_vegdFrac_YX = in_y1_vegd_YX ./ landArea_YX ;
        in_y1_bareFrac_YX = in_y1_bare_YX ./ landArea_YX ;
        in_y1_2deg_vegd_YX = sum(in_y1_2deg.maps_YXv(:,:,notBare),3) ;
        
        % Import for debugging redistribution
        if do_debug
            in_y1orig_2deg = in_y1_2deg ;
            diffO_YXv = in_y1orig_2deg.maps_YXv - in_y0orig_2deg.maps_YXv ;
        end
        
        
        % Make sure that total global loss of a land use does not exceed
        % the total global area of that land use.
        % NOTE: v10.SSP4 has Miscanthus in 2010, and some gridcells lose
        % Miscanthus from 2010 to 2011. This cannot be satisfied using my
        % 2010 maps because I have no Miscanthus. However, it shouldn't be
        % a problem as long as the grid cells PLUM specifies as losing
        % Miscanthus lose ALL of their Miscanthus.
        for i = 1:length(LUnames)
            thisLU = LUnames{i} ;
            out0_thisLU_YX = out_y0.maps_YXv(:,:,i) ;
            in0_thisLU_YX = in_y0.maps_YXv(:,:,i) ;
            in1_thisLU_YX = in_y1.maps_YXv(:,:,i) ;
            inDiff_thisLU_YX = in1_thisLU_YX - in0_thisLU_YX ;
            globLoss_thisLU = sum(inDiff_thisLU_YX(inDiff_thisLU_YX<0)) ;
            globArea_thisLU = sum(out0_thisLU_YX(:)) ;
            NETglobChg_thisLU = sum(inDiff_thisLU_YX(:)) ;
            if -globLoss_thisLU > globArea_thisLU
                disp(LUnames{i})
                fprintf('globArea_y0_thisLu = %0.4e\n',sum(in0_thisLU_YX(:))) ;
                fprintf('globArea_y1_thisLu = %0.4e\n',sum(in1_thisLU_YX(:))) ;
                fprintf('NETglobChg_thisLU  = %0.4e\n',NETglobChg_thisLU) ;
                fprintf('globLoss_thisLU    = %0.4e\n',globLoss_thisLU) ;
                fprintf('globArea_thisLU    = %0.4e\n',globArea_thisLU) ;
                if strcmp(thisLU,'Miscanthus') ...
                && thisYear-1==base_year
                    if sum(in1_thisLU_YX(in0_thisLU_YX>0)) == 0
                        warning('PLUM has Miscanthus in 2010, and some gridcells lose Miscanthus from 2010 to 2011. This cannot be satisfied using my 2010 maps because I have no Miscanthus. However, it shouldn''t be a problem because the grid cells PLUM specifies as losing Miscanthus lose ALL it.')
                    else
                        warning('PLUM has Miscanthus in 2010, and some gridcells lose Miscanthus from 2010 to 2011. This cannot be satisfied using my 2010 maps because I have no Miscanthus. But not all grid cells PLUM specifies as losing Miscanthus lose ALL of it! Difference of %0.1g km2. Ignoring.',sum(in1_thisLU_YX(in0_thisLU_YX>0)))
                    end
                    tmp = in_y0.maps_YXv(:,:,i) ;
                    tmp(inDiff_thisLU_YX<0) = in1_thisLU_YX(inDiff_thisLU_YX<0) ;
                    in_y0.maps_YXv(:,:,i) = tmp ;
                    in_y0_2deg.maps_YXv = PLUMharm_aggregate(in_y0.maps_YXv,0.5,2) ;
                    in_y0_agri_YXv = in_y0.maps_YXv(:,:,isAgri) ;
                    in_y0_2deg_agri_YXv = in_y0_2deg.maps_YXv(:,:,isAgri) ;
                elseif 0 <= NETglobChg_thisLU && NETglobChg_thisLU < 1e-3
                    warning('Total global loss of %s (%0.4e km2) exceeds its y0 area (%0.4e km2). LIKELY WILL cause infinite loop in ringRedist.\n',thisLU,-globLoss_thisLU,globArea_thisLU) ;
                else
                    warning('Total global loss of %s (%0.4e km2) exceeds its y0 area (%0.4e km2). (((May))) cause infinite loop in ringRedist.\n',thisLU,-globLoss_thisLU,globArea_thisLU) ;
                end
                pause(0.1)
            end
        end; clear *_thisLU*

        % Debugging output
        if do_debug
            PLUMharm_debugOut_deltas('iny0_to_iny1',in_y0_2deg.maps_YXv,in_y1_2deg.maps_YXv,landArea_2deg_YX,debugIJ_2deg,LUnames)
        end

        % Calculate changes in PLUM agri grids at 2 degrees
        %%% Negative indicates LOSS of thisAgri area
        agri_d_YXv = in_y1_2deg_agri_YXv - in_y0_2deg_agri_YXv ;
        nfert_d_YXv = in_y1_nfert_2deg.maps_YXv - in_y0_nfert_2deg.maps_YXv ;
        irrig_d_YXv = in_y1_irrig_2deg.maps_YXv - in_y0_irrig_2deg.maps_YXv ;

        % Apply changes to previous grid (@2deg)
        mid1_y1_2deg_agri_YXv = out_y0_2deg_agri_YXv + agri_d_YXv ;
        mid_y1_2deg_nfert_YXv = out_y0_nfert_2deg.maps_YXv + nfert_d_YXv ;
        mid_y1_2deg_irrig_YXv = out_y0_irrig_2deg.maps_YXv + irrig_d_YXv ;
        
        % Do not allow invalid management inputs. Could conceivably be
        % moved to after ringRedist, which would in some sense get at the
        % "redistribution of negative or excess management inputs," which
        % we are not currently doing.
        mid_y1_2deg_nfert_YXv(mid_y1_2deg_nfert_YXv<0) = 0 ;
        mid_y1_2deg_irrig_YXv(mid_y1_2deg_irrig_YXv<0) = 0 ;
        mid_y1_2deg_irrig_YXv(mid_y1_2deg_irrig_YXv>1) = 1 ;

        % Debugging output
        mid1_y1_2deg_ntrl_YX = out_y0_2deg.maps_YXv(:,:,strcmp(out_y0_2deg.varNames,'NATURAL')) ...
                             - sum(agri_d_YXv,3) ;
        if do_debug
            PLUMharm_debugOut_deltas('outy0_to_mid1y1', ...
                out_y0_2deg.maps_YXv, ...
                cat(3,mid1_y1_2deg_agri_YXv,mid1_y1_2deg_ntrl_YX,luh2_bare_2deg_YX), ...
                landArea_2deg_YX,debugIJ_2deg,LUnames) ;
        end

        % Get land available for conversion to agriculture,
        % as max(NATURAL-RESERVED,0)
        out_y0_2deg_ntrl_YX = out_y0_2deg.maps_YXv(:,:,strcmp(out_y0_2deg.varNames,'NATURAL')) ;
        nonResNtrl_YX = max(out_y0_2deg_ntrl_YX - resArea_2deg_YX, 0) ;

        % Check that neither crop nor pasture exceeds nonBare area;
        % check that neither crop nor pasture is negative;
        % check that nonBare area is conserved;
        % compute the total amount of crop or pasture increase / decrease that
        % is not able to be met within the 2 degree gridcells
        %%% Negative total_unmet indicates TOO MUCH LOSS of thisAgri area
        [total_unmet_agri_YXv, ...
            mid_y1_2deg_agri_YXv, mid_y1_2deg_ntrl_YX] = ...
            PLUMharm_getUnmet_cropAreaRes(...
            mid1_y1_2deg_agri_YXv, ...
            luh2_vegd_2deg_YX, ...
            resArea_2deg_YX, sum(out_y0_2deg_agri_YXv,3), debugIJ_2deg) ;

        % Check 2: Check that global area changes are (mostly) conserved
        for i = 1:Nagri
            agri_d2_YX = mid_y1_2deg_agri_YXv(:,:,i) - out_y0_2deg_agri_YXv(:,:,i) + total_unmet_agri_YXv(:,:,i) ;
            agri_d_YX = agri_d_YXv(:,:,i) ;
            agri_d_glob_1 = sum(agri_d_YX(:)) ;
            agri_d_glob_2 = sum(agri_d2_YX(:)) ;
            if abs((agri_d_glob_2-agri_d_glob_1)/agri_d_glob_1*100) > conserv_tol_pct
                error(['Global ' LUnames{i} ' area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 2)'])
            end
        end

        % Debugging output
        if do_debug
            PLUMharm_debugOut_deltas('outy0_to_midy1', ...
                out_y0_2deg.maps_YXv, ...
                cat(3,mid_y1_2deg_agri_YXv,mid_y1_2deg_ntrl_YX,luh2_bare_2deg_YX), ...
                landArea_2deg_YX,debugIJ_2deg,LUnames) ;
        end

        % Loop through every 2-degree gridcell. If a gridcell has unmet crop
        % or pasture, look for place to put this unmet amount in neighboring
        % rings, starting with gridcells that are 1 unit away, then 2, etc.
        % until all unmet has been displaced to new 2 degree cells. Track the
        % displaced crop and pasture and the # of "rings" needed for each
        % 2-degree gridcell.
        [out_y1_2deg_agri_YXv, out_y1_2deg_nfert_YXv, out_y1_2deg_irrig_YXv] = ...
            PLUMharm_ringRedist_areaCropsRes(...
            mid_y1_2deg_agri_YXv, ...
            luh2_vegd_2deg_YX, ...
            total_unmet_agri_YXv, ...
            landArea_2deg_YX, ...
            debugIJ_2deg, in_y0orig_2deg, in_y1orig_2deg, out_y0_2deg_agri_YXv, ...
            nonResNtrl_YX, LUnames_agri, ...
            mid_y1_2deg_nfert_YXv, mid_y1_2deg_irrig_YXv) ;
        
        % Debugging output
        out_y1_2deg_ntrl_YX = landArea_2deg_YX - (sum(out_y1_2deg_agri_YXv,3) + luh2_bare_2deg_YX) ;
        if do_debug
            PLUMharm_debugOut_deltas('outy0_to_outy1', ...
                out_y0_2deg.maps_YXv, ...
                cat(3,out_y1_2deg_agri_YXv,out_y1_2deg_ntrl_YX,luh2_bare_2deg_YX), ...
                landArea_2deg_YX,debugIJ_2deg,LUnames) ;
        end

        % Check 3: Check that global area changes are (mostly) conserved
        diffH_YXv = out_y1_2deg_agri_YXv - out_y0_2deg_agri_YXv ;
        for i = 1:Nagri
            agri_d_YX = agri_d_YXv(:,:,i) ;
            agri_d_glob_1 = sum(agri_d_YX(:)) ;
            agri_d3_YX = out_y1_2deg_agri_YXv(:,:,i) - out_y0_2deg_agri_YXv(:,:,i) ;
            agri_d_glob_3 = sum(agri_d3_YX(:)) ;
            if abs((agri_d_glob_3-agri_d_glob_1)/agri_d_glob_1*100) > conserv_tol_pct
                warning(['Global ' LUnames{i} ' area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (' num2str((agri_d_glob_3-agri_d_glob_1)/agri_d_glob_1*100) '%) (step 3)'])
            end
        end

        % Loop through all 2-degree gridcells and distribute the new crop and
        % pasture area changes to the half-degree gridcells within.
        %
        % If the crop or pasture change is a decrease, apply the same
        % PERCENTAGE decrease to all half-degree crop or pasture gridcells 
        % within the 2-degree cell.
        %
        % If the crop or pasture change is an increase, apply this increase to 
        % all half-degree gridcells within the 2-degree cell proportionally 
        % according to available land area. Make sure that total area within a 
        % half-degree gridcell does not exceed 1 or 0.
        out_y1_agri_YXv = ...
            PLUMharm_distDeltas_areaCrops_recursive( ...
            landArea_YX, landArea_2deg_YX, ...
            out_y0_2deg_agri_YXv, out_y1_2deg_agri_YXv, out_y0_agri_YXv, ...
            luh2_vegd_YX, luh2_vegd_YX, conserv_tol_pct) ;
        
        % Distribute management inputs from 2-degree to half-degree cells.
        out_y1_nfert_YXv = PLUMharm_distMgmt(out_y1_2deg_nfert_YXv,2,0.5) ;
        out_y1_irrig_YXv = PLUMharm_distMgmt(out_y1_2deg_irrig_YXv,2,0.5) ;

        % Ensure cells in range [0,in_y1_vegd_YX] (extreme deviations were checked
        % in loop)
        if any(isnan(out_y1_agri_YXv(:)))
            error('How did you get a NaN in out_y1_agri_YXv?')
        end
        out_y1_agri_YXv(out_y1_agri_YXv<0) = 0 ;
        luh2_vegd_YXv = repmat(luh2_vegd_YX,[1 1 Nagri]) ;
        out_y1_agri_YXv(out_y1_agri_YXv>luh2_vegd_YXv) = out_y1_agri_YXv(out_y1_agri_YXv>luh2_vegd_YXv) ;
        clear luh2_vegd_YXv

        % Check 5: Check that global area changes are (mostly) conserved
        for i = 1:Nagri
            this_d1_halfDeg_YX = in_y1_agri_YXv(:,:,i) - in_y0_agri_YXv(:,:,i) ;
            this_d5_halfDeg_YX = out_y1_agri_YXv(:,:,i) - out_y0_agri_YXv(:,:,i) ;
            this_d_glob_halfDeg_1 = sum(this_d1_halfDeg_YX(:)) ;
            this_d_glob_halfDeg_5 = sum(this_d5_halfDeg_YX(:)) ;
            if abs((this_d_glob_halfDeg_5-this_d_glob_halfDeg_1)/this_d_glob_halfDeg_1*100) > conserv_tol_pct
                warning(['Global ' LUnames{i} ' area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (' num2str((this_d_glob_halfDeg_5-this_d_glob_halfDeg_1)/this_d_glob_halfDeg_1*100) '%) (step 5)'])
            end
        end

        % Get land use areas
        out_y1_bare_YX = in_y1.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
        out_y1_ntrl_YX = landArea_YX ...
            - (sum(out_y1_agri_YXv,3) + out_y1_bare_YX) ;
        out_y1_past_YX = out_y1_agri_YXv(:,:,isAgri_isPast) ;
        out_y1_crop_YX = sum(out_y1_agri_YXv(:,:,~isAgri_isPast),3) ;
        out_y1_2deg_bare_YX = in_y1_2deg.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
        out_y1_2deg_ntrl_YX = landArea_2deg_YX ...
            - (sum(out_y1_2deg_agri_YXv,3) + out_y1_2deg_bare_YX) ;
        out_y1_2deg_past_YX = out_y1_2deg_agri_YXv(:,:,isAgri_isPast) ;
        out_y1_2deg_crop_YX = sum(out_y1_2deg_agri_YXv(:,:,~isAgri_isPast),3) ;
        
        if save_halfDeg
            % Save new LandCoverFract.txt (0.5-degree)
            unix(['mkdir -p ' PLUM_out_top num2str(thisYear)]) ;
            out_y1.varNames = {'PASTURE','CROPLAND','NATURAL','BARREN'} ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'CROPLAND')) = out_y1_crop_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'PASTURE')) = out_y1_past_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'NATURAL')) = out_y1_ntrl_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'BARREN')) = out_y1_bare_YX ;
            out_y1.maps_YXv = out_y1.maps_YXv ./ repmat(landArea_YX,[1 1 length(out_y1.varNames)]) ;
            [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
            file_out = [PLUM_out_top num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.txt'] ;
            lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                'outPrec', outPrec, ...
                'outWidth', outWidth, ...
                'delimiter', delimiter, ...
                'overwrite', overwrite, ...
                'fancy', fancy, ...
                'progress_step_pct', 20, ...
                'verbose', false) ;
            clear out_y1
            
            % Save new CropFract.txt (0.5-degree)
            out_y1.maps_YXv = out_y1_agri_YXv(:,:,~isAgri_isPast) ./ repmat(out_y1_crop_YX,[1 1 length(LPJGcrops)]) ;
            out_y1.maps_YXv(repmat(out_y1_crop_YX,[1 1 length(LPJGcrops)])==0) = 0 ;
            out_y1.varNames = LPJGcrops ;
            [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
            file_out = [PLUM_out_top num2str(thisYear) '/CropFract.base' num2str(base_year) '.txt'] ;
            lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                'outPrec', outPrec, ...
                'outWidth', outWidth, ...
                'delimiter', delimiter, ...
                'overwrite', overwrite, ...
                'fancy', fancy, ...
                'progress_step_pct', 20, ...
                'verbose', false) ;
            clear out_y1
            
            % Save new Fert.txt (0.5-degree)
            % Convert from kg/m2 to kg/ha for compatibility with original
            % PLUM style.
            out_y1.maps_YXv = 1e4*out_y1_nfert_YXv ;
            out_y1.varNames = LPJGcrops ;
            [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
            file_out = [PLUM_out_top num2str(thisYear) '/Fert.base' num2str(base_year) '.txt'] ;
            lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                'outPrec', outPrec, ...
                'outWidth', outWidth, ...
                'delimiter', delimiter, ...
                'overwrite', overwrite, ...
                'fancy', fancy, ...
                'progress_step_pct', 20, ...
                'verbose', false) ;
            clear out_y1
            
            % Save new Irrig.txt (0.5-degree)
            out_y1.maps_YXv = out_y1_irrig_YXv ;
            out_y1.varNames = LPJGcrops ;
            [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
            file_out = [PLUM_out_top num2str(thisYear) '/Irrig.base' num2str(base_year) '.txt'] ;
            lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                'outPrec', outPrec, ...
                'outWidth', outWidth, ...
                'delimiter', delimiter, ...
                'overwrite', overwrite, ...
                'fancy', fancy, ...
                'progress_step_pct', 20, ...
                'verbose', false) ;
            clear out_y1
            
        end

        if save_2deg
            % Save new LandCoverFract.txt (2-degree)
            unix(['mkdir -p ' PLUM_out_top num2str(thisYear)]) ;
            out_y1.varNames = {'PASTURE','CROPLAND','NATURAL','BARREN'} ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'CROPLAND')) = out_y1_2deg_crop_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'PASTURE')) = out_y1_2deg_past_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'NATURAL')) = out_y1_2deg_ntrl_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'BARREN')) = out_y1_2deg_bare_YX ;
            out_y1.maps_YXv = out_y1.maps_YXv ./ repmat(landArea_2deg_YX,[1 1 length(out_y1.varNames)]) ;
            [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
            file_out = [PLUM_out_top num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.2deg.txt'] ;
            lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                'outPrec', outPrec, ...
                'outWidth', outWidth, ...
                'delimiter', delimiter, ...
                'overwrite', overwrite, ...
                'fancy', fancy, ...
                'progress_step_pct', 20, ...
                'verbose', false) ;
            clear out_y1

            % Save new CropFract.txt (2-degree)
            out_y1.maps_YXv = out_y1_2deg_agri_YXv(:,:,~isAgri_isPast) ./ repmat(out_y1_2deg_crop_YX,[1 1 length(LPJGcrops)]) ;
            out_y1.maps_YXv(repmat(out_y1_2deg_crop_YX,[1 1 length(LPJGcrops)])==0) = 0 ;
            out_y1.varNames = LPJGcrops ;
            [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
            file_out = [PLUM_out_top num2str(thisYear) '/CropFract.base' num2str(base_year) '.2deg.txt'] ;
            lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                'outPrec', outPrec, ...
                'outWidth', outWidth, ...
                'delimiter', delimiter, ...
                'overwrite', overwrite, ...
                'fancy', fancy, ...
                'progress_step_pct', 20, ...
                'verbose', false) ;
            clear out_y1
            
            % Save new Fert.txt (2-degree)
            % Convert from kg/m2 to kg/ha for compatibility with original
            % PLUM style.
            out_y1.maps_YXv = 1e4*out_y1_2deg_nfert_YXv ;
            out_y1.varNames = LPJGcrops ;
            [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
            file_out = [PLUM_out_top num2str(thisYear) '/Fert.base' num2str(base_year) '.2deg.txt'] ;
            lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                'outPrec', outPrec, ...
                'outWidth', outWidth, ...
                'delimiter', delimiter, ...
                'overwrite', overwrite, ...
                'fancy', fancy, ...
                'progress_step_pct', 20, ...
                'verbose', false) ;
            clear out_y1
            
            % Save new Irrig.txt (2-degree)
            out_y1.maps_YXv = out_y1_2deg_irrig_YXv ;
            out_y1.maps_YXv(repmat(out_y1_2deg_crop_YX,[1 1 length(LPJGcrops)])==0) = 0 ;
            out_y1.varNames = LPJGcrops ;
            [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
            file_out = [PLUM_out_top num2str(thisYear) '/Irrig.base' num2str(base_year) '.2deg.txt'] ;
            lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                'outPrec', outPrec, ...
                'outWidth', outWidth, ...
                'delimiter', delimiter, ...
                'overwrite', overwrite, ...
                'fancy', fancy, ...
                'progress_step_pct', 20, ...
                'verbose', false) ;
            clear out_y1
        end

        % Prepare for next iteration
        clear out_y1_past_YX
        if y < length(years)
            in_y0 = in_y1 ;
            in_y0_2deg = in_y1_2deg ;
            in_y0_agri_YXv = in_y1_agri_YXv ;
            in_y0_2deg_agri_YXv = in_y1_2deg_agri_YXv ;
            [~,IA,IB] = intersect(out_y0.varNames,LUnames_agri,'stable') ;
            out_y0.maps_YXv(:,:,IA) = out_y1_agri_YXv(:,:,IB) ;
            clear IA IB
            out_y0.maps_YXv(:,:,strcmp(LUnames,'NATURAL')) = out_y1_ntrl_YX ;
            out_y0.maps_YXv(:,:,strcmp(LUnames,'BARREN')) = out_y1_bare_YX ;
            out_y0_2deg_agri_YXv = out_y1_2deg_agri_YXv ;
            out_y0_agri_YXv = out_y1_agri_YXv ;
            bareFrac_y0_YX = in_y1_bareFrac_YX ;
            out_y0_nfert_2deg.maps_YXv = out_y1_2deg_nfert_YXv ;
            out_y0_irrig_2deg.maps_YXv = out_y1_2deg_irrig_YXv ;
            in_y0_2deg_vegd_YX = in_y1_2deg_vegd_YX ;
            in_y0_nfert_2deg = in_y1_nfert_2deg ;
            in_y0_irrig_2deg = in_y1_irrig_2deg ;
            if do_debug
                in_y0orig_2deg = in_y1orig_2deg ;
            end
            clear *y1*
        end
        
        % Save full-precision outputs for use in restarting
        if save_2deg || save_halfDeg
            save([PLUM_out_top num2str(thisYear) 'post.base' num2str(base_year) '.mat'], ...
                '*y0*','latestPLUMin_*','-v7.3') ;
        end

        disp(['  Done (' toc_hms(toc) ').'])


    end % years loop

    disp('Done')

end






