disp('Importing reference data...')

% Make lower-left lat/lon map (for compat. with PLUM style)
lons_map_2deg = repmat(-180:2:178,[90 1]) ;
lats_map_2deg = repmat((-90:2:88)',[1 180]) ;
lons_map = repmat(-180:0.5:179.5,[360 1]) ;
lats_map = repmat((-90:0.5:89.5)',[1 720]) ;

% Read fraction of VEGETATED protected by...
%%% Terrain
PLUM_file_res_terr = '/Users/Shared/PLUM/input/protected_areas/maxcropfrac2.txt' ;
resFrac_terr = lpjgu_matlab_readTable_then2map(PLUM_file_res_terr,'verboseIfNoMat',false,'force_mat_nosave',true) ;
resFrac_terr_YX = 1 - resFrac_terr.maps_YXv ;
resFrac_terr_YX(isnan(resFrac_terr_YX)) = 0 ;
clear resFrac_terr
%%% Protected areas
PLUM_file_res_prot = '/Users/Shared/PLUM/input/protected_areas/protected_areas_with_points.txt' ;
fid = fopen(PLUM_file_res_prot, 'rt') ;
resFrac_prot = textscan(fid, ['%f' repmat(' %f',[1 720-1])], 'headerLines', 6) ;
fclose(fid) ;
resFrac_prot_YX = flipud(cell2mat(resFrac_prot)) ;
clear resFrac_prot
resFrac_prot_YX(resFrac_prot_YX==-9999) = 0 ;

% Get LUH2 land area (m2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = 1e6*transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
landArea_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear *_YXqd

% Import LUH2 base_year
luh2_file = '/Users/Shared/PLUM/input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
base = lpjgu_matlab_readTable_then2map(luh2_file, 'force_mat_save', true);%, 'verbose', true, 'verboseIfNoMat', true) ;
if ~isempty(find(base.maps_YXvy(:,:,contains(base.varNames,{'URBAN','PEATLAND'}),:)>0,1))
    error('This code is not designed to handle LUH2 inputs with any URBAN or PEATLAND area!')
end
if doHarm
    base.maps_YXv = base.maps_YXvy(:,:,~contains(base.varNames,{'URBAN','PEATLAND'}),base.yearList==base_year) ;
else
    [~,IA,~] = intersect(base.yearList,yearList_luh2) ;
    if ~isequal(IA-min(IA)+1,(1:length(yearList_luh2))')
        error('Incompatible yearList_luh2?')
    end
    base.maps_YXvy = base.maps_YXvy(:,:,~contains(base.varNames,{'URBAN','PEATLAND'}),IA) ;
    base.yearList = yearList_luh2 ;
    clear IA IB
end
base.varNames = base.varNames(~contains(base.varNames,{'URBAN','PEATLAND'})) ;

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
LUnames = base.varNames ;
Nlu = length(LUnames) ;
landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;
if ~doHarm
    Nyears_luh2 = length(yearList_luh2) ;
end

% Avoid, to high precision, sum(LU) > 1
if doHarm
    tmp_YX = sum(base.maps_YXv,3) ;
    j = 0 ;
    while any(tmp_YX(:)-1 > 3*eps)
        j = j + 1 ;
        if j > 50
            error('Possible infinite loop in "Avoid, to high precision, sum(LU) > 1"')
        end
        base.maps_YXv = base.maps_YXv ./ repmat(tmp_YX,[1 1 Nlu]) ;
        tmp_YX = sum(base.maps_YXv,3) ;
    end
    clear tmp_YX j
else
    tmp_YX1y = sum(base.maps_YXvy,3) ;
    j = 0 ;
    while any(tmp_YX1y(:)-1 > 3*eps)
        j = j + 1 ;
        if j > 50
            error('Possible infinite loop in "Avoid, to high precision, sum(LU) > 1"')
        end
        base.maps_YXvy = base.maps_YXvy ./ repmat(tmp_YX1y,[1 1 Nlu 1]) ;
        tmp_YX1y = sum(base.maps_YXvy,3) ;
    end
    clear tmp_YX j
end

% Convert from "fraction of land" to "land area (m2)"
% (This also masks where needed due to harmonization of LUH2+PLUM masks)
if doHarm
    base.maps_YXv = base.maps_YXv .* landArea_YXv ;
else
    base.maps_YXvy = base.maps_YXvy .* repmat(landArea_YXv,[1 1 1 Nyears_luh2]) ;
end

% Import base_year crop fractions
% remaps_v4 has setAside and unhandledCrops in CC3G and CC4G. In the next
% step, we will combine CC3G and CC4G into ExtraCrop. (Was called SetAside,
% but renamed to avoid confusion, considering that this is not just
% setAside but also unhandledCrops.)
base_cropf = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/remaps_v4/cropfracs.remapv4.20180214.cgFertIrr0.setaside0103.m0.txt',...
    'verboseIfNoMat',false,'force_mat_save',true) ;
% Get just base year, if needed
if doHarm && isfield(base_cropf,'yearList')
    tmp = base_cropf.maps_YXvy(:,:,:,base_cropf.yearList==base_year) ;
    base_cropf = rmfield(base_cropf,{'maps_YXvy','yearList'}) ;
    base_cropf.maps_YXv = tmp ;
    clear tmp
end
% Combine CC3G and CC4G into ExtraCrop
if any(strcmp(base_cropf.varNames,'CC3G')) && any(strcmp(base_cropf.varNames,'CC4G'))
    base_cropf.maps_YXv(:,:,strcmp(base_cropf.varNames,'CC3G')) = ...
        base_cropf.maps_YXv(:,:,strcmp(base_cropf.varNames,'CC3G')) ...
        + base_cropf.maps_YXv(:,:,strcmp(base_cropf.varNames,'CC4G')) ;
    base_cropf.maps_YXv(:,:,strcmp(base_cropf.varNames,'CC4G')) = [] ;
    base_cropf.varNames(strcmp(base_cropf.varNames,'CC4G')) = [] ;
    base_cropf.varNames{strcmp(base_cropf.varNames,'CC3G')} = 'ExtraCrop' ;
end

% Avoid, to high precision, sum(base_cropf) > 1
if doHarm
    tmp_YX = sum(base_cropf.maps_YXv,3) ;
    j = 0 ;
    while any(tmp_YX(:)-1 > 3*eps)
        j = j + 1 ;
        if j > 50
            error('Possible infinite loop in "Avoid, to high precision, sum(base_cropf) > 1"')
        end
        base_cropf.maps_YXv = base_cropf.maps_YXv ./ repmat(tmp_YX,[1 1 size(base_cropf.maps_YXv,3)]) ;
        tmp_YX = sum(base_cropf.maps_YXv,3) ;
    end
    clear tmp_YX j
end

% Get previous crop types
tmp = base_cropf.varNames ;
for c = 1:length(base_cropf.varNames)
    thisCrop = base_cropf.varNames{c} ;
    thisCropI = [thisCrop 'i'] ;
    tmp(strcmp(tmp,thisCropI)) = [] ;
end
LPJGcrops = tmp ;
Ncrops_lpjg = length(LPJGcrops) ;
Nagri = Ncrops_lpjg + 1 ;
clear tmp


if ~doHarm 
    if ~isfield(base_cropf,'yearList')
        base_cropf.yearList = base.yearList ;
        base_cropf.maps_YXvy = repmat(base_cropf.maps_YXv,[1 1 1 Nyears_luh2]) ;
        base_cropf = rmfield(base_cropf,'maps_YXv') ;
    end
end

% Get "irrigation intensity" as fraction of thisCrop that is irrigated
% (later, will interpolate to cells without any of thisCrop)
base_irrig.varNames = LPJGcrops ;
if doHarm
    base_irrig.maps_YXv = nan([size(landArea_YX) Ncrops_lpjg],'double') ;
    for c = 1:Ncrops_lpjg
        thisCrop = LPJGcrops{c} ;
        if strcmp(thisCrop,'ExtraCrop')
            base_irrig.maps_YXv(:,:,c) = zeros(size(landArea_YX)) ;
            continue
        end
        thisCropI = [thisCrop 'i'] ;
        thisCropR_YX = base_cropf.maps_YXv(:,:,strcmp(base_cropf.varNames,thisCrop)) ;
        thisCropI_YX = base_cropf.maps_YXv(:,:,strcmp(base_cropf.varNames,thisCropI)) ;
        thisCrop_YX = thisCropR_YX + thisCropI_YX ;
        tmp_YX = thisCropI_YX ./ thisCrop_YX ;
        tmp_YX(thisCrop_YX==0) = NaN ;
        base_irrig.maps_YXv(:,:,c) = tmp_YX ;
    end
else
    base_irrig.maps_YXvy = nan([size(landArea_YX) Ncrops_lpjg Nyears_luh2],'double') ;
    for c = 1:Ncrops_lpjg
        thisCrop = LPJGcrops{c} ;
        if strcmp(thisCrop,'ExtraCrop')
            base_irrig.maps_YXvy(:,:,c,:) = zeros(size(base_irrig.maps_YXvy(:,:,c,:))) ;
            continue
        end
        thisCropI = [thisCrop 'i'] ;
        thisCropR_YX = base_cropf.maps_YXvy(:,:,strcmp(base_cropf.varNames,thisCrop),:) ;
        thisCropI_YX = base_cropf.maps_YXvy(:,:,strcmp(base_cropf.varNames,thisCropI),:) ;
        thisCrop_YX = thisCropR_YX + thisCropI_YX ;
        tmp_YX = thisCropI_YX ./ thisCrop_YX ;
        tmp_YX(thisCrop_YX==0) = NaN ;
        base_irrig.maps_YXvy(:,:,c,:) = tmp_YX ;
    end
end

% Combine irrigated and rainfed, if not already done
base_cropf = PLUMharm_combineIrrig(base_cropf, LPJGcrops) ;

% Read baseline fertilization
% (later, will interpolate to cells without any of thisCrop)
base_nfert = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/remaps_v4/nfert.remapv4.20180214.cgFertIrr0.setaside0103.m0.txt',...
    'verboseIfNoMat',false,'force_mat_save',true) ;
% Get just base year, if needed
if doHarm && isfield(base_nfert,'yearList')
    tmp = base_nfert.maps_YXvy(:,:,:,base_nfert.yearList==base_year) ;
    base_nfert = rmfield(base_nfert,{'maps_YXvy','yearList'}) ;
    base_nfert.maps_YXv = tmp ;
    clear tmp
elseif ~doHarm && ~isfield(base_nfert,'yearList')
    base_nfert.yearList = base.yearList ;
    base_nfert.maps_YXvy = repmat(base_nfert.maps_YXv,[1 1 1 Nyears_luh2]) ;
    base_nfert = rmfield(base_nfert,'maps_YXv') ;
end
% Add ExtraCrop, if needed
if any(strcmp(LPJGcrops,'ExtraCrop')) && ~any(strcmp(base_nfert.varNames,'ExtraCrop'))
    base_nfert.varNames = [base_nfert.varNames {'ExtraCrop'}] ;
    if doHarm
        base_nfert.maps_YXv = cat(3,base_nfert.maps_YXv,zeros(size(landArea_YX))) ;
    else
        base_nfert.maps_YXvy = cat(3,base_nfert.maps_YXvy,zeros([size(landArea_YX) 1 Nyears_luh2])) ;
    end
end
[~,IA,IB] = intersect(LPJGcrops,base_nfert.varNames,'stable') ;
if ~isequal(IA',1:Ncrops_lpjg)
    error('Issue with crop compatibility between LPJGcrops and baseline Nfert.')
end
if doHarm
    base_nfert.maps_YXv = base_nfert.maps_YXv(:,:,IB) ;
else
    base_nfert.maps_YXvy = base_nfert.maps_YXvy(:,:,IB,:) ;
end
base_nfert.varNames = LPJGcrops ;

% Convert base to have individual crops instead of CROPLAND
base_orig = base ;
base_cropa = base_cropf ;
if doHarm
    base_cropa.maps_YXv = base_cropf.maps_YXv .* base_orig.maps_YXv(:,:,strcmp(base_orig.varNames,'CROPLAND')) ;
    clear base
    base.varNames = [base_cropa.varNames base_orig.varNames(~strcmp(base_orig.varNames,'CROPLAND'))] ;
    base.maps_YXv = cat(3,base_cropa.maps_YXv,base_orig.maps_YXv(:,:,~strcmp(base_orig.varNames,'CROPLAND'))) ;
else
    base_cropa.maps_YXvy = base_cropf.maps_YXvy .* base_orig.maps_YXvy(:,:,strcmp(base_orig.varNames,'CROPLAND'),:) ;
    clear base
    base.varNames = [base_cropa.varNames base_orig.varNames(~strcmp(base_orig.varNames,'CROPLAND'))] ;
    base.maps_YXvy = cat(3,base_cropa.maps_YXvy,base_orig.maps_YXvy(:,:,~strcmp(base_orig.varNames,'CROPLAND'),:)) ;
end
clear base_orig base_cropf base_cropa

% Get info
LUnames = base.varNames ;
maxLength_LUnames = max(cellfun(@length,LUnames)) ;
Nlu = length(LUnames) ;
notBare = ~strcmp(LUnames,'BARREN') ;
isAgri = ~strcmp(LUnames,'NATURAL') & notBare ;
LUnames_agri = LUnames(isAgri) ;
isAgri_isPast = strcmp(LUnames_agri,'PASTURE') ;
isCrop = ~strcmp(LUnames,'NATURAL') & ~strcmp(LUnames,'PASTURE') & notBare ;

% Convert NaNs to zeros for addition
if doHarm
    base.maps_YXv(isnan(base.maps_YXv)) = 0 ;
    base_nfert.maps_YXv(base.maps_YXv(:,:,isCrop)==0) = 0 ;
    base_irrig.maps_YXv(base.maps_YXv(:,:,isCrop)==0) = 0 ;
else
    base.maps_YXvy(isnan(base.maps_YXvy)) = 0 ;
    base_nfert.maps_YXvy(base.maps_YXvy(:,:,isCrop,:)==0) = 0 ;
    base_irrig.maps_YXvy(base.maps_YXvy(:,:,isCrop,:)==0) = 0 ;
end

% Get baseline fraction that is vegetated, barren
if doHarm
    base_vegd_YX = sum(base.maps_YXv(:,:,notBare),3) ;
    base_bare_YX = base.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
else
    base_vegd_YX = sum(base.maps_YXvy(:,:,notBare,1),3) ;
    base_bare_YX = base.maps_YXvy(:,:,strcmp(LUnames,'BARREN'),1) ;
end
base_vegdFrac_YX = base_vegd_YX ./ landArea_YX ;
base_bareFrac_YX = base_bare_YX ./ landArea_YX ;

% Aggregate from 0.5º to 2º
if doHarm
    base_2deg.maps_YXv = PLUMharm_aggregate(base.maps_YXv,0.5,2) ;
    base_2deg_nfert.maps_YXv = PLUMharm_aggregate_mgmt(base_nfert.maps_YXv,base.maps_YXv(:,:,isCrop),0.5,2) ;
    base_2deg_irrig.maps_YXv = PLUMharm_aggregate_mgmt(base_irrig.maps_YXv,base.maps_YXv(:,:,isCrop),0.5,2) ;
else
    base_2deg.maps_YXvy = PLUMharm_aggregate(base.maps_YXvy,0.5,2) ;
    base_2deg_nfert.maps_YXvy = PLUMharm_aggregate_mgmt(base_nfert.maps_YXvy,base.maps_YXvy(:,:,isCrop,:),0.5,2) ;
    base_2deg_irrig.maps_YXvy = PLUMharm_aggregate_mgmt(base_irrig.maps_YXvy,base.maps_YXvy(:,:,isCrop,:),0.5,2) ;
end
base_2deg.varNames = base.varNames ;
base_2deg_nfert.varNames = LPJGcrops ;
base_2deg_irrig.varNames = LPJGcrops ;
landArea_2deg_YX = PLUMharm_aggregate(landArea_YX,0.5,2) ;

% Do not allow invalid management inputs.
if doHarm
    if any(base_2deg_nfert.maps_YXv(:)<0)
        error('Negative value(s) in base_2deg_nfert.maps_YXv!')
    elseif any(base_2deg_irrig.maps_YXv(:)<0)
        error('Negative value(s) in base_2deg_irrig.maps_YXv!')
    elseif any(base_2deg_irrig.maps_YXv(:)>1+1e-6)
        error('Value(s) >1 in base_2deg_irrig.maps_YXv!')
    end
else
    if any(base_2deg_nfert.maps_YXvy(:)<0)
        error('Negative value(s) in base_2deg_nfert.maps_YXvy!')
    elseif any(base_2deg_irrig.maps_YXvy(:)<0)
        error('Negative value(s) in base_2deg_irrig.maps_YXvy!')
    elseif any(base_2deg_irrig.maps_YXvy(:)>1+1e-6)
        error('Value(s) >1 in base_2deg_irrig.maps_YXvy!')
    end
end

% Interpolate management inputs
if doHarm
    base_2deg_nfert.maps_YXv = PLUMharm_interpolateMgmt(...
        base_2deg_nfert.maps_YXv, base.maps_YXv(:,:,isCrop), landArea_2deg_YX,...
        LPJGcrops, inpaint_method) ;
    base_2deg_irrig.maps_YXv = PLUMharm_interpolateMgmt(...
        base_2deg_irrig.maps_YXv, base.maps_YXv(:,:,isCrop), landArea_2deg_YX,...
        LPJGcrops, inpaint_method) ;
end

% Get LUH2 2-deg fraction that is vegetated, barren
if doHarm
    base_2deg_vegd_YX = sum(base_2deg.maps_YXv(:,:,notBare),3) ;
    base_2deg_bare_YX = base_2deg.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
else
    base_2deg_vegd_YX = sum(base_2deg.maps_YXvy(:,:,notBare,1),3) ;
    base_2deg_bare_YX = base_2deg.maps_YXvy(:,:,strcmp(LUnames,'BARREN'),1) ;
end
base_2deg_vegdFrac_YX = base_2deg_vegd_YX ./ landArea_2deg_YX ;
base_2deg_bareFrac_YX = base_2deg_bare_YX ./ landArea_2deg_YX ;

% Get other info
if doHarm
    base_nfertTot_2deg.maps_YXv = base_2deg_nfert.maps_YXv .* base_2deg.maps_YXv(:,:,isCrop) ;
    base_irrigTot_2deg.maps_YXv = base_2deg_irrig.maps_YXv .* base_2deg.maps_YXv(:,:,isCrop) ;
else
    base_nfertTot.maps_YXvy = base_nfert.maps_YXvy .* base.maps_YXvy(:,:,isCrop,:) ;
    base_nfertTot_2deg.maps_YXvy = base_2deg_nfert.maps_YXvy .* base_2deg.maps_YXvy(:,:,isCrop,:) ;
    base_irrigTot.maps_YXvy = base_irrig.maps_YXvy .* base.maps_YXvy(:,:,isCrop,:) ;
    base_irrigTot_2deg.maps_YXvy = base_2deg_irrig.maps_YXvy .* base_2deg.maps_YXvy(:,:,isCrop,:) ;
end

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