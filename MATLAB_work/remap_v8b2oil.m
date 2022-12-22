%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Re-map area/fert data to PLUM crops, and generate extra LU file %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Sam Rabin, 2021-02-18
% For getting calibration factors.
% See previous notes for remap_v8b, with the following changes:
% - Lumps FruitAndVeg and Sugar in with Oilcrops

warning('on','all')

thisVer = 'WithFruitVegSugar_b_2oil' ;
remapVer = '8b2oil' ;
out_dir = sprintf('/Users/Shared/PLUM/input/remaps_v%s/',remapVer) ;


%% Part 1: Import data

disp('Importing...')

PLUMsetAside_frac = 0.103 ;
inpaint_method = 4 ;
force_all_rainfed = false ;

addpath(genpath(landsymm_lpjg_path()))

% Land use
lu_in = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt');%,...
%     'verboseIfNoMat',true) ;

% Fertilizer (convert kg/ha to kg/m2)
nfert_in = lpjgu_matlab_readTable_then2map('/Users/Shared/unifying_gridlist/AgGRID_nutrient_input_v1.1/AgMIP_Nfert.txt');%,...
%     'verboseIfNoMat',true) ;
nfert_in.maps_YXv = nfert_in.maps_YXv*1e-4 ;
nfert_in.maps_YXv = nfert_in.maps_YXv * 1/(1-PLUMsetAside_frac) ;   % New in remap_v5
list_cropsCombined_fert_in = nfert_in.varNames ;
NcropsCombined_fert_in = length(list_cropsCombined_fert_in) ;

% Crop fractions
croparea_in = lpjgu_matlab_readTable_then2map('/Users/sam/Geodata/MIRCA/harvested_area_grids_26crops_30mn/MIRCA.txt',...
    'verboseIfNoMat',true) ;

% Align gridlists
gridlist = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/gridlists/gridlist_62892.runAEclimOK.txt');%, 'verboseIfNoMat',true) ;
isnan_LUH2 = isnan(lu_in.maps_YXvy(:,:,1,1)) ;
if any(any(gridlist.mask_YX & isnan_LUH2))
    error('Somehow you have cells included in your gridlist that are not in the LU data.')
end
lu_in.maps_YXvy(~repmat(gridlist.mask_YX,[1 1 size(lu_in.maps_YXvy,3) size(lu_in.maps_YXvy,4)])) = NaN ;
lu_in.list_to_map = gridlist.list_to_map ;

% Process crop fractions
list_crops_frac_in = croparea_in.varNames ;
list_cropsCombined_frac_in = unique(strrep(strrep(list_crops_frac_in,'_IR',''),'_RF','')) ;
Ncrops_frac_in = length(list_crops_frac_in) ;
NcropsCombined_frac_in = length(list_cropsCombined_frac_in) ;
% Force all irrigated to rainfed, if doing so
if force_all_rainfed
    warning('GETTING RID OF IRRIGATED')
    for c = 1:NcropsCombined_frac_in
        thisCrop = list_cropsCombined_frac_in{c} ;
        thisIR = find(strcmp(croparea_in.varNames,[thisCrop '_IR'])) ;
        thisRF = find(strcmp(croparea_in.varNames,[thisCrop '_RF'])) ;
        croparea_in.maps_YXv(:,:,thisRF) = croparea_in.maps_YXv(:,:,thisRF) + croparea_in.maps_YXv(:,:,thisIR) ;
        croparea_in.maps_YXv(:,:,thisIR) = 0*croparea_in.maps_YXv(:,:,thisIR) ;
    end
end
% Calculate total area of each crop
cropareaCombined_in.maps_YXv = nan(size(croparea_in.maps_YXv,1),size(croparea_in.maps_YXv,2),NcropsCombined_frac_in) ;
for c = 1:NcropsCombined_frac_in
    thisCrop = list_cropsCombined_frac_in{c} ;
    isThisCrop = not(cellfun(@isempty,strfind(croparea_in.varNames,thisCrop))) ;
    if length(find(isThisCrop))~=2
        error('length(find(isThisCrop))~=2')
    end
    cropareaCombined_in.maps_YXv(:,:,c) = sum(croparea_in.maps_YXv(:,:,isThisCrop),3) ;
end
cropareaCombined_in.varNames = list_cropsCombined_frac_in ;

% Add manure N for year 2000 (new in remap_v5)
load('/Users/sam/Geodata/Manure_ZhangEtAl2017/zhangManure_1860to2014_agg_hd.mat') ;
manure2crop_hd_YX = manure2crop_hd_YXy(:,:,1860:2014==2000) ;
manure2crop_hd_YX(manure2crop_hd_YX<1e-6) = 0 ;
clear manure2crop_hd_YXy
if any(any(sum(cropareaCombined_in.maps_YXv,3)>0 & isnan(manure2crop_hd_YX)))
    manure2crop_hd2_YX = inpaint_nans(manure2crop_hd_YX,inpaint_method) ;
end
manure2crop_hd2_YX(isnan_LUH2) = NaN ;
nfert_in.maps_YXv = nfert_in.maps_YXv + repmat(manure2crop_hd_YX,[1 1 size(nfert_in.maps_YXv,3)]) ;

disp('Done.')

%% Part 2: Set up mapping

% Get maps
[list_cropsCombined_out, in2out_keyCombined_frac, ...
    list_ignore_frac, allVer_names, ~, ~, allVer_ignore_types] = ...
    get_remapv2_keys(thisVer) ;
NcropsCombined_out = length(list_cropsCombined_out) ;
list_crops_out = [list_cropsCombined_out strcat(list_cropsCombined_out,'i')] ;
Ncrops_out = length(list_crops_out) ;

% Make sure there are no repeats
if length(list_cropsCombined_out) ~= length(unique(list_cropsCombined_out))
    error('length(list_cropsCombined_out) ~= length(unique(list_cropsCombined_out))')
end

% Set up anonymous functions
getOi = @(x) find(strcmp(list_crops_out,x)) ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;

% Check that every crop in croparea_in is either explicitly included or ignored
for c = 1:NcropsCombined_frac_in
    thisCrop = list_cropsCombined_frac_in{c} ;
    found = 0 ;
    if any(strcmp(list_ignore_frac,thisCrop))
        found = 1 ;
    end
    for cc = 1:length(in2out_keyCombined_frac)
        found = found + length(find(strcmp(in2out_keyCombined_frac{cc},thisCrop))) ;
    end
    if found==0
        error([thisCrop ' not found in list_ignore_frac or in2out_keyCombined_frac!'])
    elseif found>1
        error([thisCrop ' found ' num2str(found) ' times in list_ignore_frac and in2out_keyCombined_frac!'])
    end
end

% Make sure nothing is left empty
if any(cellfun(@isempty,in2out_keyCombined_frac))
    error('At least one member of in2out_keyCombined_frac is empty!')
end

% Make rainfed-vs-irrigated map: frac
in2out_key_frac = {} ;
for c = 1:NcropsCombined_out
    thisCrop = list_cropsCombined_frac_in{c} ;
    in2out_key_frac{c} = strcat(in2out_keyCombined_frac{c},'_RF') ;
end
for c = 1:NcropsCombined_out
    thisCrop = list_cropsCombined_frac_in{c} ;
    in2out_key_frac{NcropsCombined_out+c} = strcat(in2out_keyCombined_frac{c},'_IR') ;
end

% Make sure that nfert croplist are all contained within frac croplist
% (not strictly necessary, but because of how following code is written)
C1 = intersect(list_cropsCombined_frac_in,list_cropsCombined_fert_in) ;
if ~isequal(sort(C1),sort(list_cropsCombined_fert_in))
    error('Nfert crops (list_cropsCombined_fert_in) are not all contained within list_cropsCombined_frac_in!')
end
clear C1

% Now get key for fert
in2out_keyCombined_fert = in2out_keyCombined_frac ;
list_crops_inFrac_notFert = setdiff(list_cropsCombined_frac_in,list_cropsCombined_fert_in) ;
luh2_file_mgmt = '/Users/sam/Geodata/LUH2/v2h/management.1850-2015.nc' ;
luh2_file_states = '/Users/sam/Geodata/LUH2/v2h/states.1850-2015.nc' ;
luh2_file_etc = '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc' ;
for c = 1:NcropsCombined_out
    thisRow = in2out_keyCombined_fert{c} ;
    thisCrop_out = list_cropsCombined_out{c} ;
    if isequal(thisRow,{'Pulses'})
        warning('Using Nfert from groundnuts+soybeans for pulses.')
        thisRow = {'GroundnutsPeanuts','Soybeans'} ;
    elseif strcmp(thisCrop_out,'DateCitGrape') || strcmp(thisCrop_out,'FruitAndVeg')
        % YEAR needs to be actual year minus 850, because units are "years
        % since 1850-01-01"
        warning('Using Nfert from LUH2 c3per for %s.', thisCrop_out)
        tmp_nfert_YX = netcdf_read_selectYears(...
            luh2_file_mgmt, 'fertl_c3per', 2000-850) ;
        tmp_carea_YX = netcdf_read_selectYears(...
            luh2_file_states, 'c3per', 2000-850) ...
            .* flip(permute(ncread(luh2_file_etc,'carea'),[2 1]),1) ;
        tmp_nfert_YX = aggregate_mgmt(tmp_nfert_YX, tmp_carea_YX, 0.25, 0.5) ;
        tmp_carea_YX = aggregate(tmp_carea_YX, 0.25, 0.5) ;
        tmp_nfert_YX(tmp_carea_YX==0) = NaN ;
        nfert_in.maps_YXv(:,:,end+1) = tmp_nfert_YX ;
        nfert_in.varNames{end+1} = thisCrop_out ;
        list_cropsCombined_fert_in{end+1} = thisCrop_out ;
        NcropsCombined_fert_in = length(list_cropsCombined_fert_in) ;
        thisRow = {thisCrop_out} ;
    else
        [C,IA] = intersect(thisRow,list_crops_inFrac_notFert) ;
        thisRow(IA) = [] ;
    end
    if isempty(thisRow)
        error('No fertilizer found for this keymap!')
    end
    in2out_keyCombined_fert{c} = thisRow ;
end

% Get ignored area
[~,I] = intersect(croparea_in.varNames,[strcat(list_ignore_frac,'_RF');strcat(list_ignore_frac,'_IR')]) ;
map_ignore_area_YX = sum(croparea_in.maps_YXv(:,:,I),3) ;


%% Part 3: Do mapping

% Crop fractions
croparea_mid.maps_YXv = nan(360,720,Ncrops_out) ;
for c = 1:Ncrops_out
    thisCrop = list_crops_out{c} ;
    thisRow = in2out_key_frac{getOi(thisCrop)} ;
    [~,IA] = intersect(croparea_in.varNames,thisRow) ;
    if length(IA)~=length(thisRow)
        error('length(IA)~=length(thisRow)')
    end
    % Straight sum
    croparea_mid.maps_YXv(:,:,c) = sum(croparea_in.maps_YXv(:,:,IA),3) ;
end
croparea_mid.varNames = list_crops_out ;

% Fertilizer
nfert_mid.maps_YXv = nan(360,720,NcropsCombined_out) ;
%%% Rainfed
for c = 1:NcropsCombined_out
    thisCrop = list_crops_out{c} ;
    if strcmp(thisCrop,'DateCitGrape') || strcmp(thisCrop,'FruitAndVeg')
        nfert_mid.maps_YXv(:,:,c) = nfert_in.maps_YXv(:,:,strcmp(nfert_in.varNames,thisCrop)) ;
    else
        thisRow = in2out_keyCombined_fert{getOi(thisCrop)} ;
        [~,IA_frac] = intersect(cropareaCombined_in.varNames,thisRow) ;
        [~,IA_fert] = intersect(nfert_in.varNames,thisRow) ;
        if length(IA_frac)~=length(thisRow)
            error('length(IA_frac)~=length(thisRow)')
        end
        if length(IA_fert)~=length(thisRow)
            error('length(IA_fert)~=length(thisRow)')
        end
        % Weighted average
        nfert_mid.maps_YXv(:,:,c) = sum(nfert_in.maps_YXv(:,:,IA_fert) .* cropareaCombined_in.maps_YXv(:,:,IA_frac) ./ repmat(sum(cropareaCombined_in.maps_YXv(:,:,IA_frac),3),[1 1 length(IA_frac)]),3) ;
    end
end
clear nfert_in
%%% Irrigated
nfert_mid.maps_YXv = cat(3,nfert_mid.maps_YXv,nfert_mid.maps_YXv) ;
nfert_mid.varNames = list_crops_out ;

% Add ExtraCrop to represent unhandled crops and setAside
% Move ignored area into ExtraCrop.
% Set up croparea_extra array
nfert0_YX = zeros(size(map_ignore_area_YX)) ;
% Assign based on pastNPP results
setaside_area_YX = sum(croparea_mid.maps_YXv * PLUMsetAside_frac,3) ;
croparea_mid.maps_YXv = croparea_mid.maps_YXv * (1-PLUMsetAside_frac) ;
croparea_extra_YX = map_ignore_area_YX + setaside_area_YX ;
% Mask pastNPP NaNs
croparea_extra_YX(isnan_LUH2) = NaN ;
nfert0_YX(isnan_LUH2) = NaN ;
% Save new croparea and nfert structures
list_cropsCombined_out = [list_cropsCombined_out {'ExtraCrop'}] ;
NcropsCombined_out = length(list_cropsCombined_out) ;
list_crops_out = [list_crops_out {'ExtraCrop'}]  ;
Ncrops_out = length(list_crops_out) ;
croparea_mid.maps_YXv = cat(3,croparea_mid.maps_YXv,croparea_extra_YX) ;
croparea_mid.varNames = list_crops_out ;
nfert_mid.maps_YXv = cat(3,nfert_mid.maps_YXv,nfert0_YX) ;
nfert_mid.varNames = list_crops_out ;

% Get fractions
cropfrac_mid.maps_YXv = croparea_mid.maps_YXv ./ repmat(sum(croparea_mid.maps_YXv,3),[1 1 Ncrops_out]) ;
% % % map_ignore_frac_YX = map_ignore_area_YX ./ (map_ignore_area_YX + sum(croparea_mid.maps_YXv,3)) ;
% % % map_include_frac_YX = 1 - map_ignore_frac_YX ;


%% Part 4: Interpolate

% Set up output arrays
cropfrac_out.varNames = list_crops_out ;
nfert_out.varNames = list_crops_out ;
cropfrac_out.maps_YXv = nan(360,720,Ncrops_out) ;
nfert_out.maps_YXv = nan(360,720,Ncrops_out) ;

% Interpolate maps
for c = 1:Ncrops_out
    thisCrop = list_crops_out{c} ;
    disp(['Interpolating ' thisCrop ' (' num2str(c) ' of ' num2str(Ncrops_out) ')...'])
    cropfrac_out.maps_YXv(:,:,c) = inpaint_nans(cropfrac_mid.maps_YXv(:,:,c), inpaint_method) ;
    if any(any(isnan(cropfrac_out.maps_YXv(:,:,c))))
        error('NaN remaining in cropfrac_out.maps_YXv(:,:,c)!')
%     elseif any(any(cropfrac_out.maps_YXv(:,:,c)<0))
%         error('Negative in cropfrac_out.maps_YXv(:,:,c)!')
    end
    thisNonIrr = find(strcmp(list_crops_out,thisCrop(1:end-1))) ;
    if length(thisNonIrr)==1 && strcmp([list_crops_out{thisNonIrr} 'i'],thisCrop)
        warning(['Nfert: Using ' list_crops_out{thisNonIrr} ' for ' list_crops_out{c}])
        nfert_out.maps_YXv(:,:,c) = nfert_out.maps_YXv(:,:,thisNonIrr) ;
    else
        nfert_out.maps_YXv(:,:,c) = inpaint_nans(nfert_mid.maps_YXv(:,:,c), inpaint_method) ;
    end
    if any(any(isnan(nfert_out.maps_YXv(:,:,c))))
        error('NaN remaining in nfert_out.maps_YXv(:,:,c)!')
%     elseif any(any(nfert_out.maps_YXv(:,:,c)<0))
%         error('Negative in nfert_out.maps_YXv(:,:,c)!')
    end
end

% % % disp('Interpolating include_frac...')
% % % map_include_frac_interpd_YX = inpaint_nans(map_include_frac_YX,inpaint_method) ;
% % % disp('Done.')

% Restrict to LUH2 gridlist
cropfrac_out.maps_YXv(repmat(isnan_LUH2,[1 1 Ncrops_out])) = NaN ;
nfert_out.maps_YXv(repmat(isnan_LUH2,[1 1 Ncrops_out])) = NaN ;
% % % map_include_frac_interpd_YX(isnan_LUH2) = NaN ;

% Set minimum of zero, if needed
if any(cropfrac_out.maps_YXv(:)<0)
    warning('Setting negative members of cropfrac_out.maps_YXv to zero.')
    cropfrac_out.maps_YXv(cropfrac_out.maps_YXv<0) = 0 ;
end
if any(nfert_out.maps_YXv(:)<0)
    warning('Setting negative members of nfert_out.maps_YXv to zero.')
    nfert_out.maps_YXv(nfert_out.maps_YXv<0) = 0 ;
end
% % % if any(map_include_frac_interpd_YX(:)<0)
% % %     warning('Setting negative members of map_include_frac_interpd_YX to zero.')
% % %     map_include_frac_interpd_YX(map_include_frac_interpd_YX<0) = 0 ;
% % % end


% Normalize cropfrac to 1, if needed
tmp = sum(cropfrac_out.maps_YXv,3) ;
if any(abs(tmp(:) - 1)>1e-6)
    cropfrac_out.maps_YXv = cropfrac_out.maps_YXv ./ repmat(tmp,[1 1 Ncrops_out]) ;
end
% % % if any(abs(map_include_frac_interpd_YX(:) - 1)>1e-6)
% % %     map_include_frac_interpd_YX(map_include_frac_interpd_YX>1) = 1 ;
% % % end

% Don't let extrapolated Nfert exceed maximum Nfert seen for this crop
for c = 1:NcropsCombined_out-1 % -1 to skip ExtraCrop
    tmp = nfert_mid.maps_YXv(:,:,c) ;
    thisMax = max(tmp(~isnan(tmp))) ;
    thisMap = nfert_out.maps_YXv(:,:,c) ;
    Nexceeded = length(find(thisMap>thisMax)) ;
    if Nexceeded>0
        warning(['Limiting nfert_out for ' list_cropsCombined_out{c} '.'])
        thisMap(thisMap>thisMax) = thisMax ;
        nfert_out.maps_YXv(:,:,c) = thisMap ;
        nfert_out.maps_YXv(:,:,c+NcropsCombined_out) = thisMap ;
    end
end
clear nfert_mid

% Add zeros for Miscanthus
list_cropsCombined_out = [list_cropsCombined_out {'Miscanthus','Miscanthusi'}] ;
NcropsCombined_out = length(list_cropsCombined_out) ;
list_crops_out = [list_crops_out {'Miscanthus','Miscanthusi'}]  ;
Ncrops_out = length(list_crops_out) ;
if force_all_rainfed
    error('Edit code to correctly add zeros for Miscanthus when force_all_rainfed.')
end
cropfrac_out.maps_YXv = cat(3,cropfrac_out.maps_YXv,zeros(size(cropfrac_out.maps_YXv(:,:,1:2)))) ;
nfert_out.maps_YXv = cat(3,nfert_out.maps_YXv,zeros(size(nfert_out.maps_YXv(:,:,1:2)))) ;
cropfrac_out.varNames = list_crops_out ;
nfert_out.varNames = list_crops_out ;
disp('Done interpolating.')


%% Part 5: Set up for save

disp('Array-ifying lu_in...')
[lu_in_array, lu_in_header_cell] = lpjgu_matlab_maps2table(lu_in,gridlist.list_to_map) ;
clear lu_in
disp('Array-ifying cropfrac_out...')
[cropfrac_out_array, cropfrac_out_header_cell] = lpjgu_matlab_maps2table(cropfrac_out,gridlist.list_to_map) ;
clear cropfrac_out
disp('Array-ifying nfert_out...')
[nfert_out_array, nfert_out_header_cell] = lpjgu_matlab_maps2table(nfert_out,gridlist.list_to_map) ;
clear nfert_out
disp('Done.')


%% Part 6: Save

%%% Options %%%%%%%%%%%%%
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(out_dir,'dir')
    mkdir(out_dir) ;
end

% Get filenames
if force_all_rainfed
    allRF_txt = '.noIrr' ;
else
    allRF_txt = '' ;
end
out_file_lu = [out_dir 'LU.remapv' remapVer '.txt'] ;
out_file_cropfrac = [out_dir 'cropfracs.remapv' remapVer '.txt'] ;
out_file_nfert = [out_dir 'nfert.remapv' remapVer '.txt'] ;

% % disp('Saving ignored crop fraction...')
% % out_file_cropfrac_ignored = strrep(out_file_cropfrac,'cropfracs','cropfracsIGNORED') ;
% % out_file_cropfrac_ignored = strrep(out_file_cropfrac_ignored,'.txt','.mat') ;
% % save(out_file_cropfrac_ignored,'ignored_LUCROParea_YX_y','-v7.3') ;

disp('Saving LU...')
check_existing_lu(thisVer, out_file_lu, allVer_names, allVer_ignore_types) ;
lpjgu_matlab_saveTable(lu_in_header_cell, lu_in_array, out_file_lu,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 1) ;

disp('Saving cropfracs...')
lpjgu_matlab_saveTable(cropfrac_out_header_cell, cropfrac_out_array, out_file_cropfrac,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20) ;

disp('Saving nfert...')
lpjgu_matlab_saveTable(nfert_out_header_cell, nfert_out_array, out_file_nfert,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20) ;


