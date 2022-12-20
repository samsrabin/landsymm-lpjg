%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Re-map area/fert data to PLUM crops, and generate extra LU file %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Sam Rabin, 2018-07-27
% This script is conceptually based on Tom Pugh's landuse_preproc_v6.m
% script.
% 
% Works the same as remap_v2, but instead of moving unhandled crop area to
% pasture, remap_v3 puts unhandled crop area into CC3G or CC4G crops. The
% partitioning depends on which has higher NPP in 2000±10 according to a 
% previous pasture-only run. Not time-evolving because other crop-specific
% outputs (cropfrac, nfert) are set at 2000 levels.
%
% If tied, assign based on latitude. Tropics (within ±22.5º of Equator)
% get C4; other areas get C3. Assigning NaN and letting interpolation
% figure it out doesn't work because you need to get from area to fraction
% before interpolation happens.
%
% Also adds zeros for Miscanthus.

% Version for crop mappings
thisVer = '20180214' ;

% Directory for pasture NPP run
pastureNPPdir = 'LPJGPLUM_1850-2010_pastureOnly_cmip5ipsl/output-2018-07-26-172215' ;

force_all_rainfed = false ;
inpaint_method = 0 ;

out_dir = '/Users/Shared/PLUM/input/remaps_v3/' ;


%% Part 1: Import data

disp('Importing...')

% Get pasture NPP file
pastureNPPdir = find_PLUM2LPJG_run(pastureNPPdir) ;
pastNPP = lpjgu_matlab_readTable_then2map([pastureNPPdir 'anpp.out'], ...
    'verboseIfNoMat',true,'force_mat_save',true) ;

% Land use
lu_in = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt');%,...
%     'verboseIfNoMat',true) ;

% Fertilizer (convert kg/ha to kg/m2)
nfert_in = lpjgu_matlab_readTable_then2map('/Users/Shared/unifying_gridlist/AgGRID_nutrient_input_v1.1/AgMIP_Nfert.txt');%,...
%     'verboseIfNoMat',true) ;
nfert_in.maps_YXv = nfert_in.maps_YXv*1e-4 ;
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
for c = 1:NcropsCombined_out
    thisRow = in2out_keyCombined_fert{c} ;
    if isequal(thisRow,{'Pulses'})
        warning('Using Nfert from groundnuts+soybeans for pulses.')
        thisRow = {'GroundnutsPeanuts','Soybeans'} ;
    else
        [~,IA] = intersect(thisRow,list_crops_inFrac_notFert) ;
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
%%% Irrigated
nfert_mid.maps_YXv = cat(3,nfert_mid.maps_YXv,nfert_mid.maps_YXv) ;
nfert_mid.varNames = list_crops_out ;

% Add CC3G and CC4G to represent unhandled crops
% Partition ignored area into CC3G and CC4G.
%%% No need to reconcile NaN masks because interpolation happens anyway
% Get pastNPP results
pastNPP_baseYear = 2000 ;
pastNPP_windowRadius = 10 ;
inThisPeriod = (pastNPP.yearList >= pastNPP_baseYear-pastNPP_windowRadius) ...
             & (pastNPP.yearList<=pastNPP_baseYear+pastNPP_windowRadius) ;
pastNPP_PC3G_YX = mean(pastNPP.maps_YXvy(:,:,strcmp(pastNPP.varNames,'PC3G'),inThisPeriod),4) ;
pastNPP_PC4G_YX = mean(pastNPP.maps_YXvy(:,:,strcmp(pastNPP.varNames,'PC4G'),inThisPeriod),4) ;
assignPC3G_YX = (pastNPP_PC3G_YX > pastNPP_PC4G_YX) ;
assignPC4G_YX = (pastNPP_PC3G_YX < pastNPP_PC4G_YX) ;
assignNone_YX = ~assignPC3G_YX & ~assignPC4G_YX ;
isnan_pastNPP_YX = isnan(pastNPP_PC3G_YX) ;
% Set default (for handling ties) based on latitude
yres = 180/size(pastNPP_PC3G_YX,1) ;
lats = ((-90+yres/2):yres:(90-yres/2))' ;
lats_YX = repmat(lats,[1 size(pastNPP_PC3G_YX,2)]) ;
istropics_YX = abs(lats_YX) <= 23.5 ;
assignPC3G_YX(~istropics_YX & assignNone_YX) = true ;
assignPC4G_YX( istropics_YX & assignNone_YX) = true ;
% Set up croparea_cc*g arrays
croparea_cc3g_YX = zeros(size(map_ignore_area_YX)) ;
croparea_cc4g_YX = zeros(size(map_ignore_area_YX)) ;
nfert0_YX = zeros(size(map_ignore_area_YX)) ;
% Assign based on pastNPP results
croparea_cc3g_YX(assignPC3G_YX) = map_ignore_area_YX(assignPC3G_YX) ;
croparea_cc4g_YX(assignPC4G_YX) = map_ignore_area_YX(assignPC4G_YX) ;
% Mask pastNPP NaNs
croparea_cc3g_YX(isnan_pastNPP_YX) = NaN ;
croparea_cc4g_YX(isnan_pastNPP_YX) = NaN ;
nfert0_YX(isnan_pastNPP_YX) = NaN ;
% Save new croparea and nfert structures
list_cropsCombined_out = [list_cropsCombined_out {'CC3G','CC4G'}] ;
NcropsCombined_out = length(list_cropsCombined_out) ;
list_crops_out = [list_crops_out {'CC3G','CC4G'}]  ;
Ncrops_out = length(list_crops_out) ;
croparea_mid.maps_YXv = cat(3,croparea_mid.maps_YXv,croparea_cc3g_YX,croparea_cc4g_YX) ;
croparea_mid.varNames = list_crops_out ;
nfert_mid.maps_YXv = cat(3,nfert_mid.maps_YXv,nfert0_YX,nfert0_YX) ;
nfert_mid.varNames = list_crops_out ;

% Get fractions
cropfrac_mid.maps_YXv = croparea_mid.maps_YXv ./ repmat(sum(croparea_mid.maps_YXv,3),[1 1 Ncrops_out]) ;
% % % map_ignore_frac_YX = map_ignore_area_YX ./ (map_ignore_area_YX + sum(croparea_mid.maps_YXv,3)) ;
% % % map_include_frac_YX = 1 - map_ignore_frac_YX ;


%% Part 4: Interpolate

cropfrac_out.maps_YXv = nan(360,720,Ncrops_out) ;
nfert_out.maps_YXv = nan(360,720,Ncrops_out) ;
for c = 1:Ncrops_out
    thisCrop = list_crops_out{c} ;
    disp(['Interpolating ' thisCrop ' (' num2str(c) ' of ' num2str(Ncrops_out) ')...'])
    cropfrac_out.maps_YXv(:,:,c) = inpaint_nans(cropfrac_mid.maps_YXv(:,:,c), inpaint_method) ;
    if c > NcropsCombined_out
        warning(['Using ' list_crops_out{c-NcropsCombined_out} ' for ' list_crops_out{c}])
        nfert_out.maps_YXv(:,:,c) = nfert_out.maps_YXv(:,:,c-NcropsCombined_out) ;
    else
        nfert_out.maps_YXv(:,:,c) = inpaint_nans(nfert_mid.maps_YXv(:,:,c), inpaint_method) ;
    end
end

cropfrac_out.varNames = list_crops_out ;
nfert_out.varNames = list_crops_out ;
% % % disp('Interpolating include_frac...')
% % % map_include_frac_interpd_YX = inpaint_nans(map_include_frac_YX,inpaint_method) ;
% % % disp('Done.')

% Restrict to LUH2 gridlist
cropfrac_out.maps_YXv(isnan(repmat(isnan_LUH2,[1 1 Ncrops_out]))) = NaN ;
nfert_out.maps_YXv(isnan(repmat(isnan_LUH2,[1 1 Ncrops_out]))) = NaN ;
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
for c = 1:NcropsCombined_out-2 % -2 to skip CC3G and CC4G
    tmp = cropfrac_mid.maps_YXv(:,:,c) ;
    thisMax = max(tmp(~isnan(tmp))) ;
    thisMap = cropfrac_out.maps_YXv(:,:,c) ;
    Nexceeded = length(find(thisMap>thisMax)) ;
    if Nexceeded>0
        warning(['Limiting nfert_out for ' list_cropsCombined_out{c} '.'])
        thisMap(thisMap>thisMax) = thisMax ;
        nfert_out.maps_YXv(:,:,c) = thisMap ;
        nfert_out.maps_YXv(:,:,c+NcropsCombined_out) = thisMap ;
    end
end

% Do not trust CC3G and CC4G anymore; reclassify into one or the other
ic3g = strcmp(cropfrac_out.varNames,'CC3G') ;
ic4g = strcmp(cropfrac_out.varNames,'CC4G') ;
cropfrac_out_cc3g_YX = cropfrac_out.maps_YXv(:,:,ic3g) ;
cropfrac_out_cc4g_YX = cropfrac_out.maps_YXv(:,:,ic4g) ;
cropfrac_out_ccXg_YX = cropfrac_out_cc3g_YX + cropfrac_out_cc4g_YX ;
cropfrac_out_cc3g_YX(assignPC3G_YX) = cropfrac_out_ccXg_YX(assignPC3G_YX) ;
cropfrac_out_cc3g_YX(~assignPC3G_YX) = 0 ;
cropfrac_out_cc4g_YX(assignPC4G_YX) = cropfrac_out_ccXg_YX(assignPC4G_YX) ;
cropfrac_out_cc4g_YX(~assignPC4G_YX) = 0 ;
cropfrac_out.maps_YXv(:,:,ic3g) = cropfrac_out_cc3g_YX ;
cropfrac_out.maps_YXv(:,:,ic4g) = cropfrac_out_cc4g_YX ;

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
disp('Array-ifying cropfrac_out...')
[cropfrac_out_array, cropfrac_out_header_cell] = lpjgu_matlab_maps2table(cropfrac_out,gridlist.list_to_map) ;
disp('Array-ifying nfert_out...')
[nfert_out_array, nfert_out_header_cell] = lpjgu_matlab_maps2table(nfert_out,gridlist.list_to_map) ;
disp('Done.')

% Get filenames
if force_all_rainfed
    allRF_txt = '.noIrr' ;
else
    allRF_txt = '' ;
end
out_file_lu = [out_dir 'LU_xtraCROPtoPAST.remapv3.' thisVer '.cgFertIrr0.m' num2str(inpaint_method) '.txt'] ;
out_file_cropfrac = [out_dir 'cropfracs.remapv3.' thisVer allRF_txt '.cgFertIrr0.m' num2str(inpaint_method) '.txt'] ;
out_file_nfert = [out_dir 'nfert.remapv3.' thisVer allRF_txt '.cgFertIrr0.m' num2str(inpaint_method) '.txt'] ;


%% Part 7: Save

%%% Options %%%%%%%%%%%%%
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
%%%%%%%%%%%%%%%%%%%%%%%%%

% % % disp('Saving ignored crop fraction...')
% % % out_file_cropfrac_ignored = strrep(out_file_cropfrac,'cropfracs','cropfracsIGNORED') ;
% % % out_file_cropfrac_ignored = strrep(out_file_cropfrac_ignored,'.txt','.mat') ;
% % % save(out_file_cropfrac_ignored,'ignored_LUCROParea_YX_y','-v7.3') ;

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

disp('Saving LU...')
check_existing_lu(thisVer, out_file_lu, allVer_names, allVer_ignore_types) ;
lpjgu_matlab_saveTable(lu_in_header_cell, lu_in_array, out_file_lu,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 1) ;





