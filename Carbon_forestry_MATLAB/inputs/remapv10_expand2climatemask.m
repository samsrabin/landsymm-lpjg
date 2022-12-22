%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Like remap_v10_g2p_garr.m, but expand gridlist to match climate mask %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sam Rabin, 2022-09-30

%%%%% Sam Rabin, 2021-02-05
% See previous notes for remap_v9_g2p, with the following changes:
% - Split Oilcrops into OilNfix and OilOther

PLUMsetAside_frac = 0.103 ;
inpaint_method = 4 ;
yearList_out = 1850:2015 ;

% Version for crop mappings
thisVer = 'WithFruitVeg_sepSugar_sepOil' ;

force_all_rainfed = false ;

remapVer = '10_g2p_isimipclimMask' ;
out_dir = sprintf('/Users/Shared/LandSyMM/inputs/LU/remaps_v%s/',remapVer) ;


%% Setup

if ~exist(out_dir,'dir')
    mkdir(out_dir) ;
end

% Start diary logging and print important info
diaryfile = sprintf('%s/matlab_log.txt', out_dir) ;
diary off
if exist(diaryfile, 'file')
    delete(diaryfile)
end
diary(diaryfile)
diary on
fprintf('PLUMsetAside_frac: %g\n', PLUMsetAside_frac)
fprintf('inpaint_method: %d\n', inpaint_method)
fprintf('yearList_out: %d-%d\n', min(yearList_out), max(yearList_out))
fprintf('thisVer: %s\n', thisVer)
fprintf('remapVer: %s\n', remapVer)
fprintf('force_all_rainfed: %d\n', force_all_rainfed)

addpath(genpath(landsymm_lpjg_path()))

Nyears_out = length(yearList_out) ;

% LUH2 input files and years they contain
file_luh2_states = '/Users/sam/Geodata/LUH2/v2h/states.1850-2015.nc' ;
yearList_luh2_states = 1850:2015 ;
if min(yearList_out) < min(yearList_luh2_states) || max(yearList_out) > max(yearList_luh2_states)
    error('yearList_out must be entirely contained within yearList_luh2_states!')
end
file_luh2_mgmts = '/Users/sam/Geodata/LUH2/v2h/management.1850-2015.nc' ;
yearList_luh2_mgmts = 1850:2015 ;
if min(yearList_out) < min(yearList_luh2_mgmts) || max(yearList_out) > max(yearList_luh2_mgmts)
    error('yearList_out must be entirely contained within yearList_luh2_mgmts!')
end

% Get info for reading input files
[~,yearIs_luh2_states,~] = intersect(yearList_luh2_states,yearList_out) ;
starts_luh2_states = [1 1 min(yearIs_luh2_states)] ;
counts_luh2_states = [Inf Inf Nyears_out] ;
[~,yearIs_luh2_mgmts,~] = intersect(yearList_luh2_mgmts,yearList_out) ;
starts_luh2_mgmts = [1 1 min(yearIs_luh2_mgmts)] ;
counts_luh2_mgmts = [Inf Inf Nyears_out] ;

% Get land use names from input; match to names in output
finfo = ncinfo(file_luh2_states) ;
list_LU_in = setdiff({finfo.Variables.Name},...
                     {'time','lon','lat','lat_bounds','lon_bounds','secma','secmb'}) ;
list_LU_out = {'NATURAL','CROPLAND','PASTURE','BARREN'} ;
Nlu_in = length(list_LU_in) ;
Nlu_out = length(list_LU_out) ;
map_LU_in2out = cell(size(list_LU_in)) ;
map_LU_in2out(contains(list_LU_in,...
    {'c3ann','c3nfx','c3per','c4ann','c4per'})) = {'CROPLAND'} ;
map_LU_in2out(contains(list_LU_in,...
    {'pastr','range'})) = {'PASTURE'} ;
map_LU_in2out(contains(list_LU_in,...
    {'primf','primn','secdf','secdn','secma','secmb'})) = {'NATURAL'} ;
map_LU_in2out(contains(list_LU_in,...
    {'urban'})) = {'BARREN'} ;
if any(isempty(map_LU_in2out))
    error('Remaining empty element(s) in map_LU_in2out!')
elseif ~isempty(setdiff(map_LU_in2out,list_LU_out))
    error('Some member of map_LU_in2out is not present in list_LU_out!')
end
warning('Should work out specification of mapping with check for duplicates on LHS.')


warning('on','all')


%% Get climate mask and gridlist

climate_file = '/Users/Shared/GGCMI/AgMIP.input/phase3/ISIMIP3/climate_land_only/climate3b/historical/UKESM1-0-LL-lpjg/ukesm1-0-ll_r1i1p1f2_w5e5_historical_pr_global_daily_1850_2014.nc4' ;
fprintf('Reading land mask from first timestep of %s', climate_file)
tmp = ncread(climate_file, 'pr', [1 1 1], [Inf Inf 1]) ;
tmp = flip(transpose(tmp), 1) ;
climate_ok_YX = ~isnan(tmp) ;
shademap(climate_ok_YX)
disp('Done reading climate file.')

% Save to gridlist
file_gridlist = sprintf('%sgridlist.remapv%s.txt', out_dir, remapVer) ;
lpjgu_saveMap2gridlist(climate_ok_YX, file_gridlist)

% Read gridlist in expected format
gridlist = lpjgu_matlab_read2geoArray(file_gridlist, ...
    'verboseIfNoMat', true, 'force_mat_save', false, 'force_mat_nosave', true) ;
Ncells = length(gridlist.list2map) ;


%% Import land uses

disp('Importing land uses...')
fprintf('  file_luh2_states: %s\n', file_luh2_states) ;

% Import cell area (km2); aggregate to half-degree
file_luh2_etc = '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc' ;
fprintf('file_luh2_etc: %s\n', file_luh2_etc) ;
larea = import_staticData(file_luh2_etc, 'icwtr', gridlist) ;
larea.garr_x = 1 - larea.garr_x ;
carea = import_staticData(file_luh2_etc, 'carea', gridlist) ;
carea_XYy = repmat(carea_XY,[1 1 Nyears_out]) ;
carea_hd_XY = coarsen_res(carea_XY,0.25,0.5) ;
carea_hd_XYy = repmat(carea_hd_XY,[1 1 Nyears_out]) ;

% Land use
out_lu = rmfield(gridlist, 'mask_YX') ;
out_lu.varNames = list_LU_out ;
out_lu.yearList = yearList_out ;
for v = 1:Nlu_in
    thisLU_in = list_LU_in{v} ;
    thisLU_out = map_LU_in2out{v} ;
    fprintf('    %s (%d of %d) to %s...\n',thisLU_in,v,Nlu_in,thisLU_out)
    i = strcmp(list_LU_out,thisLU_out) ;
    lu_in_XYy = ncread(file_luh2_states,thisLU_in,starts_luh2_states,counts_luh2_states) ;
    lu_in_XYy(isnan(lu_in_XYy)) = 0 ;
    lu_out_YXy = flipud(permute( ...
        coarsen_res(lu_in_XYy.*carea_XYy,0.25,0.5)./carea_hd_XYy, ...
        [2 1 3])) ;
    if any(any(any(lu_out_YXy-1 > 1e-6)))
        error('Some element(s) of lu_out_YXy > 1!')
    end
    clear lu_in_XYy
    % Convert to _xy
    lu_out_xy = lpjgu_YXz_to_xz(lu_out_YXy, [Ncells, size(lu_out_YXy, 3)], gridlist.list2map) ;
    clear lu_out_YXy
    
    if v==1
        out_lu.garr_xvy = zeros([Ncells length(out_lu.varNames) size(lu_out_xy, 2)]) ;
    end
    out_lu.garr_xvy(:,i,:) = out_lu.garr_xvy(:,i,:) + permute(lu_out_xy, [1 3 2]) ;
    if any(any(any(any(out_lu.garr_xvy-1 > 1e-6))))
        error('Some element(s) of out_lu.garr_xvy > 1!')
    end
    clear lu_out_xy
end
clear carea*_YXy carea_hd_XYy

disp('Finishing...')

% Add water fraction to BARREN
icwtr_YX = flipud(transpose(ncread(file_luh2_etc,'icwtr'))) ;
icwtr_YX(icwtr_YX==1) = 0 ;
icwtr_hd_YX = coarsen_res(icwtr_YX.*flipud(carea_XY'),0.25,0.5)./flipud(carea_hd_XY') ;
icwtr_hd_x = icwtr_hd_YX(gridlist.list2map) ;
v = strcmp(list_LU_out,'BARREN') ;
out_lu.garr_xvy(:,v,:) = out_lu.garr_xvy(:,v,:) ...
    + repmat(icwtr_hd_x,[1 1 Nyears_out]) ;

% Mask cells with no vegetated land according to LU dataset
bad_x = out_lu.garr_xvy(:,4,1)==1 ;
if any(bad_x)
    error('Decide how you want to handle masking of cells with no vegetated land.')
    fprintf('Removing %d cells with no vegetated land...\n', length(find(bad_x))) ;
    out_lu.garr_xvy(bad_x,:,:) = [] ;
    out_lu.list2map(bad_x) = [] ;
    out_lu.lonlats(bad_x,:) = [] ;
    gridlist.mask_YX(gridlist.list2map(bad_x)) = false ;
    gridlist.list2map(bad_x) = [] ;
    gridlist.lonlats(bad_x,:) = [] ;
    Ncells = length(gridlist.list2map) ;
end

% Force all land cells to sum to 1
lu_out_x1y = sum(out_lu.garr_xvy,2) ;
j = 0 ;
while any(any(abs(lu_out_x1y-1)>1e-6))
    j = j + 1;
    if j > 50
        error('Possible infinite loop in "Force all land cells to sum to 1".')
    end
    out_lu.garr_xvy = out_lu.garr_xvy ./ repmat(lu_out_x1y,[1 Nlu_out 1]) ;
    lu_out_x1y = sum(out_lu.garr_xvy,2) ;
end
clear lu_out_x1y

disp('Done.')

%%
tmp_x = any(any(isnan(out_lu.garr_xvy), 3), 2) ;
% tmp_x = any(any(isnan(out_lu.garr_xvy), 3), 2) & larea.garr_x>0 ;
shademap(lpjgu_vector2map(tmp_x, [360 720], gridlist.list2map));

%% Import crop fractions and process crop types

disp('Importing crop fractions...')

% Import from MIRCA.txt.maps.mat
file_cropmirca = '/Users/sam/Geodata/MIRCA/harvested_area_grids_26crops_30mn/MIRCA.txt' ;
fprintf('file_cropmirca: %s\n', file_cropmirca) ;
croparea_in = lpjgu_matlab_readTable_then2map(file_cropmirca,...
    'verboseIfNoMat',true) ;
tmp = lpjgu_YXz_to_xz(croparea_in.maps_YXv, length(gridlist.list2map), gridlist.list2map) ;
croparea_in = rmfield(croparea_in, 'maps_YXv') ;
croparea_in.garr_xv = tmp ;
clear tmp

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
        croparea_in.garr_xv(:,thisRF) = croparea_in.garr_xv(:,thisRF) + croparea_in.garr_xv(:,thisIR) ;
        croparea_in.garr_xv(:,thisIR) = 0*croparea_in.garr_xv(:,thisIR) ;
    end
end
% Check names
for c = 1:NcropsCombined_frac_in
    thisCrop = list_cropsCombined_frac_in{c} ;
    isThisCrop = not(cellfun(@isempty,strfind(croparea_in.varNames,thisCrop))) ;
    if length(find(isThisCrop))~=2
        error('length(find(isThisCrop))~=2')
    end
end

% Get mappings
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

disp('Done.')


%% Import Nfert from GGCMI phase 1
% No manure, because that wasn't included in Phase 1

disp('Importing Nfert...')
nfert_dir = '/Users/Shared/unifying_gridlist/AgGRID_nutrient_input_v1.1' ;
fprintf('nfert_dir: %s\n', nfert_dir) ;

% Map input to output crop types. Make sure we're being consistent w/r/t
% remapping of emulator output types to LPJ-GUESS types!
list_crops_out_asNfert = list_crops_out ;
for c = 1:Ncrops_out
    thisCrop_out = list_crops_out{c} ;
    if strcmp(thisCrop_out(end), 'i')
        thisCrop_out = thisCrop_out(1:end-1) ;
    end
    switch thisCrop_out
        case {'CerealsC3', 'StarchyRoots', 'FruitAndVeg', ...
                'Sugarbeet', 'OilOther', 'ExtraCrop'}
            thisCrop_in = 'wheat' ;
        case {'CerealsC4', 'Miscanthus', 'Sugarcane'}
            thisCrop_in = 'maize' ;
        case {'Oilcrops', 'Pulses', 'OilNfix'}
            thisCrop_in = 'soybean' ;
        case {'Rice'}
            thisCrop_in = 'rice' ;
        otherwise
            error('What crop from Nfert inputs should I use for %s?', thisCrop_out)
    end
    list_crops_out_asNfert{c} = thisCrop_in ;
end

mid_nfert.varNames = list_crops_out ;
mid_nfert.garr_xv = nan([Ncells Ncrops_out]) ;
list_cropsCombined_nfert_in = unique(list_crops_out_asNfert) ;
for c = 1:length(list_cropsCombined_nfert_in)
    thisCrop_in = list_cropsCombined_nfert_in{c} ;
    thisFile = sprintf('%s/agmip_%s_apprate_fill_NPK_0.5.nc4', ...
        nfert_dir, thisCrop_in) ;
    in_YX = flipud(transpose(ncread(thisFile, 'Napprate'))) ;

    % Note climate cells missing from fertilizer
    isbad_x = isnan(in_YX(gridlist.list2map)) ;
    if any(isbad_x)
        I_bad = find(isbad_x) ;
        warning('%d climate cells are missing from fertilizer for %s', ...
            length(I_bad), thisCrop_in)
    end

    v = find(strcmp(list_crops_out_asNfert, thisCrop_in)) ;
    mid_nfert.garr_xv(:,v) = repmat(in_YX(gridlist.list2map), [1 length(v)]) ;
end
mid_nfert.list2map = gridlist.list2map ;
mid_nfert.lonlats = gridlist.lonlats ;

% Get indices of irrigated crops (for troubleshooting)
ir_inds = [] ;
for c = 1:length(mid_nfert.varNames)
    thisCrop = mid_nfert.varNames{c} ;
    thisCropI = [thisCrop 'i'] ;
    ir_inds = [ir_inds find(strcmp(mid_nfert.varNames,thisCropI))] ;
end
Nirr = length(ir_inds) ;
isIr = false(size(mid_nfert.varNames)) ;
isIr(ir_inds) = true ;
isRf = ~isIr ;
rf_inds = find(isRf) ;

% Make sure that rainfed and irrigated fertilization is equal
for i = rf_inds
    thisCrop = mid_nfert.varNames{i} ;
    thisCropI = [thisCrop 'i'] ;
    if any(strcmp(mid_nfert.varNames,thisCropI))
        j = find(strcmp(mid_nfert.varNames,thisCropI)) ;
        tmpRF = mid_nfert.garr_xv(:,i) ;
        tmpIR = mid_nfert.garr_xv(:,j) ;
        tmpRF(isnan(tmpRF)) = -1 ;
        tmpIR(isnan(tmpIR)) = -1 ;
        nbad = length(find(tmpRF ~= tmpIR)) ;
        if nbad>0
            error('%s (%d,%d):\t %d\n', ...
                pad(thisCrop,max(cellfun(@length,mid_nfert.varNames))), ...
                i,j,nbad)
        end
        clear tmpRF tmpIR
    else
%             fprintf('Skipping %s.\n',thisCrop) ;
    end
    clear thisCrop thisCropI tmp*
end
disp('Rainfed and irrigated fertilization is equal.')

% Convert from kg/ha to kg/m2
mid_nfert.garr_xv = mid_nfert.garr_xv * 1e-4 ;

disp('Done.')


%% Process crop fractions

disp ('Processing crop fractions...')

% Crop fractions
croparea_mid.garr_xv = nan(Ncells,Ncrops_out) ;
for c = 1:Ncrops_out
    thisCrop = list_crops_out{c} ;
    thisRow = in2out_key_frac{getOi(thisCrop)} ;
    [~,IA] = intersect(croparea_in.varNames,thisRow) ;
    if length(IA)~=length(thisRow)
        error('length(IA)~=length(thisRow)')
    end
    % Straight sum
    croparea_mid.garr_xv(:,c) = sum(croparea_in.garr_xv(:,IA),2) ;
end
croparea_mid.varNames = list_crops_out ;

% Move ignored area (unhandled crops) and setAside area into ExtraCrop
I = contains(croparea_in.varNames,list_ignore_frac) ;
ignore_area_x = sum(croparea_in.garr_xv(:,I),2) ;
setaside_area_x = sum(croparea_mid.garr_xv * PLUMsetAside_frac,2) ;
croparea_mid.garr_xv = croparea_mid.garr_xv * (1-PLUMsetAside_frac) ;
extra_x = ignore_area_x + setaside_area_x ;
croparea_mid.garr_xv = cat(2,croparea_mid.garr_xv,extra_x) ;
list_cropsCombined_out = [list_cropsCombined_out {'ExtraCrop'}] ;
NcropsCombined_out = length(list_cropsCombined_out) ;
list_crops_out = [list_crops_out {'ExtraCrop'}]  ;
Ncrops_out = length(list_crops_out) ;

% Get fractions
cropfrac_mid.garr_xv = croparea_mid.garr_xv ./ repmat(sum(croparea_mid.garr_xv,2),[1 Ncrops_out]) ;

disp('Done.')


%% Interpolate land use inputs to climate mask

disp('Interpolating to match climate mask...')

out_cropfrac.varNames = cropfrac_mid ;
out_cropfrac.garr_xv = nan(Ncells,Ncrops_out) ;
% For checking how much this increases global area of each crop relative to
% what it would have been if, instead of interpolating, we set LUH2 crop
% area to zero.
test_y1 = 1995 ;
test_yN = 2005 ;
test_croparea_x = mean(out_lu.garr_xvy(:, strcmp(out_lu.varNames, 'CROPLAND'), ...
    yearList_out>=test_y1 & yearList_out<=test_yN),3) ...
    .* carea_YX(gridlist.list2map) ;
for c = 1:Ncrops_out
    thisCrop = list_crops_out{c} ;
    fprintf('Interpolating %s (%d of %d)...\n',thisCrop, c, Ncrops_out)
    tmp0_YX = lpjgu_xz_to_YXz(cropfrac_mid.garr_xv(:,c), size(gridlist.mask_YX), gridlist.list2map) ;
    tmp1_YX = inpaint_nans(tmp0_YX, inpaint_method) ;
    
    tmp0_YX(isnan(tmp0_YX)) = 0 ;
    area0 = sum(test_croparea_x .* tmp0_YX(gridlist.list2map)) ;
    area1 = sum(test_croparea_x .* tmp1_YX(gridlist.list2map)) ;
    areaDiff = area1 - area0 ;
    areaDiff_pct = areaDiff / area0 * 100 ;
    fprintf('    Area increased %0.2f%% (%0.1g ha)\n', ...
        areaDiff_pct, areaDiff*100) ;
    out_cropfrac.garr_xv(:,c) = tmp1_YX(gridlist.list2map) ;
    if any(isnan(out_cropfrac.garr_xv(:,c)))
        error('NaN remaining in out_cropfrac.garr_xv(:,c)!')
    end
    clear tmp*
end

% % % disp('Interpolating include_frac...')
% % % map_include_frac_interpd_YX = inpaint_nans(map_include_frac_YX,inpaint_method) ;
% % % disp('Done.')

% Set minimum of zero, if needed
if any(out_cropfrac.garr_xv(:)<0)
    warning('Setting negative members of out_cropfrac.garr_xv to zero.')
    out_cropfrac.garr_xv(out_cropfrac.garr_xv<0) = 0 ;
end
% % % if any(map_include_frac_interpd_YX(:)<0)
% % %     warning('Setting negative members of map_include_frac_interpd_YX to zero.')
% % %     map_include_frac_interpd_YX(map_include_frac_interpd_YX<0) = 0 ;
% % % end

% Normalize cropfrac to 1, if needed
tmp = sum(out_cropfrac.garr_xv,2) ;
if any(abs(tmp(:) - 1)>1e-6)
    warning('Normalizing cropfrac sums to 1')
    out_cropfrac.garr_xv = out_cropfrac.garr_xv ./ repmat(tmp,[1 Ncrops_out]) ;
end
clear tmp

disp('Done.')


%% Set up for save

% Get headers
out_lu_header_cell = [ ...
    {'Lon', 'Lat', 'Year'}, ...
    out_lu.varNames] ;
out_cropfrac_header_cell = [ ...
    {'Lon', 'Lat'}, ...
    out_cropfrac.varNames, ...
    {'Miscanthus','Miscanthusi'}] ;
out_nfert_header_cell = [ ...
    {'Lon', 'Lat'}, ...
    out_nfert.varNames, ...
    {'Miscanthus','Miscanthusi'}] ;

% Get filenames
out_file_lu = [out_dir 'LU.remapv' remapVer '.txt'] ;
out_file_cropfrac = [out_dir 'cropfracs.remapv' remapVer '.txt'] ;
out_file_nfert = [out_dir 'nfert.remapv' remapVer '.txt'] ;

% Add zeros for Miscanthus(i)
disp('Adding zeros for Miscanthus(i)...')
out_cropfrac.garr_xv = cat(2, ...
    out_cropfrac.garr_xv, ...
    zeros([Ncells 2])) ;
out_nfert.garr_xv = cat(2, ...
    out_nfert.garr_xv, ...
    zeros([Ncells 2])) ;

out_cropfrac_header_cell = [ ...
    out_cropfrac_header_cell, ...
    {'Miscanthus','Miscanthusi'}] ;
out_nfert_header_cell = [ ...
    out_nfert_header_cell, ...
    {'Miscanthus','Miscanthusi'}] ;

out_cropfrac.varNames = [ ...
    out_cropfrac.varNames, ...
    {'Miscanthus','Miscanthusi'}] ;
out_nfert.varNames = [ ...
    out_nfert.varNames, ...
    {'Miscanthus','Miscanthusi'}] ;

disp('Done.')


%% Save

%%% Options %%%%%%%%%%%%%
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
do_gzip = false ;
%%%%%%%%%%%%%%%%%%%%%%%%%

% % disp('Saving ignored crop fraction...')
% % out_file_cropfrac_ignored = strrep(out_file_cropfrac,'cropfracs','cropfracsIGNORED') ;
% % out_file_cropfrac_ignored = strrep(out_file_cropfrac_ignored,'.txt','.mat') ;
% % save(out_file_cropfrac_ignored,'ignored_LUCROParea_YX_y','-v7.3') ;

disp('Saving gridlist...')
lons4map = -179.75:0.5:179.75 ;
lats4map = -89.75:0.5:89.75 ;
lons_map = repmat(lons4map,[length(lats4map) 1]) ;
lats_map = repmat(lats4map',[1 length(lons4map)]) ;
lons_4gl = lons_map(gridlist.mask_YX) ;
lats_4gl = lats_map(gridlist.mask_YX) ;
Ncells_4gl = length(lons_4gl) ;
% Get random order for output
rng(20210106) ;
rdmsam = randsample(Ncells_4gl,Ncells_4gl) ;
% Save gridlist
outFile_gridlist = sprintf('%sgridlist.remapv%s.txt', out_dir, remapVer) ;
out_formatSpec_gridlist = '%4.2f %4.2f\n' ;
fid1_gridlist = fopen(strrep(outFile_gridlist,'\ ',' '),'w') ;
fprintf(fid1_gridlist,out_formatSpec_gridlist,[lons_4gl(rdmsam) lats_4gl(rdmsam)]') ;
fclose(fid1_gridlist) ;
fid1_gridlist = fopen(strrep(outFile_gridlist,'\ ',' '),'a+') ;
fprintf(fid1_gridlist,'%s','') ;
fclose(fid1_gridlist) ;

disp('Saving LU...')
% check_existing_lu(thisVer, out_file_lu, allVer_names, allVer_ignore_types) ;
lpjgu_matlab_saveTable(out_lu_header_cell, out_lu, out_file_lu,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 1, ...
    'gzip', do_gzip) ;
clear out_lu_array

disp('Saving cropfracs...')
lpjgu_matlab_saveTable(out_cropfrac_header_cell, out_cropfrac, out_file_cropfrac,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20, ...
    'gzip', do_gzip) ;
clear out_cropfrac_array

disp('Saving nfert...') %#ok<*UNRCH>
lpjgu_matlab_saveTable(out_nfert_header_cell, out_nfert, out_file_nfert,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20, ...
    'gzip', do_gzip) ;
clear out_nfert_array

diary off


 %% FUNCTIONS
 
function out_array = coarsen_res(in_array,in_res,out_res)

Ndims_in = length(find(size(in_array)>1)) ;
if Ndims_in>4
    error('Adapt code to work with >4 dimensions!')
end

res_ratio = out_res / in_res ;
if ~isint(res_ratio)
    error('out_res/in_res must be an integer!')
end

tmp = zeros(size(in_array,1), ...
    size(in_array,2)/res_ratio, ...
    size(in_array,3),size(in_array,4)) ;
for j = 1:res_ratio
    tmp = tmp + in_array(:,j:res_ratio:end,:,:) ;
end

out_array = zeros(size(in_array,1)/res_ratio, ...
    size(in_array,2)/res_ratio, ...
    size(in_array,3),size(in_array,4)) ;
for i = 1:res_ratio
    out_array = out_array + tmp(i:res_ratio:end,:,:,:) ;
end


 end

 
 function out_array = coarsen_res_mgmt(in_array,inArea_array,in_res,out_res)

Ndims_in = length(find(size(in_array)>1)) ;
if Ndims_in>4
    error('Adapt code to work with >4 dimensions!')
end

res_ratio = out_res / in_res ;
if ~isint(res_ratio)
    error('out_res/in_res must be an integer!')
end

% Where there is no area, set mgmt to NaN
in_array(inArea_array==0) = NaN ;

% Convert NaNs to zeros
in_array(isnan(in_array)) = 0 ;

tmp = zeros(size(in_array,1), ...
    size(in_array,2)/res_ratio, ...
    size(in_array,3),size(in_array,4)) ;
tmpArea = tmp ;
for j = 1:res_ratio
    % Area in new column
    newArea = inArea_array(:,j:res_ratio:end,:,:) ;
    % The current value of mgmt input, weighted by area totaled so far as
    % fraction of so-far plus new area
    tmpWtd = tmp .* tmpArea./(tmpArea+newArea) ;
    % The value of mgmt input in the new column, weighted by area in new
    % column as fraction of so-far plus new area
    newWtd = in_array(:,j:res_ratio:end,:,:) .* newArea./(tmpArea+newArea) ;
    % Where no so-far or new area, set old and new values to zero.
    tmpWtd(tmpArea+newArea==0) = 0 ;
    newWtd(tmpArea+newArea==0) = 0 ;
    % Update aggregated management input and area
    tmp = tmpWtd + newWtd ;
    tmpArea = tmpArea + newArea ;
end

out_array = zeros(size(in_array,1)/res_ratio, ...
    size(in_array,2)/res_ratio, ...
    size(in_array,3),size(in_array,4)) ;
outArea_array = out_array ;
for i = 1:res_ratio
    % Area in new row
    newArea = tmpArea(i:res_ratio:end,:,:,:) ;
    % The current value of mgmt input, weighted by area totaled so far as
    % fraction of so-far plus new area
    out2wtd = out_array .* outArea_array./(outArea_array+newArea) ;
    % The value of mgmt input in the new row, weighted by area in new
    % row as fraction of so-far plus new area
    new2wtd = tmp(i:res_ratio:end,:,:,:) .* newArea./(outArea_array+newArea) ;
    % Where no so-far or new area, set old and new values to zero.
    out2wtd(outArea_array+newArea==0) = 0 ;
    new2wtd(outArea_array+newArea==0) = 0 ;
    % Update aggregated management input and area
    out_array = out2wtd + new2wtd;
    outArea_array = outArea_array + newArea ;
end


% Where there is no area, set mgmt to NaN
out_array(outArea_array==0) = NaN ;


end


function S = import_staticData(file_luh2_etc, varName, gridlist)

map_XY = ncread(file_luh2_etc, varName) ;
map_YX = flipud(transpose(map_XY)) ;
map_YX = coarsen_res(map_YX,0.25,0.5) ;
S = rmfield(gridlist, 'mask_YX') ;
S.garr_x = map_YX(gridlist.mask_YX) ;

end

