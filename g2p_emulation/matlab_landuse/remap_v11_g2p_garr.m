%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Re-map area/fert data to PLUM crops, and generate extra LU file %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Sam Rabin, 2021-02-06
% See previous notes for remap_v10_g2p, with the following changes:
% - Split CerealsC3 into CerealsC3s (spring) and CerealsC3w (winter)

PLUMsetAside_frac = 0.103 ;
inpaint_method = 4 ;
yearList_out = 1850:2015 ;

% Version for crop mappings
thisVer = 'WithFruitVeg_sepSugar_sepOil_sepC3' ;

force_all_rainfed = false ;

remapVer = '11_g2p' ;
out_dir = sprintf('/Volumes/Reacher/G2P/inputs/LU/remaps_v%s/', remapVer) ;

% Output LU files will contain any cell appearing in all gridlists; output
% gridlist files will just be these with any removed that got removed
% during the process. Note that files do not need to be actual gridlist
% files; this script will intelligently convert them into gridlist structs.
gridlist_files = { ...
    '/Volumes/Reacher/G2P/inputs/gridlist/gridlist_62892.runAEclimOK.txt' ;
    '/Volumes/Reacher/GGCMI/AgMIP.input/phase3/ISIMIP3/landseamask-lpjg/gridlist_ggcmi_v1.1.gapfilled.lpjg.txt' ;
    '/Volumes/Reacher/GGCMI/AgMIP.input/phase1/processed_daily_rechunked2/gridlist_agmerra.txt' ;
    } ;

% Any gridlist listed here will serve as a further restriction on auxiliary
% gridlists that will be saved in separate directories. Note that files do 
% not need to be actual gridlist files; this script will intelligently
% convert them into gridlist structs.
restrictor_gridlist_files = { ...
    '/Users/Shared/lpj-guess/input/soil/WISE/soilmap_center_interpolated.dat' ;
    } ;


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

cd '/Users/sam/Documents/git_repos/g2p_emulation/matlab_landuse'
addpath(genpath(pwd))
addpath(genpath('/Users/sam/Documents/git_repos/g2p_emulation/matlab'))

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

% Get output gridlist(s)
for f = 1:length(gridlist_files)
    file_gridlist = gridlist_files{f} ;
    if f==1
        gridlist = read_gridlist(file_gridlist) ;
        gridlists = gridlist ;
    else
        gridlists(f) = read_gridlist(file_gridlist) ;
        
        % Expand existing gridlist
        gridlist.list2map = union(gridlist.list2map, gridlists(f).list2map, 'stable') ;
        gridlist.lonlats = union(gridlist.lonlats, gridlists(f).lonlats, 'rows', 'stable') ;
    end
end
gridlist.mask_YX(gridlist.list2map) = true ;
Ncells = length(gridlist.list2map) ;

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


%% Import land uses

disp('Importing land uses...')
fprintf('  file_luh2_states: %s\n', file_luh2_states) ;

% Import cell area (km2); aggregate to half-degree
file_luh2_etc = '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc' ;
fprintf('file_luh2_etc: %s\n', file_luh2_etc) ;
carea_XY = ncread(file_luh2_etc,'carea') ;
carea_YX = flipud(transpose(carea_XY)) ;
carea_YX = coarsen_res(carea_YX,0.25,0.5) ;
carea = rmfield(gridlist, 'mask_YX') ;
carea.garr_x = carea_YX(gridlist.mask_YX) ;
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
    ii = strcmp(list_LU_out,thisLU_out) ;
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
    out_lu.garr_xvy(:,ii,:) = out_lu.garr_xvy(:,ii,:) + permute(lu_out_xy, [1 3 2]) ;
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

% Mask cells with no vegetated land
bad_x = out_lu.garr_xvy(:,4,1)==1 ...
    | sum(out_lu.garr_xvy(:,:,1),2)==0 ...
    | any(isnan(out_lu.garr_xvy(:,:,1)),2) ;
if any(bad_x)
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

% Sanity check
if any(any(any(isnan(out_lu.garr_xvy))))
    error('NaN in out_lu.garr_xvy')
end

disp('Done.')


%% Import crop fractions and process crop types

disp('Importing crop fractions...')

% Import from MIRCA.txt.maps.mat
file_cropmirca = '/Users/sam/Geodata/MIRCA/harvested_area_grids_26crops_30mn/MIRCA.txt' ;
fprintf('file_cropmirca: %s\n', file_cropmirca) ;
croparea_in = lpjgu_matlab_readTable_then2map(file_cropmirca,...
    'verboseIfNoMat',true) ;
tmp = lpjgu_YXz_to_xz(croparea_in.maps_YXv, Ncells, gridlist.list2map) ;
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
Nsplit_cropsCombined_in = ones(NcropsCombined_frac_in, 1) ;
for c = 1:NcropsCombined_frac_in
    thisCrop = list_cropsCombined_frac_in{c} ;
    found = 0 ;
    ignored = any(strcmp(list_ignore_frac,thisCrop)) ;
    for cc = 1:length(in2out_keyCombined_frac)
        found = found + length(find(strcmp(in2out_keyCombined_frac{cc},thisCrop))) ;
    end
    if found==0
        if ~ignored
            error([thisCrop ' not found in list_ignore_frac or in2out_keyCombined_frac!'])
        end
    elseif ignored
        error('%s is in both list_ignore_frac and in2out_keyCombined_frac', thisCrop)
    elseif found>1
        warning('%s found %d times in in2out_keyCombined_frac; area will be evenly divided among output types containing it.', ...
            thisCrop, found)
    end
    Nsplit_cropsCombined_in(c) = found ;
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
        case {'CerealsC3', 'CerealsC3s', 'CerealsC3w', ...
                'StarchyRoots', 'FruitAndVeg', ...
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

out_nfert.varNames = list_crops_out ;
out_nfert.garr_xv = nan([Ncells Ncrops_out]) ;
list_cropsCombined_nfert_in = unique(list_crops_out_asNfert) ;
for c = 1:length(list_cropsCombined_nfert_in)
    thisCrop_in = list_cropsCombined_nfert_in{c} ;
    thisFile = sprintf('%s/agmip_%s_apprate_fill_NPK_0.5.nc4', ...
        nfert_dir, thisCrop_in) ;
    in_YX = flipud(transpose(ncread(thisFile, 'Napprate'))) ;

    % Trim gridlist if missing from fertilizer data
    isbad_x = isnan(in_YX(gridlist.list2map)) ;
    if any(isbad_x)
        I_bad = find(isbad_x) ;
        warning('Removing %d cells that are NaN in %s fertilizer', ...
            length(I_bad), thisCrop_in)
        out_lu.garr_xvy(I_bad,:,:) = [] ;
        out_lu.list2map(I_bad) = [] ;
        out_lu.lonlats(I_bad,:) = [] ;
        croparea_in.garr_xv(I_bad,:) = [] ;
        out_nfert.garr_xv(I_bad,:) = [] ;
        gridlist.mask_YX(gridlist.list2map(I_bad)) = false ;
        gridlist.list2map(I_bad) = [] ;
        gridlist.lonlats(I_bad,:) = [] ;
        Ncells = length(gridlist.list2map) ;
    end

    v = find(strcmp(list_crops_out_asNfert, thisCrop_in)) ;
    out_nfert.garr_xv(:,v) = repmat(in_YX(gridlist.list2map), [1 length(v)]) ;
end
out_nfert.list2map = gridlist.list2map ;
out_nfert.lonlats = gridlist.lonlats ;

% Get indices of irrigated crops (for troubleshooting)
ir_inds = [] ;
for c = 1:length(out_nfert.varNames)
    thisCrop = out_nfert.varNames{c} ;
    thisCropI = [thisCrop 'i'] ;
    ir_inds = [ir_inds find(strcmp(out_nfert.varNames,thisCropI))] ;
end
Nirr = length(ir_inds) ;
isIr = false(size(out_nfert.varNames)) ;
isIr(ir_inds) = true ;
isRf = ~isIr ;
rf_inds = find(isRf) ;

% Make sure that rainfed and irrigated fertilization is equal
for i = rf_inds
    thisCrop = out_nfert.varNames{i} ;
    thisCropI = [thisCrop 'i'] ;
    if any(strcmp(out_nfert.varNames,thisCropI))
        j = find(strcmp(out_nfert.varNames,thisCropI)) ;
        tmpRF = out_nfert.garr_xv(:,i) ;
        tmpIR = out_nfert.garr_xv(:,j) ;
        tmpRF(isnan(tmpRF)) = -1 ;
        tmpIR(isnan(tmpIR)) = -1 ;
        nbad = length(find(tmpRF ~= tmpIR)) ;
        if nbad>0
            error('%s (%d,%d):\t %d\n', ...
                pad(thisCrop,max(cellfun(@length,out_nfert.varNames))), ...
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
out_nfert.garr_xv = out_nfert.garr_xv * 1e-4 ;

% Fill NaNs with zero
out_nfert.garr_xv(isnan(out_nfert.garr_xv)) = 0 ;

% Add zeros for ExtraCrop
out_nfert.varNames{end+1} = 'ExtraCrop' ;
out_nfert.garr_xv(:,end+1) = 0 ;

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
        error('Expected %d matches of thisRow in croparea_in.varNames; found %d', ...
            length(thisRow), length(IA))
    end
    
    % Divide multiple-occurring sub-crops into even parts (e.g., split
    % Wheat evenly into CerealsC3s and CerealsC3w)
    tmp_xv = zeros(Ncells, 1) ;
    for d = 1:length(thisRow)
        inCropIrr = thisRow{d} ;
        inCrop = strrep(strrep(inCropIrr, '_RF', ''), '_IR', '') ;
        Nsplit = Nsplit_cropsCombined_in(strcmp(list_cropsCombined_frac_in, inCrop)) ;
        tmp_xv = tmp_xv + croparea_in.garr_xv(:,IA(d)) / Nsplit ;
    end
    croparea_mid.garr_xv(:,c) = tmp_xv ;
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

% Interpolate
%%% Need to do this because the LUH2 crop fraction can sometimes be
%%% positive even where MIRCA has no crop information. The alternative
%%% would be to set LUH2 crop fraction to zero there. Probably doesn't make
%%% a huge difference either way, because LUH2 crop area is usually small.
out_cropfrac.varNames = list_crops_out ;
out_cropfrac.list2map = gridlist.list2map ;
out_cropfrac.lonlats = gridlist.lonlats ;
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

disp('Saving gridlists...')
outFile_gridlist = sprintf('%sgridlist.remapv%s.txt', out_dir, remapVer) ;
save_gridlist(gridlist, outFile_gridlist)
for f = 1:length(gridlist_files)
    % Restrict gridlist
    thisGridlist = gridlists(f) ;
    thisGridlist = restrict_gridlist(thisGridlist, gridlist) ;
    % Get filename
    thisFile = dir(gridlist_files{f}) ;
    thisFile = thisFile.name ;
    thisFile = strrep(thisFile, '.txt', sprintf('.remapv%s.txt', remapVer)) ;
    outFile_thisGridlist = sprintf('%s%s', out_dir, thisFile) ;
    if strcmp(outFile_gridlist, outFile_thisGridlist)
        error('Refusing to overwrite outFile_gridlist with outFile_thisGridlist. You probably have an input gridlist named gridlist.txt, which is a problem.')
    end
    % Save
    save_gridlist(thisGridlist, outFile_thisGridlist)
    
    % Further restrict? If so, save to ancillary files.
    if ~isempty(restrictor_gridlist_files)
        for r = 1:length(restrictor_gridlist_files)
            restrictorFile = restrictor_gridlist_files{r} ;
            finfo = dir(restrictorFile) ;
            restrictor = read_gridlist(restrictorFile) ;
            thisFile_out = finfo.name ;
            if strcmp(thisFile_out, 'gridlist.txt')
                error('This will be a problem')
            end
            thisFile_out = strrep(thisFile_out, 'gridlist_', '') ;
            thisFile_out = strrep(thisFile_out, '.txt', '') ;
            thisFile_out = strrep(thisFile_out, '.dat', '') ;
            thisDir_out = thisFile_out ;
            thisDir_out = sprintf('%s%s', out_dir, thisDir_out) ;
            if ~exist(thisDir_out, 'dir')
                mkdir(thisDir_out)
            end
            outFile_thisRestrGridlist = sprintf('%s/%s', thisDir_out, thisFile) ;
            outFile_thisRestrGridlist = strrep(outFile_thisRestrGridlist, ...
                '.txt', sprintf('.%s.txt', thisFile_out)) ;
            thisGridlist_restr = restrict_gridlist(thisGridlist, restrictor) ;
            save_gridlist(thisGridlist_restr, outFile_thisRestrGridlist)
        end
    end
    
end



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

function gridlist = restrict_gridlist(gridlist, restrictor)

[toRemove, I_toRemove] = setdiff(gridlist.list2map, restrictor.list2map) ;
gridlist.list2map(I_toRemove) = [] ;
gridlist.lonlats(I_toRemove,:) = [] ;
gridlist.mask_YX(toRemove) = false ;

end


function gridlist = read_gridlist(file_gridlist)

gridlist = lpjgu_matlab_read2geoArray(file_gridlist) ;
if ~isfield(gridlist, 'mask_YX')
    tmp2.lonlats = gridlist.lonlats;
    tmp2.list2map = gridlist.list2map;
    clear gridlist
    tmp2.mask_YX = false(360,720) ;
    tmp2.mask_YX(tmp2.list2map) = true;
    gridlist = tmp2 ;
end

end


function save_gridlist(thisGridlist, outFile_gridlist)

% Set up
lons4map = -179.75:0.5:179.75 ;
lats4map = -89.75:0.5:89.75 ;
lons_map = repmat(lons4map,[length(lats4map) 1]) ;
lats_map = repmat(lats4map',[1 length(lons4map)]) ;
lons_4gl = lons_map(thisGridlist.mask_YX) ;
lats_4gl = lats_map(thisGridlist.mask_YX) ;
Ncells_4gl = length(lons_4gl) ;
% Get random order for output
rng(20210106) ;
rdmsam = randsample(Ncells_4gl,Ncells_4gl) ;
% Save gridlist

out_formatSpec_gridlist = '%4.2f %4.2f\n' ;
fid1_gridlist = fopen(strrep(outFile_gridlist,'\ ',' '),'w') ;
fprintf(fid1_gridlist,out_formatSpec_gridlist,[lons_4gl(rdmsam) lats_4gl(rdmsam)]') ;
fclose(fid1_gridlist) ;
fid1_gridlist = fopen(strrep(outFile_gridlist,'\ ',' '),'a+') ;
fprintf(fid1_gridlist,'%s','') ;
fclose(fid1_gridlist) ;

end


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

