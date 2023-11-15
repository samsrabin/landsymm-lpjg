%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create "observation"-based LU inputs for LPJ-GUESS with PLUM crops %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(landsymm_lpjg_path()))
rmpath(genpath(fullfile(landsymm_lpjg_path(), '.git')))

% PLUMharm_options.m must be somewhere on your path.
% There, specify the following variables:
%
% DIAGNOSTIC OPTIONS
% As part of remapping, we interpolate (if necessary) the land use areas and fertilizer
% maps to fill out the entire gridlist. We want to check how much area and fertilizer this
% adds relative to what we'd have if we set missing cells to 0 instead of interpolating.
% By default, this looks at the entire period specified by yearList_out. However, we can
% specify a custom period using the following variables. (The recommended values were
% chosen as the central MIRCA2000 year +/- 5 years.)
%     interp_test_y1 (optional; recommendation: 1995)
%     interp_test_yN (optional; recommendation: 2005)
%
% OUTPUT OPTIONS
%     delimiter: Recommendation: ' '.
%     do_gzip: Whether to zip up outputs from .txt. to .gz files (true/false).
%     fancy: Not sure. Recommendation: false.
%     outPrec: Precision (number of decimal places) to which output .txt files should be
%              saved. (Recommendation: 6)
%     outWidth: Width of output columns. Recommendation is 1; columns should auto-adjust.
%     overwrite: Should existing output files be overwritten? (true/false)
%     
% PATHS
%     file_gridlist_climate: (Optional.) Path of LPJ-GUESS gridlist file with the extent
%                            of climate data. This will be used to produce a diagnostic
%                            figure only.
%     file_gridlist_out: Path of LPJ-GUESS gridlist file to use as the basis for output
%                        files. Inputs will be interpolated as necessary to match. 
%                        Recommendation: '/path/to/gridlist_62892.runAEclimOK.txt'
%     geodata_dir: Directory where misc. geographic data can be found, including LUH2/,
%                  MIRCA/, AND AgGRID_nutrient_input_v1.1/.
%     lpjg_inputs_dir: Directory where LPJ-GUESS-related soil/ directory can be found.
%                      Unicluster: ~xg4606/input/lpj-guess
%     out_dir_top: A new subdirectory will be created here (if needed); that's where
%                  outputs from this script will go.
%
% REMAPPING BEHAVIOR
%     force_all_rainfed: Move all irrigated area to rainfed? Not sure why we'd want this
%                        to ever be true; maybe for troubleshooting? Recommendation: false
%     inpaint_method: Method (integer ≥0) to be used by inpaint_nans(). See help of that
%                     function for more information. Recommendation: 4.
%     lu_source: String naming the source to use for land use fractions (CROPLAND,
%                PASTURE, etc.). 
%                Options: 'HILDA+', 'LUH2'
%     PLUMsetAside_frac: PLUM parameter (crop_adj_factor) accounting for crops not
%                        represented, areas set aside and failed crops.
%     remapVer: String to be included in output filenames. E.g.,
%                   sprintf('LU.remapv%s.txt', remapVer)
%               Also determines the subdirectory within out_dir_top to which outputs will
%               be written:
%                   sprintf('remaps_v%s',remapVer)
%     thisVer: String given to get_remapv2_keys() in order to retrieve the mapping of crop
%              types between LPJ-GUESS and PLUM. E.g., 'WithFruitVeg_sepSugar_sepOil'.
%     yearList_out: Years to be saved in output files. E.g., 1850:2015.

remap_options


%% Setup

out_dir = fullfile(out_dir_top, sprintf('remaps_v%s',remapVer)) ;
if ~exist(out_dir,'dir')
    mkdir(out_dir) ;
end
out_dir_figs = fullfile(out_dir, 'figs') ;
if ~exist(out_dir_figs, 'dir')
    mkdir(out_dir_figs) ;
end

% Start diary logging and print important info
diaryfile = fullfile(out_dir, sprintf('matlab_log.txt')) ;
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

% For checking how much interpolation increases global area/nfert of each LU relative to
% what it would have been if, instead of interpolating, we set area/nfert there to zero.
if ~exist('interp_test_y1', 'var')
    interp_test_y1 = yearList_out(1) ;
end
if ~exist('interp_test_yN', 'var')
    interp_test_yN = yearList_out(end) ;
end
yearList_interp_test = interp_test_y1:interp_test_yN ;
[C, I_yrs_in_test_period] = intersect(yearList_out, yearList_interp_test) ;
if isempty(C)
    warning( ...
        ['Years provided for check of interpolation-caused increase of global LU ' ...
         'area & Nfert (%d-%d) are outside yearList_out (%d-%d); using entire ' ...
         'yearList_out' ...
         ], ...
        interp_test_y1, interp_test_yN, yearList_out(1), yearList_out(end))
    I_yrs_in_test_period = 1:Nyears_out ;
    interp_test_y1 = yearList_out(1) ;
    interp_test_yN = yearList_out(end) ;
elseif ~isequal(shiftdim(C), shiftdim(yearList_interp_test))
    warning( ...
        ['Years provided for check of interpolation-caused increase of global LU ' ...
         'area & Nfert (%d-%d) aren''t fully in yearList_out (%d-%d); using ' ...
         'intersection (%d-%d)' ...
         ], ...
        interp_test_y1, interp_test_yN, yearList_out(1), yearList_out(end), C(1), C(end))
    interp_test_y1 = C(1) ;
    interp_test_yN = C(end) ;
end
interp_test_period_str = sprintf('%d-%d', interp_test_y1, interp_test_yN) ;

warning('on','all')


%% Get output gridlist

gridlist = lpjgu_matlab_read2geoArray(file_gridlist_out, ...
    'verboseIfNoMat', false, 'force_mat_save', false, 'force_mat_nosave', true) ;
Ncells = length(gridlist.list2map) ;


%% Map what's missing from climate gridlist

if ~exist('climate_gridlist_file', 'var')
    warning('climate_gridlist_file not provided; skipping check of what''s missing')
elseif ~exist(file_gridlist_climate, 'file')
    warning('climate_gridlist_file (%s) not found; skipping check of what''s missing', ...
        file_gridlist_climate)
else
    climate_gridlist = lpjgu_matlab_read2geoArray(...
    file_gridlist_climate, ...
    'verboseIfNoMat', false, 'force_mat_save', false, 'force_mat_nosave', true) ;

    fig_title = sprintf('%d gridcells missing from ISIMIP3 climate', ...
        length(find(gridlist.mask_YX & ~climate_gridlist.mask_YX))) ;
    warning('%s (see missing_climate.csv); these will be filled by LPJ-GUESS', fig_title)
    fig_outfile = fullfile(out_dir_figs, 'cells_missing_from_isimip3_climate.png') ;
    map_YX = double(gridlist.mask_YX & ~climate_gridlist.mask_YX) ;
    map_YX(~gridlist.mask_YX) = NaN ;
    remap_shademap(map_YX, ...
        fig_title, fig_outfile) ;
    
    [C, IA] = setdiff(gridlist.list2map, climate_gridlist.list2map) ;
    missing_lonlats = gridlist.lonlats(IA, :) ;
    writematrix(missing_lonlats, fullfile(out_dir, 'missing_climate.csv'))
end


%% Import land uses

if strcmp(lu_source, 'LUH2')
    [out_lu, carea_YX] = remap_import_lu_luh2( ...
        geodata_dir, yearList_out, gridlist) ;
elseif strcmp(lu_source, 'HILDA+')
    [out_lu, carea_YX] = remap_import_lu_hildaplus( ...
        geodata_dir, yearList_out, gridlist) ;
else
    error('lu_source %s not recognized', lu_source)
end

% Note cells with no vegetated land according to LU dataset
remap_import_lu_note_noveg(out_lu)

% Force all land cells to sum to 1
% (Doing this here instead of a function avoids memory cost if this is unnecessary.)
lu_out_x1y = sum(out_lu.garr_xvy,2) ;
j = 0 ;
Nlu_out = length(out_lu.varNames) ;
while any(any(abs(lu_out_x1y-1)>1e-6))
    j = j + 1;
    if j > 50
        error('Possible infinite loop in "remap_import_lu_force_sum1()".')
    end
    out_lu.garr_xvy = out_lu.garr_xvy ./ repmat(lu_out_x1y,[1 Nlu_out 1]) ;
    lu_out_x1y = sum(out_lu.garr_xvy,2) ;
end
clear j lu_out_x1y



%% Import crop fractions and process crop types

disp('Importing crop fractions...')

% Import from MIRCA.txt.maps.mat
file_cropmirca = fullfile(geodata_dir, 'MIRCA', 'harvested_area_grids_26crops_30mn', 'MIRCA.txt') ;
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
nfert_dir = fullfile(geodata_dir, 'AgGRID_nutrient_input_v1.1') ;
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

    % Note gridlist cells missing from fertilizer
    isbad_x = isnan(in_YX(gridlist.list2map)) ;
    if any(isbad_x)
        I_bad = find(isbad_x) ;
        warning('%d gridlist cells are missing from fertilizer for %s', ...
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
mid_cropfrac.garr_xv = croparea_mid.garr_xv ./ repmat(sum(croparea_mid.garr_xv,2),[1 Ncrops_out]) ;
mid_cropfrac.varNames = list_crops_out ;
mid_cropfrac.list2map = gridlist.list2map ;
mid_cropfrac.lonlats = gridlist.lonlats ;

disp('Done.')


%% Interpolate land use inputs to gridlist

disp('Interpolating to match gridlist...')

% Make figures illustrating what will be interpolated
isbad_YX = remap_map_missing(out_lu.garr_xvy, gridlist, 'LU', out_dir_figs) ;
if any(isbad_YX(:))
    error('You also need to interpolate land use')
end
remap_map_missing(mid_cropfrac.garr_xv, gridlist, 'cropfrac', out_dir_figs) ;
remap_map_missing(mid_nfert.garr_xv, gridlist, 'nfert', out_dir_figs) ;

% Set up output structures
out_cropfrac = mid_cropfrac ;
out_cropfrac.garr_xv = nan(Ncells,Ncrops_out) ;
out_nfert = mid_cropfrac ;
out_nfert.garr_xv = nan(Ncells,Ncrops_out) ;

% For checking how much interpolation increases global area/nfert of each LU relative to
% what it would have been if, instead of interpolating, we set area/nfert there to zero.
test_croparea_x = mean(out_lu.garr_xvy(:, strcmp(out_lu.varNames, 'CROPLAND'), ...
    I_yrs_in_test_period),3) ...
    .* carea_YX(gridlist.list2map) ;

% Interpolate
for c = 1:Ncrops_out
    thisCrop = list_crops_out{c} ;
    fprintf('Interpolating %s (%d of %d)...\n',thisCrop, c, Ncrops_out)

    % cropfrac
    tmp0_YX = lpjgu_xz_to_YXz(mid_cropfrac.garr_xv(:,c), size(gridlist.mask_YX), gridlist.list2map) ;
    tmp1_YX = inpaint_nans(tmp0_YX, inpaint_method) ;
    tmp0_YX(isnan(tmp0_YX)) = 0 ;
    area0 = sum(test_croparea_x .* tmp0_YX(gridlist.list2map)) ;
    area1 = sum(test_croparea_x .* tmp1_YX(gridlist.list2map)) ;
    areaDiff = area1 - area0 ;
    areaDiff_pct = areaDiff / area0 * 100 ;
    fprintf('    Area increased %0.2f%% (%0.1g ha/yr) in %s\n', ...
        areaDiff_pct, areaDiff*100, interp_test_period_str) ;
    out_cropfrac.garr_xv(:,c) = tmp1_YX(gridlist.list2map) ;
    if any(isnan(out_cropfrac.garr_xv(:,c)))
        error('NaN remaining in out_cropfrac.garr_xv(:,c)!')
    end
    clear tmp*

    % nfert
    c2 = find(strcmp(out_nfert.varNames, thisCrop)) ;
    if length(c2) == 1
        c3 = find(strcmp(mid_nfert.varNames, thisCrop)) ;
        if length(c3) == 1
            tmp0_YX = lpjgu_xz_to_YXz(mid_nfert.garr_xv(:,c3), size(gridlist.mask_YX), gridlist.list2map) ;
        elseif isempty(c3)
            tmp0_YX = zeros(size(gridlist.mask_YX)) ;
        else
            error('%d matches found in mid_nfert.varNames for %s', length(c2), thisCrop)
        end
        tmp1_YX = inpaint_nans(tmp0_YX, inpaint_method) ;
        tmp0_YX(isnan(tmp0_YX)) = 0 ;
        % convert kg/m2 to t/km2 to t
        nfert0 = nansum(test_croparea_x*1e6 .* out_cropfrac.garr_xv(:,c2) .* tmp0_YX(gridlist.list2map)*1e-3) ;
        nfert1 = sum(test_croparea_x*1e6 .* out_cropfrac.garr_xv(:,c2) .* tmp1_YX(gridlist.list2map)*1e-3) ;
        nfertDiff = nfert1 - nfert0 ; 
        nfertDiff_pct = nfertDiff / nfert0 * 100 ;
        fprintf('    Nfert increased %0.2f%% (%0.1g t/yr) in %s\n', ...
            nfertDiff_pct, nfertDiff, interp_test_period_str) ;
        out_nfert.garr_xv(:,c2) = tmp1_YX(gridlist.list2map) ;
        if any(isnan(out_nfert.garr_xv(:,c2)))
            error('NaN remaining in out_nfert.garr_xv(:,c2)!')
        end
        clear tmp*
    else
        error('%d matches found in out_nfert.varNames for %s', length(c2), thisCrop)
    end
        
end

% % % disp('Interpolating include_frac...')
% % % map_include_frac_interpd_YX = inpaint_nans(map_include_frac_YX,inpaint_method) ;
% % % disp('Done.')

% Set minimum of zero, if needed
if any(out_cropfrac.garr_xv(:)<0)
    warning('Setting negative members of out_cropfrac.garr_xv to zero.')
    out_cropfrac.garr_xv(out_cropfrac.garr_xv<0) = 0 ;
end
if any(out_nfert.garr_xv(:)<0)
    warning('Setting negative members of out_nfert.garr_xv to zero.')
    out_nfert.garr_xv(out_nfert.garr_xv<0) = 0 ;
end

% Normalize cropfrac to 1, if needed
tmp = sum(out_cropfrac.garr_xv,2) ;
if any(abs(tmp(:) - 1)>1e-6)
    warning('Normalizing cropfrac sums to 1')
    out_cropfrac.garr_xv = out_cropfrac.garr_xv ./ repmat(tmp,[1 Ncrops_out]) ;
end
clear tmp

disp('Done.')


%% Set up for save

% Check
if any(any(isnan(out_cropfrac.garr_xv)))
    error('NaN in out_cropfrac')
elseif any(any(out_cropfrac.garr_xv < 0))
    error('Negative out_cropfrac')
elseif any(any(isnan(out_nfert.garr_xv)))
    error('NaN in out_nfert')
elseif any(any(out_nfert.garr_xv < 0))
    error('Negative out_nfert')
end


%% Get headers
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
out_file_lu = fullfile(out_dir, ['LU.remapv' remapVer '.txt']) ;
out_file_cropfrac = fullfile(out_dir, ['cropfracs.remapv' remapVer '.txt']) ;
out_file_nfert = fullfile(out_dir, ['nfert.remapv' remapVer '.txt']) ;

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
outFile_gridlist = fullfile(out_dir, sprintf('gridlist.remapv%s.txt', remapVer)) ;
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


%% Interpolate and save soil

file_soil = fullfile(lpjg_inputs_dir, 'soil', 'WISE', 'soilmap_center_interpolated.dat') ;
soil = lpjgu_matlab_read2geoArray(file_soil, ...
    'verboseIfNoMat', false, 'force_mat_save', false, 'force_mat_nosave', true) ; 

out_soil.varNames = soil.varNames ;
out_soil.lonlats = gridlist.lonlats ;
out_soil.list2map = gridlist.list2map ;
Nvars = length(soil.varNames) ;
out_soil.garr_xv = nan(Ncells, Nvars) ;
for v = 1:Nvars
    fprintf('Interpolating soil %s...\n', soil.varNames{v})
    tmp0_YX = lpjgu_xz_to_YXz(soil.garr_xv(:,v), size(gridlist.mask_YX), soil.list2map) ;
    tmp1_YX = inpaint_nans(tmp0_YX, inpaint_method) ;
    out_soil.garr_xv(:,v) = tmp1_YX(gridlist.list2map) ;
end
disp('Done.')

if any(any(isnan(out_soil.garr_xv)))
    error('NaN in out_soil.garr_xv')
elseif any(any(out_soil.garr_xv < 0))
    error('Negative in out_soil.garr_xv')
end

out_soil_header_cell = [ ...
    {'lon', 'lat'}, ...
    out_soil.varNames] ;

[~, soil_filename, soil_ext] = fileparts(file_soil) ;
out_file_soil = fullfile(out_dir, sprintf('%s.remapv%s%s', soil_filename, remapVer, soil_ext)) ;

disp('Saving soil...')
lpjgu_matlab_saveTable(out_soil_header_cell, out_soil, out_file_soil,...
    'outPrec', 3, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20, ...
    'gzip', do_gzip) ;
