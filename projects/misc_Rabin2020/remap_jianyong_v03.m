%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Re-map area/fert data to PLUM crops, and generate extra LU file %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Sam Rabin, 2019-02-02
% Uses HYDE LU, MIRCA crop fractions, AgMIP+ZhangManure fertilization.
% Cells are excluded if not in HYDE, MIRCA, or countries of interest.
% Included cells that have no value in ZhangManure, or that MIRCA/AgMIP
% did not assign any cropland, are interpolated.
%
% See previous notes for remap_jianyong_v02, with the following changes:
% - PLUMsetAside_frac set to zero.
% - Fix for limiting interpolated Nfert values.
% - Now at 5-minute resolution.

% Version for crop mappings
thisVer = 'jianyong01b' ;

% Version for output dir/files
remapVer = 'jianyong03' ;

yearList = 1850:2015 ;

PLUMsetAside_frac = 0 ;


%% Part 0: Setup

force_all_rainfed = false ;

inpaint_method = 4 ;
out_dir = sprintf('/Users/Shared/PLUM/input/remaps_%s/',remapVer) ;

warning('on','all')

Nyears = length(yearList) ;

% MIRCA crops, in correct order
cropList = {'Wheat'; 'Maize'; 'Rice'; 'Barley'; 'Rye'; 'Millet'; ...
    'Sorghum'; 'Soybeans'; 'Sunflower'; 'Potatoes'; 'Cassava'; ...
    'Sugar cane'; 'Sugar beet'; 'Oil palm'; 'Rape seed / Canola'; ...
    'Groundnuts / Peanuts'; 'Pulses'; 'Citrus'; 'Date palm'; ...
    'Grapes / Vine'; 'Cotton'; 'Cocoa'; 'Coffee'; ...
    'Others perennial'; 'Fodder grasses'; 'Others annual'} ;
Ncrops = length(cropList) ;
mirca_dir = '/Users/sam/Geodata/MIRCA/harvested_area_grids_26crops_5mn' ;

% Ancillary MIRCA crop stuff
list_crops_frac_in = [strcat(cropList,'_RF') strcat(cropList,'_IR')] ;
list_cropsCombined_frac_in = cropList ;
Ncrops_frac_in = length(list_crops_frac_in) ;
NcropsCombined_frac_in = length(list_cropsCombined_frac_in) ;

% Maybe this fixes something?
list_crops_frac_in = strrep(list_crops_frac_in,'Soybeans','Soybean') ;
list_cropsCombined_frac_in = strrep(list_cropsCombined_frac_in,'Soybeans','Soybean') ;


%% Part 1: Import data

disp('Importing ancillary data...')

addpath(genpath(landsymm_lpjg_path()))

% Import HYDE land area
inDir_hyde = '/Users/sam/Geodata/HYDE_3.2.1_lu' ;
cd(inDir_hyde)
hyde_files = dir([inDir_hyde '/*.zip']) ;
hyde_years = str2num(char(strrep({hyde_files.name}','AD_lu.zip',''))) ; %#ok<ST2NM>
maxln_cr = dlmread([inDir_hyde '/../HYDE_3.2.1/general_files/maxln_cr.asc'],' ',6,0) ;
maxln_cr = flipud(maxln_cr(:,1:4320)) ;
maxln_cr(maxln_cr==-9999) = NaN ;
garea_cr = dlmread([inDir_hyde '/../HYDE_3.2.1/general_files/garea_cr.asc'],' ',6,0) ;
garea_cr = flipud(garea_cr(:,1:4320)) ;
garea_cr(garea_cr==-9999) = NaN ;
hyde_ok_YX = maxln_cr>0 ;

% Import countries at 5-minute resolution
countries_YX = flipud(imread('/Users/Shared/PLUM/crop_calib_data/countries/ne_10m_admin_0_countries_ssrIDs_5min.tif')) ;
id_eth = 74 ;
id_ken = 121 ;
countries_ok_YX = countries_YX==id_eth | countries_YX==id_ken ;

% Import AgMIP fertilizer (for mask)
nfert_in = lpjgu_matlab_readTable_then2map('/Users/Shared/unifying_gridlist/AgGRID_nutrient_input_v1.1/AgMIP_Nfert.txt') ;
this_YX = nfert_in.maps_YXv(:,:,1) ;
tmp_5min_YX = nan(size(this_YX,1),size(countries_YX,2)) ;
for i = 1:6
    tmp_5min_YX(:,i:6:size(countries_YX,2)) = this_YX ;
end
clear this_YX
this_5min_YX = nan(size(countries_YX)) ;
for i = 1:6
    this_5min_YX(i:6:size(countries_YX,1),:) = tmp_5min_YX ;
end
clear tmp_5min_YX
nfert_ok_YX = ~isnan(this_5min_YX) ;

% Get mask
ok_YX = hyde_ok_YX & countries_ok_YX & nfert_ok_YX ;
i_ok = find(ok_YX) ;
Ncells = length(find(ok_YX)) ;
res = 5/60 ;
lat_map = repmat(transpose((-90+res/2):res:(90-res/2)),[1 size(ok_YX,2)]) ;
lon_map = repmat((-180+res/2):res:(180-res/2),[size(ok_YX,1) 1]) ;
lats = lat_map(ok_YX) ;
lons = lon_map(ok_YX) ;
[I,J] = find(ok_YX) ;
Is_ok = min(I):max(I) ;
Js_ok = min(J):max(J) ;
ok_inset_YX = ok_YX(Is_ok,Js_ok) ;
i_ok_inset = find(ok_inset_YX) ;
[I,J] = find(ok_inset_YX) ;
Is_ok_inset = min(I):max(I) ;
Js_ok_inset = min(J):max(J) ;

% Import land use (pre-transposed for later vectorization)
disp('Importing land use...')
garea_x = garea_cr(i_ok) ;
garea_yx = repmat(garea_x',[Nyears 1]) ;
hyde_crop_yx = get_HYDE_LUs( ...
    hyde_years, yearList, inDir_hyde, hyde_files, ...
    'cropland', i_ok) ./ garea_yx ;
hyde_graz_yx = get_HYDE_LUs( ...
    hyde_years, yearList, inDir_hyde, hyde_files, ...
    'grazing', i_ok) ./ garea_yx ;
hyde_bare_yx = transpose( repmat( (garea_x - maxln_cr(i_ok)) ./ garea_x, ...
                                  [1 Nyears])) ;
hyde_ntrl_yx = 1 - (hyde_crop_yx + hyde_graz_yx + hyde_bare_yx) ;

% Import crop fractions
disp('Importing crop areas...')
croparea_rf_xc = nan(Ncells,Ncrops) ;
croparea_ir_xc = nan(Ncells,Ncrops) ;
for c = 1:Ncrops
    thisCrop = cropList{c} ;
    fprintf('   %s (%d of %d)...\n', thisCrop, c, Ncrops) ;
    % Rainfed
    filename_mat = sprintf('%s/ANNUAL_AREA_HARVESTED_RFC_CROP%d_HA.mat', mirca_dir, c) ;
    if exist(filename_mat,'file')
        load(filename_mat) ;
    else
        filename_this = sprintf('%s/ANNUAL_AREA_HARVESTED_RFC_CROP%d_HA.ASC', mirca_dir, c) ;
        filename_thisGZ = [filename_this '.gz'] ;
        unix(sprintf('gunzip < %s > %s', filename_thisGZ, filename_this)) ;
        this_YX = flipud(dlmread(filename_this,' ',6,0)) ;
        save(filename_mat,'this_YX') ;
        unix(sprintf('rm %s',filename_this)) ;
    end
    croparea_rf_xc(:,c) = this_YX(i_ok) ;
    clear this_YX
    
    % Irrigated
    filename_mat = sprintf('%s/ANNUAL_AREA_HARVESTED_IRC_CROP%d_HA.mat', mirca_dir, c) ;
    if exist(filename_mat,'file')
        load(filename_mat) ;
    else
        filename_this = sprintf('%s/ANNUAL_AREA_HARVESTED_IRC_CROP%d_HA.ASC', mirca_dir, c) ;
        filename_thisGZ = [filename_this '.gz'] ;
        unix(sprintf('gunzip < %s > %s', filename_thisGZ, filename_this)) ;
        this_YX = flipud(dlmread(filename_this,' ',6,0)) ;
        save(filename_mat,'this_YX') ;
        unix(sprintf('rm %s',filename_this)) ;
    end
    croparea_ir_xc(:,c) = this_YX(i_ok) ;
    clear this_YX
end
disp('Converting to fractions...')
croparea_x = sum(croparea_rf_xc+croparea_ir_xc,2) ;
croparea_xc = repmat(croparea_x,[1 Ncrops]) ;
cropfrac_rf_xc = croparea_rf_xc ./ croparea_xc ;
cropfrac_ir_xc = croparea_ir_xc ./ croparea_xc ;
cropfrac_rf_xc(croparea_xc==0) = 0 ;
cropfrac_ir_xc(croparea_xc==0) = 0 ;
if any(isnan(cropfrac_rf_xc))
    error('NaN in cropfrac_rf_xc')
elseif any(isnan(cropfrac_ir_xc))
    error('NaN in cropfrac_ir_xc')
elseif any(isinf(cropfrac_rf_xc))
    error('Inf in cropfrac_rf_xc')
elseif any(isinf(cropfrac_ir_xc))
    error('Inf in cropfrac_ir_xc')
end
cropfrac_xc = [cropfrac_rf_xc cropfrac_ir_xc] ;

% Import fertilizer
disp('Importing fertilizer...')
% Fertilizer (convert kg/ha to kg/m2)
% nfert_in = lpjgu_matlab_readTable_then2map('/Users/Shared/unifying_gridlist/AgGRID_nutrient_input_v1.1/AgMIP_Nfert.txt') ;
nfert_in.maps_YXv = nfert_in.maps_YXv*1e-4 ;
nfert_in.maps_YXv = nfert_in.maps_YXv * 1/(1-PLUMsetAside_frac) ;
list_cropsCombined_fert_in = nfert_in.varNames ;
NcropsCombined_fert_in = length(list_cropsCombined_fert_in) ;
% Convert AgMIP-style names to pure MIRCA names
list_cropsCombined_fert_in(strcmp(list_cropsCombined_fert_in, ...
    'GroundnutsPeanuts')) = {'Groundnuts / Peanuts'} ;
list_cropsCombined_fert_in(strcmp(list_cropsCombined_fert_in, ...
    'Oilpalm')) = {'Oil palm'} ;
list_cropsCombined_fert_in(strcmp(list_cropsCombined_fert_in, ...
    'Sugarcane')) = {'Sugar cane'} ;
list_cropsCombined_fert_in(strcmp(list_cropsCombined_fert_in, ...
    'Sugarbeet')) = {'Sugar beet'} ;
list_cropsCombined_fert_in(strcmp(list_cropsCombined_fert_in, ...
    'RapeseedCanola')) = {'Rape seed / Canola'} ;
nfert_xc = nan(Ncells,NcropsCombined_fert_in) ;
for c = 1:NcropsCombined_fert_in
    thisCrop = list_cropsCombined_fert_in{c} ;
    this_YX = nfert_in.maps_YXv(:,:,c) ;
    % Kludge from 0.5 degree (30-minute) to 5 minute resolution
    tmp_5min_YX = nan(size(this_YX,1),size(countries_YX,2)) ;
    for i = 1:6
        tmp_5min_YX(:,i:6:size(countries_YX,2)) = this_YX ;
    end
    clear this_YX
    this_5min_YX = nan(size(countries_YX)) ;
    for i = 1:6
        this_5min_YX(i:6:size(countries_YX,1),:) = tmp_5min_YX ;
    end
    nfert_xc(:,c) = this_5min_YX(i_ok) ;
    if any(isnan(nfert_xc(:,c)))
        warning('%d NaN in nfert for %s',length(find(isnan(nfert_xc(:,c)))),thisCrop);
    end
    clear tmp_5min_YX this_5min_YX
end
if any(isnan(nfert_xc))
    error('NaN in nfert_xc')
elseif any(isinf(nfert_xc))
    error('Inf in nfert_xc')
end

% Add manure N for year 2000
disp('Importing manure for year 2000...')
load('/Users/sam/Geodata/Manure_ZhangEtAl2017/zhangManure_1860to2014_agg_hd.mat') ;
manure2crop_hd_YX = manure2crop_hd_YXy(:,:,1860:2014==2000) ;
manure2crop_hd_YX(manure2crop_hd_YX<1e-6) = 0 ;
clear manure2crop_hd_YXy
% Kludge from 0.5 degree (30-minute) to 5 minute resolution
tmp_5min_YX = nan(size(manure2crop_hd_YX,1),size(countries_YX,2)) ;
for i = 1:6
    tmp_5min_YX(:,i:6:size(countries_YX,2)) = manure2crop_hd_YX ;
end
clear this_YX
manure2crop_5min_YX = nan(size(countries_YX)) ;
for i = 1:6
    manure2crop_5min_YX(i:6:size(countries_YX,1),:) = tmp_5min_YX ;
end
% Interpolate
disp('Interpolating to fill gaps...')
tmp = manure2crop_5min_YX(Is_ok,Js_ok) ;
tmp = inpaint_nans(tmp,inpaint_method) ;
manure2crop_5min_YX(Is_ok,Js_ok) = tmp ;
manure2crop_5min_YX(~ok_YX) = NaN ;
if any(any(isnan(manure2crop_5min_YX) & ok_YX))
    error('NaN remaining in area of interest in manure2crop_5min_YX.')
elseif any(any(manure2crop_5min_YX<0))
    error('Negative value(s) in manure2crop_5min_YX.')
end
% Make array and add to nfert
nfert2_xc = nfert_xc + repmat(manure2crop_5min_YX(ok_YX),[1 NcropsCombined_fert_in]) ;

disp('Done.')


%% Part 2: Set up crop type mapping

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
list_cropsCombined_fert_in(strcmp(list_cropsCombined_fert_in, ...
    'Soybeans')) = {'Soybean'} ;
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
        thisRow = {'Groundnuts / Peanuts','Soybean'} ;
    else
        [~,IA] = intersect(thisRow,list_crops_inFrac_notFert) ;
        thisRow(IA) = [] ;
    end
    if isempty(thisRow)
        error('No fertilizer found for this keymap!')
    end
    in2out_keyCombined_fert{c} = thisRow ;
end

% Get ignored fraction
[~,I] = intersect(list_crops_frac_in',[strcat(list_ignore_frac,'_RF');strcat(list_ignore_frac,'_IR')]) ;
ignore_frac_x = sum(cropfrac_xc(:,I),2) ;


%% Part 3: Do mapping

disp('Mapping crop fracs...')
cropfrac_mid_xc = nan(Ncells,Ncrops_out) ;
for c = 1:Ncrops_out
    thisCrop = list_crops_out{c} ;
    thisRow = in2out_key_frac{getOi(thisCrop)} ;
    [~,IA] = intersect(list_crops_frac_in,thisRow) ;
    if length(IA)~=length(thisRow)
        error('length(IA)~=length(thisRow)')
    end
    % Straight sum
    cropfrac_mid_xc(:,c) = sum(cropfrac_xc(:,IA),2) ;
end
disp('Done.')

nfert_mid_xc = nan(Ncells,NcropsCombined_out) ;
%%% Rainfed
for c = 1:NcropsCombined_out
%     thisCrop = list_crops_out{c} ;
    thisCrop = list_cropsCombined_out{c} ;
    thisRow = in2out_keyCombined_fert{getOi(thisCrop)} ;
    [~,IA_frac] = intersect(list_cropsCombined_frac_in,thisRow) ;
    [~,IA_fert] = intersect(list_cropsCombined_fert_in,thisRow) ;
    if length(IA_frac)~=length(thisRow)
        error('length(IA_frac)~=length(thisRow)')
    end
    if length(IA_fert)~=length(thisRow)
        error('length(IA_fert)~=length(thisRow)')
    end
    % Weighted average
    nfert_mid_xc(:,c) = sum(nfert2_xc(:,IA_fert) .* cropfrac_xc(:,IA_frac) ./ repmat(sum(cropfrac_xc(:,IA_frac),2),[1 length(IA_frac)]),2) ;
end
%%% Tack on irrigated fertilization (same as rainfed)
nfert_mid_xc = cat(2,nfert_mid_xc,nfert_mid_xc) ;

% Add ExtraCrop to represent unhandled crops and setAside
% Move ignored area into ExtraCrop.
% Set up croparea_extra array
nfert0_x = zeros(size(ignore_frac_x)) ;
% Assign based on pastNPP results
setaside_area_x = sum(cropfrac_mid_xc * PLUMsetAside_frac,2) ;
cropfrac_mid_xc = cropfrac_mid_xc * (1-PLUMsetAside_frac) ;
cropfrac_extra_x = ignore_frac_x + setaside_area_x ;
% Save new croparea and nfert structures
list_cropsCombined_out = [list_cropsCombined_out {'ExtraCrop'}] ;
NcropsCombined_out = length(list_cropsCombined_out) ;
list_crops_out = [list_crops_out {'ExtraCrop'}]  ;
Ncrops_out = length(list_crops_out) ;
cropfrac_mid_xc = cat(2,cropfrac_mid_xc,cropfrac_extra_x) ;
nfert_mid_xc = cat(2,nfert_mid_xc,nfert0_x) ;


%% Part 4: Interpolate

% Interpolate cropfrac and nfert so that you have something there even if
% MIRCA/AgMIP thought there was no crop but HYDE does.

tmp = sum(cropfrac_mid_xc,2) ;
cropfrac_mid_xc(repmat(tmp,[1 size(cropfrac_mid_xc,2)])==0) = NaN ;
clear tmp
cropfrac_out_xc = nan(size(cropfrac_mid_xc)) ;
nfert_out_xc = nan(size(nfert_mid_xc)) ;
[I,J] = find(ok_YX) ;
Is = min(I):max(I) ;
Js = min(J):max(J) ;
for c = 1:Ncrops_out
    thisCrop = list_crops_out{c} ;
    disp(['Interpolating ' thisCrop ' (' num2str(c) ' of ' num2str(Ncrops_out) ')...'])
    
    % Crop fraction
    tmp = nan(size(ok_inset_YX)) ;
    tmp(i_ok_inset) = cropfrac_mid_xc(:,c) ;
    tmp = inpaint_nans(tmp,inpaint_method) ;
    cropfrac_out_xc(:,c) = tmp(i_ok_inset) ;
    if any(isnan(cropfrac_out_xc(:,c)))
        error('NaN remaining in cropfrac_out_xc(:,c)!')
    elseif any(cropfrac_out_xc(:,c)<0)
        error('Negative in cropfrac_out_xc(:,c)!')
    end
    
    thisNonIrr = find(strcmp(list_crops_out,thisCrop(1:end-1))) ;
    if length(thisNonIrr)==1 && strcmp([list_crops_out{thisNonIrr} 'i'],thisCrop)
        warning(['Nfert: Using ' list_crops_out{thisNonIrr} ' for ' list_crops_out{c}])
        if any(isnan(nfert_out_xc(:,thisNonIrr)))
            error('But %s has not been interpolated yet!', list_crops_out{thisNonIrr}) ;
        end
        nfert_out_xc(:,c) = nfert_out_xc(:,thisNonIrr) ;
    else        
        tmp = nan(size(ok_inset_YX)) ;
        tmp(i_ok_inset) = nfert_mid_xc(:,c) ;
        tmp = inpaint_nans(tmp,inpaint_method) ;
        nfert_out_xc(:,c) = tmp(i_ok_inset) ;
    end
    if any(isnan(nfert_out_xc(:,c)))
        error('NaN remaining in nfert_out_xc(:,c)!')
    elseif any(nfert_out_xc(:,c)<0)
        error('Negative in nfert_out_xc(:,c)!')
    end
end

% Set minimum of zero, if needed
if any(any(cropfrac_out_xc<0))
    warning('Setting negative members of cropfrac_out_xc to zero.')
    cropfrac_out_xc(cropfrac_out_xc<0) = 0 ;
end
if any(any(nfert_out_xc<0))
    warning('Setting negative members of nfert_out_xc to zero.')
    nfert_out_xc(nfert_out_xc<0) = 0 ;
end

% Normalize cropfrac to 1, if needed
tmp = sum(cropfrac_out_xc,2) ;
if any(abs(tmp - 1)>1e-6)
    disp('Normalizing cropfrac to 1...')
    cropfrac_out_xc = cropfrac_out_xc ./ repmat(tmp,[1 Ncrops_out]) ;
    tmp = sum(cropfrac_out_xc,2) ;
end ; clear tmp
if any(any(isnan(cropfrac_out_xc)))
    error('Normalizing cropfrac to 1 introduced NaN!')
elseif any(any(cropfrac_out_xc<0))
    error('Normalizing cropfrac to 1 introduced negative(s)!')
end

% Don't let extrapolated Nfert exceed maximum Nfert seen for this crop
for c = 1:NcropsCombined_out-1 % -1 to skip ExtraCrop
    thisMax = max(nfert_mid_xc(:,c)) ;
    thisMap = nfert_out_xc(:,c) ;
    Nexceeded = length(find(thisMap>thisMax)) ;
    if Nexceeded>0
        warning(['Limiting nfert_out for ' list_cropsCombined_out{c} '.'])
        thisMap(thisMap>thisMax) = thisMax ;
        nfert_out_xc(:,c) = thisMap ; % Rainfed
        nfert_out_xc(:,c+NcropsCombined_out) = thisMap ; % Irrigated
    end
end

disp('Done interpolating.')


%% Part 5: Set up for save

lats_yx = repmat(lats',[Nyears 1]) ;
lons_yx = repmat(lons',[Nyears 1]) ;
years_yx = repmat(yearList',[1 Ncells]) ;

disp('Array-ifying lu_in...')
lu_in_array = [lons_yx(:) lats_yx(:) years_yx(:) hyde_crop_yx(:) hyde_graz_yx(:) hyde_ntrl_yx(:) hyde_bare_yx(:)] ;
lu_in_header_cell = {'Lon','Lat','Year','CROPLAND','PASTURE','NATURAL','BARREN'} ;

disp('Array-ifying cropfrac_out...')
cropfrac_out_array = [lons lats cropfrac_out_xc] ;
cropfrac_out_header_cell = [{'Lon','Lat'} list_crops_out] ;

disp('Array-ifying nfert_out...')
nfert_out_array = [lons lats nfert_out_xc] ;
nfert_out_header_cell = [{'Lon','Lat'} list_crops_out] ;
disp('Done.')

% Get filenames
if force_all_rainfed
    allRF_txt = '.noIrr' ;
else
    allRF_txt = '' ;
end
out_file_lu = [out_dir 'LU.remap.' remapVer '.txt'] ;
out_file_cropfrac = [out_dir 'cropfracs.remap.' remapVer '.txt'] ;
out_file_nfert = [out_dir 'nfert.remap.' remapVer '.txt'] ;


%% Part 6: Save

%%% Options %%%%%%%%%%%%%
outPrec = 6 ;
outPrec_lonlat = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(out_dir,'dir')
    mkdir(out_dir) ;
end

% % disp('Saving ignored crop fraction...')
% % out_file_cropfrac_ignored = strrep(out_file_cropfrac,'cropfracs','cropfracsIGNORED') ;
% % out_file_cropfrac_ignored = strrep(out_fi le_cropfrac_ignored,'.txt','.mat') ;
% % save(out_file_cropfrac_ignored,'ignored_LUCROParea_YX_y','-v7.3') ;

disp('Saving LU...')
check_existing_lu(thisVer, out_file_lu, allVer_names, allVer_ignore_types) ;
lpjgu_matlab_saveTable(lu_in_header_cell, lu_in_array, out_file_lu,...
    'outPrec', outPrec, ...
    'outPrec_lonlat', outPrec_lonlat, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 1) ;

disp('Saving cropfracs...')
lpjgu_matlab_saveTable(cropfrac_out_header_cell, cropfrac_out_array, out_file_cropfrac,...
    'outPrec', outPrec, ...
    'outPrec_lonlat', outPrec_lonlat, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20) ;

disp('Saving nfert...')
lpjgu_matlab_saveTable(nfert_out_header_cell, nfert_out_array, out_file_nfert,...
    'outPrec', outPrec, ...
    'outPrec_lonlat', outPrec_lonlat, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20) ;




