%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Re-map area/fert data to PLUM crops, and generate extra LU file %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Sam Rabin, 2018-11-28
% See previous notes for remap_v6p2, with the following changes:
% - Including entire world
% - Not adjusting Nfert of gridcells that had more than 0 but less than
%   minimum

PLUMsetAside_frac = 0.103 ;
inpaint_method = 4 ;
yearList_out = 1850:2015 ;
% yearList_out = 1995:2005 ;

% Version for crop mappings
thisVer = '20180214' ;

force_all_rainfed = false ;

out_dir = '/Users/Shared/PLUM/input/remaps_v6p2/' ;


%% Part 0: Setup

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

% Get output gridlist
gridlist = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/gridlists/gridlist_62892.runAEclimOK.txt');%, 'verboseIfNoMat',true) ;

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


%% Part 1: Import land uses

disp('Importing land uses...')

% Import cell area (km2); aggregate to half-degree
file_luh2_etc = '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc' ;
carea_XY = ncread(file_luh2_etc,'carea') ;
carea_XYy = repmat(carea_XY,[1 1 Nyears_out]) ;
carea_hd_XY = aggregate(carea_XY,0.25,0.5) ;
carea_hd_XYy = repmat(carea_hd_XY,[1 1 Nyears_out]) ;

% Land use
lu_out_XYyv = zeros([size(carea_hd_XY) Nyears_out Nlu_out]) ;
for v = 1:Nlu_in
    thisLU_in = list_LU_in{v} ;
    thisLU_out = map_LU_in2out{v} ;
    fprintf('   %s (%d of %d) to %s...\n',thisLU_in,v,Nlu_in,thisLU_out)
    i = strcmp(list_LU_out,thisLU_out) ;
    lu_in_XYy = ncread(file_luh2_states,thisLU_in,starts_luh2_states,counts_luh2_states) ;
    lu_in_XYy(isnan(lu_in_XYy)) = 0 ;
    lu_out_XYy = aggregate(lu_in_XYy.*carea_XYy,0.25,0.5)./carea_hd_XYy ;
    if any(any(any(lu_out_XYy-1 > 1e-6)))
        error('Some element(s) of lu_out_XYy > 1!')
    end
    clear lu_in_XYy
    lu_out_XYyv(:,:,:,i) = lu_out_XYyv(:,:,:,i) + lu_out_XYy ;
    if any(any(any(any(lu_out_XYyv-1 > 1e-6))))
        error('Some element(s) of lu_out_XYyv > 1!')
    end
    clear lu_out_XYy
end
clear carea*_YXy

disp('Finishing...')

% Add water fraction to BARREN
icwtr_XY = ncread(file_luh2_etc,'icwtr') ;
icwtr_XY(icwtr_XY==1) = 0 ;
i = strcmp(list_LU_out,'BARREN') ;
lu_out_XYyv(:,:,:,i) = lu_out_XYyv(:,:,:,i) ...
    + repmat(aggregate(icwtr_XY.*carea_XY,0.25,0.5)./carea_hd_XY,[1 1 Nyears_out]) ;

% Shift to correct orientation
out_lu.varNames = list_LU_out ;
out_lu.yearList = shiftdim(yearList_out) ;
out_lu.maps_YXyv = flip(permute(lu_out_XYyv,[2 1 3 4]),1) ;
clear lu_out_XYyv

% Mask cells with no vegetated land
bad_YX = sum(out_lu.maps_YXyv(:,:,1,:),4)==0 ;
out_lu.maps_YXyv(repmat(bad_YX,[1 1 Nyears_out Nlu_out])) = NaN ;

% Force all land cells to sum to 1
lu_out_YXySum = sum(out_lu.maps_YXyv,4) ;
j = 0 ;
while any(any(any(abs(lu_out_YXySum-1)>1e-6)))
    j = j + 1;
    if j > 50
        error('Possible infinite loop in "Force all land cells to sum to 1".')
    end
    out_lu.maps_YXyv = out_lu.maps_YXyv ./ repmat(lu_out_YXySum,[1 1 1 Nlu_out]) ;
    lu_out_YXySum = sum(out_lu.maps_YXyv,4) ;
end
disp('Done.')


%% Part 2: Import crop fractions

disp('Importing crop fractions...')

% Import from MIRCA
croparea_in = lpjgu_matlab_readTable_then2map('/Users/sam/Geodata/MIRCA/harvested_area_grids_26crops_30mn/MIRCA.txt',...
    'verboseIfNoMat',true) ;

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

% Move ignored area (unhandled crops) and setAside area into ExtraCrop
I = contains(croparea_in.varNames,list_ignore_frac) ;
map_ignore_area_YX = sum(croparea_in.maps_YXv(:,:,I),3) ;
map_setaside_area_YX = sum(croparea_mid.maps_YXv * PLUMsetAside_frac,3) ;
croparea_mid.maps_YXv = croparea_mid.maps_YXv * (1-PLUMsetAside_frac) ;
map_extra_YX = map_ignore_area_YX + map_setaside_area_YX ;
croparea_mid.maps_YXv = cat(3,croparea_mid.maps_YXv,map_extra_YX) ;
list_cropsCombined_out = [list_cropsCombined_out {'ExtraCrop'}] ;
NcropsCombined_out = length(list_cropsCombined_out) ;
list_crops_out = [list_crops_out {'ExtraCrop'}]  ;
Ncrops_out = length(list_crops_out) ;

% Get fractions
cropfrac_mid.maps_YXv = croparea_mid.maps_YXv ./ repmat(sum(croparea_mid.maps_YXv,3),[1 1 Ncrops_out]) ;

% Interpolate
out_cropfrac.varNames = list_crops_out ;
out_cropfrac.yearList = shiftdim(yearList_out) ;
out_cropfrac.maps_YXv = nan(360,720,Ncrops_out) ;
for c = 1:Ncrops_out
    thisCrop = list_crops_out{c} ;
    fprintf('Interpolating %s (%d of %d)...\n',thisCrop, c, Ncrops_out)
    out_cropfrac.maps_YXv(:,:,c) = inpaint_nans(cropfrac_mid.maps_YXv(:,:,c), inpaint_method) ;
    if any(any(isnan(out_cropfrac.maps_YXv(:,:,c))))
        error('NaN remaining in out_cropfrac.maps_YXv(:,:,c)!')
    end
end

% % % disp('Interpolating include_frac...')
% % % map_include_frac_interpd_YX = inpaint_nans(map_include_frac_YX,inpaint_method) ;
% % % disp('Done.')

% Restrict to LUH2 gridlist
out_cropfrac.maps_YXv(repmat(bad_YX,[1 1 Ncrops_out])) = NaN ;
% % % map_include_frac_interpd_YX(bad_YX) = NaN ;

% Set minimum of zero, if needed
if any(out_cropfrac.maps_YXv(:)<0)
    warning('Setting negative members of out_cropfrac.maps_YXv to zero.')
    out_cropfrac.maps_YXv(out_cropfrac.maps_YXv<0) = 0 ;
end
% % % if any(map_include_frac_interpd_YX(:)<0)
% % %     warning('Setting negative members of map_include_frac_interpd_YX to zero.')
% % %     map_include_frac_interpd_YX(map_include_frac_interpd_YX<0) = 0 ;
% % % end

% Normalize cropfrac to 1, if needed
tmp = sum(out_cropfrac.maps_YXv,3) ;
if any(abs(tmp(:) - 1)>1e-6)
    out_cropfrac.maps_YXv = out_cropfrac.maps_YXv ./ repmat(tmp,[1 1 Ncrops_out]) ;
end


disp('Done.')


%% Part 3: Import Nfert and irrigation

disp('Importing Nfert and irrigation...')

% Match input and output crop names
list_crops_in = list_LU_in(strcmp(map_LU_in2out,'CROPLAND')) ;
map_crops_out2in = cell(size(list_cropsCombined_out)) ;
map_crops_out2in(contains(list_cropsCombined_out,...
    {'CerealsC3','Rice','StarchyRoots'})) = {'c3ann'} ;
map_crops_out2in(contains(list_cropsCombined_out,...
    {'CerealsC4'})) = {'c4ann'} ;
map_crops_out2in(contains(list_cropsCombined_out,...
    {'Pulses'})) = {'c3nfx'} ;
map_crops_out2in(contains(list_cropsCombined_out,...
    {'Oilcrops'})) = {'c3oil'} ;

% Get LUH2 equivalents for original MIRCA crops
luh2_equivs_in = in2out_keyCombined_frac ;
list_crops_luh2_avail = {'c3ann','c4ann','c3per','c4per','c3nfx'} ;
Ncrops_luh2_avail = length(list_crops_luh2_avail) ;
for c = 1:(NcropsCombined_out-1)
    tmp = luh2_equivs_in{c} ;
    tmp(contains(tmp,...
        {'Wheat','Barley','Rye','Sunflower','RapeseedCanola',...
        'Potatoes','Sugarbeet','Cassava','Rice'})) ...
        = {'c3ann'} ;
    tmp(contains(tmp,...
        {'Maize','Millet','Sorghum'})) ...
        = {'c4ann'} ;
    tmp(contains(tmp,...
        {'Oilpalm'})) ...
        = {'c3per'} ;
    tmp(contains(tmp,...
        {'Soybeans','GroundnutsPeanuts','Pulses'})) ...
        = {'c3nfx'} ;
    luh2_equivs_in{c} = tmp ;
end
if any(cellfun(@(x) any(~contains(x,list_crops_luh2_avail)),luh2_equivs_in))
    error('Unrecognized crop(s) when defining LUH2 equivalents to crops_in!')
end

% Get list of LUH2 crops whose managements we actually need to import
luh2_to_import = list_crops_luh2_avail ;
for c = 1:Ncrops_luh2_avail
    thisCrop = list_crops_luh2_avail{c} ;
    if ~any(cellfun(@(x) any(strcmp(x,thisCrop)),luh2_equivs_in))
        luh2_to_import(strcmp(luh2_to_import,thisCrop)) = [] ;
    end
end

% Import LUH2 originals: area, nfert, and irrig
Ncrops_luh2_in = length(luh2_to_import) ;
luh2_carea_XYyc = nan([size(carea_hd_XY) Nyears_out Ncrops_luh2_in]) ;
luh2_nfert_XYyc = nan([size(carea_hd_XY) Nyears_out Ncrops_luh2_in]) ;
luh2_irrig_XYyc = nan([size(carea_hd_XY) Nyears_out Ncrops_luh2_in]) ;
for c = 1:length(luh2_to_import)
    
    thisVar = luh2_to_import{c} ;
    fprintf('   %s: ', thisVar) ;
    
    fprintf('area... ') ;
    luh2_carea_XYy = ncread(file_luh2_states,thisVar,starts_luh2_states,counts_luh2_states) ;
    luh2_carea_XYy = luh2_carea_XYy .* carea_XYy ;
    luh2_carea_XYyc(:,:,:,c) = aggregate(luh2_carea_XYy,0.25,0.5) ;
    
    fprintf('nfert... ') ;
    thisVar = ['fertl_' luh2_to_import{c}] ;
    tmp = ncread(file_luh2_mgmts,thisVar,starts_luh2_mgmts,counts_luh2_mgmts) ;
    % Multiply by 1/(1-PLUMsetAside_frac) because LUH2 gets good global
    % total N applied without considering setAside area, and we're assuming
    % setAside area gets no N application.
    luh2_nfert_XYyc(:,:,:,c) = ...
        (1/(1-PLUMsetAside_frac)) * aggregate_mgmt(tmp,luh2_carea_XYy,0.25,0.5) ;
    clear tmp
    
    fprintf('irrig... ') ;
    thisVar = ['irrig_' luh2_to_import{c}] ;
    tmp = ncread(file_luh2_mgmts,thisVar,starts_luh2_mgmts,counts_luh2_mgmts) ;
    luh2_irrig_XYyc(:,:,:,c) = aggregate_mgmt(tmp,luh2_carea_XYy,0.25,0.5) ;
    clear tmp
    
    fprintf('\n') ;
    clear luh2_carea_XYy
end

disp('Processing...')

luh2_carea_YXyc = flip(permute(luh2_carea_XYyc,[2 1 3 4]),1) ; clear luh2_carea_XYyc
luh2_nfert_YXyc = flip(permute(luh2_nfert_XYyc,[2 1 3 4]),1) ; clear luh2_nfert_XYyc
luh2_irrig_YXyc = flip(permute(luh2_irrig_XYyc,[2 1 3 4]),1) ; clear luh2_irrig_XYyc

% Convert to output equivalents
out_nfert.varNames = list_crops_out ;
out_irrig.varNames = list_crops_out ;
out_nfert.yearList = shiftdim(yearList_out) ;
out_irrig.yearList = shiftdim(yearList_out) ;
out_nfert.maps_YXyc = nan(360,720,Nyears_out,Ncrops_out) ;
out_irrig.maps_YXyc = nan(360,720,Nyears_out,NcropsCombined_out) ;
for c = 1:NcropsCombined_out
    thisCrop = list_cropsCombined_out{c} ;
    if strcmp(thisCrop,'ExtraCrop')
        out_nfert.maps_YXyc(:,:,:,c) = zeros(360,720,Nyears_out) ;
        out_irrig.maps_YXyc(:,:,:,c) = zeros(360,720,Nyears_out) ;
    else
        these_equivs = luh2_equivs_in{c} ;
        is_match_nfert = contains(list_crops_out,thisCrop) ;
        is_match_irrig = strcmp(list_cropsCombined_out,thisCrop) ;
        Nmatch_nfert = length(find(is_match_nfert)) ;
        Nmatch_irrig = length(find(is_match_irrig)) ;
        if length(unique(these_equivs))==1
            out_nfert.maps_YXyc(:,:,:,is_match_nfert) = repmat(luh2_nfert_YXyc(:,:,:,strcmp(luh2_to_import,these_equivs{1})),[1 1 1 Nmatch_nfert]) ;
            out_irrig.maps_YXyc(:,:,:,is_match_irrig) = repmat(luh2_irrig_YXyc(:,:,:,strcmp(luh2_to_import,these_equivs{1})),[1 1 1 Nmatch_irrig]) ;
        else
            %%% Need to weight constituent crop types
            % Setup
            i_luh2 = cellfun(@(x) find(strcmp(luh2_to_import,x)),these_equivs) ;
            weightsDenom_YXyc = repmat(sum(luh2_carea_YXyc(:,:,:,i_luh2),4),[1 1 1 length(i_luh2)]) ;
            weights_YXyc = luh2_carea_YXyc(:,:,:,i_luh2) ./ weightsDenom_YXyc ;
            if any(any(any(abs(sum(weights_YXyc,4)-1)>1e-6)))
                error('Something went wrong with weighting in converting to %s.',thisCrop) ;
            end

            % Nfert
            tmp_YXy = sum(luh2_nfert_YXyc(:,:,:,i_luh2) .* weights_YXyc,4) ;
            tmp_YXy(weightsDenom_YXyc(:,:,:,1)==0) = 0 ;
            out_nfert.maps_YXyc(:,:,:,is_match_nfert) = repmat(tmp_YXy,[1 1 1 Nmatch_nfert]) ;
            clear tmp_YXy

            % Irrigation
            tmp_YXy = sum(luh2_irrig_YXyc(:,:,:,i_luh2) .* weights_YXyc,4) ;
            tmp_YXy(weightsDenom_YXyc(:,:,:,1)==0) = 0 ;
            out_irrig.maps_YXyc(:,:,:,is_match_irrig) = repmat(tmp_YXy,[1 1 1 Nmatch_irrig]) ;
            clear tmp_YXy

            clear weights_YXyc
        end
    end
end

% Irrigated crops were not filled for Nfert
for c = 1:NcropsCombined_out
    thisCrop = list_cropsCombined_out{c} ;
    thisCropI = [thisCrop 'i'] ;
    if any(strcmp(out_nfert.varNames,thisCropI))
        out_nfert.maps_YXyc(:,:,:,strcmp(out_nfert.varNames,thisCropI)) = ...
            out_nfert.maps_YXyc(:,:,:,strcmp(out_nfert.varNames,thisCrop)) ;
    end
end

out_nfert.maps_YXyc(isnan(out_nfert.maps_YXyc) & ~repmat(bad_YX,[1 1 Nyears_out Ncrops_out])) = 0 ;
out_irrig.maps_YXyc(isnan(out_irrig.maps_YXyc) & ~repmat(bad_YX,[1 1 Nyears_out NcropsCombined_out])) = 0 ;
out_nfert.maps_YXvy = permute(out_nfert.maps_YXyc,[1 2 4 3]) ;
out_nfert = rmfield(out_nfert,'maps_YXyc') ;
out_irrig.maps_YXvy = permute(out_irrig.maps_YXyc,[1 2 4 3]) ;
out_irrig = rmfield(out_irrig,'maps_YXyc') ;

% Switch to LUH2 irrigation
disp('Switching to LUH2 irrigation...')
out_cropfrac2 = out_cropfrac ;
out_cropfrac2.maps_YXvy = nan(size(out_nfert.maps_YXvy)) ;
out_cropfrac2.varNames = out_cropfrac.varNames ;
for c = 1:(NcropsCombined_out-1)
    thisCrop = list_cropsCombined_out{c} ;
    thisCropI = [thisCrop 'i'] ;
    thisCrop_YX = sum(out_cropfrac.maps_YXv(:,:,contains(list_crops_out,thisCrop),:),3) ;
    thisCrop_irrFrac_YX1y = out_irrig.maps_YXvy(:,:,c,:) ;
    stop
    out_cropfrac2.maps_YXvy(:,:,strcmp(list_crops_out,thisCrop),:)  = repmat(thisCrop_YX,[1 1 1 Nyears_out]) .* (1-thisCrop_irrFrac_YX1y) ;
    out_cropfrac2.maps_YXvy(:,:,strcmp(list_crops_out,thisCropI),:) = repmat(thisCrop_YX,[1 1 1 Nyears_out]) .*    thisCrop_irrFrac_YX1y ;
end
out_cropfrac2.maps_YXvy(:,:,end,:) = repmat(out_cropfrac.maps_YXv(:,:,end),[1 1 1 Nyears_out]) ;

% Add manure N, assuming even distribution to all crops
disp('Adding manure to Nfert...')
load('/Users/sam/Geodata/Manure_ZhangEtAl2017/zhangManure_1860to2014_agg_hd.BAD.mat') ;
yearList_manure = 1860:2014 ;
yearList_manure_missing = setdiff(yearList_out,yearList_manure) ;
if ~isempty(yearList_manure_missing)
    manure2crop_hd_YXy = cat(3, ...
        repmat(manure2crop_hd_YXy(:,:,1),[1 1 length(find(yearList_manure_missing<min(yearList_manure)))]), ...
        manure2crop_hd_YXy, ...
        repmat(manure2crop_hd_YXy(:,:,end),[1 1 length(find(yearList_manure_missing>max(yearList_manure)))]) ...
        ) ;
    yearList_manure = cat(2, ...
        yearList_manure_missing(yearList_manure_missing<min(yearList_manure)), ...
        yearList_manure, ...
        yearList_manure_missing(yearList_manure_missing>max(yearList_manure)) ...
        ) ;
end
[~,IA,~] = intersect(yearList_manure,yearList_out) ;
manure2crop_hd_YXy = manure2crop_hd_YXy(:,:,IA) ;
clear IA
manure2crop_hd_YXy(manure2crop_hd_YXy<1e-6) = 0 ;
missing_manure_YXy = out_lu.maps_YXyv(:,:,:,strcmp(list_LU_out,'CROPLAND'))>0 & isnan(manure2crop_hd_YXy) ;
tmp = permute(out_lu.maps_YXyv,[1 2 4 3]) ;
out_lu = rmfield(out_lu,'maps_YXyv') ;
out_lu.maps_YXvy = tmp ;
if any(any(any(missing_manure_YXy)))
    warning('%d cells with CROPLAND have NaN manure. Assuming 0 manure there.',length(find(missing_manure_YXy)))
end
manure2crop_hd_YXy(isnan(manure2crop_hd_YXy)) = 0 ;
out_nfert.maps_YXvy(:,:,1:end-1,:) = out_nfert.maps_YXvy(:,:,1:end-1,:) + repmat(permute(manure2crop_hd_YXy,[1 2 4 3]),[1 1 Ncrops_out-1 1]) ;

% Convert from kg/ha to kg/m2
out_nfert.maps_YXvy = out_nfert.maps_YXvy * 1e-4 ;

% Ensure no fertilization where no area
out_nfert.maps_YXvy(out_cropfrac2.maps_YXvy==0) = 0 ;

disp('Done.')


%% Part 4: Set up for save

% Array-ify
disp('Array-ifying out_lu...')
[out_lu_array, out_lu_header_cell] = lpjgu_matlab_maps2table(out_lu,gridlist.list_to_map) ;
disp('Array-ifying out_cropfrac...')
[out_cropfrac_array, out_cropfrac_header_cell] = lpjgu_matlab_maps2table(out_cropfrac2,gridlist.list_to_map) ;
disp('Array-ifying out_nfert...')
[out_nfert_array, out_nfert_header_cell] = lpjgu_matlab_maps2table(out_nfert,gridlist.list_to_map) ;

% Add zeros for Miscanthus(i)
disp('Adding zeros for Miscanthus(i)...')
out_cropfrac_array = cat(2,out_cropfrac_array,zeros(size(out_cropfrac_array,1),2)) ;
out_cropfrac_header_cell = [out_cropfrac_header_cell {'Miscanthus','Miscanthusi'}] ;
out_nfert_array = cat(2,out_nfert_array,zeros(size(out_nfert_array,1),2)) ;
out_nfert_header_cell = [out_nfert_header_cell {'Miscanthus','Miscanthusi'}] ;

% Get filenames
if force_all_rainfed
    allRF_txt = '.noIrr' ;
else
    allRF_txt = '' ;
end
out_file_lu = [out_dir 'LU.remapv6p3.txt'] ;
out_file_cropfrac = [out_dir 'cropfracs.remapv6p3.txt'] ;
out_file_nfert = [out_dir 'nfert.remapv6p3.txt'] ;

disp('Done.')


%% Part 5: Save

%%% Options %%%%%%%%%%%%%
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
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
rng(20171116) ;
rdmsam = randsample(Ncells_4gl,Ncells_4gl) ;
% Save gridlist
outFile_gridlist = [out_dir 'gridlist.txt'] ;
out_formatSpec_gridlist = '%4.2f %4.2f\n' ;
fid1_gridlist = fopen(strrep(outFile_gridlist,'\ ',' '),'w') ;
fprintf(fid1_gridlist,out_formatSpec_gridlist,[lons_4gl(rdmsam) lats_4gl(rdmsam)]') ;
fclose(fid1_gridlist) ;
fid1_gridlist = fopen(strrep(outFile_gridlist,'\ ',' '),'a+') ;
fprintf(fid1_gridlist,'%s','') ;
fclose(fid1_gridlist) ;

disp('Saving LU...')
check_existing_lu(thisVer, out_file_lu, allVer_names, allVer_ignore_types) ;
lpjgu_matlab_saveTable(out_lu_header_cell, out_lu_array, out_file_lu,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 1) ;

disp('Saving cropfracs...')
lpjgu_matlab_saveTable(out_cropfrac_header_cell, out_cropfrac_array, out_file_cropfrac,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20) ;

disp('Saving nfert...')
lpjgu_matlab_saveTable(out_nfert_header_cell, out_nfert_array, out_file_nfert,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20) ;



