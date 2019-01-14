%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Re-map area/fert data to PLUM crops, and generate extra LU file %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Sam Rabin, 2018-10-02
% See previous notes for remap_v5, with the following changes:
% - Scale fertilization and irrigated fractions according to fractional
%   changes in LUH2-equivalent crops

% Carried over from remap_v5
PLUMsetAside_frac = 0.103 ;
inpaint_method = 4 ;

% Version for crop mappings
thisVer = '20180214' ;

force_all_rainfed = false ;

out_dir = '/Users/Shared/PLUM/input/remaps_v5/' ;

warning('on','all')


%% Part 1a: Import most data

disp('Importing...')

addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;

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

disp('Done.')


%% Part 1b: Set up mapping

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


%% Part 1c: Finish import

% Scale year-to-year fert/irrig changes based on LUH2
% Ncrops_mirca_in = length(nfert_in.varNames) ;

luh2_equivs_in = in2out_keyCombined_frac ;
list_crops_luh2_avail = {'c3ann','c4ann','c3per','c4per','c3nfx'} ;
Ncrops_luh2_avail = length(list_crops_luh2_avail) ;
for c = 1:NcropsCombined_out
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

luh2_to_import = list_crops_luh2_avail ;
for c = 1:Ncrops_luh2_avail
    thisCrop = list_crops_luh2_avail{c} ;
    disp(thisCrop)
    if ~any(cellfun(@(x) any(strcmp(x,thisCrop)),luh2_equivs_in))
        luh2_to_import(strcmp(luh2_to_import,thisCrop)) = [] ;
    end
end

% Import LUH2 originals: area, nfert, and irrig
Ncrops_luh2_in = length(luh2_to_import) ;
yearList_out = 1850:2015 ;
Nyears_out = length(yearList_out) ;
luh2_file_states = '/Users/sam/Geodata/LUH2/v2h/states.1850-2015.nc' ;
luh2_file_management = '/Users/sam/Geodata/LUH2/v2h/management.1850-2015.nc' ;
carea_YX = ncread('/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc','carea') ;
carea_YXy = repmat(carea_YX,[1 1 length(1850:2015)]) ;
yearList_inFile = 1850:2015 ;
[~,yearIndices,~] = intersect(yearList_inFile,yearList_out) ;
starts = [1 1 min(yearIndices)] ;
counts = [Inf Inf Nyears_out] ;
luh2_carea_YXyc = nan(360,720,Nyears_out,Ncrops_luh2_in) ;
luh2_nfert_YXyc = nan(360,720,Nyears_out,Ncrops_luh2_in) ;
luh2_irrig_YXyc = nan(360,720,Nyears_out,Ncrops_luh2_in) ;
for c = 1:length(luh2_to_import)
    
    thisVar = luh2_to_import{c} ;
    fprintf('%s: ', thisVar) ;
    
    fprintf('area... ') ;
    luh2_carea_YXy = ncread(luh2_file_states,thisVar,starts,counts) ;
    luh2_carea_YXy = luh2_carea_YXy .* carea_YXy ;
    luh2_carea_YXyc(:,:,:,c) = flip(permute(PLUMharm_aggregate(luh2_carea_YXy,0.25,0.5),[2 1 3]),1) ;
    
    fprintf('nfert... ') ;
    thisVar = ['fertl_' luh2_to_import{c}] ;
    tmp = ncread(luh2_file_management,thisVar,starts,counts) ;
    luh2_nfert_YXyc(:,:,:,c) = flip(permute(PLUMharm_aggregate_mgmt(tmp,luh2_carea_YXy,0.25,0.5),[2 1 3]),1) ;
    clear tmp
    
    fprintf('irrig... ') ;
    thisVar = ['irrig_' luh2_to_import{c}] ;
    tmp = ncread(luh2_file_management,thisVar,starts,counts) ;
    luh2_irrig_YXyc(:,:,:,c) = flip(permute(PLUMharm_aggregate_mgmt(tmp,luh2_carea_YXy,0.25,0.5),[2 1 3]),1) ;
    clear tmp
    
    clear luh2_carea_YXy
    disp('Done.')
end
disp('Done importing. Processing...')
% Get fractional values relative to 2000
index_2000 = find(yearList_out==2000) ;
luh2_carea_YX1c = luh2_carea_YXyc(:,:,index_2000,:) ;
luh2_nfertRel_YXyc = luh2_nfert_YXyc ./ repmat(luh2_nfert_YXyc(:,:,index_2000,:),[1 1 Nyears_out 1]) ;


luh2_irrig_YXbc = luh2_irrig_YXyc(:,:,index_2000,:) ;
everGt0_YX1c = max(luh2_irrig_YXyc,[],3)>0 ;
everArea_YX1c = max(luh2_carea_YXyc,[],3)>0 ;
% If LUH2_irrig(2000)==0 but LUH2_irrig ever >0, set baseline value of each
% gridcell to the nearest value there in time.
yy = 0 ;
while any(any(any(luh2_irrig_YXbc==0 & everGt0_YX1c,4),2),1)
    yy = yy + 1 ;
    % Get upper year and lower year
    i1 = index_2000 - yy ;
    iN = index_2000 + yy ;
    % If upper xor lower year is outside year range, ignore. If both are
    % outside year range, throw error.
    if i1<1
        i1 = [] ;
    end
    if iN > Nyears_out
        iN = [] ;
    end
    indices = [i1 iN] ;
    if isempty(indices)
        error('Possible infinite loop in "while any(any(any(luh2_irrig_YXbc==0,4),2),1)"')
    end
    % Get mean of upper and lower years, ignoring either if zero
    luh2_irrig_YX2c = luh2_irrig_YXyc(:,:,[i1 iN],:) ;
    luh2_irrig_YX2c(luh2_irrig_YX2c==0) = NaN ;
    luh2_irrig_YX1c = nanmean(luh2_irrig_YX2c,3) ;
    luh2_irrig_YXbc(luh2_irrig_YXbc==0 & ~isnan(luh2_irrig_YX1c)) ...
        = luh2_irrig_YX1c(luh2_irrig_YXbc==0 & ~isnan(luh2_irrig_YX1c)) ;
    % Will repeat WHILE loop if necessary...
end

luh2_irrigRel_YXyc = luh2_irrig_YXyc ./ repmat(luh2_irrig_YXbc,[1 1 Nyears_out 1]) ;
luh2_irrigRel_YXyc(repmat(everArea_YX1c & ~everGt0_YX1c,[1 1 Nyears_out 1])) = 1 ;
if any(any(any(any(isinf(luh2_irrigRel_YXyc)))))
    error('Inf in luh2_irrigRel_YXyc!')
elseif any(any(any(any(isnan(luh2_irrigRel_YXyc)))))
    error('NaN in luh2_irrigRel_YXyc!')
end
error('Add code to reset year-2000 irrig! Maybe by setting Rel to 0?')




% Fill in missing values with global average of relative difference that
% year for that crop
luh2_nfertTot_YXyc = luh2_nfert_YXyc .* luh2_carea_YXyc ;
luh2_nfertGlob_11yc = nansum(nansum(luh2_nfertTot_YXyc,1),2) ;
luh2_nfertGlobRel_11yc = luh2_nfertGlob_11yc ./ repmat(luh2_nfertGlob_11yc(:,:,index_2000,:),[1 1 Nyears_out 1]) ;
luh2_nfertGlobRel_YXyc = repmat(luh2_nfertGlobRel_11yc,[360 720 1 1]) ;
luh2_nfertRel_YXyc(isnan(luh2_nfertRel_YXyc)) = luh2_nfertGlobRel_YXyc(isnan(luh2_nfertRel_YXyc)) ;
luh2_irrigTot_YXyc = luh2_irrig_YXyc .* luh2_carea_YXyc ;
luh2_irrigGlob_11yc = nansum(nansum(luh2_irrigTot_YXyc,1),2) ;
luh2_irrigGlobRel_11yc = luh2_irrigGlob_11yc ./ repmat(luh2_irrigGlob_11yc(:,:,index_2000,:),[1 1 Nyears_out 1]) ;
luh2_irrigGlobRel_YXyc = repmat(luh2_irrigGlobRel_11yc,[360 720 1 1]) ;
stop
luh2_irrigRel_YXyc(isnan(luh2_irrigRel_YXyc)) = luh2_irrigGlobRel_YXyc(isnan(luh2_irrigRel_YXyc)) ;
% Restrict to LUH2 gridlist
isnan_LUH2_YXyc = repmat(isnan_LUH2,[1 1 Nyears_out Ncrops_luh2_in]) ;
luh2_nfertRel_YXyc(isnan_LUH2_YXyc) = NaN ;
luh2_irrigRel_YXyc(isnan_LUH2_YXyc) = NaN ;
clear isnan_LUH2_YXyc
disp('Done.')

%% Convert to output equivalents
out_nfertRel_YXyc = nan(360,720,Nyears_out,Ncrops_out) ;
out_irrigRel_YXyc = nan(360,720,Nyears_out,Ncrops_out) ;
for c = 1:NcropsCombined_out
    thisCrop = list_cropsCombined_out{c} ;
    disp(thisCrop)
    these_equivs = luh2_equivs_in{c} ;
    is_match = contains(list_crops_out,thisCrop) ;
    Nmatch = length(find(is_match)) ;
    if length(unique(these_equivs))==1
        out_nfertRel_YXyc(:,:,:,is_match) = repmat(luh2_nfertRel_YXyc(:,:,:,strcmp(luh2_to_import,these_equivs{1})),[1 1 1 Nmatch]) ;
        out_irrigRel_YXyc(:,:,:,is_match) = repmat(luh2_irrigRel_YXyc(:,:,:,strcmp(luh2_to_import,these_equivs{1})),[1 1 1 Nmatch]) ;
    else
        i_luh2 = cellfun(@(x) find(strcmp(luh2_to_import,x)),these_equivs) ;
        weightsDenom_YXyc = repmat(sum(luh2_carea_YXyc(:,:,:,i_luh2),4),[1 1 1 length(i_luh2)]) ;
        weights_YXyc = luh2_carea_YXyc(:,:,:,i_luh2) ./ weightsDenom_YXyc ;
        if any(any(any(abs(sum(weights_YXyc,4)-1)>1e-6)))
            error('Something went wrong with weighting in converting to %s.',thisCrop) ;
        end
        tmp_YXy = sum(luh2_nfertRel_YXyc(:,:,:,i_luh2) .* weights_YXyc,4) ;
        if any(any(any(isinf(tmp_YXy))))
            error('Inf in tmp_YXy!')
        end
        luh2_nfertGlob_11yd = ...
            nansum(nansum(...
            luh2_carea_YXyc(:,:,:,i_luh2).*luh2_nfertTot_YXyc(:,:,:,i_luh2)...
            ,1),2) ;
        luh2_nfertGlob_11y = sum(luh2_nfertGlob_11yd,4) ;
        luh2_nfertGlobRel_11y = luh2_nfertGlob_11y ./ repmat(luh2_nfertGlob_11y(:,:,index_2000),[1 1 Nyears_out]) ;
        luh2_nfertGlobRel_YXy = repmat(luh2_nfertGlobRel_11y,[360 720 1]) ;
        tmp_YXy(weightsDenom_YXyc(:,:,:,1)==0) = luh2_nfertGlobRel_YXy(weightsDenom_YXyc(:,:,:,1)==0) ;
        if any(any(any(isinf(tmp_YXy))))
            error('Inf in tmp_YXy!')
        end
        out_nfertRel_YXyc(:,:,:,is_match) = repmat(tmp_YXy,[1 1 1 Nmatch]) ;
        clear tmp_YXy luh2_nfertGlob_11yd luh2_nfertGlob_11y luh2_nfertGlobRel_11y luh2_nfertGlobRel_YXy
        tmp_YXy = sum(luh2_irrigRel_YXyc(:,:,:,i_luh2) .* weights_YXyc,4) ;
        if any(any(any(isinf(tmp_YXy))))
            error('Inf in tmp_YXy!')
        end
        luh2_irrigGlob_11yd = ...
            nansum(nansum(...
            luh2_carea_YXyc(:,:,:,i_luh2).*luh2_irrigTot_YXyc(:,:,:,i_luh2)...
            ,1),2) ;
        luh2_irrigGlob_11y = sum(luh2_irrigGlob_11yd,4) ;
        luh2_irrigGlobRel_11y = luh2_irrigGlob_11y ./ repmat(luh2_irrigGlob_11y(:,:,index_2000),[1 1 Nyears_out]) ;
        luh2_irrigGlobRel_YXy = repmat(luh2_irrigGlobRel_11y,[360 720 1]) ;
        tmp_YXy(weightsDenom_YXyc(:,:,:,1)==0) = luh2_irrigGlobRel_YXy(weightsDenom_YXyc(:,:,:,1)==0) ;
        if any(any(any(isinf(tmp_YXy))))
            error('Inf in tmp_YXy!')
        end
        out_irrigRel_YXyc(:,:,:,is_match) = repmat(tmp_YXy,[1 1 1 Nmatch]) ;
        clear tmp_YXy luh2_irrigGlob_11yd luh2_irrigGlob_11y luh2_irrigGlobRel_11y luh2_irrigGlobRel_YXy
        clear weights_YXyc
    end
end


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

cropfrac_out.maps_YXv = nan(360,720,Ncrops_out) ;
nfert_out.maps_YXv = nan(360,720,Ncrops_out) ;
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

cropfrac_out.varNames = list_crops_out ;
nfert_out.varNames = list_crops_out ;

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


%% Part 5: Finish managements

% % Apply relative changes to MIRCA managements
% nfert_out2.maps_YXvy = nan([360 720 length(actual_crops_out) Nyears_out]) ;
% done_yet = false(length(actual_crops_out),1) ;
% for c = 1:Ncrops_luh2_in
%     thisCrop_luh2 = luh2_to_import{c} ;
%     disp(thisCrop_luh2)
%     isThisCrop = strcmp(luh2_equivs,thisCrop_luh2) ;
%     done_yet = done_yet | isThisCrop ;
%     Nmatched = length(find(isThisCrop)) ;
%     nfert_out2.maps_YXvy(:,:,isThisCrop,:) = ...
%         repmat(nfert_out.maps_YXv(:,:,isThisCrop), [1 1 1 Nyears_out]) ...
%         .* repmat(permute(luh2_nfertRel_YXyc(:,:,:,c), [1 2 4 3]),[1 1 Nmatched 1]) ;
%     irrig_out2.maps_YXvy(:,:,isThisCrop,:) = ...
%         repmat(irrig_out.maps_YXv(:,:,isThisCrop), [1 1 1 Nyears_out]) ...
%         .* repmat(permute(luh2_irrigRel_YXyc(:,:,:,c), [1 2 4 3]),[1 1 Nmatched 1]) ;
% 
% end
% disp('Done')

%% Add manure N for year 2000 (new in remap_v5)
load('/Users/sam/Geodata/Manure_ZhangEtAl2017/zhangManure_1860to2014_agg_hd.mat') ;
manure2crop_hd_YX = manure2crop_hd_YXy(:,:,1860:2014==2000) ;
manure2crop_hd_YX(manure2crop_hd_YX<1e-6) = 0 ;
clear manure2crop_hd_YXy
if any(any(sum(cropareaCombined_in.maps_YXv,3)>0 & isnan(manure2crop_hd_YX)))
    manure2crop_hd2_YX = inpaint_nans(manure2crop_hd_YX,inpaint_method) ;
end
manure2crop_hd2_YX(isnan_LUH2) = NaN ;
nfert_in.maps_YXv = nfert_in.maps_YXv + repmat(manure2crop_hd_YX,[1 1 size(nfert_in.maps_YXv,3)]) ;



%% Part 6: Set up for save

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
out_file_lu = [out_dir 'LU.remapv5.' thisVer '.ecFertIrr0.setaside' strrep(num2str(PLUMsetAside_frac),'.','') '.m' num2str(inpaint_method) '.txt'] ;
out_file_cropfrac = [out_dir 'cropfracs.remapv5.' thisVer allRF_txt '.ecFertIrr0.setaside' strrep(num2str(PLUMsetAside_frac),'.','') '.m' num2str(inpaint_method) '.txt'] ;
out_file_nfert = [out_dir 'nfert.remapv5.' thisVer allRF_txt '.ecFertIrr0.setaside' strrep(num2str(PLUMsetAside_frac),'.','') '.m' num2str(inpaint_method) '.txt'] ;




%% Part 7: Save

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




