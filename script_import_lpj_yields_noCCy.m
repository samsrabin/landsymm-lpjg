% Setup
crops2remove = {'CC3G','CC4G','OtHr','ExtraCrop'} ;

if strcmp(version_name,'stijn_20180119')
    if calib_ver~=11
        error('When doing stijn_20180119, you must use calib_ver=11.')
    end
    listCrops_lpj_comb = {'TeWW','TeSW','TeCo','TrRi','TrSo','GlyM','FaBe'} ;
elseif calib_ver == 11
    error('calib_ver 11 only works with stijn_20180119.')
elseif calib_ver == 12 || calib_ver == 15 || calib_ver == 16 || calib_ver == 18
    listCrops_lpj_comb = {'CerealsC3','CerealsC4','Rice','Oilcrops','StarchyRoots','Pulses'} ;
elseif calib_ver == 19
    listCrops_lpj_comb = ...
        {'CerealsC3','CerealsC4','Rice','Oilcrops',...
        'StarchyRoots','Pulses','Sugar','DateCitGrape'} ;
elseif calib_ver == 20
    listCrops_lpj_comb = ...
        {'CerealsC3','CerealsC4','Rice','Oilcrops',...
        'StarchyRoots','Pulses','Sugar','FruitAndVeg'} ;
elseif calib_ver == 13
    listCrops_lpj_comb = {'Wheat','Maize','Sorghum','Pulses','Soybeans','Rice'} ;
elseif calib_ver == 14
    listCrops_lpj_comb = {'Wheat','Maize','Sorghum','Rice'} ;
elseif calib_ver <= 10 || calib_ver == 17
    listCrops_lpj_comb = {'TeWW','TeSW','TeCo','TrRi'} ;
elseif contains(version_name,'jianyong_20190128')
    if calib_ver~=21
        error('When doing %s, you must use calib_ver= 21.', version_name)
    end
    listCrops_lpj_comb = {'FaBe','GlyM','TeCo','TeSW','TeWW','TrRi','TrSo'} ;
else
    error(['calib_ver ' num2str(calib_ver) ' not recognized for setting listCrops_lpj_comb (version_name '  ').']) ;
end

% Import land uses
disp('Importing land uses and crop fractions...')
cropfrac_lpj = lpjgu_matlab_readTable_then2map(filename_guess_cropfrac,'force_mat_save',true) ;
landuse_lpj = lpjgu_matlab_readTable_then2map(filename_guess_landuse,'force_mat_save',true) ;

% Import yield
disp('Importing simulated yield...')
if exist('filename_guess_yield', 'var')
    % LPJ-GUESS style
    if strcmp(filename_guess_yield(end-2:end),'.gz')
        filename_guess_yield = filename_guess_yield(1:end-3) ;
    end
    yield_lpj = lpjgu_matlab_readTable_then2map(filename_guess_yield,'force_mat_save',true) ;
    
    % Convert yield from kgDM/m2 to tDM/ha
    if isfield(yield_lpj,'maps_YXvy')
        yield_lpj.maps_YXvy = yield_lpj.maps_YXvy * 1e4 * 1e-3 ;
    else
        yield_lpj.maps_YXv = yield_lpj.maps_YXv * 1e4 * 1e-3 ;
    end
elseif exist('dirname_emuBL_yields', 'var')
    % GGCMI phase 2
    cropList_ggcmi = {'maize', 'rice', 'soy', 'spring_wheat', 'winter_wheat'} ;
    irrList_ggcmi = {'rf', 'ir'} ;
    NcropIrr_ggcmi = length(cropList_ggcmi)*length(irrList_ggcmi) ;
    listCrops_ggcmi = cell(NcropIrr_ggcmi, 1) ;
    yield_ggcmi.maps_YXv = nan(360, 720, NcropIrr_ggcmi) ;
    
    % Import (already in tDM/ha)
    for c = 1:length(cropList_ggcmi)
        thisCrop = cropList_ggcmi{c} ;
        thisCrop_short = e2p_get_thisCrop_short(thisCrop) ;
        thisPattern = sprintf('%s/*%s*', ...
            dirname_emuBL_yields, thisCrop) ;
        filelist = dir(thisPattern) ;
        if length(filelist) ~= 1
            error('Error finding emulated outputs for %s: expected 1, found %d', ...
                thisCrop, length(filelist))
        end
        thisFile = sprintf('%s/%s', filelist.folder, filelist.name) ;
        clear filelist
        for ii = 1:2
            v = (c-1)*2+ii ;
            listCrops_ggcmi{v} = sprintf('%s_%s', ...
                thisCrop, irrList_ggcmi{ii}) ;
            thisVar = sprintf('yield_%s_%s', irrList_ggcmi{ii}, thisCrop_short) ;
            yield_ggcmi.maps_YXv(:,:,v) = flipud(transpose(ncread(thisFile, thisVar))) ;
        end
    end
    yield_ggcmi.varNames = listCrops_ggcmi ;
    
    % Add max wheats
    I_wheats_rf = contains(yield_ggcmi.varNames, 'wheat_rf') ;
    I_wheats_ir = contains(yield_ggcmi.varNames, 'wheat_ir') ;
    maxWheat_rf_YX = max(yield_ggcmi.maps_YXv(:,:,I_wheats_rf), [], 3) ;
    maxWheat_ir_YX = max(yield_ggcmi.maps_YXv(:,:,I_wheats_ir), [], 3) ;
    cropList_ggcmi = [cropList_ggcmi {'max_wheat'}] ;
    yield_ggcmi.varNames = [yield_ggcmi.varNames ;
        {'max_wheat_rf'; 'max_wheat_ir'}] ;
    NcropIrr_ggcmi = length(yield_ggcmi.varNames) ;
    yield_ggcmi.maps_YXv = cat(3, yield_ggcmi.maps_YXv, ...
        maxWheat_rf_YX, maxWheat_ir_YX) ;
    
    % Translate crops
    I_toRemove = ~contains(cropfrac_lpj.varNames, crops2remove) ;
    varNames_out = cropfrac_lpj.varNames( ...
        ~contains(cropfrac_lpj.varNames, crops2remove)) ;
    varNames_ggcmi_lpjIrr = yield_ggcmi.varNames ;
    varNames_ggcmi_lpjIrr = regexprep(varNames_ggcmi_lpjIrr, '_rf$', '') ;
    varNames_ggcmi_lpjIrr = regexprep(varNames_ggcmi_lpjIrr, '_ir$', 'i') ;
    I_out = e2p_translate_crops_agm2out(...
        varNames_ggcmi_lpjIrr, varNames_out) ;
    [varNames_out' yield_ggcmi.varNames(I_out)]
    yield_lpj.varNames = varNames_out ;
    yield_lpj.maps_YXv = yield_ggcmi.maps_YXv(:,:,I_out) ;
else
    error('How am I supposed to import yields?')
end

disp('Processing...')

% For GGCMI, we're not comparing individual years, so find means for years
% of interest
if is_ggcmi 
    
    % yield
    if isfield(yield_lpj,'yearList') || isfield(yield_lpj,'maps_YXvy')
        error('GGCMI yields have years??')
    end
    yield_lpj.maps_YXvy = yield_lpj.maps_YXv ;
    yield_lpj = rmfield(yield_lpj, 'maps_YXv') ;
    yield_lpj.yearList = -pi ;
    
    % landuse
    if isfield(landuse_lpj,'maps_YXvy')
        landuse_lpj.maps_YXvy = nanmean(landuse_lpj.maps_YXvy(:,:,:, ...
            landuse_lpj.yearList>=year1 & landuse_lpj.yearList<=yearN), 4) ;
    else
        landuse_lpj.maps_YXvy = landuse_lpj.maps_YXv ;
        landuse_lpj = rmfield(landuse_lpj, 'maps_YXv') ;
    end
    landuse_lpj.yearList = -pi ;
    
    % cropfrac
    if isfield(cropfrac_lpj,'maps_YXvy')
        cropfrac_lpj.maps_YXvy = nanmean(cropfrac_lpj.maps_YXvy(:,:,:, ...
            cropfrac_lpj.yearList>=year1 & cropfrac_lpj.yearList<=yearN), 4) ;
    else
        cropfrac_lpj.maps_YXvy = cropfrac_lpj.maps_YXv ;
        cropfrac_lpj = rmfield(cropfrac_lpj, 'maps_YXv') ;
    end
    cropfrac_lpj.yearList = -pi ;
end

% Trim out CC3G, CC4G, ExtraCrop, and OtHr from yield_lpj
% toRemove = find(strcmp(yield_lpj.varNames,'CC3G_ic') | strcmp(yield_lpj.varNames,'CC4G_ic')) ;
% toRemove = find(strcmp(yield_lpj.varNames,'CC3G_ic') | strcmp(yield_lpj.varNames,'CC4G_ic') ...
%     | strcmp(yield_lpj.varNames,'OtHr') | strcmp(yield_lpj.varNames,'OtHri')) ;
toRemove = find(contains(yield_lpj.varNames,crops2remove)) ;
if ~isempty(toRemove)
    yield_lpj.maps_YXvy(:,:,toRemove,:) = [] ;
    yield_lpj.varNames(toRemove) = [] ;
end
% toRemove = find(strcmp(cropfrac_lpj.varNames,'CC3G_ic') | strcmp(cropfrac_lpj.varNames,'CC4G_ic') ...
%     | strcmp(cropfrac_lpj.varNames,'OtHr') | strcmp(cropfrac_lpj.varNames,'OtHri')) ;
toRemove = find(contains(cropfrac_lpj.varNames,crops2remove)) ;
if ~isempty(toRemove)
    if isfield(cropfrac_lpj,'maps_YXvy')
        cropfrac_lpj.maps_YXvy(:,:,toRemove,:) = [] ;
    else
        cropfrac_lpj.maps_YXv(:,:,toRemove) = [] ;
    end
    cropfrac_lpj.varNames(toRemove) = [] ;
end
clear toRemove

% Create cropfrac_lpj.YXvy array from YXv, if necessary
if ~is_ggcmi && ~isfield(cropfrac_lpj,'maps_YXvy')
    tmp = cropfrac_lpj.maps_YXv ;
    cropfrac_lpj = rmfield(cropfrac_lpj,'maps_YXv') ;
    cropfrac_lpj.yearList = yield_lpj.yearList ;
    cropfrac_lpj.maps_YXvy = repmat(tmp,[1 1 1 length(yield_lpj.yearList)]) ;
end

% 2014 land use is repeated for 2015
% if max(landuse_lpj.yearList==2015)
if ~is_ggcmi && max(yield_lpj.yearList)==2015
    if max(cropfrac_lpj.yearList)==2014
        cropfrac_lpj.maps_YXvy(:,:,:,size(cropfrac_lpj.maps_YXvy,4)+1) = cropfrac_lpj.maps_YXvy(:,:,:,size(cropfrac_lpj.maps_YXvy,4)) ;
        cropfrac_lpj.yearList(length(cropfrac_lpj.yearList)+1) = cropfrac_lpj.yearList(length(cropfrac_lpj.yearList)) + 1 ;
    elseif max(cropfrac_lpj.yearList) < 2015
        error('max(cropfrac_lpj.yearList)<2014')
    end
    if max(landuse_lpj.yearList)==2014
        landuse_lpj.maps_YXvy(:,:,:,size(landuse_lpj.maps_YXvy,4)+1) = landuse_lpj.maps_YXvy(:,:,:,size(landuse_lpj.maps_YXvy,4)) ;
        landuse_lpj.yearList(length(landuse_lpj.yearList)+1) = landuse_lpj.yearList(length(landuse_lpj.yearList)) + 1 ;
    elseif max(landuse_lpj.yearList) < 2015
        error('max(landuse_lpj.yearList)<2014')
    end
end

% Align years, if necessary
if ~is_ggcmi && ~isequal(yield_lpj.yearList,cropfrac_lpj.yearList)
    [~,~,IB] = intersect(yield_lpj.yearList,cropfrac_lpj.yearList) ;
    cropfrac_lpj.maps_YXvy = cropfrac_lpj.maps_YXvy(:,:,:,IB) ;
    cropfrac_lpj.yearList = cropfrac_lpj.yearList(IB) ;
end
if ~is_ggcmi && ~isequal(yield_lpj.yearList,landuse_lpj.yearList)
    [~,~,IB] = intersect(yield_lpj.yearList,landuse_lpj.yearList) ;
    landuse_lpj.maps_YXvy = landuse_lpj.maps_YXvy(:,:,:,IB) ;
    landuse_lpj.yearList = landuse_lpj.yearList(IB) ;
end
if ~is_ggcmi && (~isequal(yield_lpj.yearList,cropfrac_lpj.yearList) || ~isequal(yield_lpj.yearList,landuse_lpj.yearList))
    error('Yearlists don''t match!')
end

% Combine WW and SW, if necessary
% if strcmp(version_name,'plumTypes_assignWWSW')
if ~is_ggcmi && (any(strcmp(yield_lpj.varNames,'CerealsC3w')) || any(strcmp(yield_lpj.varNames,'CerealsC3s')))
    warning('Combining WW and SW.')
    % Rainfed
    i_yield_rf = [find(strcmp(yield_lpj.varNames,'CerealsC3w')) find(strcmp(yield_lpj.varNames,'CerealsC3s'))] ;
    i_fracs_rf = [find(strcmp(cropfrac_lpj.varNames,'CerealsC3w')) find(strcmp(cropfrac_lpj.varNames,'CerealsC3s'))] ;
    if length(i_yield_rf)~=2; error('length(i_yield_rf)~=2'); end
    if length(i_fracs_rf)~=2; error('length(i_fracs_rf)~=2'); end
    tmp_yield_rf_YXvy = yield_lpj.maps_YXvy(:,:,i_yield_rf,:) ;
    tmp_fracs_rf_YXvy = cropfrac_lpj.maps_YXvy(:,:,i_fracs_rf,:) ;
    % Irrigated
    i_yield_ir = [find(strcmp(yield_lpj.varNames,'CerealsC3wi')) find(strcmp(yield_lpj.varNames,'CerealsC3si'))] ;
    i_fracs_ir= [find(strcmp(cropfrac_lpj.varNames,'CerealsC3wi')) find(strcmp(cropfrac_lpj.varNames,'CerealsC3si'))] ;
    if length(i_yield_ir)~=2; error('length(i_yield_ir)~=2'); end
    if length(i_fracs_ir)~=2; error('length(i_fracs_ir)~=2'); end
    tmp_yield_ir_YXvy = yield_lpj.maps_YXvy(:,:,i_yield_ir,:) ;
    tmp_fracs_ir_YXvy = cropfrac_lpj.maps_YXvy(:,:,i_fracs_ir,:) ;
    % Combine
    tmp_fracs_rf_YX_y = nansum(tmp_fracs_rf_YXvy,3) ;
    tmp_yield_rf_YX_y = nansum(tmp_yield_rf_YXvy .* tmp_fracs_rf_YXvy ./ repmat(tmp_fracs_rf_YX_y,[1 1 2 1]),3) ;
    tmp_fracs_ir_YX_y = nansum(tmp_fracs_ir_YXvy,3) ;
    tmp_yield_ir_YX_y = nansum(tmp_yield_ir_YXvy .* tmp_fracs_ir_YXvy ./ repmat(tmp_fracs_ir_YX_y,[1 1 2 1]),3) ;
    % Replace
    yield_lpj.maps_YXvy(:,:,[i_yield_rf i_yield_ir],:) = [] ;
    cropfrac_lpj.maps_YXvy(:,:,[i_fracs_rf i_fracs_ir],:) = [] ;
    yield_lpj.varNames([i_yield_rf i_yield_ir]) = [] ;
    cropfrac_lpj.varNames([i_fracs_rf i_fracs_ir]) = [] ;
    yield_lpj.varNames = [yield_lpj.varNames 'CerealsC3' 'CerealsC3i'] ;
    cropfrac_lpj.varNames = [cropfrac_lpj.varNames 'CerealsC3' 'CerealsC3i'] ;
    yield_lpj.maps_YXvy = cat(3,yield_lpj.maps_YXvy,tmp_yield_rf_YX_y,tmp_yield_ir_YX_y) ;
    cropfrac_lpj.maps_YXvy = cat(3,cropfrac_lpj.maps_YXvy,tmp_fracs_rf_YX_y,tmp_fracs_ir_YX_y) ;
end

% Stijn didn't include irrigated of these?
if ~is_ggcmi && (strcmp(version_name,'stijn_20180119') || contains(version_name,'jianyong_20190128'))
    tmpRemoveList = {'FaBei','GlyMi','TrSoi'} ;
    for c = 1:length(tmpRemoveList)
        thisCrop = tmpRemoveList{c} ;
        thisIndex = find(strcmp(cropfrac_lpj.varNames,thisCrop)) ;
        if isfield(cropfrac_lpj,'maps_YXv')
            cropfrac_lpj.maps_YXv(:,:,thisIndex) = [] ;
        else
            cropfrac_lpj.maps_YXvy(:,:,thisIndex,:) = [] ;
        end
        cropfrac_lpj.varNames(thisIndex) = [] ;
    end
end

% Rearrange cropfrac_lpj so that variables are in the same order as
% yield_lpj variables
new_order = nan(length(yield_lpj.varNames),1) ;
for v = 1:length(new_order)
    this_from_yield = yield_lpj.varNames{v} ;
    new_order(v) = find(strcmp(cropfrac_lpj.varNames,this_from_yield)) ;
    clear this_from_yield
end ; clear v
cropfrac_lpj.varNames = cropfrac_lpj.varNames(new_order) ;
if ~is_ggcmi
    if ~isfield(cropfrac_lpj,'maps_YXvy')
        tmp = cropfrac_lpj.maps_YXv ;
        cropfrac_lpj = rmfield(cropfrac_lpj,'maps_YXv') ;
        cropfrac_lpj.maps_YXvy = repmat(tmp,[1 1 1 size(yield_lpj.maps_YXvy,4)]) ;
        clear tmp
        cropfrac_lpj.yearList = yield_lpj.yearList ;
    end
    cropfrac_lpj.maps_YXvy = cropfrac_lpj.maps_YXvy(:,:,new_order,:) ;
else
    if ~isfield(cropfrac_lpj,'maps_YXvy')
        cropfrac_lpj.maps_YXv = cropfrac_lpj.maps_YXv(:,:,new_order) ;
    else
        cropfrac_lpj.maps_YXvy = cropfrac_lpj.maps_YXvy(:,:,new_order,:) ;
    end
end

% % % %%%%%%%
% % % % landuse_lpj_orig = landuse_lpj ;
% % % % if exist('version_name','var') && strcmp(version_name,'2017-05-17')
% landuse_lpj.maps_YXvy = landuse_lpj.maps_YXvy(:,:,:,ismember(landuse_lpj.yearList,yield_lpj.yearList)) ;
% landuse_lpj.yearList = landuse_lpj.yearList(ismember(landuse_lpj.yearList,yield_lpj.yearList)) ;
% % % % end

% Get year info
Nyears_lpj = length(landuse_lpj.yearList) ;

% Convert cropfrac from fraction of CROPLAND to fraction of ALL LAND
if isfield(cropfrac_lpj,'maps_YXvy')
    for v = 1:size(cropfrac_lpj.maps_YXvy,3)
        cropfrac_lpj.maps_YXvy(:,:,v,:) = ...
            cropfrac_lpj.maps_YXvy(:,:,v,:) ...
            .* landuse_lpj.maps_YXvy(:,:,strcmp(landuse_lpj.varNames,'CROPLAND'),:) ;
    end
else
    ok_years = landuse_lpj.yearList>=year1 & landuse_lpj.yearList<=yearN ;
    for v = 1:size(cropfrac_lpj.maps_YXvy,3)
        cropfrac_lpj.maps_YXv(:,:,v) = ...
            cropfrac_lpj.maps_YXv(:,:,v) ...
            .* nanmean(landuse_lpj.maps_YXvy(:,:,strcmp(landuse_lpj.varNames,'CROPLAND'),ok_years),4) ;
    end
end

% Create combined irrigated+rainfed yields
Ncrops_lpj_comb = length(listCrops_lpj_comb) ;
combined_YXcy_yield = nan(size(cropfrac_lpj.maps_YXvy,1),...
                    size(cropfrac_lpj.maps_YXvy,2),...
                    Ncrops_lpj_comb,...
                    length(cropfrac_lpj.yearList)) ;
combined_YXcy_cropfrac = combined_YXcy_yield ;
for c = 1:Ncrops_lpj_comb
    thisCropR = listCrops_lpj_comb{c} ;
    thisCropI = [thisCropR 'i'] ;
%     iR_cropfrac = exact_string_in_cellarray(cropfrac_lpj.varNames,thisCropR) ;
%     iI_cropfrac = exact_string_in_cellarray(cropfrac_lpj.varNames,thisCropI) ;
    iR_cropfrac = find(strcmp(cropfrac_lpj.varNames,thisCropR)) ;
    iI_cropfrac = find(strcmp(cropfrac_lpj.varNames,thisCropI)) ;
    if strcmp(version_name,'stijn_20180119') || contains(version_name,'jianyong_20190128')
        frac_comb = cropfrac_lpj.maps_YXvy(:,:,iR_cropfrac,:) ;
        if ~isempty(iI_cropfrac)
            frac_comb = frac_comb + cropfrac_lpj.maps_YXvy(:,:,iI_cropfrac,:) ;
        else
            warning([thisCropI ' not found in cropfrac_lpj.varNames.'])
        end
    else
        frac_comb = cropfrac_lpj.maps_YXvy(:,:,exact_string_in_cellarray(cropfrac_lpj.varNames,thisCropR),:) ...
                  + cropfrac_lpj.maps_YXvy(:,:,exact_string_in_cellarray(cropfrac_lpj.varNames,thisCropI),:) ;
    end
%     iR_yield = exact_string_in_cellarray(yield_lpj.varNames,thisCropR) ;
%     iI_yield = exact_string_in_cellarray(yield_lpj.varNames,thisCropI) ;
    iR_yield = find(strcmp(yield_lpj.varNames,thisCropR)) ;
    iI_yield = find(strcmp(yield_lpj.varNames,thisCropI)) ;
    if strcmp(version_name,'stijn_20180119') || contains(version_name,'jianyong_20190128')
        if ~isempty(iI_cropfrac)
            tmp = (yield_lpj.maps_YXvy(:,:,iR_yield,:) .* cropfrac_lpj.maps_YXvy(:,:,iR_cropfrac,:) ...
                 + yield_lpj.maps_YXvy(:,:,iI_yield,:) .* cropfrac_lpj.maps_YXvy(:,:,iI_cropfrac,:)) ...
                ./ frac_comb ;
        else
            tmp = yield_lpj.maps_YXvy(:,:,iR_yield,:) .* cropfrac_lpj.maps_YXvy(:,:,iR_cropfrac,:) ...
                ./ frac_comb ;
        end
    else
        tmp = (yield_lpj.maps_YXvy(:,:,iR_yield,:) .* cropfrac_lpj.maps_YXvy(:,:,iR_cropfrac,:) ...
             + yield_lpj.maps_YXvy(:,:,iI_yield,:) .* cropfrac_lpj.maps_YXvy(:,:,iI_cropfrac,:)) ...
            ./ frac_comb ;
    end
    combined_YXcy_yield(:,:,c,:) = tmp ;
    combined_YXcy_cropfrac(:,:,c,:) = frac_comb ;
end

if isempty(find(~isnan(combined_YXcy_yield),1))
    error('isempty(find(~isnan(combined_YXcy_yield),1))')
end

if isfield(yield_lpj,'list_to_map')
    yield_lpj_comb.list_to_map = yield_lpj.list_to_map ;
end
yield_lpj_comb.varNames = listCrops_lpj_comb ;
yield_lpj_comb.maps_YXvy = combined_YXcy_yield ;
yield_lpj_comb.yearList = yield_lpj.yearList ;

cropfrac_lpj_comb.list_to_map = cropfrac_lpj.list_to_map ;
cropfrac_lpj_comb.varNames = listCrops_lpj_comb ;
cropfrac_lpj_comb.maps_YXvy = combined_YXcy_cropfrac ;
cropfrac_lpj_comb.yearList = cropfrac_lpj.yearList ;

clear combined_YXcy*

% Calculate harvested totals
croparea_lpj_YXcy_comb = cropfrac_lpj_comb.maps_YXvy .* repmat(land_area_YX,[1 1 size(cropfrac_lpj_comb.maps_YXvy,3) size(cropfrac_lpj_comb.maps_YXvy,4)]) ;
total_lpj_YXcy_comb = yield_lpj_comb.maps_YXvy .* croparea_lpj_YXcy_comb ;

disp('Done.')