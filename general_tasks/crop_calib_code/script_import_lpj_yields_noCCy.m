% Setup
crops2remove = {'CC3G','CC4G','OtHr','ExtraCrop'} ;

if strcmp(verName_calib,'stijn_20180119')
    if calib_ver~=11
        error('When doing stijn_20180119, you must use calib_ver=11.')
    end
    listCrops_lpj_comb = {'TeWW','TeSW','TeCo','TrRi','TrSo','GlyM','FaBe'} ;
elseif calib_ver < 0
    error('calib_ver %d invalid', calib_ver)
elseif calib_ver == 11
    error('calib_ver 11 only works with stijn_20180119.')
elseif calib_ver == 12 || calib_ver == 15 || calib_ver == 16 || calib_ver == 18
    listCrops_lpj_comb = {'CerealsC3','CerealsC4','Rice','Oilcrops','StarchyRoots','Pulses'} ;
elseif calib_ver == 19
    listCrops_lpj_comb = ...
        {'CerealsC3','CerealsC4','Rice','Oilcrops',...
        'StarchyRoots','Pulses','Sugar','DateCitGrape'} ;
elseif any(calib_ver == [20 23 24])
    if ~exist('remapVer', 'var')
        remapVer = '8b' ;
    end
    switch remapVer
        case {'8b'}
            listCrops_lpj_comb = ...
                {'CerealsC3','CerealsC4','Rice','Oilcrops',...
                'StarchyRoots','Pulses','Sugar','FruitAndVeg'} ;
        case {'9', '9_g2p'}
            listCrops_lpj_comb = ...
                {'CerealsC3','CerealsC4','Rice','Oilcrops',...
                'StarchyRoots','Pulses','Sugarbeet','Sugarcane','FruitAndVeg'} ;
        case {'10', '10_g2p'}
            listCrops_lpj_comb = ...
                {'CerealsC3','CerealsC4','Rice','OilNfix','OilOther',...
                'StarchyRoots','Pulses','Sugarbeet','Sugarcane','FruitAndVeg'} ;
        case {'11_g2p', '12_g2p', '13_g2p'}
            listCrops_lpj_comb = ...
                {'CerealsC3s','CerealsC3w','CerealsC4','Rice','OilNfix','OilOther',...
                'StarchyRoots','Pulses','Sugarbeet','Sugarcane','FruitAndVeg'} ;
        otherwise
            error('Are Sugar and Oilcrops split up or not?')
    end
elseif calib_ver == 13
    listCrops_lpj_comb = {'Wheat','Maize','Sorghum','Pulses','Soybeans','Rice'} ;
elseif calib_ver == 14
    listCrops_lpj_comb = {'Wheat','Maize','Sorghum','Rice'} ;
elseif calib_ver <= 10 || calib_ver == 17
    listCrops_lpj_comb = {'TeWW','TeSW','TeCo','TrRi'} ;
elseif contains(verName_calib,'jianyong_20190128')
    if calib_ver~=21
        error('When doing %s, you must use calib_ver= 21.', verName_calib)
    end
    listCrops_lpj_comb = {'FaBe','GlyM','TeCo','TeSW','TeWW','TrRi','TrSo'} ;
else
    error(['calib_ver ' num2str(calib_ver) ' not recognized for setting listCrops_lpj_comb (verName_calib '  ').']) ;
end

% Import land uses
disp('Importing land uses and crop fractions...')
landuse_lpj = lpjgu_matlab_readTable_then2map(filename_guess_landuse,'force_mat_save',true, ...
    'lon_orient', 'center', 'lat_orient', 'center', ...
    'drop_northpole', drop_northpole, 'drop_southpole', drop_southpole, ...
    'lons_centered_on_180', lons_centered_on_180) ;
cropfrac_lpj = lpjgu_matlab_readTable_then2map(filename_guess_cropfrac,'force_mat_save',true, ...
    'lon_orient', 'center', 'lat_orient', 'center', ...
    'drop_northpole', drop_northpole, 'drop_southpole', drop_southpole, ...
    'lons_centered_on_180', lons_centered_on_180) ;

% Import yield
disp('Importing simulated yield...')
if exist('filename_guess_yield', 'var')
    if exist(filename_guess_yield, 'dir')
        starting_dir = pwd ;
        cd(filename_guess_yield)
        run_segments = dir() ;
        run_segments = run_segments(3:end) ;
        clear yield_lpj
        missing_years = listYears_fao ;
        for s = 1:length(run_segments)
            if ~run_segments(s).isdir
                continue
            end
            cd(run_segments(s).name)
            fprintf('    %s...\n', run_segments(s).name)

            % Find file
            output_dirs = dir() ;
            cd(output_dirs(end).name)
            if ~exist('yield_st.out.gz', 'file')
                error('%s: yield_st.out.gz not found', run_segments(s).name)
            end

            % Read file
            tmp = lpjgu_matlab_readTable_then2map('yield_st.out.gz', ...
                'force_mat_save', true, ...
                'lon_orient', 'center', 'lat_orient', 'center', ...
                'drop_northpole', drop_northpole, 'drop_southpole', drop_southpole, ...
                'lons_centered_on_180', lons_centered_on_180) ;
            [yearList_intersect, IA] = intersect(tmp.yearList, missing_years) ;
            if isempty(yearList_intersect)
                warning('Contains no missing years; skipping.')
                clear tmp
                cd(filename_guess_yield)
                continue
            end

            % Remove any factorial experiment stands (i.e., stands with names
            % containing a digit, after char'ing digits we actually care about)
            varNames_tmp = tmp.varNames ;
            varNames_tmp = strrep(varNames_tmp, 'CerealsC3', 'CerealsCthree') ;
            varNames_tmp = strrep(varNames_tmp, 'CerealsC4', 'CerealsCfour') ;
            remove = ~cellfun(@isempty, regexp(varNames_tmp, '\d')) ;
            tmp.maps_YXvy(:,:,remove,:) = [] ;
            tmp.varNames(remove) = [] ;
            clear varNames_tmp remove

            % Save years of interest
            if ~exist('yield_lpjg', 'var')
                % Create struct
                yield_lpj = rmfield(tmp, {'maps_YXvy', 'yearList'}) ;
                yield_lpj.yearList = yearList_intersect ;
                yield_lpj.maps_YXvy = tmp.maps_YXvy(:,:,:,IA) ;
            else
                % Ensure match of existing struct
                if ~isequal(tmp.list_to_map, yield_lpj.list_to_map)
                    error('list_to_map mismatch')
                elseif ~isequal(tmp.lonlats, yield_lpj.lonlats)
                    error('lonlats mismatch')
                elseif ~isequal(tmp.lat_extent, yield_lpj.lat_extent)
                    error('lat_extent mismatch')
                elseif ~isequal(tmp.lat_orient, yield_lpj.lat_orient)
                    error('lat_orient mismatch')
                elseif ~isequal(tmp.varNames, yield_lpj.varNames)
                    error('varNames mismatch')
                end

                % Add to existing struct
                yield_lpj.maps_YXvy = cat(4, yield_lpj.maps_YXvy, tmp.maps_YXvy(:,:,:,IA)) ;
                yield_lpj.yearList = sort(cat(1, yield_lpj.yearList, yearList_intersect)) ;
            end
            clear tmp
            missing_years = setdiff(listYears_fao, yield_lpj.yearList) ;

            % All years imported? Stop reading files.
            if isequal(shiftdim(yield_lpj.yearList), shiftdim(listYears_fao))
                break
            end

            cd(filename_guess_yield)
        end

        cd(starting_dir)
        disp('Done importing simulated yield.')
    elseif exist(filename_guess_yield, 'file')
        if strcmp(filename_guess_yield(end-2:end),'.gz')
            filename_guess_yield = filename_guess_yield(1:end-3) ;
        end
        yield_lpj = lpjgu_matlab_readTable_then2map(filename_guess_yield, ...
            'force_mat_save',true, ...
            'lon_orient', 'center', 'lat_orient', 'center', ...
            'drop_northpole', drop_northpole, 'drop_southpole', drop_southpole, ...
            'lons_centered_on_180', lons_centered_on_180) ;

        % Remove any factorial experiment stands (i.e., stands with names
        % containing a digit, after char'ing digits we actually care about)
        varNames_tmp = yield_lpj.varNames ;
        varNames_tmp = strrep(varNames_tmp, 'CerealsC3', 'CerealsCthree') ;
        varNames_tmp = strrep(varNames_tmp, 'CerealsC4', 'CerealsCfour') ;
        remove = ~cellfun(@isempty, regexp(varNames_tmp, '\d')) ;
        yield_lpj.maps_YXvy(:,:,remove,:) = [] ;
        yield_lpj.varNames(remove) = [] ;
        clear varNames_tmp remove
    else
        error('filename_guess_yield not found: %s', filename_guess_yield)
    end
    
    % Check that all years were read
    if indiv_years && ~isempty(setdiff(listYears_fao, yield_lpj.yearList))
        error('yield_lpj does not contain all years in listYears_fao')
    end
    
    % Deal with negative values, used to indicate that the crop was not
    % planted this year
    if ~isempty(find(yield_lpj.maps_YXvy < 0, 1))
        yield_lpj.maps_YXvy(yield_lpj.maps_YXvy < 0) = 0 ;
    end
    
    % Convert yield from kgDM/m2 to tDM/ha
    if isfield(yield_lpj,'maps_YXvy')
        yield_lpj.maps_YXvy = yield_lpj.maps_YXvy * 1e4 * 1e-3 ;
    else
        yield_lpj.maps_YXv = yield_lpj.maps_YXv * 1e4 * 1e-3 ;
    end
        
    % Use Oilcrops as proxy for FruitAndVeg and Sugar, if needed
    yieldCrops = get_cropsCombined(yield_lpj.varNames) ;
    missingCrops = setdiff(listCrops_lpj_comb, yieldCrops) ;
    if oilcrops_proxy_fruitveg_sugar
        % Make sure we can handle this setup
        if ~isequal({'FruitAndVeg', 'Sugar'}, missingCrops)
            error('yield_lpj is missing more crops than can be handled by oilcrops_proxy_fruitveg_sugar')
        end
        % Do it
        disp('Using Oilcrops yield for FruitAndVeg and Sugar...')
        yield_lpj.maps_YXvy = cat(3, ...
            yield_lpj.maps_YXvy, ...
            yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames, 'Oilcrops'),:), ...
            yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames, 'Oilcrops'),:), ...
            yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames, 'Oilcropsi'),:), ...
            yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames, 'Oilcropsi'),:)) ;
        yield_lpj.varNames = [yield_lpj.varNames {'FruitAndVeg', 'Sugar', 'FruitAndVegi', 'Sugari'}] ;
    
    % Use TeSW and TeCo as proxy for everything, if needed
    elseif exist('remapVer', 'var') && any(strcmp(remapVer, {'12_g2p', '13_g2p'}))
        Nmissing = length(missingCrops) ;
        % Set up substitute lists
        if any(strcmp(yield_lpj.varNames, 'Oilcrops'))
            subList_oil = {'OilNfix'} ;
            subList_c3s = {'FruitAndVeg', 'OilOther', 'Pulses', 'StarchyRoots', 'Sugarbeet'} ;
        else
            subList_oil = {} ;
            subList_c3s = {'FruitAndVeg', 'OilNfix', 'OilOther', 'Pulses', 'StarchyRoots', 'Sugarbeet'} ;
        end
        subList_c4s = {'Sugarcane'} ;
        % Get substitutes
        theseSubs = cell(Nmissing, 1) ;
        if ~isempty(subList_oil)
            [~, IB] = intersect(missingCrops, subList_oil) ;
            theseSubs(IB) = {'Oilcrops'} ;
        end
        [~, IB] = intersect(missingCrops, subList_c3s) ;
        theseSubs(IB) = {'CerealsC3s'} ;
        [~, IB] = intersect(missingCrops, subList_c4s) ;
        theseSubs(IB) = {'CerealsC4'} ;
        if any(cellfun(@isempty, theseSubs))
            error('Some unknown substitute(s)')
        end
        % Add new crop names
        yield_lpj.varNames = [yield_lpj.varNames ...
            missingCrops strcat(missingCrops, 'i')] ;
        % Substitute...
        addSize = size(yield_lpj.maps_YXvy) ;
        addSize(3) = Nmissing ;
        toAdd_YXvy = nan(addSize) ;
        clear addSize
        % ...rainfed
        if ~isempty(subList_oil)
            toSub = strcmp(theseSubs, 'Oilcrops') ;
            print_missing(toSub, missingCrops, 'Oilcrops')
            toAdd_YXvy(:,:,toSub,:) = ...
                repmat(yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames, 'Oilcrops'),:), ...
                [1 1 length(find(toSub)) 1]) ;
        end
        toSub = strcmp(theseSubs, 'CerealsC3s') ;
        print_missing(toSub, missingCrops, 'CerealsC3s')
        toAdd_YXvy(:,:,toSub,:) = ...
            repmat(yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames, 'CerealsC3s'),:), ...
            [1 1 length(find(toSub)) 1]) ;
        toSub = strcmp(theseSubs, 'CerealsC4') ;
        print_missing(toSub, missingCrops, 'CerealsC4')
        toAdd_YXvy(:,:,toSub,:) = ...
            repmat(yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames, 'CerealsC4'),:), ...
            [1 1 length(find(toSub)) 1]) ;
        yield_lpj.maps_YXvy = cat(3, yield_lpj.maps_YXvy, toAdd_YXvy) ;
        % ...irrigated
        if ~isempty(subList_oil)
            toSub = strcmp(theseSubs, 'Oilcrops') ;
            toAdd_YXvy(:,:,toSub,:) = ...
                repmat(yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames, 'Oilcrops'),:), ...
                [1 1 length(find(toSub)) 1]) ;
        end
        toSub = strcmp(theseSubs, 'CerealsC3s') ;
        toAdd_YXvy(:,:,toSub,:) = ...
            repmat(yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames, 'CerealsC3si'),:), ...
            [1 1 length(find(toSub)) 1]) ;
        toSub = strcmp(theseSubs, 'CerealsC4') ;
        toAdd_YXvy(:,:,toSub,:) = ...
            repmat(yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames, 'CerealsC4i'),:), ...
            [1 1 length(find(toSub)) 1]) ;
        % Remove Oilcrops, if needed
        if ~isempty(subList_oil)
            yield_lpj.maps_YXvy(:,:,contains(yield_lpj.varNames, 'Oilcrops'),:) = [] ;
            yield_lpj.varNames(contains(yield_lpj.varNames, 'Oilcrops')) = [] ;
        end
        % ...finish
        yield_lpj.maps_YXvy = cat(3, yield_lpj.maps_YXvy, toAdd_YXvy) ;
        clear toAdd_YXvy
    end
    
    % Make sure no crops are still missing
    yieldCrops = get_cropsCombined(yield_lpj.varNames) ;
    missingCrops = setdiff(listCrops_lpj_comb, yieldCrops) ;
    if ~isempty(missingCrops)
        error('Crop(s) missing from yield output!\nMissing\n    %s\nfrom\n   %s', ...
            strjoin(missingCrops, ', '), strjoin(yieldCrops, ', '))
    end
    
elseif exist('dirname_emuBL_yields', 'var')
    % GGCMI phase 2
    cropList_ggcmi = {'maize', 'rice', 'soy', 'spring_wheat', 'winter_wheat'} ;
    irrList = {'rf', 'ir'} ;
    Nirr_ggcmi = length(irrList) ;
    NcropIrr_ggcmi = length(cropList_ggcmi)*Nirr_ggcmi ;
    listCrops_ggcmi = cell(NcropIrr_ggcmi, 1) ;
    yield_ggcmi.maps_YXv = nan(360, 720, NcropIrr_ggcmi) ;
    
    % Import (already in tDM/ha)
    for c = 1:length(cropList_ggcmi)
        thisCrop = cropList_ggcmi{c} ;
        thisCrop_short = e2p_get_thisCrop_short(thisCrop) ;
        thisPattern = sprintf('%s%s*%s*', ...
            dirname_emuBL_yields, filesep, thisCrop) ;
        filelist = dir(thisPattern) ;
        if length(filelist) ~= 1
            error('Error finding emulated outputs for %s: expected 1, found %d', ...
                thisCrop, length(filelist))
        end
        thisFile = fullfile(filelist.folder, filelist.name) ;
        clear filelist
        for ii = 1:Nirr_ggcmi
            v = (c-1)*Nirr_ggcmi+ii ;
            listCrops_ggcmi{v} = sprintf('%s_%s', ...
                thisCrop, irrList{ii}) ;
            thisVar = sprintf('yield_%s_%s', irrList{ii}, thisCrop_short) ;
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

% Get maximum of spring/winter CerealsC3, if needed
remove = ~cellfun(@isempty, regexp(yield_lpj.varNames, 'CerealsC3[sw]')) ;
if any(remove)
    I_wheats_rf = ~cellfun(@isempty, regexp(yield_lpj.varNames, 'CerealsC3[sw]$')) ;
    I_wheats_ir = ~cellfun(@isempty, regexp(yield_lpj.varNames, 'CerealsC3[sw]i$')) ;
    if isfield(yield_lpj, 'maps_YXvy')
        maxWheat_rf_YX1y = max(yield_lpj.maps_YXvy(:,:,I_wheats_rf,:), [], 3) ;
        maxWheat_ir_YX1y = max(yield_lpj.maps_YXvy(:,:,I_wheats_ir,:), [], 3) ;
        yield_lpj.maps_YXvy(:,:,remove,:) = [] ;
        yield_lpj.varNames(remove) = [] ;
        yield_lpj.maps_YXvy = cat(3, yield_lpj.maps_YXvy, ...
            maxWheat_rf_YX1y, maxWheat_ir_YX1y) ;
    else
        maxWheat_rf_YX = max(yield_lpj.maps_YXv(:,:,I_wheats_rf), [], 3) ;
        maxWheat_ir_YX = max(yield_lpj.maps_YXv(:,:,I_wheats_ir), [], 3) ;
        yield_lpj.maps_YXv(:,:,remove) = [] ;
        yield_lpj.varNames(remove) = [] ;
        yield_lpj.maps_YXv = cat(3, yield_lpj.maps_YXv, ...
            maxWheat_rf_YX, maxWheat_ir_YX) ;
    end
    yield_lpj.varNames = [yield_lpj.varNames {'CerealsC3' 'CerealsC3i'}] ;
    clear I_wheats* maxWheat_*
    remove = ~cellfun(@isempty, regexp(listCrops_lpj_comb, 'CerealsC3[sw]')) ;
    listCrops_lpj_comb(remove) = [] ;
    listCrops_lpj_comb{end+1} = 'CerealsC3' ;
end
clear remove

% Start dealing with sugar and oilcrops, if needed
listCrops_fa2o_tmp = get_FAO_info(calib_ver) ;
% Sugar
supercrop = 'Sugar' ;
combine_sugars = any(contains(listCrops_lpj_comb, ...
    {supercrop, 'Sugarbeet', 'Sugarcane'})) ...
    && length(intersect(listCrops_fa2o_tmp, ...
    {'Sugarbeet', 'Sugarcane'})) < 2 ;
if combine_sugars
    [sugrealsIrr_lpj, sugreals] = ...
        combine_subcrops_step1( ...
        supercrop, irrList, cropfrac_lpj.varNames, ...
        listCrops_lpj_comb) ;
end
% Oilcrops
supercrop = 'Oilcrops' ;
combine_oilcrops = any(contains(listCrops_lpj_comb, ...
    {'OilNfix', 'OilOther'})) ...
    && length(intersect(listCrops_fa2o_tmp, ...
    {'OilNfix', 'OilOther'})) < 2 ;
if combine_oilcrops
    [oilrealsIrr_lpj, oilreals] = ...
        combine_subcrops_step1( ...
        supercrop, irrList, cropfrac_lpj.varNames, ...
        listCrops_lpj_comb) ;
end
clear listCrops_fa2o_tmp

% Sanity check
if (isfield(yield_lpj, 'maps_YXvy') && ~isempty(find(yield_lpj.maps_YXvy < 0, 1))) ...
|| (isfield(yield_lpj, 'maps_YXv') && ~isempty(find(yield_lpj.maps_YXv < 0, 1)))
    error('Negative value(s) in yield_lpj')
end

disp('Processing...')

% For GGCMI, we're not comparing individual years, so maps_YXvy should
% just be means for years of interest (i.e., 4th dim of length 1)
if is_ggcmi 
    
    % yield
    if indiv_years
        if ~(isfield(yield_lpj,'yearList') && isfield(yield_lpj,'maps_YXvy'))
            error('indiv_years specified but yield_lpj.maps_YXvy not found')
        elseif ~isempty(setxor(yield_lpj.yearList, listYears_fao))
            [~, ~, IB] = intersect(listYears_fao, yield_lpj.yearList) ;
            yield_lpj.maps_YXvy = yield_lpj.maps_YXvy(:,:,:,IB) ;
            yield_lpj.yearList = listYears_fao ;
        end
    else
        if strcmp(model_name, 'LPJ-GUESS-sim')
            yield_lpj.maps_YXvy = nanmean(yield_lpj.maps_YXvy,4) ;
        else
            if isfield(yield_lpj,'yearList') || isfield(yield_lpj,'maps_YXvy')
                error('GGCMI yields have years??')
            end
            yield_lpj.maps_YXvy = yield_lpj.maps_YXv ;
            yield_lpj = rmfield(yield_lpj, 'maps_YXv') ;
        end
        yield_lpj.yearList = -pi ;
    end
    
    % landuse
    if indiv_years
        if ~isfield(landuse_lpj,'maps_YXvy')
            landuse_lpj.maps_YXvy = repmat(landuse_lpj.maps_YXv, [1 1 1 Nyears_fao]) ;
            landuse_lpj = rmfield(landuse_lpj, 'maps_YXv') ;
            landuse_lpj.yearList = listYears_fao ;
        elseif ~isempty(setxor(landuse_lpj.yearList, listYears_fao))
            [~, ~, IB] = intersect(listYears_fao, landuse_lpj.yearList) ;
            landuse_lpj.maps_YXvy = landuse_lpj.maps_YXvy(:,:,:,IB) ;
            landuse_lpj.yearList = listYears_fao ;
        end
    else
        if isfield(landuse_lpj,'maps_YXvy')
            landuse_lpj.maps_YXvy = nanmean(landuse_lpj.maps_YXvy(:,:,:, ...
                landuse_lpj.yearList>=year1 & landuse_lpj.yearList<=yearN), 4) ;
        else
            landuse_lpj.maps_YXvy = landuse_lpj.maps_YXv ;
            landuse_lpj = rmfield(landuse_lpj, 'maps_YXv') ;
        end
        landuse_lpj.yearList = -pi ;
    end
    
    % cropfrac
    if indiv_years
        if ~isfield(cropfrac_lpj,'maps_YXvy')
            cropfrac_lpj.maps_YXvy = repmat(cropfrac_lpj.maps_YXv, [1 1 1 Nyears_fao]) ;
            cropfrac_lpj = rmfield(cropfrac_lpj, 'maps_YXv') ;
            cropfrac_lpj.yearList = listYears_fao ;
        elseif ~isempty(setxor(cropfrac_lpj.yearList, listYears_fao))
            [~, ~, IB] = intersect(listYears_fao, cropfrac_lpj.yearList) ;
            cropfrac_lpj.maps_YXvy = cropfrac_lpj.maps_YXvy(:,:,:,IB) ;
            cropfrac_lpj.yearList = listYears_fao ;
        end
    else
        if isfield(cropfrac_lpj,'maps_YXvy')
            cropfrac_lpj.maps_YXvy = nanmean(cropfrac_lpj.maps_YXvy(:,:,:, ...
                cropfrac_lpj.yearList>=year1 & cropfrac_lpj.yearList<=yearN), 4) ;
        else
            cropfrac_lpj.maps_YXvy = cropfrac_lpj.maps_YXv ;
            cropfrac_lpj = rmfield(cropfrac_lpj, 'maps_YXv') ;
        end
        cropfrac_lpj.yearList = -pi ;
    end
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
if ~isempty(intersect(cropfrac_lpj.varNames,{'CerealsC3w','CerealsC3s'}))
    warning('Combining WW and SW.')
    if is_ggcmi % You've already combined yields, above
        % Rainfed
        i_fracs_rf = [find(strcmp(cropfrac_lpj.varNames,'CerealsC3w')) find(strcmp(cropfrac_lpj.varNames,'CerealsC3s'))] ;
        if length(i_fracs_rf)~=2; error('length(i_fracs_rf)~=2'); end
        tmp_fracs_rf_YXvy = cropfrac_lpj.maps_YXvy(:,:,i_fracs_rf,:) ;
        % Irrigated
        i_fracs_ir= [find(strcmp(cropfrac_lpj.varNames,'CerealsC3wi')) find(strcmp(cropfrac_lpj.varNames,'CerealsC3si'))] ;
        if length(i_fracs_ir)~=2; error('length(i_fracs_ir)~=2'); end
        tmp_fracs_ir_YXvy = cropfrac_lpj.maps_YXvy(:,:,i_fracs_ir,:) ;
        % Combine
        tmp_fracs_rf_YX_y = nansum(tmp_fracs_rf_YXvy,3) ;
        tmp_fracs_ir_YX_y = nansum(tmp_fracs_ir_YXvy,3) ;
        % Replace
        cropfrac_lpj.maps_YXvy(:,:,[i_fracs_rf i_fracs_ir],:) = [] ;
        cropfrac_lpj.varNames([i_fracs_rf i_fracs_ir]) = [] ;
        cropfrac_lpj.varNames = [cropfrac_lpj.varNames 'CerealsC3' 'CerealsC3i'] ;
        cropfrac_lpj.maps_YXvy = cat(3,cropfrac_lpj.maps_YXvy,tmp_fracs_rf_YX_y,tmp_fracs_ir_YX_y) ;
    else
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
end

% Stijn didn't include irrigated of these?
if ~is_ggcmi && (strcmp(verName_calib,'stijn_20180119') || contains(verName_calib,'jianyong_20190128'))
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

% Remove Miscanthus area
if ~any(strcmp(listCrops_lpj_comb, 'Miscanthus'))
    if any(contains(cropfrac_lpj.varNames, 'Miscanthus'))
        if isfield(yield_lpj, 'maps_YXv')
            if any(any(cropfrac_lpj.maps_YXv(:,:,contains(cropfrac_lpj.varNames, 'Miscanthus')) > 1e-6))
                error('Large amounts of Miscanthus being removed (max %g)', ...
                    max(max(cropfrac_lpj.maps_YXv(:,:,contains(cropfrac_lpj.varNames, 'Miscanthus'))))) ;
            end
            cropfrac_lpj.maps_YXv(:,:,contains(cropfrac_lpj.varNames, 'Miscanthus')) = [] ;
        end
        if isfield(yield_lpj, 'maps_YXvy')
            if any(any(any(cropfrac_lpj.maps_YXvy(:,:,contains(cropfrac_lpj.varNames, 'Miscanthus'),:) > 1e-6)))
                error('Large amounts of Miscanthus being removed (max %g)', ...
                    max(max(max(cropfrac_lpj.maps_YXvy(:,:,contains(cropfrac_lpj.varNames, 'Miscanthus'),:))))) ;
            end
            cropfrac_lpj.maps_YXvy(:,:,contains(cropfrac_lpj.varNames, 'Miscanthus'),:) = [] ;
        end
        cropfrac_lpj.varNames(contains(cropfrac_lpj.varNames, 'Miscanthus')) = [] ;
    end
    if any(contains(yield_lpj.varNames, 'Miscanthus'))
        if isfield(yield_lpj, 'maps_YXv')
            yield_lpj.maps_YXv(:,:,contains(yield_lpj.varNames, 'Miscanthus')) = [] ;
        end
        if isfield(yield_lpj, 'maps_YXvy')
            yield_lpj.maps_YXvy(:,:,contains(yield_lpj.varNames, 'Miscanthus'),:) = [] ;
        end
        yield_lpj.varNames(contains(yield_lpj.varNames, 'Miscanthus')) = [] ;
    end
end

% Make sure cropfrac_lpj and yield_lpj contain the same variables
[C,IA,IB] = setxor(cropfrac_lpj.varNames, yield_lpj.varNames) ;
if ~isempty(C)
    error('Crops differ between cropfrac_lpj and yield_lpj.\nMissing from yield_lpj:\n    %s\nMissing from cropfrac_lpj:\n    %s', ...
        strjoin(cropfrac_lpj.varNames(IA), ', '), strjoin(yield_lpj.varNames(IB), ', '))
end

% Rearrange cropfrac_lpj so that variables are in the same order as
% yield_lpj variables
if ~isequal(sort(listCrops_lpj_comb), get_cropsCombined(yield_lpj.varNames))
    error('Crop list mismatch') ;
end
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

if ~isequal(size(yield_lpj.maps_YXvy), size(cropfrac_lpj.maps_YXvy))
    error('Size mismatch between yield_lpj.maps_YXvy and cropfrac_lpj.maps_YXvy')
end

% Deal with missing simulated cells
tmp_croparea_YXvy = cropfrac_lpj.maps_YXvy .* repmat(land_area_YX, ...
    [1 1 size(cropfrac_lpj.maps_YXvy,3) size(cropfrac_lpj.maps_YXvy,4)]) ;
% NaN in simulation, but 0 crop area anyway
if do_remove_area_dueto_NaNsim && find(isnan(yield_lpj.maps_YXvy) & tmp_croparea_YXvy==0, 1)
    yield_lpj.maps_YXvy(isnan(yield_lpj.maps_YXvy) & tmp_croparea_YXvy==0) = 0 ;
end
% NaN in simulation but there is some crop area
removed_area_dueto_NaNsim = do_remove_area_dueto_NaNsim && find( ...
    isnan(yield_lpj.maps_YXvy) & tmp_croparea_YXvy>0, 1) ;
if removed_area_dueto_NaNsim
    croparea_lpj_removed = rmfield(cropfrac_lpj, 'maps_YXvy') ;
    croparea_lpj_removed.maps_YXvy = zeros(size(cropfrac_lpj.maps_YXvy)) ;
    for c = 1:length(yield_lpj.varNames)
        thisCrop = yield_lpj.varNames{c} ;
        tmp_yield_YX1y = yield_lpj.maps_YXvy(:,:,c,:) ;
        tmp_croparea_YX1y = tmp_croparea_YXvy(:,:,c,:) ;
        I_bad = find(isnan(tmp_yield_YX1y) & tmp_croparea_YX1y>0) ;
        if ~isempty(I_bad)
            % Track how much area is being removed
            removed_YX1y = zeros(size(tmp_croparea_YX1y)) ;
            removed_YX1y(I_bad) = tmp_croparea_YX1y(I_bad) ;
            croparea_lpj_removed.maps_YXvy(:,:,c,:) = removed_YX1y ;
            clear removed_YX1y
        end
        clear *_YX1y I_bad
    end
    % Remove it
    cropfrac_lpj.maps_YXvy(croparea_lpj_removed.maps_YXvy > 0) = 0 ;
    yield_lpj.maps_YXvy(croparea_lpj_removed.maps_YXvy > 0) = 0 ;
else
    croparea_lpj_removed = [] ;
end
clear tmp_croparea_YXvy

% Combine Sugarbeet and Sugarcane, if needed
if combine_sugars
    [cropfrac_lpj, yield_lpj, croparea_lpj_removed, ...
    listCrops_lpj_comb] = ...
    combine_subcrops_step2( ...
    'Sugar', sugreals, sugrealsIrr_lpj, listCrops_lpj_comb, ...
    cropfrac_lpj, yield_lpj, croparea_lpj_removed) ;
end

% Combine OilNfix and OilOther, if needed
if combine_oilcrops
    [cropfrac_lpj, yield_lpj, croparea_lpj_removed, ...
    listCrops_lpj_comb] = ...
    combine_subcrops_step2( ...
    'Oilcrops', oilreals, oilrealsIrr_lpj, listCrops_lpj_comb, ...
    cropfrac_lpj, yield_lpj, croparea_lpj_removed) ;
end

% Create combined irrigated+rainfed yields
Ncrops_lpj_comb = length(listCrops_lpj_comb) ;
combined_YXcy_yield = nan(size(cropfrac_lpj.maps_YXvy,1),...
                    size(cropfrac_lpj.maps_YXvy,2),...
                    Ncrops_lpj_comb,...
                    length(cropfrac_lpj.yearList)) ;
combined_YXcy_cropfrac = combined_YXcy_yield ;
if removed_area_dueto_NaNsim
    cropareaRemoved_lpj_YXcy_comb = zeros(size(combined_YXcy_yield)) ;
else
    cropareaRemoved_lpj_YXcy_comb = [] ;
end
for c = 1:Ncrops_lpj_comb
    
    % Find indices
    thisCropR = listCrops_lpj_comb{c} ;
    thisCropI = [thisCropR 'i'] ;
    iR_cropfrac = find(strcmp(cropfrac_lpj.varNames,thisCropR)) ;
    if ~isequal(iR_cropfrac, exact_string_in_cellarray(cropfrac_lpj.varNames,thisCropR))
        % exact_string_in_cellarray(cropfrac_lpj.varNames,thisCropR)
        % was previously used in getting frac_comb under "else" below
        error('???')
    end
    iI_cropfrac = find(strcmp(cropfrac_lpj.varNames,thisCropI)) ;
    if ~isequal(iI_cropfrac, exact_string_in_cellarray(cropfrac_lpj.varNames,thisCropI))
        % exact_string_in_cellarray(cropfrac_lpj.varNames,thisCropI)
        % was previously used in getting frac_comb under "else" below
        error('???')
    end
    
    % Get rainfed and irrigated weights
    if strcmp(verName_calib,'stijn_20180119') || contains(verName_calib,'jianyong_20190128')
        frac_comb = cropfrac_lpj.maps_YXvy(:,:,iR_cropfrac,:) ;
        if ~isempty(iI_cropfrac)
            frac_comb = frac_comb + cropfrac_lpj.maps_YXvy(:,:,iI_cropfrac,:) ;
        else
            warning([thisCropI ' not found in cropfrac_lpj.varNames.'])
        end
    else
        frac_comb = cropfrac_lpj.maps_YXvy(:,:,iR_cropfrac,:) ...
                  + cropfrac_lpj.maps_YXvy(:,:,iI_cropfrac,:) ;
    end
    weights_rf_YXvy = cropfrac_lpj.maps_YXvy(:,:,iR_cropfrac,:) ./ frac_comb ;
    weights_ir_YXvy = cropfrac_lpj.maps_YXvy(:,:,iI_cropfrac,:) ./ frac_comb ;
    
    % Combine rainfed and irrigated
    iR_yield = find(strcmp(yield_lpj.varNames,thisCropR)) ;
    iI_yield = find(strcmp(yield_lpj.varNames,thisCropI)) ;
    if isempty(iI_cropfrac) && (strcmp(verName_calib,'stijn_20180119') || contains(verName_calib,'jianyong_20190128'))
        tmp = yield_lpj.maps_YXvy(:,:,iR_yield,:) .* weights_rf_YXvy ;
    else
        tmp = yield_lpj.maps_YXvy(:,:,iR_yield,:) .* weights_rf_YXvy ...
            + yield_lpj.maps_YXvy(:,:,iI_yield,:) .* weights_ir_YXvy ;
    end
    clear weights_*_YXvy
    if removed_area_dueto_NaNsim
        tmp(frac_comb==0) = 0 ;
    end
    combined_YXcy_yield(:,:,c,:) = tmp ;
    clear tmp
    combined_YXcy_cropfrac(:,:,c,:) = frac_comb ;
    clear frac_comb
    
    % Check whether this process set any positive yields to NaN
    check_nanified_yield(yield_lpj.maps_YXvy(:,:,iR_yield,:), combined_YXcy_yield(:,:,c,:), thisCropR) ;
    check_nanified_yield(yield_lpj.maps_YXvy(:,:,iI_yield,:), combined_YXcy_yield(:,:,c,:), thisCropI) ;
    
    if removed_area_dueto_NaNsim
        cropareaRemoved_lpj_YXcy_comb(:,:,c,:) = ...
            sum(croparea_lpj_removed.maps_YXvy(:,:,[iR_cropfrac iI_cropfrac],:),3) ;
    end
    clear tmp
end

if isempty(find(~isnan(combined_YXcy_yield),1))
    error('isempty(find(~isnan(combined_YXcy_yield),1))')
end

% Set up "combined" arrays
if isfield(yield_lpj,'list_to_map')
    yield_lpj_comb.list_to_map = yield_lpj.list_to_map ;
end
yield_lpj_comb.yearList = yield_lpj.yearList ;
% clear yield_lpj
yield_lpj_comb.varNames = listCrops_lpj_comb ;
yield_lpj_comb.maps_YXvy = combined_YXcy_yield ;
% clear combined_YXcy_yield
cropfrac_lpj_comb.list_to_map = cropfrac_lpj.list_to_map ;
cropfrac_lpj_comb.yearList = cropfrac_lpj.yearList ;
% clear cropfrac_lpj
cropfrac_lpj_comb.varNames = listCrops_lpj_comb ;
cropfrac_lpj_comb.maps_YXvy = combined_YXcy_cropfrac ;
% clear combined_YXcy_cropfrac

% Calculate harvested totals
croparea_lpj_YXcy_comb = cropfrac_lpj_comb.maps_YXvy .* repmat(land_area_YX,[1 1 size(cropfrac_lpj_comb.maps_YXvy,3) size(cropfrac_lpj_comb.maps_YXvy,4)]) ;
total_lpj_YXcy_comb = yield_lpj_comb.maps_YXvy .* croparea_lpj_YXcy_comb ;

% Warn about removals, if necessary
if removed_area_dueto_NaNsim
    tmp_varNames = pad(strcat(listCrops_lpj_comb,':')) ;
    for c = 1:Ncrops_lpj_comb
        thisCrop = tmp_varNames{c} ;
        cropareaRemoved_lpj_YX1y_comb = cropareaRemoved_lpj_YXcy_comb(:,:,c,:) ;
        croparea_lpj_YX1y_comb = croparea_lpj_YXcy_comb(:,:,c,:) ...
            + cropareaRemoved_lpj_YX1y_comb;
        croparea_removed_ha = sum(cropareaRemoved_lpj_YX1y_comb(~isnan(cropareaRemoved_lpj_YX1y_comb))) ;
        croparea_removed_pct = 100 * croparea_removed_ha ...
            / sum(croparea_lpj_YX1y_comb(~isnan(croparea_lpj_YX1y_comb))) ;

        warning(['NaN sim yield -> 0 area: %s %0.1f%% area ' ...
            '(%0.2g ha) '], ...
            thisCrop, croparea_removed_pct, croparea_removed_ha) ;
        clear area_removed pct_removed croparea_lpj_YX1y_comb
    end
    clear tmp_varNames
end

disp('Done.')


%% FUNCTIONS

function list_out = get_cropsCombined(list_in)

list_out = list_in(cellfun(@isempty, regexp(list_in, 'i$'))) ;
list_out = setdiff(list_out, {'CC3G_ic', 'CC4G_ic', 'ExtraCrop', 'Miscanthus'}) ;
        
end


function [wt_sub1_YX, wt_sub2_YX] = get_subcrop_wts(sub1_YX, sub2_YX, supreals)

total_YX = sub1_YX + sub2_YX ;
total_YX(total_YX==0) = 1 ;

wt_sub1_YX = sub1_YX ./ total_YX ;
wt_sub2_YX = sub2_YX ./ total_YX ;

if any(wt_sub1_YX(~isnan(wt_sub1_YX)) + wt_sub2_YX(~isnan(wt_sub2_YX)) - 1 > eps) 
    error('%s + %s weights don''t sum to 1 within eps', ...
        supreals{1}, supreals{2})
end

end


function [suprealsIrr_lpj, supreals] = ...
    combine_subcrops_step1( ...
    supercrop, irrList, cropfrac_lpj_varNames, ...
    listCrops_lpj_comb)

Nirr_ggcmi = length(irrList) ;

% Construct types
if strcmp(supercrop, 'Sugar')
    supreals = {'Sugarbeet', 'Sugarcane'} ;
    supfakes = {'spring_wheat', 'maize'} ;
elseif strcmp(supercrop, 'Oilcrops')
    supreals = {'OilNfix', 'OilOther'} ;
    supfakes = {'soybean', 'spring_wheat'} ;
else
    error('Error trying to construct subcrop types: supercrop %s not recognized', ...
        supercrop)
end
Nsub = length(supreals) ;
suprealsIrr_ggcmi = cell(Nsub*Nirr_ggcmi, 1) ;
for c = 1:Nsub
    thisReal = supreals{c} ;
    thisFake = supfakes{c} ;
    for ii = 1:Nirr_ggcmi
        thisIrr = irrList{ii} ;
        thisRealIrr = sprintf('%s_%s', thisReal, thisIrr) ;
        thisFakeIrr = sprintf('%s_%s', thisFake, thisIrr) ;
        v = (c-1)*Nirr_ggcmi + ii ;
        suprealsIrr_ggcmi{v} = thisRealIrr ;
    end
end
suprealsIrr_lpj = suprealsIrr_ggcmi' ;
suprealsIrr_lpj = strrep(suprealsIrr_lpj, '_rf', '') ;
suprealsIrr_lpj = strrep(suprealsIrr_lpj, '_ir', 'i') ;

if length(intersect(cropfrac_lpj_varNames, suprealsIrr_lpj)) == length(suprealsIrr_lpj)
    if length(intersect(listCrops_lpj_comb, supreals)) ~= length(supreals)
        error('%s supreals not found in listCrops_lpj_comb', supercrop)
    end
else
    error('Do I need to split area map of %s into %s and %s?', ...
        supercrop, supreals{1}, supreals{2})
end


end


function [cropfrac_lpj, yield_lpj, croparea_lpj_removed, ...
    listCrops_lpj_comb] = ...
    combine_subcrops_step2( ...
    supercrop, supreals, suprealsIrr_lpj, listCrops_lpj_comb, ...
    cropfrac_lpj, yield_lpj, croparea_lpj_removed)

supreal1 = supreals{1} ;
supreal1i = [supreals{1} 'i'] ;
supreal2 = supreals{2} ;
supreal2i = [supreals{2} 'i'] ;

removed_area_dueto_NaNsim = ~isempty(croparea_lpj_removed) ;

% Remove subcrops
[~, ~, I_sup] = intersect(suprealsIrr_lpj, cropfrac_lpj.varNames, 'stable') ;
sup1_frac_YXvy = cropfrac_lpj.maps_YXvy(:,:,I_sup,:) ;
[wt_sub1_rf_YX, wt_sub2_rf_YX] = get_subcrop_wts( ...
    sup1_frac_YXvy(:,:,strcmp(suprealsIrr_lpj, supreal1)), ...
    sup1_frac_YXvy(:,:,strcmp(suprealsIrr_lpj, supreal2)), ...
    supreals) ;
[wt_sub1_ir_YX, wt_sub2_ir_YX] = get_subcrop_wts( ...
    sup1_frac_YXvy(:,:,strcmp(suprealsIrr_lpj, supreal1i)), ...
    sup1_frac_YXvy(:,:,strcmp(suprealsIrr_lpj, supreal2i)), ...
    supreals) ;
sup1_yield_YXvy = yield_lpj.maps_YXvy(:,:,I_sup,:) ;
if removed_area_dueto_NaNsim
    sup1_areaRemoved_YXvy = croparea_lpj_removed.maps_YXvy(:,:,I_sup,:) ;
end
cropfrac_lpj.maps_YXvy(:,:,I_sup,:) = [] ;
cropfrac_lpj.varNames(I_sup) = [] ;
yield_lpj.maps_YXvy(:,:,I_sup,:) = [] ;
yield_lpj.varNames(I_sup) = [] ;
if removed_area_dueto_NaNsim
    croparea_lpj_removed.maps_YXvy(:,:,I_sup,:) = [] ;
    croparea_lpj_removed.varNames(I_sup) = [] ;
end
listCrops_lpj_comb(contains(listCrops_lpj_comb, supreals)) = [] ;

% Combine fractions/areas
[~, I_rf] = intersect(suprealsIrr_lpj, supreals) ;
[~, I_ir] = intersect(suprealsIrr_lpj, strcat(supreals, 'i')) ;
sup2_frac_YXvy = cat(3, ...
    sum(sup1_frac_YXvy(:,:,I_rf,:),3), ...
    sum(sup1_frac_YXvy(:,:,I_ir,:),3)) ;
if removed_area_dueto_NaNsim
    sup2_areaRemoved_YXvy = cat(3, ...
        sum(sup1_areaRemoved_YXvy(:,:,I_rf,:),3), ...
        sum(sup1_areaRemoved_YXvy(:,:,I_ir,:),3)) ;
end

% Combine yields
sup2_yield_YXvy = cat(3, ...
    sup1_yield_YXvy(:,:,strcmp(suprealsIrr_lpj, supreal1)) .* wt_sub1_rf_YX ...
    + sup1_yield_YXvy(:,:,strcmp(suprealsIrr_lpj, supreal2)) .* wt_sub2_rf_YX, ...
    sup1_yield_YXvy(:,:,strcmp(suprealsIrr_lpj, supreal1i)) .* wt_sub1_ir_YX ...
    + sup1_yield_YXvy(:,:,strcmp(suprealsIrr_lpj, supreal2i)) .* wt_sub2_ir_YX) ;

% Add back
to_add = {supercrop, [supercrop 'i']} ;
cropfrac_lpj.varNames = [cropfrac_lpj.varNames to_add] ;
cropfrac_lpj.maps_YXvy = cat(3, ...
    cropfrac_lpj.maps_YXvy, sup2_frac_YXvy) ;
yield_lpj.varNames = [yield_lpj.varNames to_add] ;
yield_lpj.maps_YXvy = cat(3, ...
    yield_lpj.maps_YXvy, sup2_yield_YXvy) ;
if removed_area_dueto_NaNsim
    croparea_lpj_removed.varNames = [croparea_lpj_removed.varNames to_add] ;
    croparea_lpj_removed.maps_YXvy = cat(3, ...
        croparea_lpj_removed.maps_YXvy, sup2_areaRemoved_YXvy) ;
end
listCrops_lpj_comb = [listCrops_lpj_comb supercrop] ;


end


function print_missing(toSub, missingCrops, thisSub)

toSub_I = find(toSub) ;
toSub_txt = '' ;
for c = 1:length(toSub_I)
    C = toSub_I(c) ;
    toSub_txt = [toSub_txt missingCrops{C}] ;
    if c < length(toSub_I)
        toSub_txt = [toSub_txt ', '] ;
    end
end

fprintf('Using %s yield for: %s\n', thisSub, toSub_txt) ;

end


function check_nanified_yield(yield1, yield2, thisCrop)

isbad = yield1 > 0 & isnan(yield2) ;
Nbad = length(find(isbad)) ;
if Nbad > 0
    worst = max(yield1(isbad)) ;
    warning('%s: %d positive yields (max %0.3f) NaN-ified after combining rf+ir', ...
        thisCrop, Nbad, worst)
end

end


