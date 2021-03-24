%% Setup

% Define reference file names
if strcmp(thisSystem, 'ssr_mac')
    crop_calib_data_dir = '/Users/Shared/PLUM/crop_calib_data' ;
    hotspot_tif = '/Users/sam/Geodata/BiodiversityHotspotsRevisited_ConservationInternational_2004/data/hotspots_revisited_2004.outerlimit.tif' ;
    hotspot_shp = '/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/LPJGP_paper02_Sam/hotspots_clipByGridlist.shp' ;
    ecoid_file = '/Users/sam/Geodata/General/WWF terrestrial ecosystems/wwf_terr_ecos_UnpackClip.halfDeg.ECO_ID.tif' ;
    cslf_shp = '/Users/sam/Geodata/General/WWF terrestrial ecosystems/wwf_terr_ecos_UnpackClip.halfDeg.CSLF.shp' ;
    fpu_file = '/Users/Shared/PLUM/food_production_units/FPU.asc' ;
    pop_file = '/Users/Shared/PLUM/other_plum_data/SSP/ssp.csv' ;
else
    crop_calib_data_dir = '~/PLUM/crop_calib_data' ;
    other_data_dir = '~/PLUM/other_data' ;
    if ~exist(other_data_dir, 'dir')
        error('other_data_dir not found: %s', other_data_dir)
    end
    hotspot_tif = sprintf('%s/hotspots_revisited_2004.outerlimit.tif', other_data_dir) ;
    hotspot_shp = sprintf('%s/hotspots_clipByGridlist.shp', other_data_dir) ;
    ecoid_file = sprintf('%s/wwf_terr_ecos_UnpackClip.halfDeg.ECO_ID.tif', other_data_dir) ;
    cslf_shp = sprintf('%s/wwf_terr_ecos_UnpackClip.halfDeg.CSLF.shp', other_data_dir) ;
    fpu_file = sprintf('%s/FPU.asc', other_data_dir) ;
    pop_file = sprintf('%s/ssp.csv', other_data_dir) ;
end
if ~exist(crop_calib_data_dir, 'dir')
    error('crop_calib_data_dir not found: %s', crop_calib_data_dir)
end
landarea_file = sprintf('%s/input_data/staticData_quarterdeg.nc', plumharm_repo_path) ;
continents_shp = sprintf('%s/input_data/continents_from_countries.shp', plumharm_repo_path) ;
fao_file = sprintf('%s/fao/FAOdata_1961-2010_calibVer16_Production.mat', ...
            crop_calib_data_dir) ;

% Define raster reference object and missing value
map_size = [360 720] ;
R = georasterref('RasterSize', map_size, ...
    'RasterInterpretation', 'cells', ...
    'ColumnsStartFrom', 'north', ...
    'LatitudeLimits', [-90 90], ...
    'LongitudeLimits', [-180 180]) ;
gtif_missing = -1e20 ;

% Define function to calculate sem
sem_ssr = @(data,yrs) std(data,find(size(data)==length(yrs))) / sqrt(length(yrs)) ;

% Define special characters
en_dash = char(8211) ;
plusminus = char(177) ;

include_fao = true ;
if strcmp(thisVer,'20180424agmip7') || strcmp(thisVer,'20180424agmip7_asPLUMout2011-2015')
    warning('not including fao data!')
    include_fao = false ;
end

define_runDirs_etc

% Deal with years
if max(yearList_baseline) >= min(yearList_future)
    error('Baseline and future yearLists overlap!')
elseif max(yearList_baseline) ~= min(yearList_future)-1
    error('Baseline and future yearLists are offset!')
end
Nyears_bl = length(yearList_baseline) ;
Nyears_fu = length(yearList_future) ;

% Output directory
if strcmp(thisSystem, 'ssr_mac')
    outDir_base = addslashifneeded(['~/Documents/Dropbox/2016_KIT/LandSyMM/LPJGP_paper02_Sam/'...
        'figures_' thisVer '_SI']) ;
else
    outDir_base = addslashifneeded(sprintf('~/PLUM/outputs/%s/figs', thisVer)) ;
end
outDir_maps = addslashifneeded([outDir_base 'maps']) ;
outDir_gtif = addslashifneeded([outDir_base 'gtif']) ;
outDir_ts = [outDir_base 'TS'] ;
if rebase
    outDir_ts = addslashifneeded([outDir_ts '_rebased']) ;
else
    outDir_ts = addslashifneeded(outDir_ts) ;
end

if ~exist(outDir_ts,'dir')
    mkdir(outDir_ts)
end
if ~exist(outDir_maps,'dir')
    mkdir(outDir_maps)
end
if ~exist(outDir_gtif,'dir')
    mkdir(outDir_gtif)
end

% Conversion factors
cf_kg2Mt = 1e-3*1e-6 ;
cf_t2kg = 1e3 ;   % For FAO only
cf_ha2m2 = 1e4 ; % For FAO only
cf_kgPm2_to_tonsPha = 1e-3*1e4 ;
cf_kg2Tg = 1e-9 ;
cf_kg2Pg = 1e-12 ;
cf_m3_to_km3 = (1e-3)^3 ;
cf_kcalEcal = 1e-15 ;

% Compare run list with list of directories
Nruns = length(runList) ;
if Nruns ~= length(runDirs)
    error('Length mismatch between runList and runDirs!')
end

% Get complete paths to baseline and future runs
baselineDir = find_PLUM2LPJG_run(baselineDir) ;
for r = 1:Nruns
    runDirs{r} = find_PLUM2LPJG_run(runDirs{r}) ;
end

stdLegend = ['Baseline',runList] ;
stdLegend_plusFAO = ['Baseline','FAO',runList] ;

% When doing 2x3 baseline vs. SSP maps, what subplots do subplots go into?
ssp_plot_index = [5 3 6 2] ;

% Import FAO data
%    Area harvested: ha
%    Production:     metric tons
%    Yield:          Hg/ha
if include_fao 
    if strcmp(thisVer,'20180424agmip7') ...
    || strcmp(thisVer,'20180424agmip7_asPLUMout2011-2015') ...
    || strcmp(thisVer,'v3s1_v20180426')...
    || strcmp(thisVer,'v3s1_v20180426_asPLUM2011')...
    || strcmp(thisVer,'v4s1_v20180426') ...
    || strcmp(thisVer,'v4s1_v20180426_asPLUMout2011') ...
    || strcmp(thisVer,'v4s1_v20180426_asLUH2_2010') ...
    || strcmp(thisVer,'v6s1_v20180703') ...
    || strcmp(thisVer,'v10s1_v20180801') ...
    || contains(thisVer,'harm2') ...
    || contains(thisVer,'harm3')
        fao = load(fao_file) ;
    else
        error('thisVer not recognized while loading FAO data')
    end
end

% Get plot titles
runList_titles = runList ;
if contains(thisVer, '_attr')
    tmp = strsplit(runList{1}, '-') ;
    tmp_ssp = tmp{1}(end) ;
    tmp_rcp = tmp{2} ;
    ii = find(strcmp(runList, sprintf('SSP%s-%s', tmp_ssp, tmp_rcp))) ;
    if length(ii)~=1
        error('length(ii)~=1')
    end
    runList_titles{ii} = sprintf('%s (s%slum_r%sclico2)', ...
        runList{ii}, tmp_ssp, tmp_rcp) ;
    ii = find(strcmp(runList, 'constLU')) ;
    runList_titles{ii} = sprintf('Const. LU (r%sclico2)', ...
        tmp_rcp) ;
    ii = find(strcmp(runList, 'constClimCO2')) ;
    runList_titles{ii} = sprintf('Const. Clim./CO2 (s%slum)', ...
        tmp_ssp) ;
    ii = find(strcmp(runList, 'constCO2')) ;
    runList_titles{ii} = sprintf('Const. CO2 (s%slum_r%scli)', ...
        tmp_ssp, tmp_rcp) ;
end
runList_titles = strrep(runList_titles, '_', '\_') ;


%% Import baseline
disp('Importing baseline...')

% Baseline run
ts_tmp = load([baselineDir 'timeseries.mat']) ;
bl_ts_fields = fieldnames(ts_tmp) ;
if ~isempty(ignored_crops)
    ignoreText = 'Removing ' ;
    for c = 1:length(ignored_crops)
        if c == length(ignored_crops)
            ignoreText = [ignoreText 'and '] ;
        end
        ignoreText = [ignoreText ignored_crops{c}] ;
        if c < length(ignored_crops)
            ignoreText = [ignoreText ', '] ;
        end
    end
    ignoreText = [ignoreText ' results.'] ;
    warning(ignoreText)
    clear ignoreText
    for c = 1:length(ignored_crops)
        thisI = find(endsWith(bl_ts_fields,ignored_crops{c})) ;
        if isempty(thisI)
            warning(['(Actually ' ignored_crops{c} ' not present in results.)'])
        end
        bl_ts_fields(thisI) = [] ;
    end
end

%
% Parse into original vs. before extra cropland to pasture
isb4xtracrop2past = ~cellfun(@isempty,regexpi(bl_ts_fields,'^[a-z]+0_')) ...
                  | ~cellfun(@isempty,regexpi(bl_ts_fields,'^LUarea_ts_[a-z])+0') ...
                  ) ;
bl_ts_fields0 = bl_ts_fields(isb4xtracrop2past) ;
bl_ts_fields = bl_ts_fields(~isb4xtracrop2past) ;
for f = 1:length(bl_ts_fields)
    thisField_in = bl_ts_fields{f} ;
    thisName_out = ['ts_' strrep(thisField_in,'_ts','') '_bl'] ;
    eval([thisName_out ' = ts_tmp.' thisField_in ' ;']) ;
end
if ~isempty(bl_ts_fields0)
    for f = 1:length(bl_ts_fields0)
        thisField_in = bl_ts_fields0{f} ;
        thisName_out = ['ts_' strrep(thisField_in,'_ts','') '_bl'] ;
        eval([thisName_out ' = ts_tmp.' thisField_in ' ;']) ;
    end
end

% Get CFT names
CFTnames = strrep({bl_ts_fields{getCellsWithString(bl_ts_fields,'croparea')}}','croparea_ts_','') ;
disp('CFTnames:') ; disp(CFTnames)

Ncrops = length(CFTnames) ;

% Check whether rainfed and irrigated crops are separated
combine_rf_ir = false ;
for c = 1:Ncrops
    thisCrop = CFTnames{c} ;
    thisCropi = [thisCrop 'i'] ;
    if strcmp(CFTnames,thisCropi)
        warning('Combining rainfed and irrigated outputs...')
        combine_rf_ir = true ;
        break
    end
end

% Add "irrigated" to "rainfed"
if combine_rf_ir
    error('Finish code to do this for baseline run!')
%     for c = 1:Ncrops_lpjg
%         thisCrop = CFTnames{c} ;
%         thisCropi = [thisCrop 'i'] ;
%         eval(['ts_croparea_' thisCrop '_bl = ts_croparea_' thisCrop '_bl + ts_tmp.croparea_ts_' thisCrop 'i ;']) ;
%         eval(['ts_gsirrig_' thisCrop '_bl = ts_gsirrig_' thisCrop '_bl + ts_tmp.gsirrig_ts_' thisCrop 'i ;']) ;
%         eval(['ts_yield_' thisCrop '_bl = ts_yield_' thisCrop '_bl + ts_tmp.yield_ts_' thisCrop 'i ;']) ;
%     end
end

clear ts_tmp

% Adjust FAO data
if include_fao
    if Ncrops>0
        if any(strcmp(fao.listCrops_fa2o,'Wheat'))
            fao.listCrops_fa2o{strcmp(fao.listCrops_fa2o,'Wheat')} = 'CerealsC3' ;
        end
        if any(strcmp(fao.listCrops_fa2o,'Maize'))
            fao.listCrops_fa2o{strcmp(fao.listCrops_fa2o,'Maize')} = 'CerealsC4' ;
        end
        if any(strcmp(fao.listCrops_fa2o,'Starchy roots'))
            fao.listCrops_fa2o{strcmp(fao.listCrops_fa2o,'Starchy roots')} = 'StarchyRoots' ;
        end
        if any(contains(CFTnames,'_plum'))
            fao.listCrops_fa2o = strcat(fao.listCrops_fa2o,'_plum') ;
        end
        for c = 1:Ncrops
            thisCrop = CFTnames{c} ;
            % % %         thisCrop = strrep(CFTnames{c},'_plum','') ;
            if isfield(fao.listCrops_fa2o,'a') || isfield(fao.listCrops_fa2o,'p')
                if ~(isfield(fao.listCrops_fa2o,'a') && isfield(fao.listCrops_fa2o,'p'))
                    error('One but not both present of fao.listCrops_fa2o.a and .b')
                end
            elseif length(find(strcmp(fao.listCrops_fa2o,thisCrop)))==1
                eval(['ts_cropprod_' thisCrop '_fao = single(cf_t2kg*squeeze(nansum(fao.tmp_total_fa2_Ccy(:,strcmp(fao.listCrops_fa2o,thisCrop),:),1))) ;']) ;
                eval(['ts_croparea_' thisCrop '_fao = single(cf_ha2m2*squeeze(nansum(fao.tmp_croparea_fa2_Ccy(:,strcmp(fao.listCrops_fa2o,thisCrop),:),1))) ;']) ;
            else
                warning(['Assuming zeros for FAO data for ' thisCrop])
                eval(['ts_cropprod_' thisCrop '_fao = zeros(length(fao.tmp_fao_yearList),1,''single'') ;']) ;
                eval(['ts_croparea_' thisCrop '_fao = zeros(length(fao.tmp_fao_yearList),1,''single'') ;']) ;
            end
        end
    end
end

firstdec_tmp = load([baselineDir 'first_decade.garr.mat']) ;
bl_map_fields = fieldnames(firstdec_tmp) ;
for f = 1:length(bl_map_fields)
    thisField = bl_map_fields{f} ;
    eval(['list2map_test = firstdec_tmp.' thisField '.list2map ;']) ;
    if ~exist('list2map', 'var')
        list2map = list2map_test ;
    elseif ~isequal(list2map, list2map_test)
        error('Mismatch in list2map!')
    end; clear list2map_test
    eval(['garr_' thisField ' = renameStructField(firstdec_tmp.' thisField ',''garr_xvy'',''garr_xvyB'') ;']) ;
    firstdec_tmp = rmfield(firstdec_tmp, thisField) ;
end
clear firstdec_tmp
lastdec_tmp = load([baselineDir 'last_decade.garr.mat']) ;
bl_map_fields = fieldnames(lastdec_tmp) ;
for f = 1:length(bl_map_fields)
    thisField = bl_map_fields{f} ;
    eval(['list2map_test = lastdec_tmp.' thisField '.list2map ;']) ;
    if ~exist('list2map', 'var')
        list2map = list2map_test ;
    elseif ~isequal(list2map, list2map_test)
        error('Mismatch in list2map!')
    end; clear list2map_test
    eval(['garr_' thisField ' = renameStructField(lastdec_tmp.' thisField ',''garr_xvy'',''garr_xvyB'') ;']) ;
    lastdec_tmp = rmfield(lastdec_tmp, thisField) ;
end
clear lastdec_tmp
last30yrs_tmp = load([baselineDir 'last_30yrs.garr.mat']) ;
bl_map_fields = fieldnames(last30yrs_tmp) ;
for f = 1:length(bl_map_fields)
    thisField = bl_map_fields{f} ;
    eval(['list2map_test = last30yrs_tmp.' thisField '.list2map ;']) ;
    if ~exist('list2map', 'var')
        list2map = list2map_test ;
    elseif ~isequal(list2map, list2map_test)
        error('Mismatch in list2map!')
    end; clear list2map_test
    eval(['garr_' thisField ' = renameStructField(last30yrs_tmp.' thisField ',''garr_xvs'',''garr_xvsB'') ;']) ;
    last30yrs_tmp = rmfield(last30yrs_tmp, thisField) ;
end
clear last30yrs_tmp

% Get variable names
tmp = whos('ts_*_bl') ;
vars_ts_bl = {tmp.name}' ;
clear tmp
% Parse into original vs. before extra cropland to pasture
isb4xtracrop2past = ~cellfun(@isempty,regexpi(vars_ts_bl,'^ts_[a-z]+0_')) ...
                  | ~cellfun(@isempty,regexpi(vars_ts_bl,'^ts_LUarea_[a-z]+0_')) ;
vars_ts_bl0 = vars_ts_bl(isb4xtracrop2past) ;
vars_ts_bl = vars_ts_bl(~isb4xtracrop2past) ;

disp('Done importing baseline.')


%% Future runs

if combine_rf_ir
    error('Write code to do this for future runs!')
end

disp('Importing future...')

% Future runs
Nyears_2trim = 0 ;
for r = 1:Nruns
    disp(['Importing run ' num2str(r) ' of ' num2str(Nruns) '...']) ;
    ts_tmp = load([runDirs{r} 'timeseries.mat']) ;
    yr_ts_fields = fieldnames(ts_tmp) ;
    yr_ts_fields(endsWith(yr_ts_fields,ignored_crops)) = [] ;
    is_expYields = ~cellfun(@isempty,regexp(yr_ts_fields,'^cropprodExp_')) ;
    have_expYields = any(is_expYields) ;
    firstdec_tmp = load([runDirs{r} 'first_decade.garr.mat']) ;
    lastdec_tmp = load([runDirs{r} 'last_decade.garr.mat']) ;
    last30yrs_tmp = load([runDirs{r} 'last_30yrs.garr.mat']) ;
    if r == 1
        if length(ts_tmp.LUarea_ts_bare) > Nyears_fu
            Nyears_2trim = length(ts_tmp.LUarea_ts_bare) - Nyears_fu ;
            warning(['Trimming first ' num2str(Nyears_2trim) ' years of future runs.'])
        elseif length(ts_tmp.LUarea_ts_bare) < Nyears_fu
            error('length(ts_tmp.LUarea_ts_bare) < Nyears_fu!')
        end
        y1 = Nyears_2trim + 1 ;
        
        % Set up empty arrays
        for v = 1:length(vars_ts_bl)
            thisVar_in = vars_ts_bl{v} ;
            thisVar_out = strrep(thisVar_in,'_bl','_yr') ;
            eval([thisVar_out ' = nan(Nyears_fu,Nruns,''single'') ;']) ;
            clear thisVar*
        end ; clear v
        if have_expYields
            expYield_names = yr_ts_fields(is_expYields) ;
            for c = 1:length(expYield_names)
                thisVar_in = expYield_names{c} ;
                thisVar_out = strcat(strrep(thisVar_in,'cropprodExp_ts','ts_cropprodExp'),'_yr') ;
                eval([thisVar_out ' = nan(Nyears_fu,Nruns,''single'') ;']) ;
            end
        end
        
        % Get garrs that were read in baseline
        tmp = whos('garr_*') ;
        vars_garrs_bl = {tmp.name}' ;
        clear tmp
        for v = 1:length(vars_garrs_bl)
            thisVar_in = vars_garrs_bl{v} ;
            if contains(thisVar_in,'garr_LU0')
                continue
            end
            thisVar_out = strrep(thisVar_in,'garr_','') ;
            
            % If not present in future run, remove from analysis
            if ~(isfield(firstdec_tmp,thisVar_out) || isfield(lastdec_tmp,thisVar_out) || isfield(last30yrs_tmp,thisVar_out))
                warning([thisVar_out ' not found in future run ' num2str(r) '. Removing.']) ;
                clear(thisVar_in)
                continue
            end
            
            if contains(thisVar_in,'_d1')
                eval([thisVar_in '.garr_xvyr = nan([size(firstdec_tmp.' thisVar_out '.garr_xvy) Nruns],''single'') ;']) ;
            elseif contains(thisVar_in,'_d9')
                eval([thisVar_in '.garr_xvyr = nan([size(lastdec_tmp.' thisVar_out '.garr_xvy) Nruns],''single'') ;']) ;
            elseif contains(thisVar_in,'_last30')
                eval([thisVar_in '.garr_xvsr = nan([size(last30yrs_tmp.' thisVar_out '.garr_xvs) Nruns],''single'') ;']) ;
            else
                error('How did this happen?')
            end
            clear thisVar*
        end ; clear v
        
        % Get vars_garrs_bl again, now that you've cleared variables not
        % present in future run(s)
        tmp = whos('garr_*') ;
        vars_garrs_bl = {tmp.name}' ;
        clear tmp
    end % if r==1
    
    remove = false(size(bl_ts_fields)) ;
    for f = 1:length(bl_ts_fields)
        thisField_in = bl_ts_fields{f} ;
        thisName_out = ['ts_' strrep(thisField_in,'_ts','') '_yr'] ;
        % If not present in future run, remove from analysis
        if contains(thisField_in,'crop0') || contains(thisField_in,'past0')
            continue
        elseif ~isfield(ts_tmp,thisField_in)
            warning([thisField_in ' not found in future run ' num2str(r) '. Removing.']) ;
            thisName_in = ['ts_' strrep(thisField_in,'_ts','') '_bl'] ;
            clear(thisName_in)
            remove(f) = true ;
            if exist(thisName_out,'var')
                clear(thisName_out)
            end
            continue
        end
        thisName_out = ['ts_' strrep(thisField_in,'_ts','') '_yr'] ;
        eval([thisName_out '(:,r) = ts_tmp.' thisField_in '(y1:end) ;']) ;
    end ; clear f
    
    % Update bl_ts_fields now that you've removed variables not present in
    % future run(s)
    bl_ts_fields(remove) = [] ;
    
    % Expected yields
    if have_expYields
        expYield_names = yr_ts_fields(is_expYields) ;
        for c = 1:length(expYield_names)
            thisVar_in = expYield_names{c} ;
            thisVar_out = strcat(strrep(thisVar_in,'cropprodExp_ts','ts_cropprodExp'),'_yr') ;
            eval([thisVar_out '(:,r) = ts_tmp.' thisVar_in '(y1:end) ;']) ;
        end
    end
    clear ts_tmp
    
    % Maps
    vars_garrs_bl_toRemove = false(size(vars_garrs_bl)) ;
    for v = 1:length(vars_garrs_bl)
        thisVar_in = vars_garrs_bl{v} ;
        if contains(thisVar_in,'garr_LU0')
            continue
        end
        thisVar_out = strrep(thisVar_in,'garr_','') ;
        if contains(thisVar_in,'d1')
            if ~isfield(firstdec_tmp,thisVar_out)
                warning([thisVar_in ' not found in future run ' num2str(r) '. Removing.']) ;
                vars_garrs_bl_toRemove(v) = true ;
                eval(['clear ' thisVar_in]) ;
                continue
            end
            eval(['list2map_test = firstdec_tmp.' thisVar_out '.list2map ;']) ;
            if ~exist('list2map', 'var')
                list2map = list2map_test ;
            elseif ~isequal(list2map, list2map_test)
                error('Mismatch in list2map!')
            end; clear list2map_test
            eval(['isequal_varNames = isequal(' thisVar_in '.varNames, firstdec_tmp.' thisVar_out '.varNames) ; ']) ;
            if ~isequal_varNames
                eval(['isequal_varNames_afterSort = isequal(sort(' thisVar_in '.varNames), sort(firstdec_tmp.' thisVar_out '.varNames)) ; ']) ;
                if ~isequal_varNames_afterSort
                    error(['~isequal_varNames (' thisVar_in '). Not fixable.'])
                end
                warning(['~isequal_varNames (' thisVar_in '). Fixing.'])
                eval(['[~,~,IB] = intersect(' thisVar_in '.varNames,firstdec_tmp.' thisVar_out '.varNames,''stable'') ;']) ;
                eval([thisVar_in '.garr_xvyr(:,:,:,r) = firstdec_tmp.' thisVar_out '.garr_xvy(:,IB,:) ;']) ;
            else
                eval([thisVar_in '.garr_xvyr(:,:,:,r) = firstdec_tmp.' thisVar_out '.garr_xvy ;']) ;
            end
        elseif contains(thisVar_in,'d9')
            if ~isfield(lastdec_tmp,thisVar_out)
                warning([thisVar_in ' not found in future run ' num2str(r) '. Removing.']) ;
                vars_garrs_bl_toRemove(v) = true ;
                eval(['clear ' thisVar_in]) ;
                continue
            end
            eval(['list2map_test = lastdec_tmp.' thisVar_out '.list2map ;']) ;
            if ~exist('list2map', 'var')
                list2map = list2map_test ;
            elseif ~isequal(list2map, list2map_test)
                error('Mismatch in list2map!')
            end; clear list2map_test
            eval(['isequal_varNames = isequal(' thisVar_in '.varNames, lastdec_tmp.' thisVar_out '.varNames) ; ']) ;
            if ~isequal_varNames
                eval(['isequal_varNames_afterSort = isequal(sort(' thisVar_in '.varNames), sort(lastdec_tmp.' thisVar_out '.varNames)) ; ']) ;
                if ~isequal_varNames_afterSort
                    error(['~isequal_varNames (' thisVar_in '). Not fixable.'])
                end
                warning(['~isequal_varNames (' thisVar_in '). Fixing.'])
                eval(['[~,~,IB] = intersect(' thisVar_in '.varNames,lastdec_tmp.' thisVar_out '.varNames,''stable'') ;']) ;
                eval([thisVar_in '.garr_xvyr(:,:,:,r) = lastdec_tmp.' thisVar_out '.garr_xvy(:,IB,:) ;']) ;
            else
                eval([thisVar_in '.garr_xvyr(:,:,:,r) = lastdec_tmp.' thisVar_out '.garr_xvy ;']) ;
            end
        elseif contains(thisVar_in,'last30')
            if ~isfield(last30yrs_tmp,thisVar_out)
                warning([thisVar_in ' not found in future run ' num2str(r) '. Removing.']) ;
                vars_garrs_bl_toRemove(v) = true ;
                eval(['clear ' thisVar_in]) ;
                continue
            end
            eval(['list2map_test = last30yrs_tmp.' thisVar_out '.list2map ;']) ;
            if ~exist('list2map', 'var')
                list2map = list2map_test ;
            elseif ~isequal(list2map, list2map_test)
                error('Mismatch in list2map!')
            end; clear list2map_test
            eval(['isequal_varNames = isequal(' thisVar_in '.varNames, last30yrs_tmp.' thisVar_out '.varNames) ; ']) ;
            if ~isequal_varNames
                eval(['isequal_varNames_afterSort = isequal(sort(' thisVar_in '.varNames), sort(last30yrs_tmp.' thisVar_out '.varNames)) ; ']) ;
                if ~isequal_varNames_afterSort
                    error(['~isequal_varNames (' thisVar_in '). Not fixable.'])
                end
                warning(['~isequal_varNames (' thisVar_in '). Fixing.'])
                eval(['[~,~,IB] = intersect(' thisVar_in '.varNames,last30yrs_tmp.' thisVar_out '.varNames,''stable'') ;']) ;
                eval([thisVar_in '.garr_xvsr(:,:,:,r) = last30yrs_tmp.' thisVar_out '.garr_xvs(:,IB,:) ;']) ;
            else
                eval([thisVar_in '.garr_xvsr(:,:,:,r) = last30yrs_tmp.' thisVar_out '.garr_xvs ;']) ;
            end
            eval(['isequal_statList = isequal(' thisVar_in '.statList, last30yrs_tmp.' thisVar_out '.statList) ; ']) ;
            if ~isequal_statList
                eval(['isequal_varNames_afterSort = isequal(sort(' thisVar_in '.statList), sort(last30yrs_tmp.' thisVar_out '.statList)) ; ']) ;
                if ~isequal_varNames_afterSort
                    error(['~isequal_statList (' thisVar_in '). Not fixable.'])
                end
                warning(['~isequal_statList (' thisVar_in '). Fixing.'])
                eval(['[~,~,IB] = intersect(' thisVar_in '.statList,last30yrs_tmp.' thisVar_out '.statList,''stable'') ;']) ;
                eval([thisVar_in '.garr_xvsr(:,:,:,r) = last30yrs_tmp.' thisVar_out '.garr_xvs(:,:,IB) ;']) ;
            else
                eval([thisVar_in '.garr_xvsr(:,:,:,r) = last30yrs_tmp.' thisVar_out '.garr_xvs ;']) ;
            end
        else
            error('How did this happen?')
        end
    end
    
    vars_garrs_bl(vars_garrs_bl_toRemove) = [] ;
    clear vars_garrs_bl_toRemove
    
    if isfield(firstdec_tmp,'expyield_d1')
        if r == 1
            garr_yieldExp_d1 = renameStructField(firstdec_tmp.expyield_d1,'garr_xvy','garr_xvyr') ;
        else
            garr_yieldExp_d1.garr_xvyr(:,:,:,r) = firstdec_tmp.expyield_d1.maps_YXvy ;
        end
    end
    if isfield(lastdec_tmp,'expyield_d9')
        if r == 1
            garr_yieldExp_d9 = renameStructField(lastdec_tmp.expyield_d9,'garr_xvy','garr_xvyr') ;
        else
            garr_yieldExp_d9.garr_xvy(:,:,:,r) = lastdec_tmp.expyield_d9.garr_xvy ;
        end
    end
    if isfield(last30yrs_tmp,'expyield_last30')
        if r == 1
            garr_yieldExp_last30 = renameStructField(last30yrs_tmp.expyield_last30,'garr_xvs','garr_xvsr') ;
        else
            garr_yieldExp_last30.garr_xvsr(:,:,:,:,r) = last30yrs_tmp.expyield_last30.garr_xvs ;
        end
    end

    clear firstdec_tmp
    clear lastdec_tmp
    clear last30yrs_tmp

end ; clear r
%%
disp('Processing...')
nanmask = lpjgu_vector2map(zeros(size(list2map)), map_size, list2map) ;
nanmask(isnan(nanmask)) = 1 ;
nanmask = logical(nanmask) ;

% Import land area (km2)
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = gcel_area_YXqd(:,1:2:1440) + gcel_area_YXqd(:,2:2:1440) ;
gcel_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
gcel_area_x = gcel_area_YX(list2map) ;
gcel_area_YX(nanmask) = NaN ;
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
land_area_x = land_area_YX(list2map) ;
land_area_unmasked_YX = land_area_YX ;
land_area_unmasked_weights_YX = land_area_unmasked_YX ./ nansum(nansum(land_area_unmasked_YX)) ;
land_area_unmasked_weights_x = land_area_unmasked_weights_YX(list2map) ;
land_area_YX(nanmask) = NaN ;
%%% Convert to m2
land_area_YX = land_area_YX*1e6 ;
land_area_x = land_area_x*1e6 ;
gcel_area_YX = gcel_area_YX*1e6 ;
gcel_area_x = gcel_area_x*1e6 ;
clear tmp land_frac_YXqd land_area_YXqd

% area_xBH: "Baseline" as last year of Historical run
% area_xBFr: "Baseline" as first year of Future runs
% diff_xrH: Difference from End-Historical to End-Future
% diff_xrF: Difference from Begin-Future to End-Future
tmp_list_4 = {'ntrl','bare','crop','past'} ;
tmp_list_full = {'NATURAL','BARREN','CROPLAND','PASTURE'} ;
for i = 1:length(tmp_list_4)
    this_4 = tmp_list_4{i} ;
    this_full = tmp_list_full{i} ;
    eval([this_4 '_area_xBH = gcel_area_x .* garr_LU_d9.garr_xvyB(:,strcmp(garr_LU_d9.varNames,this_full),end) ;']) ;
    eval([this_4 '_area_xBFr = gcel_area_x .* squeeze(garr_LU_d1.garr_xvyr(:,strcmp(garr_LU_d1.varNames,this_full),1,:)) ;']) ;
    eval([this_4 '_area_xr = repmat(gcel_area_x,[1 Nruns]) .* squeeze(garr_LU_d9.garr_xvyr(:,strcmp(garr_LU_d9.varNames,this_full),end,:)) ;']) ;
    eval([this_4 '_diff_xrH = ' this_4 '_area_xr - repmat(' this_4 '_area_xBH,[1 Nruns]) ;']) ;
    eval([this_4 '_diff_xrF = ' this_4 '_area_xr - ' this_4 '_area_xBFr ;']) ;
end

if exist('garr_LU0_d9','var')
    warning('Not tested after maps->garr conversion. Problems here might be a result of that.')
    crop0_area_xBH = gcel_area_x .* garr_LU0_d9.garr_xvyB(:,strcmp(garr_LU_d9.varNames,'CROPLAND'),end) ;
    crop0_diff_xrH = crop_area_xr - repmat(crop0_area_xBH,[1 Nruns]) ;
    past0_area_xBH = gcel_area_x .* garr_LU0_d9.garr_xvyB(:,strcmp(garr_LU_d9.varNames,'PASTURE'),end) ;
    past0_diff_xrH = past_area_xr - repmat(past0_area_xBH,[1 Nruns]) ;
end

agri_area_xBH = crop_area_xBH + past_area_xBH ;
agri_area_xBFr = crop_area_xBFr + past_area_xBFr ;
agri_area_xr = crop_area_xr + past_area_xr ;
agri_diff_xrH = crop_diff_xrH + past_diff_xrH ;
agri_diff_xrF = crop_diff_xrF + past_diff_xrF ;

% Remove ignored crops
if ~isempty(ignored_crops)
    garr_cropfracs_d1.garr_xvyB(:,contains(garr_cropfracs_d1.varNames,ignored_crops),:) = [] ;
    garr_cropfracs_d1.garr_xvyr(:,contains(garr_cropfracs_d1.varNames,ignored_crops),:,:) = [] ;
    garr_cropfracs_d1.varNames(contains(garr_cropfracs_d1.varNames,ignored_crops)) = [] ;
    garr_cropfracs_d9.garr_xvyB(:,contains(garr_cropfracs_d9.varNames,ignored_crops),:) = [] ;
    garr_cropfracs_d9.garr_xvyr(:,contains(garr_cropfracs_d9.varNames,ignored_crops),:,:) = [] ;
    garr_cropfracs_d9.varNames(contains(garr_cropfracs_d9.varNames,ignored_crops)) = [] ;
    garr_yield_d1.garr_xvyB(:,contains(garr_yield_d1.varNames,ignored_crops),:) = [] ;
    garr_yield_d1.garr_xvyr(:,contains(garr_yield_d1.varNames,ignored_crops),:,:) = [] ;
    garr_yield_d1.varNames(contains(garr_yield_d1.varNames,ignored_crops)) = [] ;
    garr_yield_d9.garr_xvyB(:,contains(garr_yield_d9.varNames,ignored_crops),:) = [] ;
    garr_yield_d9.garr_xvyr(:,contains(garr_yield_d9.varNames,ignored_crops),:,:) = [] ;
    garr_yield_d9.varNames(contains(garr_yield_d9.varNames,ignored_crops)) = [] ;
    if exist('garr_yieldExp_d1','var')
        % (Baseline has no expected yield) garr_yieldExp_d1.garr_xvyB(:,contains(garr_yieldExp_d1.varNames,ignored_crops),:) = [] ;
        garr_yieldExp_d1.garr_xvyr(:,contains(garr_yieldExp_d1.varNames,ignored_crops),:,:) = [] ;
        garr_yieldExp_d1.varNames(contains(garr_yieldExp_d1.varNames,ignored_crops)) = [] ;
    end
    if exist('garr_yieldExp_d9','var')
        % (Baseline has no expected yield) garr_yieldExp_d9.garr_xvyB(:,contains(garr_yieldExp_d9.varNames,ignored_crops),:) = [] ;
        garr_yieldExp_d9.garr_xvyr(:,contains(garr_yieldExp_d9.varNames,ignored_crops),:,:) = [] ;
        garr_yieldExp_d9.varNames(contains(garr_yieldExp_d9.varNames,ignored_crops)) = [] ;
    end
    if exist('garr_Nfert_d1','var')
        garr_Nfert_d1.garr_xvyB(:,contains(garr_Nfert_d1.varNames,ignored_crops),:) = [] ;
        garr_Nfert_d1.garr_xvyr(:,contains(garr_Nfert_d1.varNames,ignored_crops),:,:) = [] ;
        garr_Nfert_d1.varNames(contains(garr_Nfert_d1.varNames,ignored_crops)) = [] ;
    end
    if exist('garr_Nfert_d9','var')
        garr_Nfert_d9.garr_xvyB(:,contains(garr_Nfert_d9.varNames,ignored_crops),:) = [] ;
        garr_Nfert_d9.garr_xvyr(:,contains(garr_Nfert_d9.varNames,ignored_crops),:,:) = [] ;
        garr_Nfert_d9.varNames(contains(garr_Nfert_d9.varNames,ignored_crops)) = [] ;
    end
end

if ~isequal(garr_cropfracs_d1.varNames,garr_cropfracs_d9.varNames)
    error('~isequal(garr_cropfracs_d1.varNames,garr_cropfracs_d9.varNames)')
else
    CFTnames_garr = garr_cropfracs_d1.varNames ;
end
garr_cropareas_d1 = garr_cropfracs_d1 ;
garr_cropareas_d1.garr_xvyB = garr_cropfracs_d1.garr_xvyB .* repmat(gcel_area_x,[1 Ncrops size(garr_cropfracs_d1.garr_xvyB,3)]) .* repmat(garr_LU_d1.garr_xvyB(:,strcmp(garr_LU_d1.varNames,'CROPLAND'),:),[1 Ncrops 1]) ;
garr_cropareas_d1.garr_xvyr = garr_cropfracs_d1.garr_xvyr .* repmat(gcel_area_x,[1 Ncrops size(garr_cropfracs_d1.garr_xvyB,3) Nruns]) .* repmat(garr_LU_d1.garr_xvyr(:,strcmp(garr_LU_d1.varNames,'CROPLAND'),:,:),[1 Ncrops 1 1]) ;
garr_cropareas_d9 = garr_cropfracs_d9 ;
garr_cropareas_d9.garr_xvyB = garr_cropfracs_d9.garr_xvyB .* repmat(gcel_area_x,[1 Ncrops size(garr_cropfracs_d9.garr_xvyB,3)]) .* repmat(garr_LU_d9.garr_xvyB(:,strcmp(garr_LU_d9.varNames,'CROPLAND'),:),[1 Ncrops 1]) ;
garr_cropareas_d9.garr_xvyr = garr_cropfracs_d9.garr_xvyr .* repmat(gcel_area_x,[1 Ncrops size(garr_cropfracs_d9.garr_xvyB,3) Nruns]) .* repmat(garr_LU_d9.garr_xvyr(:,strcmp(garr_LU_d9.varNames,'CROPLAND'),:,:),[1 Ncrops 1 1]) ;

garr_cropareas_xvBH =          garr_cropareas_d9.garr_xvyB(:,:,end) ;
garr_cropareas_xvBFr = squeeze(garr_cropareas_d1.garr_xvyr(:,:,1,  :)) ;
garr_cropareas_xvr   = squeeze(garr_cropareas_d9.garr_xvyr(:,:,end,:)) ;
garr_cropareasDiffs_xvrH = garr_cropareas_xvr - repmat(garr_cropareas_xvBH,[1 1 Nruns]) ;
garr_cropareasDiffs_xvrF = garr_cropareas_xvr - garr_cropareas_xvBFr ;

disp('Done importing future.')


%% Perform secondary calculations

disp('Performing secondary calculations...')

% Adjust for technological improvement, if doing so
if do_adjYieldTech
    
    % Garrs
    garr_yield_d1 = adj_yield_tech(garr_yield_d1, yearList_baseline(end-9:end), yearList_future(1:10)) ;
    garr_yield_d9 = adj_yield_tech(garr_yield_d9, yearList_baseline(end-9:end), yearList_future(end-9:end)) ;
    
    % Time series
    tmp = whos('ts_croppro*') ;
    tmp_name = {tmp.name}' ;
    for i = 1:length(tmp_name)
        thisName_cropprod = tmp_name{i} ;
        
        % Do not apply adjustment to FAO data or PLUM-expected numbers
        if strcmp(thisName_cropprod(end-2:end),'fao') || contains(thisName_cropprod,'cropproExp') || contains(thisName_cropprod,'cropprodExp')
            continue
        end
        
        % Get correct yearList depending on baseline vs. future data
        if strcmp(thisName_cropprod(end-1:end),'bl')
            yearList_tmp = yearList_baseline ;
        elseif strcmp(thisName_cropprod(end-1:end),'yr')
            yearList_tmp = yearList_future ;
        else
            error('Problem parsing variable name (%s) for tech. improvement adjustment.', thisName_cropprod)
        end
        
        % Do the adjustment
        thisCmd = sprintf('%s = adj_yield_tech(%s, yearList_tmp) ;', ...
            thisName_cropprod, thisName_cropprod) ;
        eval(thisCmd) ;
        
        clear yearList_tmp thisCmd
    end ; clear i
    clear tmp tmp_name
    
end

% Combine evaporation and transpiration
if exist('ts_aevap_bl','var') && exist('ts_aaet_bl','var')
    ts_aevapaaet_bl = ts_aevap_bl + ts_aaet_bl ;
end
if exist('ts_aevap_yr','var') && exist('ts_aaet_yr','var')
    ts_aevapaaet_yr = ts_aevap_yr + ts_aaet_yr ;
end

% Combine gaseous and liquid N fluxes
if exist('ts_nflux_flux_bl','var') && exist('ts_nflux_leach_bl','var')
    ts_nloss_bl = ts_nflux_flux_bl + ts_nflux_leach_bl ;
end
if exist('ts_nflux_flux_yr','var') && exist('ts_nflux_leach_yr','var')
    ts_nloss_yr = ts_nflux_flux_yr + ts_nflux_leach_yr ;
end

% Calculate actual YIELD (production per area, kg/m2)
tmp = whos('ts_croppro*') ;
tmp_name = {tmp.name}' ;
for i = 1:length(tmp_name)
    thisName_cropprod = tmp_name{i} ;
    if contains(thisName_cropprod,'cropproExp')
        thisname_yield = strrep(thisName_cropprod,'cropproExp','yielExp') ;
        thisname_croparea = strrep(thisName_cropprod,'cropproExp','croparea') ;
    elseif contains(thisName_cropprod,'cropprodExp')
        thisname_yield = strrep(thisName_cropprod,'cropprodExp','yieldExp') ;
        thisname_croparea = strrep(thisName_cropprod,'cropprodExp','croparea') ;
    else
        thisname_yield = strrep(thisName_cropprod,'cropprod','yield') ;
        thisname_croparea = strrep(thisName_cropprod,'cropprod','croparea') ;
    end
    eval([thisname_yield ' = ' thisName_cropprod ' ./ ' thisname_croparea ' ;']) ;
    % Try with pre-extracrop-to-pasture variable
    thisname_croparea = strrep(thisName_cropprod,'cropprod','croparea0') ;
    thisname_yield = strrep(thisName_cropprod,'cropprod','yield0') ;
    if exist(thisname_croparea,'var') && exist(thisname_yield,'var')
        eval([thisname_yield ' = ' thisName_cropprod ' ./ ' thisname_croparea ' ;']) ;
    end
    clear thisN* thisn
end ; clear i
clear tmp*

% Calculate calorie production from LPJ-GUESS
tmp = whos('ts_cropprod_*') ;
tmp_name = {tmp.name}' ;
ts_kcal_bl = zeros(size(Nyears_bl,1)) ;
ts_kcal_fao = zeros(size(Nyears_bl,1)) ;
ts_kcal_yr = zeros(size(Nyears_bl,Nruns)) ;
for i = 1:length(tmp_name)
    thisCrop = strrep(strrep(strrep(strrep(tmp_name{i},'ts_cropprod_',''),'_bl',''),'_yr',''),'_fao','') ;
    if strcmp(thisCrop,'Miscanthus')
        continue
    end
    thisSuffix = strrep(tmp_name{i},['ts_cropprod_' thisCrop '_'],'') ;
    kcal_per_g = get_kcalDensity2(thisCrop) ;
    kcal_per_kg = 1e3 * kcal_per_g ;
    eval(['ts_kcal_' thisSuffix ' = ts_kcal_' thisSuffix ' + kcal_per_kg * eval(tmp_name{i}) ;']) ;
end ; clear i

% Calculate calorie production expected by PLUM
tmp = whos('ts_yieldExp_*') ;
if ~isempty(tmp)
    tmp_name = {tmp.name}' ;
    ts_kcalExp_yr = zeros(size(Nyears_bl,Nruns)) ;
    ts_kcalExp2_yr = zeros(size(Nyears_bl,Nruns)) ;
    for i = 1:length(tmp_name)
        thisCrop = strrep(strrep(tmp_name{i},'ts_yieldExp_',''),'_yr','') ;
        if strcmp(thisCrop,'Miscanthus')
            continue
        end
        thisSuffix = strrep(tmp_name{i},['ts_yieldExp_' thisCrop '_'],'') ;
        if ~strcmp(thisSuffix,'yr')
            error('~strcmp(thisSuffix,''yr'') ???')
        end
        kcal_per_g = get_kcalDensity2(thisCrop) ;
        kcal_per_kg = 1e3 * kcal_per_g ;
        eval(['ts_kcalExp2_' thisSuffix ' = ts_kcalExp2_' thisSuffix ' + kcal_per_kg * (' tmp_name{i} '.* ts_croparea_' thisCrop '_yr) ;']) ;
    end ; clear i
end ; clear tmp

% Peak monthly runoff maps
if exist('maps_mon_runoff_d1','var')
    maps_pk_runoff_d1 = maps_mon_runoff_d1 ;
    maps_pk_runoff_d1.maps_xvyB = max(maps_pk_runoff_d1.maps_xvyB,[],2) ;
    maps_pk_runoff_d1.maps_xvyr = max(maps_pk_runoff_d1.maps_xvyr,[],2) ;
    maps_pk_runoff_d1.varNames = {'Max'} ;
    maps_pk_runoff_d9 = maps_mon_runoff_d9 ;
    maps_pk_runoff_d9.maps_xvyB = max(maps_pk_runoff_d9.maps_xvyB,[],2) ;
    maps_pk_runoff_d9.maps_xvyr = max(maps_pk_runoff_d9.maps_xvyr,[],2) ;
    maps_pk_runoff_d9.varNames = {'Max'} ;
    
    % Change in 5th and 95th percentile of peak monthly runoff
    % After Asadieh & Krakauer (2017)
    maps_pk_runoff_d1.maps_p05_xB = prctile(maps_pk_runoff_d1.maps_xvyB,5,3) ;
    maps_pk_runoff_d1.maps_p95_xB = prctile(maps_pk_runoff_d1.maps_xvyB,95,3) ;
    maps_pk_runoff_d9.maps_p05_xB = prctile(maps_pk_runoff_d9.maps_xvyB,5,3) ;
    maps_pk_runoff_d9.maps_p95_xB = prctile(maps_pk_runoff_d9.maps_xvyB,95,3) ;
    maps_pk_runoff_d1.maps_p05_xr = squeeze(prctile(maps_pk_runoff_d1.maps_xvyr,5,3)) ;
    maps_pk_runoff_d1.maps_p95_xr = squeeze(prctile(maps_pk_runoff_d1.maps_xvyr,95,3)) ;
    maps_pk_runoff_d9.maps_p05_xr = squeeze(prctile(maps_pk_runoff_d9.maps_xvyr,5,3)) ;
    maps_pk_runoff_d9.maps_p95_xr = squeeze(prctile(maps_pk_runoff_d9.maps_xvyr,95,3)) ;
    below_thresh_x = mean(maps_awater_d1.maps_xvyB(:,strcmp(maps_awater_d1.varNames,'Runoff'),:),3)/365 < 0.01 ;
    maps_pk_runoff_d1.maps_p05_xB(below_thresh_x) = NaN ;
    maps_pk_runoff_d1.maps_p95_xB(below_thresh_x) = NaN ;
    maps_pk_runoff_d9.maps_p05_xB(below_thresh_x) = NaN ;
    maps_pk_runoff_d9.maps_p95_xB(below_thresh_x) = NaN ;
    maps_pk_runoff_d1.maps_p05_xr(repmat(below_thresh_x, [1 Nruns])) = NaN ;
    maps_pk_runoff_d1.maps_p95_xr(repmat(below_thresh_x, [1 Nruns])) = NaN ;
    maps_pk_runoff_d9.maps_p05_xr(repmat(below_thresh_x, [1 Nruns])) = NaN ;
    maps_pk_runoff_d9.maps_p95_xr(repmat(below_thresh_x, [1 Nruns])) = NaN ;

end

% Total crop production
tmp = whos('ts_cropprod_*') ;
tmp_name = {tmp.name}' ;
ts_cropprod_bl = zeros(size(Nyears_bl,1)) ;
ts_cropprod_fao = zeros(size(Nyears_bl,1)) ;
ts_cropprod_yr = zeros(size(Nyears_bl,Nruns)) ;
for i = 1:length(tmp_name)
    thisCrop = strrep(strrep(strrep(strrep(tmp_name{i},'ts_cropprod_',''),'_bl',''),'_yr',''),'_fao','') ;
    if strcmp(thisCrop,'Miscanthus')
        continue
    end
    thisSuffix = strrep(tmp_name{i},['ts_cropprod_' thisCrop '_'],'') ;
    eval(sprintf('ts_cropprod_%s = ts_cropprod_%s + %s ;', ...
        thisSuffix, thisSuffix, tmp_name{i})) ;
end ; clear i

% Non-bare land area
bare_area_xmean = mean(cat(2,bare_area_xBH,bare_area_xr),2) ;
bare_frac_xmean = bare_area_xmean ./ gcel_area_x ;
vegd_frac_xmean = 1 - bare_frac_xmean ;
vegd_area_xmean = vegd_frac_xmean .* gcel_area_x ;

disp('Done performing secondary calculations.')


%% Import biodiversity hotspots

disp('Importing biodiversity hotspots...')
hotspot_YX = flipud(imread(hotspot_tif)) ;
hotspot_YX(nanmask) = NaN ;
hotspot_YX = 1==hotspot_YX ;
hotspot_area_YX = hotspot_YX.*gcel_area_YX ;
hotspot_x = hotspot_YX(list2map) ;
hotspot_area_x = hotspot_area_YX(list2map) ;

disp('Importing Congolian swamp and lowland forests...')
ecoid_YX = flipud(imread(ecoid_file)) ;
cslf_YX = false(size(ecoid_YX)) ;
cslf_IDs = [ ...
    30129 ; ... Western Congolian swamp forests
    30126 ; ... Northwest Congolian lowland forests
    30124 ; ... Northeast Congolian lowland forests
    30110 ; ... Eastern Congolian swamp forests
    30104 ; ... Central Congolian lowland forests
    ] ;
for j = 1:length(cslf_IDs)
    cslf_YX(ecoid_YX==cslf_IDs(j)) = true ;
end
clear ecoid_YX
cslf_x = cslf_YX(list2map) ;
hotspotCSLF_area_YX = (hotspot_YX | cslf_YX).*gcel_area_YX ;
hotspotCSLF_area_x = hotspotCSLF_area_YX(list2map) ;


%% Import food production units and basins

disp('Importing FPUs...')

fpu_YX = flipud(dlmread(fpu_file,' ',6,0)) ;
fpu_YX(fpu_YX==-9999) = NaN ;
fpu_x = fpu_YX(list2map) ;

% Combine some FPUs to create the Amazon and Nile basins
basins_YX = fpu_YX ;
basin_groups = { ...
    [70 72 74 146 149 150 203] ; % Nile
    [9 10 12 13] ; % Amazon
    } ;
for b = 1:length(basin_groups)
    thisGroup = basin_groups{b} ;
    for ii = 2:length(thisGroup)
        thisFPU = thisGroup(ii) ;
        basins_YX(basins_YX==thisFPU) = thisGroup(1) ;
    end ; clear ii
    clear thisGroup thisFPU
end ; clear b basin_groups

% Mask
basins_YX(nanmask) = NaN ;
basins_x = basins_YX(list2map) ;

% Get basin numbers
basin_list = unique(basins_YX(~isnan(basins_YX))) ;
Nbasins = length(basin_list) ;


%% Import population

disp('Importing population...')

pop = readtable(pop_file) ;
pop(:,strcmp(pop.Properties.VariableNames,'Model')) = [] ;
pop = pop(contains(pop.Scenario, 'v9_130325'),:) ;
pop.Scenario = strrep(pop.Scenario,'_v9_130325','') ;
[~,~,sortCols] = intersect({'Scenario','Country','Year'},pop.Properties.VariableNames,'stable') ;
pop = sortrows(pop,sortCols) ; % Needed to produce _ycr dimensioned array

yearList_pop_orig = unique(pop.Year) ;
Nyears_pop_orig = length(yearList_pop_orig) ;
yearList_pop = min(yearList_pop_orig):max(yearList_pop_orig) ;
Nyears_pop = length(yearList_pop) ;
countryList_pop = unique(pop.Country) ;
Ncountries_pop = length(countryList_pop) ;
runList_pop = unique(pop.Scenario) ;
Nruns_pop = length(runList_pop) ;

% Get population array, interpolating and converting millions of people to people
pop_ycr_preInterp = reshape(pop.population*1e6,[Nyears_pop_orig Ncountries_pop Nruns_pop]) ;
pop_ycr_interpd = nan(Nyears_pop, Ncountries_pop, Nruns_pop) ;
for c = 1:Ncountries_pop
    for r = 1:Nruns_pop
        pop_ycr_interpd(:,c,r) = interp1(yearList_pop_orig, pop_ycr_preInterp(:,c,r), yearList_pop) ;
    end
end

clear pop_ycr_preInterp

% Rearrange to correspond to runDirs_plum
rearranged = zeros(size(runDirs_plum)) ;
for r = 1:Nruns
    thisSSP = regexp(runDirs_plum{r},'SSP.','match') ;
    i = find(strcmp(runList_pop,thisSSP)) ;
    if isempty(i)
        error('thisSSP (%s) not found in runList_pop', thisSSP)
    elseif length(i) > 1
        error('thisSSP (%s) found more than once in runList_pop', thisSSP)
    end
    rearranged(r) = i ;
    clear i
end; clear r
pop_ycr = pop_ycr_interpd(:,:,rearranged) ;
clear pop_ycr_interpd rearranged

% Get global array
pop_yr = squeeze(sum(pop_ycr,2)) ;


%% Import global demand (Mt to kg, kcal)

disp('Importing global demand...')

warning('off','MATLAB:table:ModifiedAndSavedVarnames')
for r = 1:Nruns
    thisDir = runDirs_plum{r} ;
    thisFile = sprintf('%s/demand.txt', thisDir) ;
    thisTable = readtable(thisFile) ;
    [~,~,sortCols] = intersect({'Commodity','Year'},thisTable.Properties.VariableNames, 'stable') ;
    thisTable = sortrows(thisTable,sortCols) ; % Needed to produce _yvr dimensioned array
    if r==1
        commods = unique(thisTable.Commodity) ; % Sorts alphabetically
        Ncommods = length(commods) ;
        yearList_PLUMout = unique(thisTable.Year) ;
        Nyears_PLUMout = length(yearList_PLUMout) ;
        ts_commodDemand_yvr = nan(Nyears_PLUMout, Ncommods, Nruns) ;
    end
    ts_commodDemand_yvr(:,:,r) = reshape(thisTable.Amount_Mt_,[Nyears_PLUMout Ncommods]) ;
    clear thisDir thisFile thisTable sortCols
end
warning('on','MATLAB:table:ModifiedAndSavedVarnames')

% Add total crop and livestock demand
commods_livestock = {'ruminants','monogastrics'} ;
[~, i_crop] = setdiff(commods, commods_livestock) ;
ts_commodDemand_yvr(:,end+1,:) = sum(ts_commodDemand_yvr(:,i_crop,:),2) ;
commods{end+1} = 'crops' ;
[~, i_livestock] = intersect(commods, commods_livestock) ;
ts_commodDemand_yvr(:,end+1,:) = sum(ts_commodDemand_yvr(:,i_livestock,:),2) ;
commods{end+1} = 'livestock' ;
Ncommods = length(commods) ;

% Convert Mt to kg
ts_commodDemand_yvr = ts_commodDemand_yvr * 1e6*1e3 ;

% Get calories (crops only)
ts_commodDemand_kcal_yvr = nan(Nyears_PLUMout, Ncommods, Nruns) ;
for ii = 1:length(i_crop)
    thisCrop = commods{i_crop(ii)} ;
    kcal_per_g = get_kcalDensity2(thisCrop) ;
    kcal_per_kg = 1e3 * kcal_per_g ;
    ts_commodDemand_kcal_yvr(:,ii,:) = kcal_per_kg*ts_commodDemand_yvr(:,ii,:) ;
end
ts_commodDemand_kcal_yvr(:,strcmp(commods,'crops'),:) = nansum(ts_commodDemand_kcal_yvr,2) ;

% Get per-capita demand (kg/person)
if ~isequal(shiftdim(yearList_pop), shiftdim(yearList_PLUMout))
    error('This code assumes population and demand have identical yearLists!')
end
pop_yvr = repmat(permute(pop_yr,[1 3 2]), [1 size(ts_commodDemand_yvr,2) 1]) ;
ts_commodDemandPC_yvr = ts_commodDemand_yvr ./ pop_yvr ;
clear pop_yvr


%% Import LUH1 


disp('Done.')

