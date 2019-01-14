%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLUM2LPJG figures: Multiple runs: %%%
%%% PLUM-style native %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% thisVer = '20180424agmip7' ;
% thisVer = '20180424agmip7_asPLUMout2011-2015' ;
thisVer = 'v3s1_v20180426' ;

do_save = true ;
rebase = false ;
pngres = 150 ;

       
%% Setup

include_fao = true ;
if strcmp(thisVer,'20180424agmip7') || strcmp(thisVer,'20180424agmip7_asPLUMout2011-2015')
    warning('not including fao data!')
    include_fao = false ;
end

addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work/')
addpath(genpath('~/Documents/Dropbox/Dissertation/MATLAB work'))

trimFirstFuture = 0 ;
if strcmp(thisVer,'20180424agmip7')
    runList = {'SSP1/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runDirs = {
        'PLUM2LPJG_SSP1_RCP45_v3s1/output-2018-04-23-145614' ;
        'PLUM2LPJG_SSP3_RCP60_v3s1/output-2018-04-23-145614' ;
        'PLUM2LPJG_SSP4_RCP60_v3s1/output-2018-04-23-145614' ;
        'PLUM2LPJG_SSP5_RCP85_v3s1/output-2018-04-23-145614' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_PLUM6xtra/matlab_merge_20180424050342' ;
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20180424agmip7_asPLUMout2011-2015')
    warning('SOMETHING WENT WRONG WITH SSP3; IGNORING')
%     runList = {'SSP1/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runList = {'SSP1/RCP4.5','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runDirs = {
        'PLUM2LPJG_SSP1_RCP45_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-134335' ;
%         'PLUM2LPJG_SSP3_RCP60_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-131116' ;
        'PLUM2LPJG_SSP4_RCP60_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-134415' ;
        'PLUM2LPJG_SSP5_RCP85_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-132213' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_PLUM6xtra/matlab_merge_20180424050342' ;
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'v3s1_v20180426')
    runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    runDirs = {
        'PLUM2LPJG_SSP1_RCP45_v3s1_v20180426/output-2018-05-01-164708' ;
        'PLUM2LPJG_SSP3_RCP60_v3s1_v20180426/output-2018-04-30-125214' ;
        'PLUM2LPJG_SSP4_RCP60_v3s1_v20180426/output-2018-04-30-125218' ;
        'PLUM2LPJG_SSP5_RCP85_v3s1_v20180426/output-2018-05-01-024615' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180502104840' ;
    yearList_baseline = 1850:2010 ;
else
    error(['thisVer (' thisVer ') not recognized!'])
end

% Deal with years
if max(yearList_baseline) >= min(yearList_future)
    error('Baseline and future yearLists overlap!')
elseif max(yearList_baseline) ~= min(yearList_future)-1
    error('Baseline and future yearLists are offset!')
end
Nyears_bl = length(yearList_baseline) ;
Nyears_fu = length(yearList_future) ;

% Output directory
outDir_maps = addslashifneeded(['/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/'...
                           'figures_' thisVer '/maps']) ;
outDir_ts = ['/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/'...
                           'figures_' thisVer '/TS'] ;
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

% Conversion factors
cf_kg2Mt = 1e-3*1e-6 ;
cf_t2kg = 1e3 ;
cf_ha2Mkm2 = 1e-8 ;

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
    || strcmp(thisVer,'v3s1_v20180426')
        fao = load('/Users/Shared/PLUM/crop_calib_data/fao/FAOdata_1961-2010_calibVer16.mat') ;
    else
        error('thisVer not recognized: loading FAO data')
    end
end


%% Import baseline

disp('Importing baseline...')

% Baseline run
ts_tmp = load([baselineDir 'timeseries.mat']) ;
bl_ts_fields = fieldnames(ts_tmp) ;
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
                eval(['ts_cropprod_' thisCrop '_fao = cf_t2kg*squeeze(nansum(fao.tmp_total_fa2_Ccy(:,strcmp(fao.listCrops_fa2o,thisCrop),:),1)) ;']) ;
                eval(['ts_croparea_' thisCrop '_fao = cf_ha2Mkm2*squeeze(nansum(fao.tmp_croparea_fa2_Ccy(:,strcmp(fao.listCrops_fa2o,thisCrop),:),1)) ;']) ;
            else
                warning(['Assuming zeros for FAO data for ' thisCrop])
                eval(['ts_cropprod_' thisCrop '_fao = zeros(length(fao.tmp_fao_yearList),1) ;']) ;
                eval(['ts_croparea_' thisCrop '_fao = zeros(length(fao.tmp_fao_yearList),1) ;']) ;
            end
        end
    end
end

firstdec_tmp = load([baselineDir 'first_decade.mat']) ;
bl_map_fields = fieldnames(firstdec_tmp) ;
for f = 1:length(bl_map_fields)
    thisField = bl_map_fields{f} ;
    eval(['maps_' thisField ' = renameStructField(firstdec_tmp.' thisField ',''maps_YXvy'',''maps_YXvyB'') ;']) ;
end
lastdec_tmp = load([baselineDir 'last_decade.mat']) ;
bl_map_fields = fieldnames(lastdec_tmp) ;
for f = 1:length(bl_map_fields)
    thisField = bl_map_fields{f} ;
    eval(['maps_' thisField ' = renameStructField(lastdec_tmp.' thisField ',''maps_YXvy'',''maps_YXvyB'') ;']) ;
end

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
    have_expYields = false ;
    if exist([runDirs{r} 'timeseries_PLUMexp.mat'],'file')
        have_expYields = true ;
        tsExpYield_tmp = load([runDirs{r} 'timeseries_PLUMexp.mat']) ;
    end
    firstdec_tmp = load([runDirs{r} 'first_decade.mat']) ;
    lastdec_tmp = load([runDirs{r} 'last_decade.mat']) ;
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
            eval([thisVar_out ' = nan(Nyears_fu,Nruns) ;']) ;
            clear thisVar*
        end ; clear v
        
        % Get maps that were read in baseline
        tmp = whos('maps_*') ;
        vars_maps_bl = {tmp.name}' ;
        clear tmp
        for v = 1:length(vars_maps_bl)
            thisVar_in = vars_maps_bl{v} ;
            if contains(thisVar_in,'maps_LU0')
                continue
            end
            thisVar_out = strrep(thisVar_in,'maps_','') ;
            if contains(thisVar_in,'_d1')
                eval([thisVar_in '.maps_YXvyr = nan([size(firstdec_tmp.' thisVar_out '.maps_YXvy) Nruns]) ;']) ;
            elseif contains(thisVar_in,'_d9')
                eval([thisVar_in '.maps_YXvyr = nan([size(lastdec_tmp.' thisVar_out '.maps_YXvy) Nruns]) ;']) ;
            else
                error('How did this happen?')
            end
            clear thisVar*
        end ; clear v
    end
    
    for f = 1:length(bl_ts_fields)
        thisField_in = bl_ts_fields{f} ;
        if contains(thisField_in,'crop0') || contains(thisField_in,'past0')
            continue
        end
        thisName_out = ['ts_' strrep(thisField_in,'_ts','') '_yr'] ;
        eval([thisName_out '(:,r) = ts_tmp.' thisField_in '(y1:end) ;']) ;
    end ; clear f
    
    % Expected yields
    if have_expYields
        tsExpYield_fields = fieldnames(tsExpYield_tmp) ;
        for f = 1:length(tsExpYield_fields)
            thisField_in = tsExpYield_fields{f} ;
            thisName_out = strrep(['ts_' strrep(thisField_in,'_ts','') '_yr'],'cropprod','cropproExp') ;
            if r==1
                eval([thisName_out ' = nan(Nyears_fu,Nruns) ;']) ;
            end
            eval([thisName_out '(:,r) = tsExpYield_tmp.' thisField_in '(y1:end) ;']) ;
        end ; clear f
    end
    clear ts_tmp
    
    % Maps
    for v = 1:length(vars_maps_bl)
        thisVar_in = vars_maps_bl{v} ;
        if contains(thisVar_in,'maps_LU0')
            continue
        end
        thisVar_out = strrep(thisVar_in,'maps_','') ;
        if contains(thisVar_in,'d1')
            eval([thisVar_in '.maps_YXvyr(:,:,:,:,r) = firstdec_tmp.' thisVar_out '.maps_YXvy ;']) ;
        elseif contains(thisVar_in,'d9')
            eval([thisVar_in '.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.' thisVar_out '.maps_YXvy ;']) ;
        else
            error('How did this happen?')
        end
    end
end ; clear r

disp('Processing...')
nanmask = isnan(maps_LU_d1.maps_YXvyB(:,:,1,1)) ;

% Import land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
land_area_YX(nanmask) = NaN ;
land_area_YX_m2 = land_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd

% area_YXBH: "Baseline" as last year of Historical run
% area_YXBFr: "Baseline" as first year of Future runs
% diff_YXrH: Difference from End-Historical to End-Future
% diff_YXrF: Difference from Begin-Future to End-Future
tmp_list_4 = {'ntrl','bare','crop','past'} ;
tmp_list_full = {'NATURAL','BARREN','CROPLAND','PASTURE'} ;
for i = 1:length(tmp_list_4)
    this_4 = tmp_list_4{i} ;
    this_full = tmp_list_full{i} ;
    eval([this_4 '_area_YXBH = land_area_YX .* maps_LU_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,this_full),end) ;']) ;
    eval([this_4 '_area_YXBFr = land_area_YX .* squeeze(maps_LU_d1.maps_YXvyr(:,:,strcmp(maps_LU_d1.varNames,this_full),1,:)) ;']) ;
    eval([this_4 '_area_YXr = repmat(land_area_YX,[1 1 Nruns]) .* squeeze(maps_LU_d9.maps_YXvyr(:,:,strcmp(maps_LU_d9.varNames,this_full),end,:)) ;']) ;
    eval([this_4 '_diff_YXrH = ' this_4 '_area_YXr - repmat(' this_4 '_area_YXBH,[1 1 Nruns]) ;']) ;
    eval([this_4 '_diff_YXrF = ' this_4 '_area_YXr - ' this_4 '_area_YXBFr ;']) ;
end
if false
% ntrl_area_YXBH = land_area_YX .* maps_LU_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'NATURAL'),end) ;
% ntrl_area_YXBFr = land_area_YX .* squeeze(maps_LU_d1.maps_YXvyr(:,:,strcmp(maps_LU_d1.varNames,'NATURAL'),1,:)) ;
% ntrl_area_YXr = repmat(land_area_YX,[1 1 Nruns]) .* squeeze(maps_LU_d9.maps_YXvyr(:,:,strcmp(maps_LU_d9.varNames,'NATURAL'),end,:)) ;
% ntrl_diff_YXrH = ntrl_area_YXr - repmat(ntrl_area_YXBH,[1 1 Nruns]) ;
% ntrl_diff_YXrF = ntrl_area_YXr - ntrl_area_YXBFr ;
% 
% bare_area_YXBH = land_area_YX .* maps_LU_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'BARREN'),end) ;
% bare_area_YXBFr = land_area_YX .* squeeze(maps_LU_d1.maps_YXvyr(:,:,strcmp(maps_LU_d1.varNames,'BARREN'),1,:)) ;
% bare_area_YXr = repmat(land_area_YX,[1 1 Nruns]) .* squeeze(maps_LU_d9.maps_YXvyr(:,:,strcmp(maps_LU_d9.varNames,'BARREN'),end,:)) ;
% bare_diff_YXrH = bare_area_YXr - repmat(bare_area_YXBH,[1 1 Nruns]) ;
% bare_diff_YXrF = bare_area_YXr - bare_area_YXBFr ;
% 
% crop_area_YXBH = land_area_YX .* maps_LU_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'CROPLAND'),end) ;
% crop_area_YXBFr = land_area_YX .* squeeze(maps_LU_d1.maps_YXvyr(:,:,strcmp(maps_LU_d1.varNames,'CROPLAND'),1,:)) ;
% crop_area_YXr = repmat(land_area_YX,[1 1 Nruns]) .* squeeze(maps_LU_d9.maps_YXvyr(:,:,strcmp(maps_LU_d9.varNames,'CROPLAND'),end,:)) ;
% crop_diff_YXrH = crop_area_YXr - repmat(crop_area_YXBH,[1 1 Nruns]) ;
% crop_diff_YXrF = crop_area_YXr - crop_area_YXBFr ;
% past_area_YXBH = land_area_YX .* maps_LU_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'PASTURE'),end) ;
% past_area_YXBFr = land_area_YX .* squeeze(maps_LU_d1.maps_YXvyr(:,:,strcmp(maps_LU_d1.varNames,'PASTURE'),1,:)) ;
% past_area_YXr = repmat(land_area_YX,[1 1 Nruns]) .* squeeze(maps_LU_d9.maps_YXvyr(:,:,strcmp(maps_LU_d9.varNames,'PASTURE'),end,:)) ;
% past_diff_YXrH = past_area_YXr - repmat(past_area_YXBH,[1 1 Nruns]) ;
% past_diff_YXrF = past_area_YXr - past_area_YXBFr ;
end

if exist('maps_LU0_d9','var')
    crop0_area_YXBH = land_area_YX .* maps_LU0_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'CROPLAND'),end) ;
    crop0_diff_YXrH = crop_area_YXr - repmat(crop0_area_YXBH,[1 1 Nruns]) ;
end
if exist('maps_LU0_d9','var')
    past0_area_YXBH = land_area_YX .* maps_LU0_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'PASTURE'),end) ;
    past0_diff_YXrH = past_area_YXr - repmat(past0_area_YXBH,[1 1 Nruns]) ;
end

agri_area_YXBH = crop_area_YXBH + past_area_YXBH ;
agri_area_YXBFr = crop_area_YXBFr + past_area_YXBFr ;
agri_area_YXr = crop_area_YXr + past_area_YXr ;
agri_diff_YXrH = crop_diff_YXrH + past_diff_YXrH ;
agri_diff_YXrF = crop_diff_YXrF + past_diff_YXrF ;

% Check that order of maps is correct
for v = 1:length(maps_LU_d1.varNames)
    thisLU = maps_LU_d1.varNames{v} ;
    if strcmp(thisLU,'NATURAL')
        thisLU_short = 'ntrl' ;
    elseif strcmp(thisLU,'CROPLAND')
        thisLU_short = 'crop' ;
    elseif strcmp(thisLU,'PASTURE')
        thisLU_short = 'past' ;
    elseif strcmp(thisLU,'BARREN')
        thisLU_short = 'bare' ;
    end
    eval(['fromTS = ts_LUarea_' thisLU_short '_bl(end) ;']) ;
    eval(['fromMap = 1e-6*nansum(nansum(' thisLU_short '_area_YXBH)) ;']) ;
    if abs(fromTS-fromMap) > 1e-5
        error('LU variables are in wrong order!')
    end
end
for v = 1:length(maps_LU_d9.varNames)
    thisLU = maps_LU_d9.varNames{v} ;
    if strcmp(thisLU,'NATURAL')
        thisLU_short = 'ntrl' ;
    elseif strcmp(thisLU,'CROPLAND')
        thisLU_short = 'crop' ;
    elseif strcmp(thisLU,'PASTURE')
        thisLU_short = 'past' ;
    elseif strcmp(thisLU,'BARREN')
        thisLU_short = 'bare' ;
    end
    eval(['fromTS = ts_LUarea_' thisLU_short '_bl(end) ;']) ;
    eval(['fromMap = 1e-6*nansum(nansum(' thisLU_short '_area_YXBH)) ;']) ;
    if abs(fromTS-fromMap) > 1e-5
        error('LU variables are in wrong order!')
    end
end

if ~isequal(maps_cropfracs_d1.varNames,maps_cropfracs_d9.varNames)
    error('~isequal(maps_cropfracs_d1.varNames,maps_cropfracs_d9.varNames)')
else
    CFTnames_maps = maps_cropfracs_d1.varNames ;
end
maps_cropareas_d1 = maps_cropfracs_d1 ;
maps_cropareas_d1.maps_YXvyB = maps_cropfracs_d1.maps_YXvyB .* repmat(land_area_YX,[1 1 Ncrops size(maps_cropfracs_d1.maps_YXvyB,4)]) .* repmat(maps_LU_d1.maps_YXvyB(:,:,strcmp(maps_LU_d1.varNames,'CROPLAND'),:),[1 1 Ncrops 1]) ;
maps_cropareas_d1.maps_YXvyr = maps_cropfracs_d1.maps_YXvyr .* repmat(land_area_YX,[1 1 Ncrops size(maps_cropfracs_d1.maps_YXvyB,4) Nruns]) .* repmat(maps_LU_d1.maps_YXvyr(:,:,strcmp(maps_LU_d1.varNames,'CROPLAND'),:,:),[1 1 Ncrops 1 1]) ;
maps_cropareas_d9 = maps_cropfracs_d9 ;
maps_cropareas_d9.maps_YXvyB = maps_cropfracs_d9.maps_YXvyB .* repmat(land_area_YX,[1 1 Ncrops size(maps_cropfracs_d9.maps_YXvyB,4)]) .* repmat(maps_LU_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'CROPLAND'),:),[1 1 Ncrops 1]) ;
maps_cropareas_d9.maps_YXvyr = maps_cropfracs_d9.maps_YXvyr .* repmat(land_area_YX,[1 1 Ncrops size(maps_cropfracs_d9.maps_YXvyB,4) Nruns]) .* repmat(maps_LU_d9.maps_YXvyr(:,:,strcmp(maps_LU_d9.varNames,'CROPLAND'),:,:),[1 1 Ncrops 1 1]) ;

maps_cropareas_YXvBH = maps_cropareas_d9.maps_YXvyB(:,:,:,end) ;
maps_cropareas_YXvBFr = squeeze(maps_cropareas_d1.maps_YXvyr(:,:,:,1,:)) ;
maps_cropareas_YXvr = squeeze(maps_cropareas_d9.maps_YXvyr(:,:,:,end,:)) ;
maps_cropareasDiffs_YXvrH = maps_cropareas_YXvr - repmat(maps_cropareas_YXvBH,[1 1 1 Nruns]) ;
maps_cropareasDiffs_YXvrF = maps_cropareas_YXvr - maps_cropareas_YXvBFr ;

% Calculate actual YIELD (production per area)
kg_2_t = 1e-3 ;
Mkm2_2_ha = 1e6*1e2 ;
tmp = whos('ts_croppro*') ;
tmp_name = {tmp.name}' ;
for i = 1:length(tmp_name)
    thisName_cropprod = tmp_name{i} ;
    is_expYield = contains(thisName_cropprod,'cropproExp') ;
    if is_expYield
        thisname_yield = strrep(thisName_cropprod,'cropproExp','yielExp') ;
        thisname_croparea = strrep(thisName_cropprod,'cropproExp','croparea') ;
    else
        thisname_yield = strrep(thisName_cropprod,'cropprod','yield') ;
        thisname_croparea = strrep(thisName_cropprod,'cropprod','croparea') ;
    end
    eval([thisname_yield ' = (' thisName_cropprod ' * kg_2_t) ./ (' thisname_croparea '*Mkm2_2_ha) ;']) ;
    % Try with pre-extracrop-to-pasture variable
    thisname_croparea = strrep(thisName_cropprod,'cropprod','croparea0') ;
    thisname_yield = strrep(thisName_cropprod,'cropprod','yield0') ;
    if exist(thisname_croparea,'var') && exist(thisname_yield,'var')
        eval([thisname_yield ' = (' thisName_cropprod ' * kg_2_t) ./ (' thisname_croparea '*Mkm2_2_ha) ;']) ;
    end
    clear this*
end ; clear i
clear tmp*

disp('Done importing future.')


%% Map changes in LU area: End-Historical to End-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = figurePos ;
nx = 3 ;
ny = 2 ;
colorBarLoc = 'SouthOutside' ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;

% Natural area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_ntrl_area_YXBH = ntrl_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(ntrl_area_YXBH(:,:))) ;
    else
        this_ntrl_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_ntrl_area_YXBH, ntrl_diff_YXrH(:,:,r), total_bl, ...
        thisY1, thisYN, '"Natural"', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
    close
end

% Cropland area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_crop_area_YXBH = crop_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(crop_area_YXBH(:,:))) ;
    else
        this_crop_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_crop_area_YXBH, crop_diff_YXrH(:,:,r), total_bl, ...
        thisY1, thisYN, 'Cropland', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
%
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
    close
end

% "Pasture" area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_past_area_YXBH = past_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(past_area_YXBH(:,:))) ;
    else
        this_past_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_past_area_YXBH, past_diff_YXrH(:,:,r), total_bl, ...
        thisY1, thisYN, '"Pasture"', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
    close
end

% Agricultural area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_agri_area_YXBH = agri_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(agri_area_YXBH(:,:))) ;
    else
        this_agri_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_agri_area_YXBH, agri_diff_YXrH(:,:,r), total_bl, ...
        thisY1, thisYN, 'Agricultural', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_agri.png'],['-r' num2str(pngres)])
    close
end


%% Map changes in LU area: End-Historical to End-Future (crop0 and past0)

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = figurePos ;
nx = 3 ;
ny = 2 ;
colorBarLoc = 'SouthOutside' ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;

% Cropland area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_crop0_area_YXBH = crop0_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(crop0_area_YXBH(:,:))) ;
    else
        this_crop0_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_crop0_area_YXBH, crop0_diff_YXrH(:,:,r), total_bl, ...
        thisY1, thisYN, 'Cropland0', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop0.png'],['-r' num2str(pngres)])
    close
end

% "Pasture" area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_past0_area_YXBH = past0_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(past0_area_YXBH(:,:))) ;
    else
        this_past0_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_past0_area_YXBH, past0_diff_YXrH(:,:,r), total_bl, ...
        thisY1, thisYN, '"Pasture0"', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past0.png'],['-r' num2str(pngres)])
    close
end


%% Map changes in LU area: End-Historical to Begin-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = figurePos ;
nx = 3 ;
ny = 2 ;
colorBarLoc = 'SouthOutside' ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(1) ;

% Natural area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_ntrl_area_YXBH = ntrl_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(ntrl_area_YXBH(:,:))) ;
    else
        this_ntrl_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_ntrl_area_YXBH, ntrl_area_YXBFr(:,:,r) - ntrl_area_YXBH, total_bl, ...
        thisY1, thisYN, '"Natural"', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
    close
end

% Cropland area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_crop_area_YXBH = crop_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(crop_area_YXBH(:,:))) ;
    else
        this_crop_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_crop_area_YXBH, crop_area_YXBFr(:,:,r) - crop_area_YXBH, total_bl, ...
        thisY1, thisYN, 'Cropland', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
    close
end

% "Pasture" area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_past_area_YXBH = past_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(past_area_YXBH(:,:))) ;
    else
        this_past_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_past_area_YXBH, past_area_YXBFr(:,:,r) - past_area_YXBH, total_bl, ...
        thisY1, thisYN, '"Pasture"', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
    close
end

% Agricultural area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_agri_area_YXBH = agri_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(agri_area_YXBH(:,:))) ;
    else
        this_agri_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_agri_area_YXBH, agri_area_YXBFr(:,:,r) - agri_area_YXBH, total_bl, ...
        thisY1, thisYN, 'Agricultural', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_agri.png'],['-r' num2str(pngres)])
    close
end


%% Map changes in LU area: End-Historical to Begin-Future (crop0 and past0)

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = figurePos ;
nx = 3 ;
ny = 2 ;
colorBarLoc = 'SouthOutside' ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(1) ;

% Cropland area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_crop_area_YXBH = crop0_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(crop0_area_YXBH(:,:))) ;
    else
        this_crop_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_crop_area_YXBH, crop_area_YXBFr(:,:,r) - crop0_area_YXBH, total_bl, ...
        thisY1, thisYN, 'Cropland0', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop0.png'],['-r' num2str(pngres)])
    close
end

% "Pasture" area
figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if r==1
        this_past_area_YXBH = past0_area_YXBH ;
        total_bl = 1e-6*nansum(nansum(past0_area_YXBH(:,:))) ;
    else
        this_past_area_YXBH = [] ;
    end
    i2 = ssp_plot_index(r) ;
    [h, caxis_max] = make_LUdiff_fig(...
        this_past_area_YXBH, past_area_YXBFr(:,:,r) - past0_area_YXBH, total_bl, ...
        thisY1, thisYN, '"Pasture0"', '', runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past0.png'],['-r' num2str(pngres)])
    close
end


%% Map changes in LU area: Begin-Future to End-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.05 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = [1 33 935 772] ;
nx = 2 ;
ny = 4 ;
colorBarLoc = 'EastOutside' ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_future(1) ;
thisYN = yearList_future(end) ;

% Natural area
figure('Color','w','Position',thisPos) ;
for r = 1:Nruns
    i1 = (r-1)*2+1  ;
    i2 = (r-1)*2+2  ;
    [~, ~] = make_LUdiff_fig(...
        ntrl_area_YXBFr(:,:,r), ntrl_diff_YXrF(:,:,r), [], ...
        thisY1, thisYN, '"Natural"', runList{r}, runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, i1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
    close
end

% Cropland area
figure('Color','w','Position',thisPos) ;
for r = 1:Nruns
    i1 = (r-1)*2+1  ;
    i2 = (r-1)*2+2  ;
    [~, ~] = make_LUdiff_fig(...
        crop_area_YXBFr(:,:,r), crop_diff_YXrF(:,:,r), [], ...
        thisY1, thisYN, 'Cropland', runList{r}, runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, i1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
    close
end

% "Pasture" area
figure('Color','w','Position',thisPos) ;
for r = 1:Nruns
    i1 = (r-1)*2+1  ;
    i2 = (r-1)*2+2  ;
    [~, ~] = make_LUdiff_fig(...
        past_area_YXBFr(:,:,r), past_diff_YXrF(:,:,r), [], ...
        thisY1, thisYN, '"Pasture"', runList{r}, runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, i1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
    close
end

% Agricultural area
figure('Color','w','Position',thisPos) ;
for r = 1:Nruns
    i1 = (r-1)*2+1  ;
    i2 = (r-1)*2+2  ;
    [~, ~] = make_LUdiff_fig(...
        agri_area_YXBFr(:,:,r), agri_diff_YXrF(:,:,r), [], ...
        thisY1, thisYN, 'Agricultural', runList{r}, runList{r}, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, i1, i2, colorBarLoc) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
end
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_agri.png'],['-r' num2str(pngres)])
    close
end


%% Map changes in each crop area: End-Historical to End-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = figurePos ;
nx = 3 ;
ny = 2 ;
colorBarLoc = 'SouthOutside' ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;

for c = 1:Ncrops
    figure('Color','w','Position',thisPos) ;
    gcas = {} ;
    colorlim = 0 ;
    for r = 1:Nruns
        if r==1
            thiscrop_area_YXBH = squeeze(maps_cropareas_YXvBH(:,:,c)) ;
            total_bl = 1e-6*nansum(nansum(thiscrop_area_YXBH)) ;
        else
            thiscrop_area_YXBH = [] ;
        end
        i2 = ssp_plot_index(r) ;
        [h, caxis_max] = make_LUdiff_fig(...
            thiscrop_area_YXBH, maps_cropareasDiffs_YXvrH(:,:,c,r), total_bl, ...
            thisY1, thisYN, CFTnames_maps{c}, '', runList{r}, ...
            spacing, fontSize, textX, textY_1, textY_2, ...
            nx, ny, 1, i2, colorBarLoc) ;
        caxis([-max(abs(caxis)) max(abs(caxis))])
        gcas = [gcas h] ;
        colorlim = max(colorlim,caxis_max) ;
    end
    for r = 1:Nruns
        caxis(gcas(r),[-colorlim colorlim])
    end
    if do_save
        export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_maps{c} '.png'],['-r' num2str(pngres)])
        close
    end
end


%% Map changes in each crop area: End-Historical to Begin-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = figurePos ;
nx = 3 ;
ny = 2 ;
colorBarLoc = 'SouthOutside' ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(1) ;

for c = 1:Ncrops
    figure('Color','w','Position',thisPos) ;
    gcas = {} ;
    colorlim = 0 ;
    for r = 1:Nruns
        if r==1
            thiscrop_area_YXBH = squeeze(maps_cropareas_YXvBH(:,:,c)) ;
            total_bl = 1e-6*nansum(nansum(thiscrop_area_YXBH)) ;
        else
            thiscrop_area_YXBH = [] ;
        end
        i2 = ssp_plot_index(r) ;
        [h, caxis_max] = make_LUdiff_fig(...
            thiscrop_area_YXBH, maps_cropareas_YXvBFr(:,:,c,r)-squeeze(maps_cropareas_YXvBH(:,:,c)), total_bl, ...
            thisY1, thisYN, CFTnames_maps{c}, '', runList{r}, ...
            spacing, fontSize, textX, textY_1, textY_2, ...
            nx, ny, 1, i2, colorBarLoc) ;
        caxis([-max(abs(caxis)) max(abs(caxis))])
        gcas = [gcas h] ;
        colorlim = max(colorlim,caxis_max) ;
    end
    for r = 1:Nruns
        caxis(gcas(r),[-colorlim colorlim])
    end
    if do_save
        export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_maps{c} '.png'],['-r' num2str(pngres)])
        close
    end
end


%% Map changes in each crop area: Begin-Future to End-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.05 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = [1 33 935 772] ;
nx = 2 ;
ny = 4 ;
colorBarLoc = 'EastOutside' ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_future(1) ;
thisYN = yearList_future(end) ;

% Natural area
for c = 1:Ncrops
    figure('Color','w','Position',thisPos) ;
    for r = 1:Nruns
        i1 = (r-1)*2+1  ;
        i2 = (r-1)*2+2  ;
        [~, ~] = make_LUdiff_fig(...
            maps_cropareas_YXvBFr(:,:,c,r), maps_cropareasDiffs_YXvrF(:,:,c,r), [], ...
            thisY1, thisYN, CFTnames_maps{c}, runList{r}, runList{r}, ...
            spacing, fontSize, textX, textY_1, textY_2, ...
            nx, ny, i1, i2, colorBarLoc) ;
        caxis([-max(abs(caxis)) max(abs(caxis))])
    end
    if do_save
        export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_maps{c} '.png'],['-r' num2str(pngres)])
        close
    end
end


%% Map changes in BD hotspot area: CI

% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edgecolor = 0.6*ones(3,1) ;
latlim = [-60,80];
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
cbarOrient = 'SouthOutside' ;
lineWidth = 0.25 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(version('-release'),'2014b')
    
    % Biodiversity hotspots
    hotspot_shp = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/hotspots_clipByGridlist.shp' ;
    hotspot_YX = dlmread('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/hotspots_raster.txt',...
        ' ', 6, 0) ;
    hotspot_YX(hotspot_YX==-9999) = NaN ;
    hotspot_YX(:,721) = [] ;
    hotspot_YX = flipud(hotspot_YX) ;
    hotspot_YX(nanmask) = NaN ;
    hotspot_YX(~nanmask & isnan(hotspot_YX)) = 0 ;
    hotspot_area_YX = hotspot_YX.*land_area_YX ;
    hotspot_area_YXB = hotspot_area_YX .* maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'NATURAL'),end) ;
    % hotspot_area_YXB = hotspot_area_YX ;
    hotspot_area_YXr = repmat(hotspot_area_YX,[1 1 Nruns]) .* squeeze(maps_LU.maps_YXvyr(:,:,strcmp(maps_LU.varNames,'NATURAL'),end,:)) ;
    
    hotspot_diff_YXr = hotspot_area_YXr - repmat(hotspot_area_YXB,[1 1 Nruns]) ;
    
    map_hotspot_diffs(...
        hotspot_area_YXB, hotspot_diff_YXr, hotspot_YX, hotspot_shp, ...
        spacing, latlim, edgecolor, cbarOrient, fontSize, ...
        textX, textY_1, textY_2, ssp_plot_index, lineWidth, ...
        yearList_baseline, yearList_future, runList)
    
    if do_save
        export_fig([outDir_maps 'areaDiff_BDhotspots_CI.png'],['-r' num2str(pngres)])
        close
    end
    
else
    warning('Skipping hotspots')
end


%% Map changes in BD hotspot area: glob200

% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edgecolor = 0.6*ones(3,1) ;
latlim = [-60,80];
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
cbarOrient = 'SouthOutside' ;
lineWidth = 0.25 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(version('-release'),'2014b')

    % Biodiversity hotspots
    hotspot_YX = imread('/Users/sam/Geodata/global200ecoregions/g200_terr_raster0.5deg.tif') ;
    hotspot_shp = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work/g200_terr_from0.5raster.shp' ;
    hotspot_YX = flipud(hotspot_YX) ;
    hotspot_area_YX = hotspot_YX.*land_area_YX ;
    hotspot_area_YXB = hotspot_area_YX .* maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'NATURAL'),end) ;
    % hotspot_area_YXB = hotspot_area_YX ;
    hotspot_area_YXr = repmat(hotspot_area_YX,[1 1 Nruns]) .* squeeze(maps_LU.maps_YXvyr(:,:,strcmp(maps_LU.varNames,'NATURAL'),end,:)) ;
    
    hotspot_diff_YXr = hotspot_area_YXr - repmat(hotspot_area_YXB,[1 1 Nruns]) ;
    
    map_hotspot_diffs(...
        hotspot_area_YXB, hotspot_diff_YXr, hotspot_YX, hotspot_shp, ...
        spacing, latlim, edgecolor, cbarOrient, fontSize, ...
        textX, textY_1, textY_2, ssp_plot_index, lineWidth, ...
        yearList_baseline, yearList_future, runList)
    
    if do_save
        export_fig([outDir_maps 'areaDiff_BDhotspots_glob200.png'],['-r' num2str(pngres)])
        close
    end

else
    warning('Skipping hotspots')
end


%% Map changes from 2010 to 2011: LU


%% Plot timeseries: Land uses

rebase = false ;

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
blYears = yearList_baseline ;
% blYears = 1960:yearList_baseline(end) ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_LUarea_crop_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_crop_bl, ts_LUarea_crop_yr, ignYrs, yearList_future) ;
[tmp_ts_LUarea_past_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_past_bl, ts_LUarea_past_yr, ignYrs, yearList_future) ;
[tmp_ts_LUarea_ntrl_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_ntrl_bl, ts_LUarea_ntrl_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
tmp = ts_LUarea_past_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
hold on
tmp = ts_LUarea_crop_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,movmean(tmp,Nsmth),'--k','LineWidth',lineWidth)
plot(yearList_future,tmp_ts_LUarea_past_yr,'-','LineWidth',lineWidth)
set(gca,'ColorOrderIndex',1) ;
plot(yearList_future,tmp_ts_LUarea_crop_yr,'--','LineWidth',lineWidth)
if exist('ts_LUarea_past0_bl','var')
    tmp = ts_LUarea_past0_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
    plot(blYears,movmean(tmp,Nsmth),'-','LineWidth',lineWidth,'Color',0.75*ones(1,3))
end
if exist('ts_LUarea_crop0_bl','var')
    tmp = ts_LUarea_crop0_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
    plot(blYears,movmean(tmp,Nsmth),'--','LineWidth',lineWidth,'Color',0.75*ones(1,3))
end
hold off
% legend([strcat(stdLegend,', pasture') strcat(stdLegend,', cropland')], ...
%        'Location','NorthEastOutside') ;
legend('Pasture','Cropland', ...
       'Location','SouthEast') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('Million km^2')
title(['Agricultural area' title_suffix])

subplot_tight(1,2,2,spacing)
tmp = ts_LUarea_ntrl_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_LUarea_ntrl_yr,'LineWidth',lineWidth)
hold off
legend(stdLegend, ...
       'Location','SouthWest') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('Million km^2')
title(['Natural area' title_suffix])

if do_save
    export_fig([outDir_ts 'landUse' file_suffix '.pdf'])
    close
end


%% Plot timeseries: N fertilizer on each crop

% Options %%%%%%%%%
thisVar = 'nflux_fert' ;
units = 'Mt N' ;
title_prefix = 'N applied' ;
conv_fact = 1e-6 ;   % tons to Mt
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
thisLegend = stdLegend ;

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, {}, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;
if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end

% Again, but with fertilizer according to crop0
thisVar = 'nflux0_fert' ;
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, {}, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;
if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries: N fertilizer PER HECTARE on each crop

% Options %%%%%%%%%
thisVar = 'nflux_fert' ;
units = 'tons ha^{-1}' ;
title_prefix = 'N applied per hectare' ;
conv_fact = 1 ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
thisLegend = stdLegend ;

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cell_bl_nflux_fert = cell_bl ;
cell_yr_nflux_fert = cell_yr ;
clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, 'croparea', 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, 'croparea', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cell_bl_croparea = cell_bl ;
cell_yr_croparea = cell_yr ;
clear cell_bl cell_yr
% Get tons/ha
for i = 1:length(cell_bl_nflux_fert)
    cell_bl{i} = (cell_bl_nflux_fert{i}) ./ (1e6*1e2*cell_bl_croparea{i}) ;
    cell_yr{i} = (cell_yr_nflux_fert{i}) ./ (1e6*1e2*cell_yr_croparea{i}) ;
end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, {}, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;
if do_save
    export_fig([outDir_ts 'nfertPerHa' file_suffix '.pdf'])
    close
end

%% Plot timeseries: Production of each crop

% Options %%%%%%%%%
thisVar = 'cropprod' ;
units = 'Mt DM' ;
title_prefix = 'Production' ;
conv_fact = cf_kg2Mt ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
thisLegend = stdLegend_plusFAO ;

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries: Production of each crop, EXPECTED

% Options %%%%%%%%%
thisVar = 'cropproExp' ;
expected = true ;
units = 'Mt DM' ;
title_prefix = 'Production (exp.)' ;
conv_fact = cf_kg2Mt ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
thisLegend = stdLegend_plusFAO ;
    
clear cell_bl cell_yr
if strcmp(thisVar,'cropproExp')
    cmds = get_cell_forPlot(whos, 'cropprod', 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
else
    cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if strcmp(thisVar,'cropproExp')
    cmds = get_cell_forPlot(whos, 'cropprod', 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
else
    cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries: Area of each crop

% Options %%%%%%%%%
thisVar = 'croparea' ;
units = 'Million km^2' ;
title_prefix = 'Area' ;
conv_fact = 1 ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

if include_fao
    thisLegend = stdLegend_plusFAO ;
else
    thisLegend = stdLegend ;
end

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if include_fao
    cmds = get_cell_forPlot(whos, thisVar, 'fao', CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
else
    cell_fao = {} ;
    fao.tmp_fao_yearList = [] ;
end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    CFTnames, units, title_prefix, thisLegend, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries: Yield of each crop

% Options %%%%%%%%%
thisVar = 'yield' ;
units = 't ha^{-1}' ;
title_prefix = 'Yield' ;
conv_fact = 1 ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
if include_fao
    thisLegend = stdLegend_plusFAO ;
else
    thisLegend = stdLegend ;
end

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if include_fao
    cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
else
    cell_fao = {} ;
    fao.tmp_fao_yearList = [] ;
end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries: Yield of each crop, plum EXPECTED

% Options %%%%%%%%%
thisVar = 'yielExp' ;
units = 't ha^{-1}' ;
title_prefix = 'Yield (exp.)' ;
whichCFTs = 'plum' ;
conv_fact = 1 ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

if strcmp(whichCFTs,'lpjg')
    theseCFTnames = CFTnames ;
    thisLegend = stdLegend ;
elseif strcmp(whichCFTs,'plum')
    theseCFTnames = CFTnames ;
    thisLegend = stdLegend_plusFAO ;
else
    error(['whichCFTs (' whichCFTs ') not recognized!'])
end

clear cell_bl cell_yr
if strcmp(thisVar,'yielExp')
    cmds = get_cell_forPlot(whos, 'yield', 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
else
    cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if strcmp(whichCFTs,'plum')
    if strcmp(thisVar,'yielExp')
        cmds = get_cell_forPlot(whos, 'yield', 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    else
        cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    end
else
    cell_fao = {} ;
end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;

if do_save
    export_fig([outDir_ts thisVar '_' whichCFTs file_suffix '.pdf'])
    close
end


%% Plot timeseries: total area (FOR TROUBLESHOOTING)

% figure ;
% plot(yearList_baseline,movmean(ts_LUarea_crop_bl+ts_LUarea_past_bl+ts_LUarea_ntrl_bl+ts_LUarea_bare_bl,Nsmth),'-k','LineWidth',lineWidth)
% hold on
% plot(yearList_future,ts_LUarea_crop_yr+ts_LUarea_past_yr+ts_LUarea_ntrl_yr+ts_LUarea_bare_yr,'LineWidth',lineWidth)
% hold off
% legend(stdLegend, ...
%        'Location','NorthEastOutside') ;
% set(gca,'FontSize',fontSize)
% xlabel('Year')
% ylabel('Million km^2')
% title('All land area')


%% Plot timeseries: N loss

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
ignYrs = 0 ;
Nsmth = 5 ;
spacing = [0.1 0.1] ;   % [vert, horz]
%%%%%%%%%%%%%%%%%%%

% ts_nloss_bl = ts_nflux_flux_bl + ts_nflux_harvest_bl + ts_nflux_leach_bl + ts_nflux_LUch_bl ;
% ts_nloss_yr = ts_nflux_flux_yr + ts_nflux_harvest_yr + ts_nflux_leach_yr + ts_nflux_LUch_yr ;
ts_nloss_bl = ts_nflux_flux_bl + ts_nflux_leach_bl ;
ts_nloss_yr = ts_nflux_flux_yr + ts_nflux_leach_yr ;

[tmp_ts_nflux_flux_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_nflux_flux_bl, ts_nflux_flux_yr, ignYrs, yearList_future) ;
[tmp_ts_nflux_leach_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_nflux_leach_bl, ts_nflux_leach_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
plot(yearList_baseline,movmean(ts_nflux_flux_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_nflux_flux_yr,'LineWidth',lineWidth)
hold off
% legend(stdLegend, ...
%        'Location','NorthWest') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('TgN')
% title(['N loss: Gaseous' title_suffix])
title(['N loss: Gaseous'])

subplot_tight(1,2,2,spacing)
plot(yearList_baseline,movmean(ts_nflux_leach_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_nflux_leach_yr,'LineWidth',lineWidth)
hold off
% legend(stdLegend, ...
%        'Location','NorthWest') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('TgN')
% title(['N loss: Dissolved' title_suffix])
title(['N loss: Dissolved'])

clear tmp_*

if do_save
    export_fig([outDir_ts 'Nloss' file_suffix '.pdf'])
    close
end


%% Plot timeseries: N fertilizer

% % Options %%%%%%%%%
% lineWidth = 2 ;
% fontSize = 14 ;
% ignYrs = 0 ;
% Nsmth = 1 ;
% %%%%%%%%%%%%%%%%%%%
% 
% [tmp_ts_nflux_fert_yr, title_suffix, file_suffix] = ...
%     rebase_future2baseline(rebase, Nsmth, ts_nflux_fert_bl, ts_nflux_fert_yr, ignYrs, yearList_future) ;
% 
% figure('Position',figurePos,'Color','w') ;
% plot(yearList_baseline,movmean(-ts_nflux_fert_bl,Nsmth),'-k','LineWidth',lineWidth)
% hold on
% plot(yearList_future,-tmp_ts_nflux_fert_yr,'LineWidth',lineWidth)
% hold off
% legend(stdLegend, ...
%        'Location','NorthEastOutside') ;
% set(gca,'FontSize',fontSize)
% xlabel('Year')
% ylabel('TgN')
% title(['N fertilization' title_suffix])
% 
% clear tmp_*
% 
% if do_save
%     export_fig([outDir_ts 'Nfert' file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: Irrigation

% % Options %%%%%%%%%
% lineWidth = 2 ;
% fontSize = 14 ;
% ignYrs = 0 ;
% Nsmth = 5 ;
% %%%%%%%%%%%%%%%%%%%
% 
% [tmp_ts_irrig_yr, title_suffix, file_suffix] = ...
%     rebase_future2baseline(rebase, Nsmth, ts_irrig_bl, ts_irrig_yr, ignYrs, yearList_future) ;
% 
% figure('Position',figurePos,'Color','w') ;
% plot(yearList_baseline,movmean(ts_irrig_bl,Nsmth),'-k','LineWidth',lineWidth)
% hold on
% plot(yearList_future,tmp_ts_irrig_yr,'LineWidth',lineWidth)
% hold off
% legend(stdLegend, ...
%        'Location','NorthEastOutside') ;
% set(gca,'FontSize',fontSize)
% xlabel('Year')
% ylabel('1000 km^3')
% title(['Irrigation' title_suffix])
% 
% clear tmp_*
% 
% if do_save
%     export_fig([outDir_ts 'irrigation' file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: N fertilizer, irrigation

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
ignYrs = 0 ;
Nsmth = 1 ;
spacing = [0.1 0.1] ;   % vert, horiz
blYears = yearList_baseline ;
% blYears = 1960:yearList_baseline(end) ;
forPres = true ;
perArea = true ;
%%%%%%%%%%%%%%%%%%%

this_ts_nflux_fert_bl = ts_nflux_fert_bl ; % TgN
this_ts_nflux_fert_yr = ts_nflux_fert_yr ; % TgN
units_nflux_fert = 'TgN' ;
this_ts_irrig_bl = ts_irrig_bl ; % 1000km3
this_ts_irrig_yr = ts_irrig_yr ; % 1000km3
units_irrig = '1000 km^3' ;
if perArea
    this_ts_nflux_fert_bl = (ts_nflux_fert_bl * 1e9) ... % Tg to kg
                         ./ (ts_LUarea_crop_bl * 1e6*1e2) ; % Mkm2 to ha
    this_ts_nflux_fert_yr = (ts_nflux_fert_yr * 1e9) ... % Tg to kg
                         ./ (ts_LUarea_crop_yr * 1e6*1e2) ; % Mkm2 to ha
    units_nflux_fert = 'kg/ha' ;
end
if perArea
    this_ts_irrig_bl = (ts_irrig_bl * 1e3*1e6) ... % 1000km3 to L
                         ./ (ts_LUarea_crop_bl * 1e6*1e2) ; % Mkm2 to ha
    this_ts_irrig_yr = (ts_irrig_yr * 1e3*1e6) ... % Tg to kg
                         ./ (ts_LUarea_crop_yr * 1e6*1e2) ; % Mkm2 to ha
    units_irrig = 'L/ha' ;
end


[tmp_this_ts_nflux_fert_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, this_ts_nflux_fert_bl, this_ts_nflux_fert_yr, ignYrs, yearList_future) ;

if forPres
    title_suffix = '' ;
    file_suffix = [file_suffix '_4pres'] ;
end

figure('Position',figurePos,'Color','w') ;
subplot_tight(1,2,1,spacing)
tmp = this_ts_nflux_fert_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,movmean(-tmp,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,-tmp_this_ts_nflux_fert_yr,'LineWidth',lineWidth)
hold off
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel(units_nflux_fert)
title(['N fertilization' title_suffix])

clear tmp_*


% Plot timeseries: Irrigation
% Options %%%%%%%%%
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_this_ts_irrig_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, this_ts_irrig_bl, this_ts_irrig_yr, ignYrs, yearList_future) ;

if forPres
    title_suffix = '' ;
    file_suffix = [file_suffix '_4pres'] ;
end

subplot_tight(1,2,2,spacing)
tmp = this_ts_irrig_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_this_ts_irrig_yr,'LineWidth',lineWidth)
hold off
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel(units_irrig)
title(['Irrigation' title_suffix])

clear tmp_*

if do_save
    export_fig([outDir_ts 'mgmt_inputs' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Albedo

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_albedo1_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_albedo1_bl, ts_albedo1_yr, ignYrs, yearList_future) ;
[tmp_ts_albedo7_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_albedo7_bl, ts_albedo7_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;
plot(yearList_baseline,movmean(ts_albedo1_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_albedo1_yr,'LineWidth',lineWidth)
plot(yearList_baseline,movmean(ts_albedo7_bl,Nsmth),'--k','LineWidth',lineWidth)
ax = gca;
ax.ColorOrderIndex = 1;
plot(yearList_future,tmp_ts_albedo7_yr,'--','LineWidth',lineWidth)
hold off
legend([strcat(stdLegend,', Jan.') strcat(stdLegend,', Jul.')],...
       'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('Albedo')
title(['Albedo' title_suffix])

clear tmp_*

if do_save
    export_fig([outDir_ts 'albedo' file_suffix '.pdf'])
    close
end


%% Plot timeseries: C pools

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_cpool_VegC_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_cpool_VegC_bl, ts_cpool_VegC_yr, ignYrs, yearList_future) ;
[tmp_ts_cpool_Total_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_cpool_Total_bl, ts_cpool_Total_yr, ignYrs, yearList_future) ;


figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
plot(yearList_baseline,movmean(ts_cpool_VegC_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_cpool_VegC_yr,'LineWidth',lineWidth)
hold off
% legend(stdLegend, ...
%        'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('PgC')
title(['Vegetation C' title_suffix])

subplot_tight(1,2,2,spacing)
plot(yearList_baseline,movmean(ts_cpool_Total_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_cpool_Total_yr,'LineWidth',lineWidth)
hold off
% legend(stdLegend, ...
%        'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('PgC')
title(['Total C' title_suffix])

clear tmp_*

if do_save
    export_fig([outDir_ts 'cpools' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Monoterpene

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_amon_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_amon_bl, ts_amon_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;
plot(yearList_baseline,movmean(ts_amon_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_amon_yr,'LineWidth',lineWidth)
hold off
legend(stdLegend, ...
       'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('TgC')
title(['Monoterpene emissions' title_suffix])

clear tmp_*

if do_save
    export_fig([outDir_ts 'bvoc_mon' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Isoprene

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_aiso_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_aiso_bl, ts_aiso_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;
plot(yearList_baseline,movmean(ts_aiso_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_aiso_yr,'LineWidth',lineWidth)
hold off
legend(stdLegend, ...
       'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('TgC')
title(['Isoprene emissions' title_suffix])

clear tmp_*

if do_save
    export_fig([outDir_ts 'bvoc_iso' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Evapotranspiration

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 5 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_aevap_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_aevap_bl, ts_aevap_yr, ignYrs, yearList_future) ;
[tmp_ts_aaet_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_aaet_bl, ts_aaet_yr, ignYrs, yearList_future) ;

tmp_ts_aevapaaet_bl = ts_aevap_bl + ts_aaet_bl ;
tmp_ts_aevapaaet_yr = tmp_ts_aevap_yr + tmp_ts_aaet_yr ;

figure('Position',figurePos,'Color','w') ;
plot(yearList_baseline,movmean(tmp_ts_aevapaaet_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_aevapaaet_yr,'LineWidth',lineWidth)
hold off
legend(stdLegend, ...
       'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('1000 km^3')
title(['Evapotranspiration' title_suffix])

clear tmp_*

if do_save
    export_fig([outDir_ts 'evapotranspiration' file_suffix '.pdf'])
    close
end



%% Plot timeseries: Runoff

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 5 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_tot_runoff_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_tot_runoff_bl, ts_tot_runoff_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;
plot(yearList_baseline,movmean(ts_tot_runoff_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_tot_runoff_yr,'LineWidth',lineWidth)
hold off
legend(stdLegend, ...
       'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('1000 km^3')
title(['Runoff' title_suffix])

clear tmp_*

if do_save
    export_fig([outDir_ts 'runoff' file_suffix '.pdf'])
    close
end


%% Map end-of-run mean cropfracs

% Options %%%%%%%%%
title_prefix = 'Frac. crop' ;
whichCFTs = 'lpjg' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
y2include = 65:350 ;
%%%%%%%%%%%%%%%%%%%

Ny = 2 ;
Nmaps = length(runList) + 1 ;
if strcmp(whichCFTs,'lpjg')
    Ncrops = Ncrops_lpjg ;
    theseCrops = CFTnames ;
    Nx = 3 ;
    theseMaps = maps_cropfracs ;
    figure_position = figurePos ;
elseif strcmp(whichCFTs,'plum')
    Nx = 4 ;
    theseMaps = maps_cropfracs_plum7 ;
    figure_position = figurePos ;
else
    error(['whichCFTs (' whichCFTs ') not recognized!']) ;
end


for c = 1:Ncrops
    figure('Position',figure_position,'Color','w') ;
    thisCrop = theseCrops{c} ;
    for p = 1:Nmaps
        subplot_tight(Ny,Nx,p,spacing)
        if p==1
            tmp = theseMaps.maps_YXvyB(:,:,strcmp(theseMaps.varNames,thisCrop),:) ;
            tmp(maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'CROPLAND'),:)==0) = 0 ;
            tmp = mean(tmp,4) ;
            pcolor(tmp(y2include,:)) ;
            thisTitle = 'Baseline' ;
        else
            pcolor(mean(theseMaps.maps_YXvyr(y2include,:,strcmp(theseMaps.varNames,thisCrop),:,p-1),4)) ;
            thisTitle = runList{p-1} ;
        end
        shading flat ; axis equal tight off
        caxis([0 1]) ; colorbar('SouthOutside')
        title([thisCrop ' (' thisTitle ')'])
        set(gca,'FontSize',fontSize) ;
    end
    
    if do_save
        export_fig([outDir_maps 'maps_cropfracs_' thisCrop '.png'],['-r' num2str(pngres)])
        close
    end
end

clear theseMaps Ncrops Nmaps Nx Ny







