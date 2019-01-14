%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLUM2LPJG figures: Multiple runs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% thisVer = '20171129' ;
% thisVer = '20171201' ;
% thisVer = '20171205_SSP' ;
% thisVer = '20171205_SSP_irrBL' ;
thisVer = '20171212_SSP_irrBL' ;
% thisVer = '20171212_SSP_irrBL_constLU' ;
% thisVer = '20171212_SSP_irrBL_constLUco2' ;
% thisVer = '20171227' ;   % Implied _irrBL from here on
% thisVer = '20180108a' ;
% thisVer = '20180108b' ;
% thisVer = '20180111' ;
% thisVer = '20180115' ;
% thisVer = '20180125' ;

do_save = false ;
rebase = false ;

       
%% Setup

addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work/')
addpath(genpath('~/Documents/Dropbox/Dissertation/MATLAB work'))

trimFirstFuture = 0 ;
if strcmp(thisVer,'20171129')
    runList = {'RCP2.6','RCP4.5','RCP6.0','RCP8.5'} ;
    runDirs = {'PLUM2LPJG_26_s1_2011-2100/output-2017-11-25-102553', ...
               'PLUM2LPJG_45_s1_2011-2100/output-2017-11-25-110709', ...
               'PLUM2LPJG_60_s1_2011-2100/output-2017-11-25-104811', ...
               'PLUM2LPJG_85_s1_2011-2100/output-2017-11-25-112417'} ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJG_26_s1_1850-2010/output-2017-11-24-175840' ;
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20171201')
    runList = {'RCP2.6','RCP4.5','RCP6.0','RCP8.5'} ;
    runDir_base = addslashifneeded('/Users/Shared/PLUM/trunk_runs') ;
    runDirs = {'PLUM2LPJG_26_s1_2011-2100/output-2017-11-30-125032', ...
               'PLUM2LPJG_45_s1_2011-2100/output-2017-11-30-165834', ...
               'PLUM2LPJG_60_s1_2011-2100/output-2017-11-30-164245', ...
               'PLUM2LPJG_85_s1_2011-2100/output-2017-11-30-162354'} ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJG_26_s1_1850-2010/output-2017-11-29-163216' ;
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20171205_SSP')
%     runList = {'SSP1/RCP4.5','SSP2/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runList = {'SSP1/RCP4.5','SSP2/RCP4.5','SSP3/RCP6.0','SSP5/RCP8.5'} ;
    runDirs = {'PLUM2LPJG_SSP1_RCP45_mdnPLUM/output-2017-12-05-093016', ...
               'PLUM2LPJG_SSP2_RCP45_mdnPLUM/output-2017-12-05-093222', ...
               'PLUM2LPJG_SSP3_RCP60_mdnPLUM/output-2017-12-05-113149', ...
               'PLUM2LPJG_SSP4_RCP60_mdnPLUM/output-2017-12-05-121202', ...
               'PLUM2LPJG_SSP5_RCP85_mdnPLUM/output-2017-12-05-122918' ...
               } ;
    yearList_future = 2011:2100 ;
%     baselineDir = 'PLUM2LPJG_26_s1_1850-2010/output-2017-11-29-163216' ;
    baselineDir = 'PLUM2LPJG_26_s1_1850-2010/output-2017-12-02-203502' ; warning('Not using actual baselineDir!')
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20171205_SSP_irrBL')
%     runList = {'SSP1/RCP4.5','SSP2/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runList = {'SSP1/RCP4.5','SSP2/RCP4.5','SSP3/RCP6.0','SSP5/RCP8.5'} ;
    runDirs = {'PLUM2LPJG_SSP1_RCP45_mdnPLUM/output-2017-12-05-093016', ...
               'PLUM2LPJG_SSP2_RCP45_mdnPLUM/output-2017-12-05-093222', ...
               'PLUM2LPJG_SSP3_RCP60_mdnPLUM/output-2017-12-05-113149', ...
               ...'PLUM2LPJG_SSP4_RCP60_mdnPLUM/output-2017-12-05-121202', ...
               'PLUM2LPJG_SSP5_RCP85_mdnPLUM/output-2017-12-05-122918' ...
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJGblIrr_26_s1_1850-2010/output-2017-12-06-001347' ; warning('Not using actual baselineDir!')
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20171212_SSP_irrBL')
%     runList = {'SSP1/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runList = {'SSP1 (R45)','SSP3 (R60)','SSP4 (R60)','SSP5 (R85)'} ;
    runDirs = {'PLUM2LPJGblIrr_SSP1_RCP45_mdnPLUM/output-2017-12-07-095646', ...
               'PLUM2LPJGblIrr_SSP3_RCP60_mdnPLUM/output-2017-12-07-084408', ...
               'PLUM2LPJGblIrr_SSP4_RCP60_mdnPLUM/output-2017-12-07-144106', ...
               'PLUM2LPJGblIrr_SSP5_RCP85_mdnPLUM/output-2017-12-07-090531' ...
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJGblIrr_26_s1_1850-2010/output-2017-12-06-001347' ;
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20171212_SSP_irrBL_constLU')
    runList = {'SSP1_constLU/RCP4.5','SSP3_constLU/RCP6.0','SSP4_constLU/RCP6.0','SSP5_constLU/RCP8.5'} ;
    runDirs = {'PLUM2LPJGblIrr_45_s1_2011-2100_constSSP1LU/output-2017-12-10-232848', ...
               'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP3LU/output-2017-12-14-125023', ...
               'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP4LU/output-2017-12-15-015520', ...
               'PLUM2LPJGblIrr_85_s1_2011-2100_constSSP5LU/output-2017-12-15-011713' ...
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJGblIrr_26_s1_1850-2010/output-2017-12-06-001347' ;
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20171212_SSP_irrBL_constLUco2')
    runList = {'SSP1_constLU/RCP4.5_constCO2','SSP3_constLU/RCP6.0_constCO2','SSP4_constLU/RCP6.0_constCO2','SSP5_constLU/RCP8.5_constCO2'} ;
    runDirs = {'PLUM2LPJGblIrr_45_s1_2011-2100_constSSP1LU_constCO2/output-2017-12-11-115548', ...
               'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP3LU_constCO2/output-2017-12-11-180709', ...
               'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP4LU_constCO2/output-2017-12-11-203339', ...
               'PLUM2LPJGblIrr_85_s1_2011-2100_constSSP5LU_constCO2/output-2017-12-11-205323' ...
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJGblIrr_26_s1_1850-2010/output-2017-12-06-001347' ;
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20171227')
    runList = {'SSP1/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runDirs = {'PLUM2LPJGblIrr_SSP1v2_RCP45/output-2017-12-25-140114', ...
               'PLUM2LPJGblIrr_SSP3v2_RCP60/output-2017-12-25-133608', ...
               'PLUM2LPJGblIrr_SSP4v2_RCP60/output-2017-12-25-134011', ...
               'PLUM2LPJGblIrr_SSP5v2_RCP85/output-2017-12-25-135831' ...
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJGblIrr_26_s1_1850-2010/output-2017-12-06-001347' ;
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20180108a')
    runList = {'SSP1/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runDirs = {'PLUM2LPJGblIrr_SSP1v2_RCP45/output-2017-12-25-140114', ...
               'PLUM2LPJGblIrr_SSP3v2_RCP60/output-2017-12-25-133608', ...
               'PLUM2LPJGblIrr_SSP4v2_RCP60/output-2017-12-25-134011', ...
               'PLUM2LPJGblIrr_SSP5v2_RCP85/output-2017-12-25-135831' ...
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJGblIrrNoFG_26_s1_1850-2010_moreProcs/matlab_merge_20180108114616' ; warning('Not using actual baselineDir!')
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20180108b')
    runList = {'SSP1/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runDirs = {'PLUM2LPJGblIrr_SSP1v2_RCP45/output-2017-12-25-140114', ...
               'PLUM2LPJGblIrr_SSP3v2_RCP60/output-2017-12-25-133608', ...
               'PLUM2LPJGblIrr_SSP4v2_RCP60/output-2017-12-25-134011', ...
               'PLUM2LPJGblIrr_SSP5v2_RCP85/output-2017-12-25-135831' ...
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJGblIrrNoFGnoRye_26_s1_1850-2010_moreProcs/matlab_merge_20180108120140' ; warning('Not using actual baselineDir!')
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20180111')
    runList = {'SSP1/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runDirs = {'PLUM2LPJGblIrr_SSP1v2_RCP45.20180109/output-2018-01-11-081646' ;
               'PLUM2LPJGblIrr_SSP3v2_RCP60.20180109/output-2018-01-11-023306' ;
               'PLUM2LPJGblIrr_SSP4v2_RCP60.20180109/output-2018-01-11-023748' ;
               'PLUM2LPJGblIrr_SSP5v2_RCP85.20180109/output-2018-01-11-021841' ;
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJGblIrrNoFGnoRye_26_s1_1850-2010_moreProcs/matlab_merge_20180108120140' ; warning('Not using actual baselineDir!')
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20180115')
    runList = {'SSP1/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runDirs = {'PLUM2LPJGblIrr_SSP1v2_RCP45.20180109/output-2018-01-13-230513' ;
               'PLUM2LPJGblIrr_SSP3v2_RCP60.20180109/output-2018-01-13-231353' ;
               'PLUM2LPJGblIrr_SSP4v2_RCP60.20180109/output-2018-01-13-231930' ;
               'PLUM2LPJGblIrr_SSP5v2_RCP85.20180109/output-2018-01-13-225751' ;
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'PLUM2LPJGblIrrNoFGnoRye_26_s1_1850-2010_moreProcs/matlab_merge_20180108120140' ; warning('Not using actual baselineDir!')
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20180125')
    runList = {'SSP1/RCP4.5','SSP3/RCP6.0','SSP4/RCP6.0','SSP5/RCP8.5'} ;
    runDirs = {'PLUM2LPJG.PLUM7.SSP1v2.RCP45/output-2018-02-02-054242' ;
               'PLUM2LPJG.PLUM7.SSP3v2.RCP60/output-2018-01-25-090934' ;
               'PLUM2LPJG.PLUM7.SSP4v2.RCP60/output-2018-01-25-090916' ;
               'PLUM2LPJG.PLUM7.SSP5v2.RCP85/output-2018-01-24-182830' ;
               } ;
    yearList_future = 2011:2100 ;
%     baselineDir = 'PLUM2LPJG_PLUM7_1850-2010/matlab_merge_20180124084436' ; warning('Not using actual baselineDir!')
    baselineDir = 'PLUM2LPJG_PLUM7_1850-2010/matlab_merge_20180422140716' ; warning('Not using actual baselineDir!')
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
if strcmp(thisVer,'20180108b')
    fao = load('/Users/Shared/PLUM/crop_calib_data/fao/FAOdata_1961-2010_calibVer10.mat') ;
else
    fao = load('/Users/Shared/PLUM/crop_calib_data/fao/FAOdata_1961-2010_calibVer9.mat') ;
end


%% Import baseline

% Baseline run
ts_tmp = load([baselineDir 'timeseries.mat']) ;
bl_ts_fields = fieldnames(ts_tmp) ;
for f = 1:length(bl_ts_fields)
    thisField_in = bl_ts_fields{f} ;
    thisName_out = ['ts_' strrep(thisField_in,'_ts','') '_bl'] ;
    eval([thisName_out ' = ts_tmp.' thisField_in ' ;']) ;
end

% Get CFT names
CFTnames = strrep({bl_ts_fields{getCellsWithString(bl_ts_fields,'croparea')}}','croparea_ts_','') ;
if any(not(cellfun(@isempty,strfind(CFTnames,'TeWW'))))
    is_lpjg = cellfun(@length,CFTnames)==4 & ~strcmp(CFTnames,'Rice') ;
    CFTnames_lpjg = {CFTnames{is_lpjg}}' ;
    CFTnames_plum = sort({CFTnames{~is_lpjg}}') ;
elseif any(not(cellfun(@isempty,strfind(CFTnames,'_lpjg'))))
    CFTnames_lpjg = {CFTnames{not(cellfun(@isempty,strfind(CFTnames,'_lpjg')))}}' ;
    CFTnames_plum = {CFTnames{not(cellfun(@isempty,strfind(CFTnames,'_plum')))}}' ;
    if length(CFTnames_lpjg)+length(CFTnames_plum) ~= length(CFTnames)
        error('length(CFTnames_lpjg)+length(CFTnames_plum) ~= length(CFTnames))')
    end
else
    tmp = {'CerealsC3','CerealsC3s','CerealsC3w','CerealsC4','Miscanthus','Oilcrops','Pulses','Rice','StarchyRoots'} ;
    error('This isn''t going to work out well. You''ve combined LPJG and PLUM crops.')
    [C,IA,IB] = intersect(CFTnames,tmp) ;
    if length(CFTnames)==length(tmp) && isequal(sort(IA)',1:length(CFTnames))
        CFTnames_lpjg = CFTnames ;
        CFTnames_lpjg(strcmp(CFTnames_lpjg,'CerealsC3')) = [] ;
        CFTnames_plum = CFTnames ;
        CFTnames_plum(strcmp(CFTnames_plum,'CerealsC3s')) = [] ;
        CFTnames_plum(strcmp(CFTnames_plum,'CerealsC3w')) = [] ;
    else
        error('I don''t know how to parse this set of CFTnames!')
    end
    clear tmp
end
if isempty(CFTnames_lpjg)
    warning('isempty(CFTnames_lpjg)')
elseif isempty(CFTnames_plum)
    warning('isempty(CFTnames_plum)')
end
disp('CFTnames_lpjg:') ; disp(CFTnames_lpjg)
disp('CFTnames_plum:') ; disp(CFTnames_plum)

Ncrops_lpjg = length(CFTnames_lpjg) ;
Ncrops_plum = length(CFTnames_plum) ;

% Check whether rainfed and irrigated crops are separated
combine_rf_ir = false ;
for c = 1:Ncrops_lpjg
    thisCrop = CFTnames_lpjg{c} ;
    thisCropi = [thisCrop 'i'] ;
    if strcmp(CFTnames_lpjg,thisCropi)
        warning('Combining rainfed and irrigated outputs...')
        combine_rf_ir = true ;
        break
    end
end

% Add "irrigated" to "rainfed"
if combine_rf_ir
    error('Finish code to do this for baseline run!')
%     for c = 1:Ncrops_lpjg
%         thisCrop = CFTnames_lpjg{c} ;
%         thisCropi = [thisCrop 'i'] ;
%         eval(['ts_croparea_' thisCrop '_bl = ts_croparea_' thisCrop '_bl + ts_tmp.croparea_ts_' thisCrop 'i ;']) ;
%         eval(['ts_gsirrig_' thisCrop '_bl = ts_gsirrig_' thisCrop '_bl + ts_tmp.gsirrig_ts_' thisCrop 'i ;']) ;
%         eval(['ts_yield_' thisCrop '_bl = ts_yield_' thisCrop '_bl + ts_tmp.yield_ts_' thisCrop 'i ;']) ;
%     end
end

clear ts_tmp

% Adjust FAO data
if Ncrops_plum>0
    if any(strcmp(fao.listCrops_fa2o,'Wheat'))
        fao.listCrops_fa2o{strcmp(fao.listCrops_fa2o,'Wheat')} = 'CerealsC3' ;
    end
    if any(strcmp(fao.listCrops_fa2o,'Maize'))
        fao.listCrops_fa2o{strcmp(fao.listCrops_fa2o,'Maize')} = 'CerealsC4' ;
    end
    if any(strcmp(fao.listCrops_fa2o,'Starchy roots'))
        fao.listCrops_fa2o{strcmp(fao.listCrops_fa2o,'Starchy roots')} = 'StarchyRoots' ;
    end
    if any(not(isempty(strfind(CFTnames_plum,'_plum'))))
        fao.listCrops_fa2o = strcat(fao.listCrops_fa2o,'_plum') ;
    end
    for c = 1:Ncrops_plum
        thisCrop = CFTnames_plum{c} ;
% % %         thisCrop = strrep(CFTnames_plum{c},'_plum','') ;
        if length(find(strcmp(fao.listCrops_fa2o,thisCrop)))==1
            eval(['ts_cropprod_' thisCrop '_fao = cf_t2kg*squeeze(nansum(fao.tmp_total_fa2_Ccy(:,strcmp(fao.listCrops_fa2o,thisCrop),:),1)) ;']) ;
            eval(['ts_croparea_' thisCrop '_fao = cf_ha2Mkm2*squeeze(nansum(fao.tmp_croparea_fa2_Ccy(:,strcmp(fao.listCrops_fa2o,thisCrop),:),1)) ;']) ;
        else
            warning(['Assuming zeros for FAO data for ' thisCrop])
            eval(['ts_cropprod_' thisCrop '_fao = zeros(length(fao.tmp_fao_yearList),1) ;']) ;
            eval(['ts_croparea_' thisCrop '_fao = zeros(length(fao.tmp_fao_yearList),1) ;']) ;
        end
    end
end

lastdec_tmp = load([baselineDir 'last_decade.mat']) ;
maps_LU.maps_YXvyB = lastdec_tmp.LU.maps_YXvy ;
maps_cropfracs.maps_YXvyB = lastdec_tmp.cropfracs.maps_YXvy ;
if isfield(lastdec_tmp,'cropfracs_plum7')
    maps_cropfracs_plum7.maps_YXvyB = lastdec_tmp.cropfracs_plum7.maps_YXvy ;
end
% clear lastdec_tmp

% Get variable names
tmp = whos('ts_*_bl') ;
vars_ts_bl = {tmp.name}' ;
clear tmp


%% Future runs

if combine_rf_ir
    error('Write code to do this for future runs!')
end

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
        maps_all.list_to_map = lastdec_tmp.LU.list_to_map ;
        maps_all.yearList = lastdec_tmp.LU.yearList ;
        
        % Get maps that were read in baseline
        tmp = whos('maps_*') ;
        vars_maps_bl = {tmp.name}' ;
        vars_maps_bl(strcmp(vars_maps_bl,'maps_all')) = [] ;
        clear tmp
        for v = 1:length(vars_maps_bl)
            thisVar_in = vars_maps_bl{v} ;
            thisVar_out = strrep(thisVar_in,'maps_','') ;
            eval([thisVar_in '.maps_YXvyr = nan([size(lastdec_tmp.' thisVar_out '.maps_YXvy) Nruns]) ;']) ;
            eval([thisVar_in '.varNames = lastdec_tmp.' thisVar_out '.varNames ;']) ;
            clear thisVar*
        end ; clear v
%         maps_LU.maps_YXvyr = nan([size(lastdec_tmp.LU.maps_YXvy) Nruns]) ;
%         maps_LU.varNames = lastdec_tmp.LU.varNames ;
%         maps_cropfracs.maps_YXvyr = nan([size(lastdec_tmp.cropfracs.maps_YXvy) Nruns]) ;
%         maps_cropfracs.varNames = lastdec_tmp.cropfracs.varNames ;
%         maps_cropfracs_plum7.maps_YXvyr = nan([size(lastdec_tmp.cropfracs_plum7.maps_YXvy) Nruns]) ;
%         maps_cropfracs_plum7.varNames = lastdec_tmp.cropfracs_plum7.varNames ;
    end
    
    for f = 1:length(bl_ts_fields)
        thisField_in = bl_ts_fields{f} ;
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
        thisVar_out = strrep(thisVar_in,'maps_','') ;
        eval([thisVar_in '.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.' thisVar_out '.maps_YXvy ;']) ;
    end
%     maps_LU.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.LU.maps_YXvy ;
end ; clear r
disp('Done with import.')

disp('Processing...')
nanmask = isnan(maps_LU.maps_YXvyB(:,:,1,1)) ;

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

ntrl_area_YXB = land_area_YX .* maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'NATURAL'),end) ;
ntrl_area_YXr = repmat(land_area_YX,[1 1 Nruns]) .* squeeze(maps_LU.maps_YXvyr(:,:,strcmp(maps_LU.varNames,'NATURAL'),end,:)) ;
ntrl_diff_YXr = ntrl_area_YXr - repmat(ntrl_area_YXB,[1 1 Nruns]) ;

bare_area_YXB = land_area_YX .* maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'BARREN'),end) ;
bare_area_YXr = repmat(land_area_YX,[1 1 Nruns]) .* squeeze(maps_LU.maps_YXvyr(:,:,strcmp(maps_LU.varNames,'BARREN'),end,:)) ;
bare_diff_YXr = bare_area_YXr - repmat(bare_area_YXB,[1 1 Nruns]) ;

crop_area_YXB = land_area_YX .* maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'CROPLAND'),end) ;
crop_area_YXr = repmat(land_area_YX,[1 1 Nruns]) .* squeeze(maps_LU.maps_YXvyr(:,:,strcmp(maps_LU.varNames,'CROPLAND'),end,:)) ;
crop_diff_YXr = crop_area_YXr - repmat(crop_area_YXB,[1 1 Nruns]) ;

past_area_YXB = land_area_YX .* maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'PASTURE'),end) ;
past_area_YXr = repmat(land_area_YX,[1 1 Nruns]) .* squeeze(maps_LU.maps_YXvyr(:,:,strcmp(maps_LU.varNames,'PASTURE'),end,:)) ;
past_diff_YXr = past_area_YXr - repmat(past_area_YXB,[1 1 Nruns]) ;

agri_area_YXB = crop_area_YXB + past_area_YXB ;
agri_area_YXr = crop_area_YXr + past_area_YXr ;
agri_diff_YXr = crop_diff_YXr + past_diff_YXr ;

% Check that order of maps is correct
for v = 1:length(maps_LU.varNames)
    thisLU = maps_LU.varNames{v} ;
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
    eval(['fromMap = 1e-6*nansum(nansum(' thisLU_short '_area_YXB)) ;']) ;
    if abs(fromTS-fromMap) > 1e-9
        error('LU variables are in wrong order!')
    end
end

% Calculate actual YIELD (production per area)
kg_2_t = 1e-3 ;
Mkm2_2_ha = 1e6*1e2 ;
tmp = whos('ts_croppro*') ;
tmp_name = {tmp.name}' ;
for i = 1:length(tmp_name)
    thisName_cropprod = tmp_name{i} ;
    is_expYield = ~isempty(strfind(thisName_cropprod,'cropproExp')) ;
    if is_expYield
        thisname_yield = strrep(thisName_cropprod,'cropproExp','yielExp') ;
        thisname_croparea = strrep(thisName_cropprod,'cropproExp','croparea') ;
    else
        thisname_yield = strrep(thisName_cropprod,'cropprod','yield') ;
        thisname_croparea = strrep(thisName_cropprod,'cropprod','croparea') ;
    end
    eval([thisname_yield ' = (' thisName_cropprod ' * kg_2_t) ./ (' thisname_croparea '*Mkm2_2_ha) ;']) ;
    clear this*
end ; clear i
clear tmp*

disp('Done.')


%% Map changes in LU area

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
%%%%%%%%%%%%%%%%%%%

% Natural area
figure('Color','w','Position',figurePos) ;
subplot_tight(2,3,1,spacing) ;
pcolor(ntrl_area_YXB(69:end,:)) ; shading flat ; axis equal tight off
total_bl = 1e-6*nansum(nansum(ntrl_area_YXB)) ;
text(textX,textY_1,[num2str(round(total_bl,1)) ' Mkm^2'],'FontSize',fontSize) ;
% colormap(gca,brewermap(64,'YlOrBr'))
colorbar('SouthOutside')
set(gca,'XTick',[],'YTick',[])
set(gca,'FontSize',fontSize)
title(['"Natural" area, ' num2str(yearList_baseline(end)) ' (km^2)'])
colorlim = 0 ;
for r = 1:Nruns
    subplot_tight(2,3,ssp_plot_index(r),spacing) ;
    pcolor(ntrl_diff_YXr(69:end,:,r)) ; shading flat ; axis equal tight off
    colorbar('SouthOutside')
    colormap(gca,brighten(brewermap(64,'RdBu_ssr'),-0.3))
    set(gca,'XTick',[],'YTick',[])
    set(gca,'FontSize',fontSize)
    title(['\Delta "Natural" area, ' num2str(yearList_future(end)) ': ' runList{r} ' (km^2)'])
    total_yr = 1e-6*nansum(nansum(ntrl_diff_YXr(:,:,r))) ;
    text(textX,textY_1,[num2str(round(total_yr,1)) ' Mkm^2'],'FontSize',fontSize) ;
    pctDiff = round(100*total_yr./total_bl,1) ;
    pctDiff_str = num2str(pctDiff) ;
    if pctDiff>0
        pctDiff_str = ['+' pctDiff_str] ;
    end
    text(textX,textY_2,['(' pctDiff_str ' %)'],'FontSize',fontSize) ;
    % Set up for changing color axis limits
    gcas{r} = gca ;
    colorlim = max(colorlim, max(abs(caxis))) ;
end
for r = 1:Nruns
    caxis(gcas{r},[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_LU_ntrl.png'],'-r300')
    close
end

% Cropland area
figure('Color','w','Position',figurePos) ;
subplot_tight(2,3,1,spacing) ;
pcolor(crop_area_YXB(69:end,:)) ; shading flat ; axis equal tight off
total_bl = 1e-6*nansum(nansum(crop_area_YXB)) ;
text(textX,textY_1,[num2str(round(total_bl,1)) ' Mkm^2'],'FontSize',fontSize) ;
% colormap(gca,brewermap(64,'YlGn'))
colorbar('SouthOutside')
set(gca,'XTick',[],'YTick',[])
set(gca,'FontSize',fontSize)
title(['Cropland area, ' num2str(yearList_baseline(end)) ' (km^2)'])
colorlim = 0 ;
for r = 1:Nruns
    subplot_tight(2,3,ssp_plot_index(r),spacing) ;
    pcolor(crop_diff_YXr(69:end,:,r)) ; shading flat ; axis equal tight off
    colorbar('SouthOutside')
    colormap(gca,brighten(brewermap(64,'RdBu_ssr'),-0.3))
    set(gca,'XTick',[],'YTick',[])
    set(gca,'FontSize',fontSize)
    title(['\Delta Cropland area, ' num2str(yearList_future(end)) ': ' runList{r} ' (km^2)'])
    total_yr = 1e-6*nansum(nansum(crop_diff_YXr(:,:,r))) ;
    text(textX,textY_1,[num2str(round(total_yr,1)) ' Mkm^2'],'FontSize',fontSize) ;
    pctDiff = round(100*total_yr./total_bl,1) ;
    pctDiff_str = num2str(pctDiff) ;
    if pctDiff>0
        pctDiff_str = ['+' pctDiff_str] ;
    end
    text(textX,textY_2,['(' pctDiff_str ' %)'],'FontSize',fontSize) ;
    % Set up for changing color axis limits
    gcas{r} = gca ;
    colorlim = max(colorlim, max(abs(caxis))) ;
end
for r = 1:Nruns
    caxis(gcas{r},[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_LU_crop.png'],'-r300')
    close
end

% Pasture area
figure('Color','w','Position',figurePos) ;
subplot_tight(2,3,1,spacing) ;
pcolor(past_area_YXB(69:end,:)) ; shading flat ; axis equal tight off
total_bl = 1e-6*nansum(nansum(past_area_YXB)) ;
text(textX,textY_1,[num2str(round(total_bl,1)) ' Mkm^2'],'FontSize',fontSize) ;
% colormap(gca,brewermap(64,'YlGn'))
colorbar('SouthOutside')
set(gca,'XTick',[],'YTick',[])
set(gca,'FontSize',fontSize)
title(['"Pasture" area, ' num2str(yearList_baseline(end)) ' (km^2)'])
colorlim = 0 ;
for r = 1:Nruns
    subplot_tight(2,3,ssp_plot_index(r),spacing) ;
    pcolor(past_diff_YXr(69:end,:,r)) ; shading flat ; axis equal tight off
    colorbar('SouthOutside')
    colormap(gca,brighten(brewermap(64,'RdBu_ssr'),-0.3))
    set(gca,'XTick',[],'YTick',[])
    set(gca,'FontSize',fontSize)
    title(['\Delta "Pasture" area, ' num2str(yearList_future(end)) ': ' runList{r} ' (km^2)'])
    total_yr = 1e-6*nansum(nansum(past_diff_YXr(:,:,r))) ;
    text(textX,textY_1,[num2str(round(total_yr,1)) ' Mkm^2'],'FontSize',fontSize) ;
    pctDiff = round(100*total_yr./total_bl,1) ;
    pctDiff_str = num2str(pctDiff) ;
    if pctDiff>0
        pctDiff_str = ['+' pctDiff_str] ;
    end
    text(textX,textY_2,['(' pctDiff_str ' %)'],'FontSize',fontSize) ;
    % Set up for changing color axis limits
    gcas{r} = gca ;
    colorlim = max(colorlim, max(abs(caxis))) ;
end
for r = 1:Nruns
    caxis(gcas{r},[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_LU_past.png'],'-r300')
    close
end

% Agricultural area
figure('Color','w','Position',figurePos) ;
subplot_tight(2,3,1,spacing) ;
pcolor(agri_area_YXB(69:end,:)) ; shading flat ; axis equal tight off
total_bl = 1e-6*nansum(nansum(agri_area_YXB)) ;
text(textX,textY_1,[num2str(round(total_bl,1)) ' Mkm^2'],'FontSize',fontSize) ;
% colormap(gca,brewermap(64,'YlGn'))
colorbar('SouthOutside')
set(gca,'XTick',[],'YTick',[])
set(gca,'FontSize',fontSize)
title(['Agricultural area, ' num2str(yearList_baseline(end)) ' (km^2)'])
colorlim = 0 ;
for r = 1:Nruns
    subplot_tight(2,3,ssp_plot_index(r),spacing) ;
    pcolor(agri_diff_YXr(69:end,:,r)) ; shading flat ; axis equal tight off
    colorbar('SouthOutside')
    colormap(gca,brighten(brewermap(64,'RdBu_ssr'),-0.3))
    set(gca,'XTick',[],'YTick',[])
    set(gca,'FontSize',fontSize)
    title(['\Delta Agricultural area, ' num2str(yearList_future(end)) ': ' runList{r} ' (km^2)'])
    total_yr = 1e-6*nansum(nansum(agri_diff_YXr(:,:,r))) ;
    text(textX,textY_1,[num2str(round(total_yr,1)) ' Mkm^2'],'FontSize',fontSize) ;
    pctDiff = round(100*total_yr./total_bl,1) ;
    pctDiff_str = num2str(pctDiff) ;
    if pctDiff>0
        pctDiff_str = ['+' pctDiff_str] ;
    end
    text(textX,textY_2,['(' pctDiff_str ' %)'],'FontSize',fontSize) ;
    % Set up for changing color axis limits
    gcas{r} = gca ;
    colorlim = max(colorlim, max(abs(caxis))) ;
end
for r = 1:Nruns
    caxis(gcas{r},[-colorlim colorlim])
end
if do_save
    export_fig([outDir_maps 'areaDiff_LU_agri.png'],'-r300')
    close
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
        export_fig([outDir_maps 'areaDiff_BDhotspots_CI.png'],'-r300')
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
        export_fig([outDir_maps 'areaDiff_BDhotspots_glob200.png'],'-r300')
        close
    end

else
    warning('Skipping hotspots')
end


%% Plot timeseries: Land uses

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
% blYears = yearList_baseline ;
blYears = 1960:yearList_baseline(end) ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_LUarea_crop_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_crop_bl, ts_LUarea_crop_yr, ignYrs, yearList_future) ;
[tmp_ts_LUarea_past_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_past_bl, ts_LUarea_past_yr, ignYrs, yearList_future) ;
[tmp_ts_LUarea_ntrl_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_past_bl, ts_LUarea_ntrl_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
tmp = ts_LUarea_past_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_LUarea_past_yr,'-','LineWidth',lineWidth)
tmp = ts_LUarea_crop_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,movmean(tmp,Nsmth),'--k','LineWidth',lineWidth)
set(gca,'ColorOrderIndex',1) ;
plot(yearList_future,tmp_ts_LUarea_crop_yr,'--','LineWidth',lineWidth)
hold off
% legend([strcat(stdLegend,', pasture') strcat(stdLegend,', cropland')], ...
%        'Location','NorthEastOutside') ;
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
% legend(stdLegend, ...
%        'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('Million km^2')
title(['Natural area' title_suffix])

if do_save
    export_fig([outDir_ts 'landUse' file_suffix '.pdf'])
    close
end



%% Plot timeseries: total area

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
Nsmth = 5 ;
spacing = [0.1 0.1] ;   % vert, horiz
% blYears = yearList_baseline ;
blYears = 1960:yearList_baseline(end) ;
forPres = true ;
rebase = true ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_nflux_fert_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_nflux_fert_bl, ts_nflux_fert_yr, ignYrs, yearList_future) ;

if forPres
    title_suffix = '' ;
    file_suffix = [file_suffix '_4pres'] ;
end

figure('Position',figurePos,'Color','w') ;
subplot_tight(1,2,1,spacing)
tmp = ts_nflux_fert_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,movmean(-tmp,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,-tmp_ts_nflux_fert_yr,'LineWidth',lineWidth)
hold off
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('TgN')
title(['N fertilization' title_suffix])

clear tmp_*


% Plot timeseries: Irrigation
% Options %%%%%%%%%
% ignYrs = 0 ;
% Nsmth = 5 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_irrig_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_irrig_bl, ts_irrig_yr, ignYrs, yearList_future) ;

if forPres
    title_suffix = '' ;
    file_suffix = [file_suffix '_4pres'] ;
end

subplot_tight(1,2,2,spacing)
tmp = ts_irrig_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,tmp_ts_irrig_yr,'LineWidth',lineWidth)
hold off
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('1000 km^3')
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
    theseCrops = CFTnames_lpjg ;
    Nx = 3 ;
    theseMaps = maps_cropfracs ;
    figure_position = figurePos ;
elseif strcmp(whichCFTs,'plum')
    Ncrops = Ncrops_plum ;
    theseCrops = CFTnames_plum ;
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
        export_fig([outDir_maps 'maps_cropfracs_' thisCrop '.png'],'-r150')
        close
    end
end

clear theseMaps Ncrops Nmaps Nx Ny


%% Plot timeseries: Production of each crop, lpjg

% Options %%%%%%%%%
thisVar = 'cropprod' ;
units = 'Mt DM' ;
title_prefix = 'Production' ;
whichCFTs = 'lpjg' ;
conv_fact = cf_kg2Mt ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

if strcmp(whichCFTs,'lpjg')
    theseCFTnames = CFTnames_lpjg ;
    thisLegend = stdLegend ;
elseif strcmp(whichCFTs,'plum')
    theseCFTnames = CFTnames_plum ;
    thisLegend = stdLegend_plusFAO ;
else
    error(['whichCFTs (' whichCFTs ') not recognized!'])
end

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if strcmp(whichCFTs,'plum')
    cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
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


%% Plot timeseries: Production of each crop, plum

% Options %%%%%%%%%
thisVar = 'cropprod' ;
units = 'Mt DM' ;
title_prefix = 'Production' ;
whichCFTs = 'plum' ;
conv_fact = cf_kg2Mt ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

if strcmp(whichCFTs,'lpjg')
    theseCFTnames = CFTnames_lpjg ;
    thisLegend = stdLegend ;
elseif strcmp(whichCFTs,'plum')
    theseCFTnames = CFTnames_plum ;
    thisLegend = stdLegend_plusFAO ;
else
    error(['whichCFTs (' whichCFTs ') not recognized!'])
end

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if strcmp(whichCFTs,'plum')
    cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
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


%% Plot timeseries: Production of each crop, plum EXPECTED

% Options %%%%%%%%%
thisVar = 'cropproExp' ;
expected = true ;
units = 'Mt DM' ;
title_prefix = 'Production (exp.)' ;
whichCFTs = 'plum' ;
conv_fact = cf_kg2Mt ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

if strcmp(whichCFTs,'lpjg')
    theseCFTnames = CFTnames_lpjg ;
    thisLegend = stdLegend ;
elseif strcmp(whichCFTs,'plum')
    theseCFTnames = CFTnames_plum ;
    thisLegend = stdLegend_plusFAO ;
else
    error(['whichCFTs (' whichCFTs ') not recognized!'])
end

clear cell_bl cell_yr
if strcmp(thisVar,'cropproExp')
    cmds = get_cell_forPlot(whos, 'cropprod', 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
else
    cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if strcmp(whichCFTs,'plum')
    if strcmp(thisVar,'cropproExp')
        cmds = get_cell_forPlot(whos, 'cropprod', 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
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


%% Plot timeseries: Area of each crop, lpjg

% Options %%%%%%%%%
thisVar = 'croparea' ;
units = 'Million km^2' ;
title_prefix = 'Area' ;
whichCFTs = 'lpjg' ;
conv_fact = 1 ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

if strcmp(whichCFTs,'lpjg')
    theseCFTnames = CFTnames_lpjg ;
    thisLegend = stdLegend ;
elseif strcmp(whichCFTs,'plum')
    theseCFTnames = CFTnames_plum ;
    thisLegend = stdLegend_plusFAO ;
else
    error(['whichCFTs (' whichCFTs ') not recognized!'])
end

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if strcmp(whichCFTs,'plum')
    cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
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


%% Plot timeseries: Area of each crop, plum

% Options %%%%%%%%%
thisVar = 'croparea' ;
units = 'Million km^2' ;
title_prefix = 'Area' ;
whichCFTs = 'plum' ;
conv_fact = 1 ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

if strcmp(whichCFTs,'lpjg')
    theseCFTnames = CFTnames_lpjg ;
    thisLegend = stdLegend ;
elseif strcmp(whichCFTs,'plum')
    theseCFTnames = CFTnames_plum ;
    thisLegend = stdLegend_plusFAO ;
else
    error(['whichCFTs (' whichCFTs ') not recognized!'])
end

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if strcmp(whichCFTs,'plum')
    cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
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


%% Plot timeseries: Yield of each crop, lpjg

% Options %%%%%%%%%
thisVar = 'yield' ;
units = 't ha^{-1}' ;
title_prefix = 'Yield' ;
whichCFTs = 'lpjg' ;
conv_fact = 1 ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

if strcmp(whichCFTs,'lpjg')
    theseCFTnames = CFTnames_lpjg ;
    thisLegend = stdLegend ;
elseif strcmp(whichCFTs,'plum')
    theseCFTnames = CFTnames_plum ;
    thisLegend = stdLegend_plusFAO ;
else
    error(['whichCFTs (' whichCFTs ') not recognized!'])
end

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if strcmp(whichCFTs,'plum')
    cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
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


%% Plot timeseries: Yield of each crop, plum

% Options %%%%%%%%%
thisVar = 'yield' ;
units = 't ha^{-1}' ;
title_prefix = 'Yield' ;
whichCFTs = 'plum' ;
conv_fact = 1 ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 5 ;
%%%%%%%%%%%%%%%%%%%

if strcmp(whichCFTs,'lpjg')
    theseCFTnames = CFTnames_lpjg ;
    thisLegend = stdLegend ;
elseif strcmp(whichCFTs,'plum')
    theseCFTnames = CFTnames_plum ;
    thisLegend = stdLegend_plusFAO ;
else
    error(['whichCFTs (' whichCFTs ') not recognized!'])
end

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if strcmp(whichCFTs,'plum')
    cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
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
    theseCFTnames = CFTnames_lpjg ;
    thisLegend = stdLegend ;
elseif strcmp(whichCFTs,'plum')
    theseCFTnames = CFTnames_plum ;
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






