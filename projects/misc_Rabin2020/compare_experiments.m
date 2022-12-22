%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLUM2LPJG figures: Compare experiments %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% thisVer = '20171212_irrBL_SSP_vs_constLU' ;
% thisVer = '20171212_irrBL_SSP_vs_constLUco2' ;
% thisVer = '20171212_irrBL_constLU_vs_constLUco2' ;
% thisVer = '20180424_agmip7_SSP_vs_constLU' ;
thisVer = '20180606_SSP_vs_constPLUM2011' ;

do_save = false ;
rebase = false ;


%% Setup

if strcmp(thisVer,'20171212_irrBL_SSP_vs_constLU')
    runList = {'SSP1 (R45)','SSP3 (R60)','SSP4 (R60)','SSP5 (R85)',...
               'SSP1constLU (R45)','SSP3constLU (R60)','SSP4constLU (R60)','SSP5constLU (R85)'} ;
    runDir_base = addslashifneeded('/Users/Shared/PLUM/trunk_runs') ;
    runDirs = {[runDir_base 'PLUM2LPJGblIrr_SSP1_RCP45_mdnPLUM/output-2017-12-07-095646'], ...
               [runDir_base 'PLUM2LPJGblIrr_SSP3_RCP60_mdnPLUM/output-2017-12-07-084408'], ...
               [runDir_base 'PLUM2LPJGblIrr_SSP4_RCP60_mdnPLUM/output-2017-12-07-144106'], ...
               [runDir_base 'PLUM2LPJGblIrr_SSP5_RCP85_mdnPLUM/output-2017-12-07-090531'], ...
               [runDir_base 'PLUM2LPJGblIrr_45_s1_2011-2100_constSSP1LU/output-2017-12-10-232848'], ...
               [runDir_base 'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP3LU/output-2017-12-14-125023'], ...
               [runDir_base 'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP4LU/output-2017-12-15-015520'], ...
               [runDir_base 'PLUM2LPJGblIrr_85_s1_2011-2100_constSSP5LU/output-2017-12-15-011713'] ...
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = addslashifneeded([runDir_base 'PLUM2LPJGblIrr_26_s1_1850-2010/output-2017-12-06-001347']) ;
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20171212_irrBL_SSP_vs_constLUco2')
    runList = {'SSP1 (R45)','SSP3 (R60)','SSP4 (R60)','SSP5 (R85)',...
               'SSP1constLU (R45constCO2)','SSP3constLU (R60constCO2)','SSP4constLU (R60constCO2)','SSP5constLU (R85constCO2)'} ;
    runDir_base = addslashifneeded('/Users/Shared/PLUM/trunk_runs') ;
    runDirs = {[runDir_base 'PLUM2LPJGblIrr_SSP1_RCP45_mdnPLUM/output-2017-12-07-095646'], ...
               [runDir_base 'PLUM2LPJGblIrr_SSP3_RCP60_mdnPLUM/output-2017-12-07-084408'], ...
               [runDir_base 'PLUM2LPJGblIrr_SSP4_RCP60_mdnPLUM/output-2017-12-07-144106'], ...
               [runDir_base 'PLUM2LPJGblIrr_SSP5_RCP85_mdnPLUM/output-2017-12-07-090531'], ...
               [runDir_base 'PLUM2LPJGblIrr_45_s1_2011-2100_constSSP1LU_constCO2/output-2017-12-11-115548'], ...
               [runDir_base 'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP3LU_constCO2/output-2017-12-11-180709'], ...
               [runDir_base 'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP4LU_constCO2/output-2017-12-11-203339'], ...
               [runDir_base 'PLUM2LPJGblIrr_85_s1_2011-2100_constSSP5LU_constCO2/output-2017-12-11-205323'] ...
               } ;
elseif strcmp(thisVer,'20171212_irrBL_constLU_vs_constLUco2')
    runList = {'SSP1 (R45)','SSP3 (R60)','SSP4 (R60)','SSP5 (R85)',...
               'SSP1constLU (R45constCO2)','SSP3constLU (R60constCO2)','SSP4constLU (R60constCO2)','SSP5constLU (R85constCO2)'} ;
    runDir_base = addslashifneeded('/Users/Shared/PLUM/trunk_runs') ;
    runDirs = {[runDir_base 'PLUM2LPJGblIrr_45_s1_2011-2100_constSSP1LU/output-2017-12-10-232848'], ...
               [runDir_base 'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP3LU/output-2017-12-14-125023'], ...
               [runDir_base 'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP4LU/output-2017-12-15-015520'], ...
               [runDir_base 'PLUM2LPJGblIrr_85_s1_2011-2100_constSSP5LU/output-2017-12-15-011713'], ...
               [runDir_base 'PLUM2LPJGblIrr_45_s1_2011-2100_constSSP1LU_constCO2/output-2017-12-11-115548'], ...
               [runDir_base 'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP3LU_constCO2/output-2017-12-11-180709'], ...
               [runDir_base 'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP4LU_constCO2/output-2017-12-11-203339'], ...
               [runDir_base 'PLUM2LPJGblIrr_85_s1_2011-2100_constSSP5LU_constCO2/output-2017-12-11-205323'] ...
               } ;
    yearList_future = 2011:2100 ;
elseif strcmp(thisVer,'20180424_agmip7_SSP_vs_constLU')
        runList = {'SSP1 (R45)','SSP4 (R60)','SSP5 (R85)',...
               'SSP1constLU (R45)','SSP4constLU (R60)','SSP5constLU (R85)'} ;
    runDir_base = addslashifneeded('/Users/Shared/PLUM/trunk_runs') ;
    runDirs = {[runDir_base 'PLUM2LPJG_SSP1_RCP45_v3s1/output-2018-04-23-145614'], ...
               [runDir_base 'PLUM2LPJG_SSP4_RCP60_v3s1/output-2018-04-23-145614'], ...
               [runDir_base 'PLUM2LPJG_SSP5_RCP85_v3s1/output-2018-04-23-145614'], ...
               [runDir_base 'PLUM2LPJG_SSP1_RCP45_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-134335'], ...
               [runDir_base 'PLUM2LPJG_SSP4_RCP60_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-134415'], ...
               [runDir_base 'PLUM2LPJG_SSP5_RCP85_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-132213'], ...
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = addslashifneeded([runDir_base 'LPJGPLUM_1850-2010_PLUM6xtra/matlab_merge_20180424050342']) ;
    yearList_baseline = 1850:2010 ;
elseif strcmp(thisVer,'20180606_SSP_vs_constPLUM2011')
        runList = {'SSP1 (R45)','SSP3 (R60)','SSP4 (R60)','SSP5 (R85)',...
               'SSP1constPLUM2011 (R45)','SSP3constPLUM2011 (R60)','SSP4constPLUM2011 (R60)','SSP5constPLUM2011 (R85)'} ;
    runDir_base = addslashifneeded('/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs') ;
    runDirs = {find_PLUM2LPJG_run('PLUM2LPJG_SSP1_RCP45_v4s1_v20180426/output-2018-05-17-205051'), ...
               find_PLUM2LPJG_run('PLUM2LPJG_SSP3_RCP60_v4s1_v20180426/output-2018-05-17-213445'), ...
               find_PLUM2LPJG_run('PLUM2LPJG_SSP4_RCP60_v4s1_v20180426/output-2018-05-18-035126'), ...
               find_PLUM2LPJG_run('PLUM2LPJG_SSP5_RCP85_v4s1_v20180426/output-2018-05-18-055638'), ...
               find_PLUM2LPJG_run('PLUM2LPJG_SSP1_RCP45_v4s1_v20180426_asPLUMout2011/output-2018-05-22-003904'), ...
               find_PLUM2LPJG_run('PLUM2LPJG_SSP3_RCP60_v4s1_v20180426_asPLUMout2011/output-2018-06-02-140439'), ...
               find_PLUM2LPJG_run('PLUM2LPJG_SSP4_RCP60_v4s1_v20180426_asPLUMout2011/output-2018-06-02-140613'), ...
               find_PLUM2LPJG_run('PLUM2LPJG_SSP5_RCP85_v4s1_v20180426_asPLUMout2011/output-2018-05-21-022133'), ...
               } ;
    yearList_future = 2011:2100 ;
    baselineDir = addslashifneeded(find_PLUM2LPJG_run('LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180605160616')) ;
    yearList_baseline = 1850:2010 ;

else
    error(['thisVer (' thisVer ') not recognized!'])
end

% Output directory
outDir_ts = addslashifneeded(...
    ['/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/'...
     'figures_' thisVer '/diffTS']) ;

if ~exist(outDir_ts,'dir')
    mkdir(outDir_ts)
end

% Check runs
Nruns = length(runList) ;
if Nruns ~= length(runDirs)
    error('Length mismatch between runList and runDirs!')
end
for r = 1:Nruns
    if ~exist(runDirs{r},'dir')
        error(['runDir(' num2str(r) ') not found!'])
    end
    runDirs{r} = addslashifneeded(runDirs{r}) ;
end ; clear r

stdLegend = ['Baseline',runList] ;

addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work/')


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
CFTnames = strrep({bl_ts_fields{getCellsWithString(bl_ts_fields,'croparea_')}}','croparea_ts_','') ;
if any(contains(CFTnames,'TeWW'))
    is_lpjg = cellfun(@length,CFTnames)==4 & ~strcmp(CFTnames,'Rice') ;
    CFTnames_lpjg = {CFTnames{is_lpjg}}' ;
    CFTnames_plum = sort({CFTnames{~is_lpjg}}') ;
elseif any(contains(CFTnames,'_lpjg'))
    CFTnames_lpjg = {CFTnames{not(cellfun(@isempty,strfind(CFTnames,'_lpjg')))}}' ;
    CFTnames_plum = {CFTnames{not(cellfun(@isempty,strfind(CFTnames,'_plum')))}}' ;
    if length(CFTnames_lpjg)+length(CFTnames_plum) ~= length(CFTnames)
        error('length(CFTnames_lpjg)+length(CFTnames_plum) ~= length(CFTnames))')
    end
else
    tmp = {'CerealsC3','CerealsC3s','CerealsC3w','CerealsC4','Miscanthus','Oilcrops','Pulses','Rice','StarchyRoots'} ;
    error('This isn''t going to work out well. Have you combined LPJG and PLUM crops?')
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



%% Import

% Future runs
for r = 1:Nruns
    ts_tmp = load([runDirs{r} 'timeseries.mat']) ;
    lastdec_tmp = load([runDirs{r} 'last_decade.mat']) ;
    if r == 1
        Nyears = length(ts_tmp.LUarea_ts_bare) ;
        if Nyears ~= length(yearList_future)
            error('Nyears ~= length(yearList_future)')
        end
%         ts_aaet_yr = nan(Nyears,Nruns) ;
%         ts_aevap_yr = nan(Nyears,Nruns) ;
%         ts_aiso_yr = nan(Nyears,Nruns) ;
%         ts_amon_yr = nan(Nyears,Nruns) ;
%         ts_cpool_VegC_yr = nan(Nyears,Nruns) ;
%         ts_cpool_LitterSoilC_yr = nan(Nyears,Nruns) ;
%         ts_cpool_HarvSlowC_yr = nan(Nyears,Nruns) ;
%         ts_cpool_total_yr = nan(Nyears,Nruns) ;
        ts_nflux_flux_yr = nan(Nyears,Nruns) ;
        ts_nflux_leach_yr = nan(Nyears,Nruns) ;
        ts_nflux_harvest_yr = nan(Nyears,Nruns) ;
        ts_nflux_LUch_yr = nan(Nyears,Nruns) ;
        ts_nfert_yr = nan(Nyears,Nruns) ;
%         ts_totRunoff_yr = nan(Nyears,Nruns) ;
%         ts_albedo1_yr = nan(Nyears,Nruns) ;
%         ts_albedo7_yr = nan(Nyears,Nruns) ;
        ts_LUarea_ntrl_yr = nan(Nyears,Nruns) ;
        ts_LUarea_bare_yr = nan(Nyears,Nruns) ;
        ts_LUarea_crop_yr = nan(Nyears,Nruns) ;
        ts_LUarea_past_yr = nan(Nyears,Nruns) ;
        ts_croparea_TeWW_yr = nan(Nyears,Nruns) ;
        ts_croparea_TeSW_yr = nan(Nyears,Nruns) ;
        ts_croparea_TeCo_yr = nan(Nyears,Nruns) ;
        ts_croparea_TrRi_yr = nan(Nyears,Nruns) ;
        ts_irrig_yr = nan(Nyears,Nruns) ;
        ts_gsirrig_TeWW_yr = nan(Nyears,Nruns) ;
        ts_gsirrig_TeSW_yr = nan(Nyears,Nruns) ;
        ts_gsirrig_TeCo_yr = nan(Nyears,Nruns) ;
        ts_gsirrig_TrRi_yr = nan(Nyears,Nruns) ;
        ts_yield_TeWW_yr = nan(Nyears,Nruns) ;
        ts_yield_TeSW_yr = nan(Nyears,Nruns) ;
        ts_yield_TeCo_yr = nan(Nyears,Nruns) ;
        ts_yield_TrRi_yr = nan(Nyears,Nruns) ;
%         maps_all.list_to_map = lastdec_tmp.cpool.list_to_map ;
%         maps_all.yearList = lastdec_tmp.cpool.yearList ;
%         maps_cpool.maps_YXvyr = nan([size(lastdec_tmp.cpool.maps_YXvy) Nruns]) ;
%         maps_aaet.maps_YXvyr = nan([size(lastdec_tmp.aaet.maps_YXvy) Nruns]) ;
% %         maps_aevap.maps_YXvyr = nan([size(lastdec_tmp.aevap.maps_YXvy) Nruns]) ;
%         maps_tot_runoff.maps_YXvyr = nan([size(lastdec_tmp.tot_runoff.maps_YXvy) Nruns]) ;
%         maps_albedo.maps_YXvyr = nan([size(lastdec_tmp.albedo.maps_YXvy) Nruns]) ;
%         maps_aiso.maps_YXvyr = nan([size(lastdec_tmp.aiso.maps_YXvy) Nruns]) ;
%         maps_amon.maps_YXvyr = nan([size(lastdec_tmp.amon.maps_YXvy) Nruns]) ;
%         maps_nflux.maps_YXvyr = nan([size(lastdec_tmp.nflux.maps_YXvy) Nruns]) ;
%         maps_LU.maps_YXvyr = nan([size(lastdec_tmp.LU.maps_YXvy) Nruns]) ;
%         maps_cropfracs.maps_YXvyr = nan([size(lastdec_tmp.cropfracs.maps_YXvy) Nruns]) ;
%         maps_irrig.maps_YXvyr = nan([size(lastdec_tmp.irrig.maps_YXvy) Nruns]) ;
%         maps_gsirrig.maps_YXvyr = nan([size(lastdec_tmp.gsirrig.maps_YXvy) Nruns]) ;
%         maps_cpool.varNames = lastdec_tmp.cpool.varNames ;
%         maps_aaet.varNames = lastdec_tmp.aaet.varNames ;
% %         maps_aevap.varNames = lastdec_tmp.aevap.varNames ;
%         maps_tot_runoff.varNames = lastdec_tmp.tot_runoff.varNames ;
%         maps_albedo.varNames = lastdec_tmp.albedo.varNames ;
%         maps_aiso.varNames = lastdec_tmp.aiso.varNames ;
%         maps_amon.varNames = lastdec_tmp.amon.varNames ;
%         maps_nflux.varNames = lastdec_tmp.nflux.varNames ;
%         maps_LU.varNames = lastdec_tmp.LU.varNames ;
%         maps_cropfracs.varNames = lastdec_tmp.cropfracs.varNames ;
%         maps_irrig.varNames = lastdec_tmp.irrig.varNames ;
%         maps_gsirrig.varNames = lastdec_tmp.gsirrig.varNames ;
%         maps_yield.varNames = lastdec_tmp.yield.varNames ;
    end
%     ts_aaet_yr(:,r) = ts_tmp.aaet_ts ;
%     ts_aevap_yr(:,r) = ts_tmp.aevap_ts ;
%     ts_aiso_yr(:,r) = ts_tmp.aiso_ts ;
%     ts_amon_yr(:,r) = ts_tmp.amon_ts ;
%     ts_cpool_VegC_yr(:,r) = ts_tmp.cpool_ts_VegC ;
%     ts_cpool_LitterSoilC_yr(:,r) = ts_tmp.cpool_ts_LitterSoilC ;
%     ts_cpool_HarvSlowC_yr(:,r) = ts_tmp.cpool_ts_HarvSlowC ;
%     ts_cpool_total_yr(:,r) = ts_tmp.cpool_ts_Total ;
    ts_nflux_flux_yr(:,r) = ts_tmp.nflux_ts_flux ;
    ts_nflux_leach_yr(:,r) = ts_tmp.nflux_ts_leach ;
    ts_nflux_harvest_yr(:,r) = ts_tmp.nflux_ts_harvest ;
    ts_nflux_LUch_yr(:,r) = ts_tmp.nflux_ts_LU_ch ;
    ts_nfert_yr(:,r) = ts_tmp.nflux_ts_fert ;
%     ts_totRunoff_yr(:,r) = ts_tmp.tot_runoff_ts ;
%     ts_albedo1_yr(:,r) = ts_tmp.albedo_ts.jan ;
%     ts_albedo7_yr(:,r) = ts_tmp.albedo_ts.jul ;
    ts_LUarea_ntrl_yr(:,r) = ts_tmp.LUarea_ts_ntrl ;
    ts_LUarea_bare_yr(:,r) = ts_tmp.LUarea_ts_bare ;
    ts_LUarea_crop_yr(:,r) = ts_tmp.LUarea_ts_crop ;
    ts_LUarea_past_yr(:,r) = ts_tmp.LUarea_ts_past ;
    ts_croparea_TeWW_yr(:,r) = ts_tmp.croparea_ts_TeWW_lpjg ;
    ts_croparea_TeSW_yr(:,r) = ts_tmp.croparea_ts_TeSW_lpjg ;
    ts_croparea_TeCo_yr(:,r) = ts_tmp.croparea_ts_TeCo_lpjg ;
    ts_croparea_TrRi_yr(:,r) = ts_tmp.croparea_ts_TrRi_lpjg ;
    if isfield(ts_tmp,'croparea_ts_TeWWi')
        ts_croparea_TeWW_yr(:,r) = ts_croparea_TeWW_yr(:,r) + ts_tmp.croparea_ts_TeWWi_lpjg ;
        ts_croparea_TeSW_yr(:,r) = ts_croparea_TeSW_yr(:,r) + ts_tmp.croparea_ts_TeSWi_lpjg ;
        ts_croparea_TeCo_yr(:,r) = ts_croparea_TeCo_yr(:,r) + ts_tmp.croparea_ts_TeCoi_lpjg ;
        ts_croparea_TrRi_yr(:,r) = ts_croparea_TrRi_yr(:,r) + ts_tmp.croparea_ts_TrRii_lpjg ;
    end
    ts_irrig_yr(:,r) = ts_tmp.irrig_ts ;
%     ts_gsirrig_TeWW_yr(:,r) = ts_tmp.gsirrig_ts_TeWW ;
%     ts_gsirrig_TeSW_yr(:,r) = ts_tmp.gsirrig_ts_TeSW ;
%     ts_gsirrig_TeCo_yr(:,r) = ts_tmp.gsirrig_ts_TeCo ;
%     ts_gsirrig_TrRi_yr(:,r) = ts_tmp.gsirrig_ts_TrRi ;
%     if isfield(ts_tmp,'gsirrig_ts_TeWWi')
%         ts_gsirrig_TeWW_yr(:,r) = ts_gsirrig_TeWW_yr(:,r) + ts_tmp.gsirrig_ts_TeWWi ;
%         ts_gsirrig_TeSW_yr(:,r) = ts_gsirrig_TeSW_yr(:,r) + ts_tmp.gsirrig_ts_TeSWi ;
%         ts_gsirrig_TeCo_yr(:,r) = ts_gsirrig_TeCo_yr(:,r) + ts_tmp.gsirrig_ts_TeCoi ;
%         ts_gsirrig_TrRi_yr(:,r) = ts_gsirrig_TrRi_yr(:,r) + ts_tmp.gsirrig_ts_TrRii ;
%     end
    ts_yield_TeWW_yr(:,r) = ts_tmp.cropprod_ts_TeWW_lpjg ;
    ts_yield_TeSW_yr(:,r) = ts_tmp.cropprod_ts_TeSW_lpjg ;
    ts_yield_TeCo_yr(:,r) = ts_tmp.cropprod_ts_TeCo_lpjg ;
    ts_yield_TrRi_yr(:,r) = ts_tmp.cropprod_ts_TrRi_lpjg ;
    if isfield(ts_tmp,'yield_ts_TeWWi')
        ts_yield_TeWW_yr(:,r) = ts_yield_TeWW_yr(:,r) + ts_tmp.cropprod_ts_TeWWi_lpjg ;
        ts_yield_TeSW_yr(:,r) = ts_yield_TeSW_yr(:,r) + ts_tmp.cropprod_ts_TeSWi_lpjg ;
        ts_yield_TeCo_yr(:,r) = ts_yield_TeCo_yr(:,r) + ts_tmp.cropprod_ts_TeCoi_lpjg ;
        ts_yield_TrRi_yr(:,r) = ts_yield_TrRi_yr(:,r) + ts_tmp.cropprod_ts_TrRii_lpjg ;
    end
    clear ts_tmp
%     maps_cpool.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.cpool.maps_YXvy ;
%     maps_aaet.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.aaet.maps_YXvy ;
%     maps_aevap.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.aevap.maps_YXvy ;
%     maps_tot_runoff.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.tot_runoff.maps_YXvy ;
%     maps_albedo.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.albedo.maps_YXvy ;
%     maps_aiso.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.aiso.maps_YXvy ;
%     maps_amon.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.amon.maps_YXvy ;
%     maps_nflux.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.nflux.maps_YXvy ;
%     maps_LU.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.LU.maps_YXvy ;
%     maps_cropfracs.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.cropfracs.maps_YXvy ;
%     maps_irrig.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.irrig.maps_YXvy ;
%     maps_gsirrig.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.gsirrig.maps_YXvy ;
%     maps_yield.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.yield.maps_YXvy ;
%     clear lastdec_tmp
end ; clear r


%% Calculate time series differences

indices1 = 1:Nruns/2 ;
indices2 = (Nruns/2+1):Nruns ;

tsNames = {...
'ts_LUarea_bare_yr', 'ts_LUarea_crop_yr', 'ts_LUarea_ntrl_yr', 'ts_LUarea_past_yr', ...
...'ts_aaet_yr', 'ts_aevap_yr', 'ts_totRunoff_yr', ...
...'ts_aiso_yr', 'ts_amon_yr', ...
...'ts_albedo1_yr', 'ts_albedo7_yr', ...
...'ts_cpool_HarvSlowC_yr', 'ts_cpool_LitterSoilC_yr', 'ts_cpool_VegC_yr', 'ts_cpool_total_yr', ...
'ts_croparea_TeCo_yr', 'ts_croparea_TeSW_yr', 'ts_croparea_TeWW_yr', 'ts_croparea_TrRi_yr', ...
...'ts_gsirrig_TeCo_yr', 'ts_gsirrig_TeSW_yr', 'ts_gsirrig_TeWW_yr', 'ts_gsirrig_TrRi_yr', ...
'ts_irrig_yr', 'ts_nfert_yr', ...
'ts_nflux_LUch_yr', 'ts_nflux_flux_yr', 'ts_nflux_harvest_yr', 'ts_nflux_leach_yr', ...
'ts_yield_TeCo_yr', 'ts_yield_TeSW_yr', 'ts_yield_TeWW_yr', 'ts_yield_TrRi_yr'} ;

for t = 1:length(tsNames)
    thisName = tsNames{t} ;
    eval([thisName 'D = ' thisName '(:,indices1) - ' thisName '(:,indices2) ;']) ;
end


%% Plot timeseries: yield

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.08] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 0 ;
rebase = false ;
%%%%%%%%%%%%%%%%%%%

figure('Position',figurePos,'Color','w') ;

for c = 1:length(CFTnames)
    
    thisCrop = strrep(CFTnames{c},'_lpjg','') ;
    
    eval(['ts_yield_thisCrop_bl = ts_cropprod_' thisCrop '_lpjg_bl*1e-3 ;']) ; % kg to tons
    eval(['ts_yield_thisCrop_yr = ts_yield_' thisCrop '_yr*1e-3 ;']) ; % kg to tons
    eval(['ts_croparea_thisCrop_bl = ts_croparea_' thisCrop '_lpjg_bl*1e6*1e2 ;']) ; % Mkm2 to ha
    eval(['ts_croparea_thisCrop_yr = ts_croparea_' thisCrop '_yr*1e6*1e2 ;']) ; % Mkm2 to ha
    
    ts_YIELD_thisCrop_bl = ts_yield_thisCrop_bl ./ ts_croparea_thisCrop_bl ;
    ts_YIELD_thisCrop_yr = ts_yield_thisCrop_yr ./ ts_croparea_thisCrop_yr ;
    
    [tmp_ts_YIELD_thisCrop_yr, title_suffix, file_suffix] = ...
        rebase_future2baseline(rebase, Nsmth, ts_YIELD_thisCrop_bl, ts_YIELD_thisCrop_yr, ignYrs, yearList_future) ;
    
    subplot_tight(2,2,c,spacing)
    plot(yearList_future,tmp_ts_YIELD_thisCrop_yr(:,indices1),'LineWidth',lineWidth) ;
    hold on
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList_future,tmp_ts_YIELD_thisCrop_yr(:,indices2),'--','LineWidth',lineWidth) ;
    hold off
    xlims = get(gca,'XLim') ;
    ylims = get(gca,'YLim') ;
    if min(ylims)<0 && max(ylims)>0
        hold on
        plot([xlims(1) xlims(2)],[0 0],'--k','LineWidth',lineWidth)
        hold off
    end
    set(gca,'FontSize',fontSize,'XLim',xlims,'YLim',ylims)
    xlabel('Year')
    ylabel('tons ha^{-1}')
    title(['Yield: ' thisCrop])
%     legend(runList,'Location','northwest')
end

if do_save
    export_fig([outDir_ts 'yields_compare.pdf'])
    close
end


%% Plot timeseries: cropprod

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 10 ;
rebase = true ;
%%%%%%%%%%%%%%%%%%%

ts_yield_thisCrop_bl = ts_cropprod_TeWW_lpjg_bl ;
ts_yield_thisCrop_yr = ts_yield_TeWW_yr ;

[tmp_ts_YIELD_thisCrop_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_yield_thisCrop_bl, ts_yield_thisCrop_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

plot(yearList_future,tmp_ts_YIELD_thisCrop_yr(:,indices1),'LineWidth',lineWidth) ;
hold on
set(gca,'ColorOrderIndex',1) ;
plot(yearList_future,tmp_ts_YIELD_thisCrop_yr(:,indices2),'--','LineWidth',lineWidth) ;
hold off
xlims = get(gca,'XLim') ;
ylims = get(gca,'YLim') ;
if min(ylims)<0 && max(ylims)>0
    hold on
    plot([xlims(1) xlims(2)],[0 0],'--k','LineWidth',lineWidth)
    hold off
end
set(gca,'FontSize',fontSize,'XLim',xlims,'YLim',ylims)
xlabel('Year')
ylabel('PgC')
title('Vegetation C')
legend(runList)


if do_save
    export_fig([outDir_ts 'cpools_compare.pdf'])
    close
end



%% Plot timeseries: C pools v2

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%


figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
plot(yearList_future,ts_cpool_VegC_yr(:,indices1),'LineWidth',lineWidth) ;
hold on
set(gca,'ColorOrderIndex',1) ;
plot(yearList_future,ts_cpool_VegC_yr(:,indices2),'--','LineWidth',lineWidth) ;
hold off
xlims = get(gca,'XLim') ;
ylims = get(gca,'YLim') ;
if min(ylims)<0 && max(ylims)>0
    hold on
    plot([xlims(1) xlims(2)],[0 0],'--k','LineWidth',lineWidth)
    hold off
end
set(gca,'FontSize',fontSize,'XLim',xlims,'YLim',ylims)
xlabel('Year')
ylabel('PgC')
title('Vegetation C')

subplot_tight(1,2,2,spacing)
plot(yearList_future,ts_cpool_total_yr(:,indices1),'LineWidth',lineWidth)
hold on
set(gca,'ColorOrderIndex',1) ;
plot(yearList_future,ts_cpool_total_yr(:,indices2),'--','LineWidth',lineWidth) ;
hold off
xlims = get(gca,'XLim') ;
ylims = get(gca,'YLim') ;
if min(ylims)<0 && max(ylims)>0
    hold on
    plot([xlims(1) xlims(2)],[0 0],'--k','LineWidth',lineWidth)
    hold off
end
set(gca,'FontSize',fontSize,'XLim',xlims,'YLim',ylims)
xlabel('Year')
ylabel('PgC')
title('Total C')

clear tmp_*

if do_save
    export_fig([outDir_ts 'cpools_compare.pdf'])
    close
end


%% Plot timeseries differences: C pools v2

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%


figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
plot(yearList_future,ts_cpool_VegC_yrD,'LineWidth',lineWidth) ;
xlims = get(gca,'XLim') ;
ylims = get(gca,'YLim') ;
if min(ylims)<0 && max(ylims)>0
    hold on
    plot([xlims(1) xlims(2)],[0 0],'--k','LineWidth',lineWidth)
    hold off
end
set(gca,'FontSize',fontSize,'XLim',xlims,'YLim',ylims)
xlabel('Year')
ylabel('PgC')
title('Vegetation C')

subplot_tight(1,2,2,spacing)
plot(yearList_future,ts_cpool_total_yrD,'LineWidth',lineWidth)
if min(ylims)<0 && max(ylims)>0
    hold on
    plot([xlims(1) xlims(2)],[0 0],'--k','LineWidth',lineWidth)
    hold off
end
set(gca,'FontSize',fontSize,'XLim',xlims,'YLim',ylims)
xlabel('Year')
ylabel('PgC')
title('Total C')

clear tmp_*

if do_save
    export_fig([outDir_ts 'cpools.pdf'])
    close
end


%% Plot timeseries differences: C pools v2, percentage

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%


figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
plot(yearList_future,100*(ts_cpool_VegC_yr(:,indices1)-ts_cpool_VegC_yr(:,indices2))./ts_cpool_VegC_yr(:,indices2),'LineWidth',lineWidth) ;
xlims = get(gca,'XLim') ;
ylims = get(gca,'YLim') ;
if min(ylims)<0 && max(ylims)>0
    hold on
    plot([xlims(1) xlims(2)],[0 0],'--k','LineWidth',lineWidth)
    hold off
end
set(gca,'FontSize',fontSize,'XLim',xlims,'YLim',ylims)
xlabel('Year')
ylabel('Pct. difference')
title('Vegetation C')

subplot_tight(1,2,2,spacing)
plot(yearList_future,100*(ts_cpool_total_yr(:,indices1)-ts_cpool_total_yr(:,indices2))./ts_cpool_total_yr(:,indices2),'LineWidth',lineWidth) ;
if min(ylims)<0 && max(ylims)>0
    hold on
    plot([xlims(1) xlims(2)],[0 0],'--k','LineWidth',lineWidth)
    hold off
end
set(gca,'FontSize',fontSize,'XLim',xlims,'YLim',ylims)
xlabel('Year')
ylabel('Pct. difference')
title('Total C')

clear tmp_*

if do_save
    export_fig([outDir_ts 'cpools_pctDiff.pdf'])
    close
end









