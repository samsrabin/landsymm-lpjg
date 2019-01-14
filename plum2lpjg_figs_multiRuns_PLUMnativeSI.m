%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLUM2LPJG figures: Multiple runs: %%%
%%% PLUM-style native %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % thisVer = '20180424agmip7' ;
% % thisVer = '20180424agmip7_asPLUMout2011-2015' ;
% thisVer = 'v4s1_v20180426' ;
% thisVer = 'v4s1_v20180426_asPLUMout2011' ;
% thisVer = 'v4s1_v20180426_asLUH2_2010' ;
% thisVer = 'v6s1_v20180703' ;
% thisVer = 'v10s1_v20180801' ;
% thisVer = 'v10s1_v20180801' ;
% thisVer = 'harm2' ;
% thisVer = 'harm2_constLU' ;
% thisVer = 'harm2.1' ;
% thisVer = 'harm2.1_constClimCO2' ;
% thisVer = 'harm2.1_constLU' ;
% thisVer = 'harm2.1_S1R4.5_attr' ;
% thisVer = 'harm2.1_S3R6.0_attr' ;
% thisVer = 'harm2.1_S4R6.0_attr' ;
% thisVer = 'harm2.1_S5R8.5_attr' ;
% thisVer = 'harm2.2' ;
% thisVer = 'harm2.3' ;
% thisVer = 'harm2.3_constClimCO2' ;
% thisVer = 'harm2.3_constLU' ;
% thisVer = 'harm2.3_S1R4.5_attr' ;
% thisVer = 'harm2.3_S3R6.0_attr' ;
% thisVer = 'harm2.3_S4R6.0_attr' ;
% thisVer = 'harm2.3_S5R8.5_attr' ;
thisVer = 'harm2.3_constClim' ;

unhCropFrac = 0 ; % Set to zero for previous behavior. v10 = 0.177

% ignored_crops = {'CC3G','CC4G'} ;
% ignored_crops = {'CC3G','CC4G','Miscanthus'} ;
ignored_crops = {'CC3G','CC4G','ExtraCrop'} ;

do_save = true ;
rebase = false ;
pngres = 150 ;

do_caps = -1 ;

       
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
    skip3rdColor = false ;
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
    skip3rdColor = false ;
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
    skip3rdColor = false ;
elseif strcmp(thisVer,'v3s1_v20180426_asPLUM2011')
%     runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    runList = {'S1/R4.5','S4/R6.0','S5/R8.5'} ;
    runDirs = {
    'PLUM2LPJG_SSP1_RCP45_v3s1_constLUmgmt_asPLUMout2011/output-2018-05-03-060629' ;
%     'PLUM2LPJG_SSP3_RCP60_v3s1_constLUmgmt_asPLUMout2011/output-2018-05-03-053433' ;
    'PLUM2LPJG_SSP4_RCP60_v3s1_constLUmgmt_asPLUMout2011/output-2018-05-03-074359' ;
    'PLUM2LPJG_SSP5_RCP85_v3s1_constLUmgmt_asPLUMout2011/output-2018-05-03-080058' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180502104840' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'v4s1_v20180426')
    runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    runDirs = {
    'PLUM2LPJG_SSP1_RCP45_v4s1_v20180426/output-2018-05-17-205051' ;
    'PLUM2LPJG_SSP3_RCP60_v4s1_v20180426/output-2018-05-17-213445' ;
    'PLUM2LPJG_SSP4_RCP60_v4s1_v20180426/output-2018-05-18-035126' ;
%     'PLUM2LPJG_SSP5_RCP85_v4s1_v20180426/output-2018-05-18-055638' ;
    'PLUM2LPJG_SSP5_RCP85_v4s1_v20180426/output-2018-06-12-021143' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180605160616' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'v4s1_v20180426_asPLUMout2011')
    runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    runDirs = {
    'PLUM2LPJG_SSP1_RCP45_v4s1_v20180426_asPLUMout2011/output-2018-05-22-003904' ;
    'PLUM2LPJG_SSP3_RCP60_v4s1_v20180426_asPLUMout2011/output-2018-06-02-140439' ;
    'PLUM2LPJG_SSP4_RCP60_v4s1_v20180426_asPLUMout2011/output-2018-06-02-140613' ;
%     'PLUM2LPJG_SSP5_RCP85_v4s1_v20180426_asPLUMout2011/output-2018-05-21-022133' ;
    'PLUM2LPJG_SSP5_RCP85_v4s1_v20180426_asPLUMout2011/output-2018-06-12-021136' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180605160616' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'v4s1_v20180426_asLUH2_2010')
    runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    % Intentionally repeated SSP3, because only thing that matters in asLUH2
    % is RCP
    runDirs = {
    'PLUM2LPJG_SSP1_RCP45_v3s1_constLUmgmt_asLUH2_2010/output-2018-06-13-213037' ;
    'PLUM2LPJG_SSP3_RCP60_v3s1_constLUmgmt_asLUH2_2010/output-2018-06-14-002019' ;
    'PLUM2LPJG_SSP3_RCP60_v3s1_constLUmgmt_asLUH2_2010/output-2018-06-14-002019' ;
    'PLUM2LPJG_SSP5_RCP85_v3s1_constLUmgmt_asLUH2_2010/output-2018-06-14-002646' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180605160616' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'v6s1_v20180703')
    runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    runDirs = {
    'PLUM2LPJG_SSP1_RCP45_v6s1_v20180426/output-2018-07-02-151649' ;
    'PLUM2LPJG_SSP3_RCP60_v6s1_v20180426/output-2018-07-02-170609' ;
    'PLUM2LPJG_SSP4_RCP60_v6s1_v20180426/output-2018-07-03-020625' ;
    'PLUM2LPJG_SSP5_RCP85_v6s1_v20180426/output-2018-07-03-020605' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180605160616' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'v10s1_v20180801')
    runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    runDirs = {
        'PLUM2LPJG_SSP1_RCP45_v10s1_v20180730/output-2018-08-01-050407' ;
        'PLUM2LPJG_SSP3_RCP60_v10s1_v20180730/output-2018-08-01-135841' ;
        'PLUM2LPJG_SSP4_RCP60_v10s1_v20180730/output-2018-08-01-163815' ;
        'PLUM2LPJG_SSP5_RCP85_v10s1_v20180730/output-2018-08-01-164201' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_PLUM6xtraMisc_cg/matlab_merge_20180801105737' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2')
    runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP1_RCP45/output-2018-10-27-195714' ;
        'LPJGPLUM_2011-2100_harm2_SSP3_RCP60/output-2018-10-29-133454' ;
        'LPJGPLUM_2011-2100_harm2_SSP4_RCP60/output-2018-10-29-133454' ;
        'LPJGPLUM_2011-2100_harm2_SSP5_RCP85/output-2018-10-29-135226' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6/output-2018-10-27-073916' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2_constLU')
    runList = {'RCP4.5','RCP6.0','RCP8.5'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_constLU_RCP45/output-2018-10-30-124833' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP60/output-2018-10-30-220353' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP85/output-2018-10-30-214312' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6/output-2018-10-27-073916' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = true ;
elseif strcmp(thisVer,'harm2.1')
    runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP1_RCP45/output-2018-11-05-085615' ;
        'LPJGPLUM_2011-2100_harm2_SSP3_RCP60/output-2018-11-05-071233' ;
        'LPJGPLUM_2011-2100_harm2_SSP4_RCP60/output-2018-11-05-071814' ;
        'LPJGPLUM_2011-2100_harm2_SSP5_RCP85/output-2018-11-05-110842' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6/output-2018-11-03-234931' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.1_constClimCO2')
    runList = {'SSP1','SSP3','SSP4','SSP5'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP1_constClimCO2/output-2018-11-03-233941' ;
        'LPJGPLUM_2011-2100_harm2_SSP3_constClimCO2/output-2018-11-03-221433' ;
        'LPJGPLUM_2011-2100_harm2_SSP4_constClimCO2/output-2018-11-03-221734' ;
        'LPJGPLUM_2011-2100_harm2_SSP5_constClimCO2/output-2018-11-03-233913' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6/output-2018-11-03-234931' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.1_constLU')
    runList = {'RCP4.5','RCP6.0','RCP8.5'} ;
    runColNames = {'RCP45','RCP60','RCP85'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_constLU_RCP45/output-2018-11-05-003042' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP60/output-2018-11-05-003344' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP85/output-2018-11-05-000411' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6/output-2018-11-03-234931' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = true ;
elseif strcmp(thisVer,'harm2.1_S1R4.5_attr')
    runList = {'S1/R4.5','constLU','constClimCO2'} ;
    runColNames = {'Full','constLU','constClimCO2'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP1_RCP45/output-2018-11-05-085615' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP45/output-2018-11-05-003042' ;
        'LPJGPLUM_2011-2100_harm2_SSP1_constClimCO2/output-2018-11-03-233941' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6/output-2018-11-03-234931' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.1_S3R6.0_attr')
    runList = {'S3/R6.0','constLU','constClimCO2'} ;
    runColNames = {'Full','constLU','constClimCO2'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP3_RCP60/output-2018-11-05-071233' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP60/output-2018-11-05-003344' ;
        'LPJGPLUM_2011-2100_harm2_SSP3_constClimCO2/output-2018-11-03-221433' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6/output-2018-11-03-234931' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.1_S4R6.0_attr')
    runList = {'S4/R6.0','constLU','constClimCO2'} ;
    runColNames = {'Full','constLU','constClimCO2'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP4_RCP60/output-2018-11-05-071814' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP60/output-2018-11-05-003344' ;
        'LPJGPLUM_2011-2100_harm2_SSP4_constClimCO2/output-2018-11-03-221734' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6/output-2018-11-03-234931' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.1_S5R8.5_attr')
    runList = {'S5/R8.5','constLU','constClimCO2'} ;
    runColNames = {'Full','constLU','constClimCO2'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP5_RCP85/output-2018-11-05-110842' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP85/output-2018-11-05-000411' ;
        'LPJGPLUM_2011-2100_harm2_SSP5_constClimCO2/output-2018-11-03-233913' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6/output-2018-11-03-234931' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.2')
    runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP1_RCP45/output-2018-12-11-000445' ;
        'LPJGPLUM_2011-2100_harm2_SSP3_RCP60/output-2018-12-10-221610' ;
        'LPJGPLUM_2011-2100_harm2_SSP4_RCP60/output-2018-12-10-221802' ;
        'LPJGPLUM_2011-2100_harm2_SSP5_RCP85/output-2018-12-10-235151' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p3/output-2018-12-01-082125' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.3')
    runList = {'S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP1_RCP45/output-2018-12-11-000445' ;
        'LPJGPLUM_2011-2100_harm2_SSP3_RCP60/output-2018-12-10-221610' ;
        'LPJGPLUM_2011-2100_harm2_SSP4_RCP60/output-2018-12-10-221802' ;
        'LPJGPLUM_2011-2100_harm2_SSP5_RCP85/output-2018-12-10-235151' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p3/output-2018-12-09-071305' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.3_constClimCO2')
    runList = {'SSP1','SSP3','SSP4','SSP5'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP1_constClimCO2/output-2018-12-10-234221' ;
        'LPJGPLUM_2011-2100_harm2_SSP3_constClimCO2/output-2018-12-10-220740' ;
        'LPJGPLUM_2011-2100_harm2_SSP4_constClimCO2/output-2018-12-10-220855' ;
        'LPJGPLUM_2011-2100_harm2_SSP5_constClimCO2/output-2018-12-10-233802' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p3/output-2018-12-09-071305' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.3_constLU')
    runList = {'RCP4.5','RCP6.0','RCP8.5'} ;
    runColNames = {'RCP45','RCP60','RCP85'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_constLU_RCP45/output-2018-12-11-032126' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP60/output-2018-12-11-031352' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP85/output-2018-12-11-025144' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p3/output-2018-12-09-071305' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = true ;
elseif strcmp(thisVer,'harm2.3_S1R4.5_attr')
    runList = {'S1/R4.5','constLU','constClimCO2'} ;
    runColNames = {'Full','constLU','constClimCO2'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP1_RCP45/output-2018-12-11-000445' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP45/output-2018-12-11-032126' ;
        'LPJGPLUM_2011-2100_harm2_SSP1_constClimCO2/output-2018-12-10-234221' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p3/output-2018-12-09-071305' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.3_S3R6.0_attr')
    runList = {'S3/R6.0','constLU','constClimCO2'} ;
    runColNames = {'Full','constLU','constClimCO2'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP3_RCP60/output-2018-12-10-221610' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP60/output-2018-12-11-031352' ;
        'LPJGPLUM_2011-2100_harm2_SSP3_constClimCO2/output-2018-12-10-220740' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p3/output-2018-12-09-071305' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.3_S4R6.0_attr')
    runList = {'S4/R6.0','constLU','constClimCO2'} ;
    runColNames = {'Full','constLU','constClimCO2'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP4_RCP60/output-2018-12-10-221802' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP60/output-2018-12-11-031352' ;
        'LPJGPLUM_2011-2100_harm2_SSP4_constClimCO2/output-2018-12-10-220855' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p3/output-2018-12-09-071305' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.3_S5R8.5_attr')
    runList = {'S5/R8.5','constLU','constClimCO2'} ;
    runColNames = {'Full','constLU','constClimCO2'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP5_RCP85/output-2018-12-10-235151' ;
        'LPJGPLUM_2011-2100_harm2_constLU_RCP85/output-2018-12-11-025144' ;
        'LPJGPLUM_2011-2100_harm2_SSP5_constClimCO2/output-2018-12-10-233802' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p3/output-2018-12-09-071305' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm2.3_constClim')
    runList = {'S1/R4.5co2','S3/R6.0co2','S4/R6.0co2','S5/R8.5co2'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm2_SSP1_constClim/output-2018-12-30-225152' ;
        'LPJGPLUM_2011-2100_harm2_SSP3_constClim/output-2018-12-31-045830' ;
        'LPJGPLUM_2011-2100_harm2_SSP4_constClim/output-2018-12-31-045506' ;
        'LPJGPLUM_2011-2100_harm2_SSP5_constClim/output-2019-01-07-080321' ;
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p3/output-2018-12-09-071305' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
else
    error(['thisVer (' thisVer ') not recognized!'])
end

if ~exist('runColNames','var')
    runColNames = {'SSP1','SSP3','SSP4','SSP5'} ;
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
outDir_base = addslashifneeded(['/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/'...
                           'figures_' thisVer '_SI']) ;
outDir_maps = addslashifneeded([outDir_base 'maps']) ;
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

% Conversion factors
cf_kg2Mt = 1e-3*1e-6 ;
cf_t2kg = 1e3 ;   % For FAO only
cf_ha2m2 = 1e4 ; % For FAO only
cf_kgPm2_to_tonsPha = 1e-3*1e4 ;
cf_kg2Tg = 1e-9 ;
cf_kg2Pg = 1e-12 ;
cf_m3_to_km3 = (1e-3)^3 ;
cf_kcal2Ecal = 1e-15 ;

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
    || contains(thisVer,'harm2')
        fao = load('/Users/Shared/PLUM/crop_calib_data/fao/FAOdata_1961-2010_calibVer16_Production.mat') ;
    else
        error('thisVer not recognized while loading FAO data')
    end
end


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
    yr_ts_fields = fieldnames(ts_tmp) ;
    yr_ts_fields(endsWith(yr_ts_fields,ignored_crops)) = [] ;
    is_expYields = ~cellfun(@isempty,regexp(yr_ts_fields,'^cropprodExp_')) ;
    have_expYields = any(is_expYields) ;
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
            
            % If not present in future run, remove from analysis
            if ~(isfield(firstdec_tmp,thisVar_out) || isfield(lastdec_tmp,thisVar_out))
                warning([thisVar_out ' not found in future run ' num2str(r) '. Removing.']) ;
                clear(thisVar_in)
                continue
            end
            
            if contains(thisVar_in,'_d1')
                eval([thisVar_in '.maps_YXvyr = nan([size(firstdec_tmp.' thisVar_out '.maps_YXvy) Nruns],''single'') ;']) ;
            elseif contains(thisVar_in,'_d9')
                eval([thisVar_in '.maps_YXvyr = nan([size(lastdec_tmp.' thisVar_out '.maps_YXvy) Nruns],''single'') ;']) ;
            else
                error('How did this happen?')
            end
            clear thisVar*
        end ; clear v
        
        % Get vars_maps_bl again, now that you've cleared variables not
        % present in future run(s)
        tmp = whos('maps_*') ;
        vars_maps_bl = {tmp.name}' ;
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
    vars_maps_bl_toRemove = false(size(vars_maps_bl)) ;
    for v = 1:length(vars_maps_bl)
        thisVar_in = vars_maps_bl{v} ;
        if contains(thisVar_in,'maps_LU0')
            continue
        end
        thisVar_out = strrep(thisVar_in,'maps_','') ;
        if contains(thisVar_in,'d1')
            if ~isfield(firstdec_tmp,thisVar_out)
                warning([thisVar_in ' not found in future run ' num2str(r) '. Removing.']) ;
                vars_maps_bl_toRemove(v) = true ;
                eval(['clear ' thisVar_in]) ;
                continue
            end
            eval(['isequal_varNames = isequal(' thisVar_in '.varNames, firstdec_tmp.' thisVar_out '.varNames) ; ']) ;
            if ~isequal_varNames
                eval(['isequal_varNames_afterSort = isequal(sort(' thisVar_in '.varNames), sort(firstdec_tmp.' thisVar_out '.varNames)) ; ']) ;
                if ~isequal_varNames_afterSort
                    error(['~isequal_varNames (' thisVar_in '). Not fixable.'])
                end
                warning(['~isequal_varNames (' thisVar_in '). Fixing.'])
                eval(['[~,~,IB] = intersect(' thisVar_in '.varNames,firstdec_tmp.' thisVar_out '.varNames,''stable'') ;']) ;
                eval([thisVar_in '.maps_YXvyr(:,:,:,:,r) = firstdec_tmp.' thisVar_out '.maps_YXvy(:,:,IB,:) ;']) ;
            else
                eval([thisVar_in '.maps_YXvyr(:,:,:,:,r) = firstdec_tmp.' thisVar_out '.maps_YXvy ;']) ;
            end
        elseif contains(thisVar_in,'d9')
            if ~isfield(lastdec_tmp,thisVar_out)
                warning([thisVar_in ' not found in future run ' num2str(r) '. Removing.']) ;
                vars_maps_bl_toRemove(v) = true ;
                eval(['clear ' thisVar_in]) ;
                continue
            end
            eval(['isequal_varNames = isequal(' thisVar_in '.varNames, lastdec_tmp.' thisVar_out '.varNames) ; ']) ;
            if ~isequal_varNames
                eval(['isequal_varNames_afterSort = isequal(sort(' thisVar_in '.varNames), sort(lastdec_tmp.' thisVar_out '.varNames)) ; ']) ;
                if ~isequal_varNames_afterSort
                    error(['~isequal_varNames (' thisVar_in '). Not fixable.'])
                end
                warning(['~isequal_varNames (' thisVar_in '). Fixing.'])
                eval(['[~,~,IB] = intersect(' thisVar_in '.varNames,lastdec_tmp.' thisVar_out '.varNames,''stable'') ;']) ;
                eval([thisVar_in '.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.' thisVar_out '.maps_YXvy(:,:,IB,:) ;']) ;
            else
                eval([thisVar_in '.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.' thisVar_out '.maps_YXvy ;']) ;
            end
        else
            error('How did this happen?')
        end
    end
    
    vars_maps_bl(vars_maps_bl_toRemove) = [] ;
    clear vars_maps_bl_toRemove
    
    if isfield(firstdec_tmp,'expyield_d1')
        if r == 1
            maps_yieldExp_d1 = renameStructField(firstdec_tmp.expyield_d1,'maps_YXvy','maps_YXvyr') ;
        else
            maps_yieldExp_d1.maps_YXvyr(:,:,:,:,r) = firstdec_tmp.expyield_d1.maps_YXvy ;
        end
    end
    if isfield(lastdec_tmp,'expyield_d9')
        if r == 1
            maps_yieldExp_d9 = renameStructField(lastdec_tmp.expyield_d9,'maps_YXvy','maps_YXvyr') ;
        else
            maps_yieldExp_d9.maps_YXvyr(:,:,:,:,r) = lastdec_tmp.expyield_d9.maps_YXvy ;
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
tmp = gcel_area_YXqd(:,1:2:1440) + gcel_area_YXqd(:,2:2:1440) ;
gcel_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
gcel_area_YX(nanmask) = NaN ;
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
land_area_YX(nanmask) = NaN ;
%%% Convert to m2
land_area_YX = land_area_YX*1e6 ;
gcel_area_YX = gcel_area_YX*1e6 ;
clear tmp land_frac_YXqd land_area_YXqd

% area_YXBH: "Baseline" as last year of Historical run
% area_YXBFr: "Baseline" as first year of Future runs
% diff_YXrH: Difference from End-Historical to End-Future
% diff_YXrF: Difference from Begin-Future to End-Future
tmp_list_4 = {'ntrl','bare','crop','past'} ;
tmp_list_full = {'NATURAL','BARREN','CROPLAND','PASTURE'} ;
for i = 1:length(tmp_list_4)
    this_4 = tmp_list_4{i} ;
    this_full = tmp_list_full{i} ;
    eval([this_4 '_area_YXBH = gcel_area_YX .* maps_LU_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,this_full),end) ;']) ;
    eval([this_4 '_area_YXBFr = gcel_area_YX .* squeeze(maps_LU_d1.maps_YXvyr(:,:,strcmp(maps_LU_d1.varNames,this_full),1,:)) ;']) ;
    eval([this_4 '_area_YXr = repmat(gcel_area_YX,[1 1 Nruns]) .* squeeze(maps_LU_d9.maps_YXvyr(:,:,strcmp(maps_LU_d9.varNames,this_full),end,:)) ;']) ;
    eval([this_4 '_diff_YXrH = ' this_4 '_area_YXr - repmat(' this_4 '_area_YXBH,[1 1 Nruns]) ;']) ;
    eval([this_4 '_diff_YXrF = ' this_4 '_area_YXr - ' this_4 '_area_YXBFr ;']) ;
end

if exist('maps_LU0_d9','var')
    crop0_area_YXBH = gcel_area_YX .* maps_LU0_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'CROPLAND'),end) ;
    crop0_diff_YXrH = crop_area_YXr - repmat(crop0_area_YXBH,[1 1 Nruns]) ;
end
if exist('maps_LU0_d9','var')
    past0_area_YXBH = gcel_area_YX .* maps_LU0_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'PASTURE'),end) ;
    past0_diff_YXrH = past_area_YXr - repmat(past0_area_YXBH,[1 1 Nruns]) ;
end

agri_area_YXBH = crop_area_YXBH + past_area_YXBH ;
agri_area_YXBFr = crop_area_YXBFr + past_area_YXBFr ;
agri_area_YXr = crop_area_YXr + past_area_YXr ;
agri_diff_YXrH = crop_diff_YXrH + past_diff_YXrH ;
agri_diff_YXrF = crop_diff_YXrF + past_diff_YXrF ;

% Remove ignored crops
if ~isempty(ignored_crops)
    maps_cropfracs_d1.maps_YXvyB(:,:,contains(maps_cropfracs_d1.varNames,ignored_crops),:) = [] ;
    maps_cropfracs_d1.maps_YXvyr(:,:,contains(maps_cropfracs_d1.varNames,ignored_crops),:,:) = [] ;
    maps_cropfracs_d1.varNames(contains(maps_cropfracs_d1.varNames,ignored_crops)) = [] ;
    maps_cropfracs_d9.maps_YXvyB(:,:,contains(maps_cropfracs_d9.varNames,ignored_crops),:) = [] ;
    maps_cropfracs_d9.maps_YXvyr(:,:,contains(maps_cropfracs_d9.varNames,ignored_crops),:,:) = [] ;
    maps_cropfracs_d9.varNames(contains(maps_cropfracs_d9.varNames,ignored_crops)) = [] ;
    maps_yield_d1.maps_YXvyB(:,:,contains(maps_yield_d1.varNames,ignored_crops),:) = [] ;
    maps_yield_d1.maps_YXvyr(:,:,contains(maps_yield_d1.varNames,ignored_crops),:,:) = [] ;
    maps_yield_d1.varNames(contains(maps_yield_d1.varNames,ignored_crops)) = [] ;
    maps_yield_d9.maps_YXvyB(:,:,contains(maps_yield_d9.varNames,ignored_crops),:) = [] ;
    maps_yield_d9.maps_YXvyr(:,:,contains(maps_yield_d9.varNames,ignored_crops),:,:) = [] ;
    maps_yield_d9.varNames(contains(maps_yield_d9.varNames,ignored_crops)) = [] ;
    if exist('maps_yieldExp_d1','var')
        % (Baseline has no expected yield) maps_yieldExp_d1.maps_YXvyB(:,:,contains(maps_yieldExp_d1.varNames,ignored_crops),:) = [] ;
        maps_yieldExp_d1.maps_YXvyr(:,:,contains(maps_yieldExp_d1.varNames,ignored_crops),:,:) = [] ;
        maps_yieldExp_d1.varNames(contains(maps_yieldExp_d1.varNames,ignored_crops)) = [] ;
    end
    if exist('maps_yieldExp_d9','var')
        % (Baseline has no expected yield) maps_yieldExp_d9.maps_YXvyB(:,:,contains(maps_yieldExp_d9.varNames,ignored_crops),:) = [] ;
        maps_yieldExp_d9.maps_YXvyr(:,:,contains(maps_yieldExp_d9.varNames,ignored_crops),:,:) = [] ;
        maps_yieldExp_d9.varNames(contains(maps_yieldExp_d9.varNames,ignored_crops)) = [] ;
    end
    if exist('maps_Nfert_d1','var')
        maps_Nfert_d1.maps_YXvyB(:,:,contains(maps_Nfert_d1.varNames,ignored_crops),:) = [] ;
        maps_Nfert_d1.maps_YXvyr(:,:,contains(maps_Nfert_d1.varNames,ignored_crops),:,:) = [] ;
        maps_Nfert_d1.varNames(contains(maps_Nfert_d1.varNames,ignored_crops)) = [] ;
    end
    if exist('maps_Nfert_d9','var')
        maps_Nfert_d9.maps_YXvyB(:,:,contains(maps_Nfert_d9.varNames,ignored_crops),:) = [] ;
        maps_Nfert_d9.maps_YXvyr(:,:,contains(maps_Nfert_d9.varNames,ignored_crops),:,:) = [] ;
        maps_Nfert_d9.varNames(contains(maps_Nfert_d9.varNames,ignored_crops)) = [] ;
    end
end

if ~isequal(maps_cropfracs_d1.varNames,maps_cropfracs_d9.varNames)
    error('~isequal(maps_cropfracs_d1.varNames,maps_cropfracs_d9.varNames)')
else
    CFTnames_maps = maps_cropfracs_d1.varNames ;
end
maps_cropareas_d1 = maps_cropfracs_d1 ;
maps_cropareas_d1.maps_YXvyB = maps_cropfracs_d1.maps_YXvyB .* repmat(gcel_area_YX,[1 1 Ncrops size(maps_cropfracs_d1.maps_YXvyB,4)]) .* repmat(maps_LU_d1.maps_YXvyB(:,:,strcmp(maps_LU_d1.varNames,'CROPLAND'),:),[1 1 Ncrops 1]) ;
maps_cropareas_d1.maps_YXvyr = maps_cropfracs_d1.maps_YXvyr .* repmat(gcel_area_YX,[1 1 Ncrops size(maps_cropfracs_d1.maps_YXvyB,4) Nruns]) .* repmat(maps_LU_d1.maps_YXvyr(:,:,strcmp(maps_LU_d1.varNames,'CROPLAND'),:,:),[1 1 Ncrops 1 1]) ;
maps_cropareas_d9 = maps_cropfracs_d9 ;
maps_cropareas_d9.maps_YXvyB = maps_cropfracs_d9.maps_YXvyB .* repmat(gcel_area_YX,[1 1 Ncrops size(maps_cropfracs_d9.maps_YXvyB,4)]) .* repmat(maps_LU_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'CROPLAND'),:),[1 1 Ncrops 1]) ;
maps_cropareas_d9.maps_YXvyr = maps_cropfracs_d9.maps_YXvyr .* repmat(gcel_area_YX,[1 1 Ncrops size(maps_cropfracs_d9.maps_YXvyB,4) Nruns]) .* repmat(maps_LU_d9.maps_YXvyr(:,:,strcmp(maps_LU_d9.varNames,'CROPLAND'),:,:),[1 1 Ncrops 1 1]) ;

maps_cropareas_YXvBH =          maps_cropareas_d9.maps_YXvyB(:,:,:,end) ;
maps_cropareas_YXvBFr = squeeze(maps_cropareas_d1.maps_YXvyr(:,:,:,1,  :)) ;
maps_cropareas_YXvr   = squeeze(maps_cropareas_d9.maps_YXvyr(:,:,:,end,:)) ;
maps_cropareasDiffs_YXvrH = maps_cropareas_YXvr - repmat(maps_cropareas_YXvBH,[1 1 1 Nruns]) ;
maps_cropareasDiffs_YXvrF = maps_cropareas_YXvr - maps_cropareas_YXvBFr ;

disp('Done importing future.')


%% Calculate secondary variables

disp('Calculating secondary variables...')

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
ts_kcal2_bl = zeros(size(Nyears_bl,1)) ;
ts_kcal2_fao = zeros(size(Nyears_bl,1)) ;
ts_kcal2_yr = zeros(size(Nyears_bl,Nruns)) ;
for i = 1:length(tmp_name)
    thisCrop = strrep(strrep(strrep(strrep(tmp_name{i},'ts_cropprod_',''),'_bl',''),'_yr',''),'_fao','') ;
    if strcmp(thisCrop,'Miscanthus')
        continue
    end
    thisSuffix = strrep(tmp_name{i},['ts_cropprod_' thisCrop '_'],'') ;
    kcal_per_g = get_kcalDensity(thisCrop) ;
    kcal_per_kg = 1e3 * kcal_per_g ;
    eval(['ts_kcal_' thisSuffix ' = ts_kcal_' thisSuffix ' + kcal_per_kg * eval(tmp_name{i}) ;']) ;
    kcal2_per_g = get_kcalDensity2(thisCrop) ;
    kcal2_per_kg = 1e3 * kcal2_per_g ;
    eval(['ts_kcal2_' thisSuffix ' = ts_kcal2_' thisSuffix ' + kcal2_per_kg * eval(tmp_name{i}) ;']) ;
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
        kcal_per_g = get_kcalDensity(thisCrop) ;
        kcal_per_kg = 1e3 * kcal_per_g ;
        eval(['ts_kcalExp_' thisSuffix ' = ts_kcalExp_' thisSuffix ' + kcal_per_kg * (' tmp_name{i} '.* ts_croparea_' thisCrop '_yr) ;']) ;
        kcal2_per_g = get_kcalDensity2(thisCrop) ;
        kcal2_per_kg = 1e3 * kcal2_per_g ;
        eval(['ts_kcalExp2_' thisSuffix ' = ts_kcalExp2_' thisSuffix ' + kcal2_per_kg * (' tmp_name{i} '.* ts_croparea_' thisCrop '_yr) ;']) ;
    end ; clear i
end ; clear tmp

% Peak monthly runoff maps
if exist('maps_mon_runoff_d1','var')
    maps_pk_runoff_d1 = maps_mon_runoff_d1 ;
    maps_pk_runoff_d1.maps_YXvyB = max(maps_pk_runoff_d1.maps_YXvyB,[],3) ;
    maps_pk_runoff_d1.maps_YXvyr = max(maps_pk_runoff_d1.maps_YXvyr,[],3) ;
    maps_pk_runoff_d1.varNames = {'Max'} ;
    maps_pk_runoff_d9 = maps_mon_runoff_d9 ;
    maps_pk_runoff_d9.maps_YXvyB = max(maps_pk_runoff_d9.maps_YXvyB,[],3) ;
    maps_pk_runoff_d9.maps_YXvyr = max(maps_pk_runoff_d9.maps_YXvyr,[],3) ;
    maps_pk_runoff_d9.varNames = {'Max'} ;
end

disp('Done calculating secondary variables.')


%% Save outputs

if false
    disp('Saving outputs...')
    save(['/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/fig_script_outputs/' thisVer '.mat'], ...
        'agri_*', 'bare_*', 'bl_*', 'crop*', 'firstdec_tmp', 'lastdec_tmp', ...
        'land_area_YX', 'gcel_area_YX', 'maps_*', 'nanmask', 'ntrl_*', 'past*', ...
        'ts_*_bl','ts_*_yr')
    disp('Done.')
end


%% Table after Krause et al. (2017) Table 2

disp('Making table...')

years_endh = 2000:2009 ;
years_begf = 2011:2020 ;
years_endf = 2090:2099 ;

% Name, code, conversion factor, formatSpec mean, formatSpec SEM
rowInfo = {'Vegetation C (GtC)', 'cpool_VegC', cf_kg2Pg, '%d', '%d' ;
           'Soil and litter C (GtC)', 'cpool_LitterSoilC', cf_kg2Pg, '%d', '%d' ;
           'Product C (GtC)', 'cpool_HarvSlowC', cf_kg2Pg, '%0.1f', '%0.1f' ;
           'Total C (GtC)', 'cpool_Total', cf_kg2Pg, '%d', '%d' ;
           'January albedo', 'albedo1', 1, '%0.3f', '%0.3f' ;
           'July albedo', 'albedo7', 1, '%0.3f', '%0.3f' ;
           'Evapotranspiration (1000 km^3)', 'aevapaaet', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
           'Runoff (1000 km^3)', 'tot_runoff', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
           'Peak monthly runoff (1000 km^3)', 'pkrunoff', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
           'Crop production (Ecal)', 'kcal', cf_kcal2Ecal, '%0.1f', '%0.1f' ;
           'N loss (TgN)', 'nloss', cf_kg2Tg, '%0.1f', '%0.1f' ;
           'N loss: Gaseous (TgN)', 'nflux_flux', cf_kg2Tg, '%0.1f', '%0.1f' ;
           'N loss: Dissolved (TgN)', 'nflux_leach', cf_kg2Tg, '%0.1f', '%0.1f' ;
           'Isoprene emissions (TgC)', 'aiso', cf_kg2Tg, '%0.1f', '%0.1f' ;
           'Monoterpene emissions (TgC)', 'amon', cf_kg2Tg, '%0.1f', '%0.1f' ;
           } ;

Nvars = size(rowInfo,1) ;
mean_endh = nan(Nvars,1) ;
mean_begf = nan(Nvars,Nruns) ;
mean_endf = nan(Nvars,Nruns) ;
sem_endh = nan(Nvars,1) ;
sem_begf = nan(Nvars,Nruns) ;
sem_endf = nan(Nvars,Nruns) ;
string_endh = cell(Nvars,1) ;
string_begf = cell(Nvars,Nruns) ;
string_endf = cell(Nvars,Nruns) ;
for c = 1:Nvars
    
    % Get values
    thisVar = rowInfo{c,2} ;
    thisConv = rowInfo{c,3} ;
    mean_endh(c) = thisConv*eval(['mean(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
    sem_endh(c) = thisConv*eval(['std(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
    mean_begf(c,:) = thisConv*eval(['mean(ts_' thisVar '_yr(yearList_future>=min(years_begf) & yearList_future<=max(years_begf),:))']) ;
    if strcmp(thisVar,'nloss')
        x=1;
    end
    sem_begf(c,:) = thisConv*eval(['std(ts_' thisVar '_yr(yearList_future>=min(years_begf) & yearList_future<=max(years_begf),:))']) ;
    mean_endf(c,:) = thisConv*eval(['mean(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
    sem_endf(c,:) = thisConv*eval(['std(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
    
    % Turn into strings
    if strcmp(rowInfo{c,4},'%d')
        thisMean = round(mean_endh(c)) ;
    else
        thisMean = mean_endh(c) ;
    end
    if strcmp(rowInfo{c,4},'%d')
        thisSD = round(sem_endh(c)) ;
    else
        thisSD = sem_endh(c) ;
    end
    string_endh{c} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean thisSD]) ;
    for r = 1:Nruns
        if strcmp(rowInfo{c,4},'%d')
            thisMean_begf = round(mean_begf(c,r)) ;
            thisMean_endf = round(mean_endf(c,r)) ;
        else
            thisMean_begf = mean_begf(c,r) ;
            thisMean_endf = mean_endf(c,r) ;
        end
        if strcmp(rowInfo{c,4},'%d')
            thisSD_begf = round(sem_begf(c,r)) ;
            thisSD_endf = round(sem_endf(c,r)) ;
        else
            thisSD_begf = sem_begf(c,r) ;
            thisSD_endf = sem_endf(c,r) ;
        end
        string_begf{c,r} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean_begf thisSD_begf]) ;
        string_endf{c,r} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean_endf thisSD_endf]) ;
    end

end

table_out = table(collate_empties(rowInfo(:,1)),...
                  collate_empties(string_endh)) ;
for r = 1:Nruns
    table_out = [table_out collate_twocells(string_begf(:,r),string_endf(:,r))] ;
end
table_out.Properties.VariableNames = [{'Ecosystem_function','Baseline'} runColNames] ;

if do_save
    writetable(table_out,[outDir_base 'summary_table.xlsx'],'Sheet',1) ;
end
disp('Done making table.')


%% Map differences from beginning to end of future: Isoprene emissions 

% Options %%%%%%%%%
filename_base = [outDir_maps 'bvoc_iso'] ;
title_text = 'isoprene emissions' ;
sumvars = 'Total' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1e-3 ;   % mgC/m2 to gC/m2
units_map = 'gC m^{-2} yr ^{-1}' ;
conv_fact_total = 1e-6*cf_kg2Tg ;   % mgC to TgC
units_total = 'TgC yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_aiso_d1, maps_aiso_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: Monoterpene emissions

% Options %%%%%%%%%
filename_base = [outDir_maps 'bvoc_mon'] ;
title_text = 'monoterpene emissions' ;
sumvars = 'Total' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1e-3 ;   % mgC/m2 to gC/m2
units_map = 'gC m^{-2} yr ^{-1}' ;
conv_fact_total = 1e-6*cf_kg2Tg ;   % mgC to TgC
units_total = 'TgC yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_amon_d1, maps_amon_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: Total N loss

% Options %%%%%%%%%
filename_base = [outDir_maps 'nloss'] ;
title_text = 'total N loss' ;
sumvars = {'flux','leach'} ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % kgN/m2
units_map = 'kgN m^{-2} yr ^{-1}' ;
conv_fact_total = cf_kg2Tg ;   % kgN to TgN
units_total = 'TgN yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_nflux_d1, maps_nflux_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: Gaseous N loss

% Options %%%%%%%%%
filename_base = [outDir_maps 'nloss_gas'] ;
title_text = 'gaseous N loss' ;
sumvars = 'flux' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % kgN/m2
units_map = 'kgN m^{-2} yr ^{-1}' ;
conv_fact_total = cf_kg2Tg ;   % kgN to TgN
units_total = 'TgN yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_nflux_d1, maps_nflux_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: Dissolved N loss

% Options %%%%%%%%%
filename_base = [outDir_maps 'nloss_liq'] ;
title_text = 'dissolved N loss' ;
sumvars = 'leach' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % kgN/m2
units_map = 'kgN m^{-2} yr ^{-1}' ;
conv_fact_total = cf_kg2Tg ;   % kgN to TgN
units_total = 'TgN yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_nflux_d1, maps_nflux_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: January albedo

% warning('Skipping January albedo: pass per-year non-bare area as weighting!')

% Options %%%%%%%%%
filename_base = [outDir_maps 'albedo_jan'] ;
title_text = 'January albedo' ;
sumvars = 'January' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % unitless
units_map = '' ;
conv_fact_total = [] ;   % unitless
units_total = '' ;
this_land_area_map = land_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_albedo_d1, maps_albedo_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: July albedo

% warning('Skipping July albedo: pass per-year non-bare area as weighting!')

% Options %%%%%%%%%
filename_base = [outDir_maps 'albedo_jul'] ;
title_text = 'July albedo' ;
sumvars = 'July' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % unitless
units_map = '' ;
conv_fact_total = [] ;   % unitless
units_total = '' ;
this_land_area_map = land_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_albedo_d1, maps_albedo_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: Vegetation C

% Options %%%%%%%%%
filename_base = [outDir_maps 'cpool_veg'] ;
title_text = 'vegetation C' ;
sumvars = 'VegC' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1e-6*1e4 ;   % kgC/m2 to tonsC/ha
units_map = 'tons C ha^{-1}' ;
conv_fact_total = cf_kg2Pg ;   % kgC to PgC
units_total = 'PgC' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_cpool_d1, maps_cpool_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: Total C

% Options %%%%%%%%%
filename_base = [outDir_maps 'cpool_tot'] ;
title_text = 'total C' ;
sumvars = 'Total' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1e-6*1e4 ;   % kgC/m2 to tonsC/ha
units_map = 'tons C ha^{-1}' ;
conv_fact_total = cf_kg2Pg ;   % kgC to PgC
units_total = 'PgC' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_cpool_d1, maps_cpool_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: Evapotranspiration

% Options %%%%%%%%%
filename_base = [outDir_maps 'water_evapotransp'] ;
title_text = 'evapotranspiration' ;
sumvars = {'Evap','Transp'} ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % mm
units_map = 'mm yr^{-1}' ;
conv_fact_total = 1e-3*cf_m3_to_km3*1e-3 ;   % mm to 1000 km3
units_total = '1000 km^3 yr^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_awater_d1, maps_awater_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: Annual runoff

% Options %%%%%%%%%
filename_base = [outDir_maps 'water_runoff'] ;
title_text = 'annual runoff' ;
sumvars = 'Runoff' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % mm
units_map = 'mm yr^{-1}' ;
conv_fact_total = 1e-3*cf_m3_to_km3*1e-3 ;   % mm to 1000 km3
units_total = '1000 km^3 yr^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_awater_d1, maps_awater_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from beginning to end of future: Peak runoff

% Options %%%%%%%%%
filename_base = [outDir_maps 'water_runoff_peak'] ;
title_text = 'peak monthly runoff' ;
sumvars = 'Max' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % mm
units_map = 'mm mon^{-1}' ;
conv_fact_total = [] ;   % Left empty to avoid calculation of total
units_total = '' ;
this_land_area_map = [] ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs(do_save, maps_pk_runoff_d1, maps_pk_runoff_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


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
conv_fact_map = 1e-6 ;   % m2 to km2
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
units_map = 'km^2' ;
units_total = 'Mkm^2' ;
only1bl = true ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;

% Natural area
make_LUdiff_fig_v2(...
    ntrl_area_YXBH, ntrl_diff_YXrH, ...
    thisY1, thisYN, '"Natural"', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
    close
end

% Cropland area
make_LUdiff_fig_v2(...
    crop_area_YXBH, crop_diff_YXrH, ...
    thisY1, thisYN, 'Cropland', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
    close
end

% "Pasture" area
make_LUdiff_fig_v2(...
    past_area_YXBH, past_diff_YXrH, ...
    thisY1, thisYN, '"Pasture"', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
    close
end

% Agricultural area
make_LUdiff_fig_v2(...
    agri_area_YXBH, agri_diff_YXrH, ...
    thisY1, thisYN, 'Agricultural', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_agri.png'],['-r' num2str(pngres)])
    close
end

% Cropland0 area
if exist('crop0_area_YXBH','var')
    make_LUdiff_fig_v2(...
        crop0_area_YXBH, crop0_diff_YXrH, ...
        thisY1, thisYN, 'Cropland0', runList, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
        Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
    if do_save
        export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop0.png'],['-r' num2str(pngres)])
        close
    end
end

% "Pasture0" area
if exist('past0_area_YXBH','var')
    make_LUdiff_fig_v2(...
        past0_area_YXBH, past0_diff_YXrH, ...
        thisY1, thisYN, '"Pasture0"', runList, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
        Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
    if do_save
        export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past0.png'],['-r' num2str(pngres)])
        close
    end
end


%% Map changes in LU area: End-Historical to Begin-Future

% % Options %%%%%%%%%
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% textX = 25 ;
% textY_1 = 50 ;
% textY_2 = 20 ;
% thisPos = figurePos ;
% nx = 3 ;
% ny = 2 ;
% colorBarLoc = 'SouthOutside' ;
% conv_fact_map = 1e-6 ;   % m2 to km2
% conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
% units_map = 'km^2' ;
% units_total = 'Mkm^2' ;
% only1bl = true ;
% %%%%%%%%%%%%%%%%%%%
% 
% thisY1 = yearList_baseline(end) ;
% thisYN = yearList_future(1) ;
% 
% % Natural area
% make_LUdiff_fig_v2(...
%     ntrl_area_YXBH, ntrl_area_YXBFr - repmat(ntrl_area_YXBH,[1 1 Nruns]), ...
%     thisY1, thisYN, '"Natural"', runList, ...
%     spacing, fontSize, textX, textY_1, textY_2, ...
%     nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%     Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
% if do_save
%     export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
%     close
% end
% 
% % Cropland area
% make_LUdiff_fig_v2(...
%     crop_area_YXBH, crop_area_YXBFr - repmat(crop_area_YXBH,[1 1 Nruns]), ...
%     thisY1, thisYN, 'Cropland', runList, ...
%     spacing, fontSize, textX, textY_1, textY_2, ...
%     nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%     Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
% if do_save
%     export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
%     close
% end
% 
% % "Pasture" area
% make_LUdiff_fig_v2(...
%     past_area_YXBH, past_area_YXBFr - repmat(past_area_YXBH,[1 1 Nruns]), ...
%     thisY1, thisYN, '"Pasture"', runList, ...
%     spacing, fontSize, textX, textY_1, textY_2, ...
%     nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%     Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
% if do_save
%     export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
%     close
% end
% 
% % Agricultural area
% make_LUdiff_fig_v2(...
%     agri_area_YXBH, agri_area_YXBFr - repmat(agri_area_YXBH,[1 1 Nruns]), ...
%     thisY1, thisYN, 'Agricultural', runList, ...
%     spacing, fontSize, textX, textY_1, textY_2, ...
%     nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%     Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
% if do_save
%     export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_agri.png'],['-r' num2str(pngres)])
%     close
% end
% 
% % Cropland0 area
% if exist('crop0_area_YXBH','var')
%     make_LUdiff_fig_v2(...
%         crop0_area_YXBH, crop_area_YXBFr - repmat(crop0_area_YXBH,[1 1 Nruns]), ...
%         thisY1, thisYN, 'Cropland0', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop0.png'],['-r' num2str(pngres)])
%         close
%     end
% end
% 
% % "Pasture0" area
% if exist('past0_area_YXBH','var')
%     make_LUdiff_fig_v2(...
%         past0_area_YXBH, past_area_YXBFr - repmat(past0_area_YXBH,[1 1 Nruns]), ...
%         thisY1, thisYN, '"Pasture0"', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past0.png'],['-r' num2str(pngres)])
%         close
%     end
% end


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
conv_fact_map = 1e-6 ;   % m2 to km2
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
units_map = 'km^2' ;
units_total = 'Mkm^2' ;
only1bl = false ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_future(1) ;
thisYN = yearList_future(end) ;

% Natural area
make_LUdiff_fig_v2(...
    ntrl_area_YXBFr, ntrl_diff_YXrF, ...
    thisY1, thisYN, '"Natural"', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
    close
end

% Cropland area
make_LUdiff_fig_v2(...
    crop_area_YXBFr, crop_diff_YXrF, ...
    thisY1, thisYN, 'Cropland', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
    close
end

% "Pasture" area
make_LUdiff_fig_v2(...
    past_area_YXBFr, past_diff_YXrF, ...
    thisY1, thisYN, '"Pasture"', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
    close
end

% Agricultural area
make_LUdiff_fig_v2(...
    agri_area_YXBFr, agri_diff_YXrF, ...
    thisY1, thisYN, 'Agricultural', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
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
conv_fact_map = 1e-6 ;   % m2 to km2
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
units_map = 'km^2' ;
units_total = 'Mkm^2' ;
only1bl = true ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;
for c = 1:Ncrops
    make_LUdiff_fig_v2(...
        maps_cropareas_YXvBH(:,:,c), squeeze(maps_cropareasDiffs_YXvrH(:,:,c,:)), ...
        thisY1, thisYN, CFTnames_maps{c}, runList, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
        Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
    if do_save
        export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_maps{c} '.png'],['-r' num2str(pngres)])
        close
    end
end


%% Map changes in each crop area: End-Historical to Begin-Future

% % Options %%%%%%%%%
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% textX = 25 ;
% textY_1 = 50 ;
% textY_2 = 20 ;
% thisPos = figurePos ;
% nx = 3 ;
% ny = 2 ;
% colorBarLoc = 'SouthOutside' ;
% conv_fact_map = 1e-6 ;   % m2 to km2
% conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
% units_map = 'km^2' ;
% units_total = 'Mkm^2' ;
% only1bl = true ;
% %%%%%%%%%%%%%%%%%%%
% 
% thisY1 = yearList_baseline(end) ;
% thisYN = yearList_future(1) ;
% 
% for c = 1:Ncrops
%     make_LUdiff_fig_v2(...
%         maps_cropareas_YXvBH(:,:,c), squeeze(maps_cropareas_YXvBFr(:,:,c,:))-maps_cropareas_YXvBH(:,:,c), ...
%         thisY1, thisYN, CFTnames_maps{c}, runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_maps{c} '.png'],['-r' num2str(pngres)])
%         close
%     end
% end


%% Map changes in each crop area: Begin-Future to End-Future

% % Options %%%%%%%%%
% fontSize = 14 ;
% spacing = [0.05 0.05] ;   % [vert, horz]
% textX = 25 ;
% textY_1 = 50 ;
% textY_2 = 20 ;
% thisPos = [1 33 935 772] ;
% nx = 2 ;
% ny = 4 ;
% colorBarLoc = 'EastOutside' ;
% conv_fact_map = 1e-6 ;   % m2 to km2
% conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
% units_map = 'km^2' ;
% units_total = 'Mkm^2' ;
% only1bl = false ;
% %%%%%%%%%%%%%%%%%%%
% 
% thisY1 = yearList_future(1) ;
% thisYN = yearList_future(end) ;
% 
% for c = 1:Ncrops
%     make_LUdiff_fig_v2(...
%         squeeze(maps_cropareas_YXvBFr(:,:,c,:)), squeeze(maps_cropareasDiffs_YXvrF(:,:,c,:)), ...
%         thisY1, thisYN, CFTnames_maps{c}, runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_maps{c} '.png'],['-r' num2str(pngres)])
%         close
%     end
% end


%% Plot timeseries: Land uses

rebase = false ;

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 18 ;
spacing = [0.15 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
blYears = yearList_baseline ;
% blYears = 1960:yearList_baseline(end) ;
conv_fact = 1e-6*1e-6 ;   % m2 to Mkm2
units = 'Million km^2' ;
% thisPos = figurePos ;
thisPos = [1 383 1440 422] ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_LUarea_crop_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_crop_bl, ts_LUarea_crop_yr, ignYrs, yearList_future) ;
[tmp_ts_LUarea_past_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_past_bl, ts_LUarea_past_yr, ignYrs, yearList_future) ;
[tmp_ts_LUarea_ntrl_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_ntrl_bl, ts_LUarea_ntrl_yr, ignYrs, yearList_future) ;

figure('Position',thisPos,'Color','w') ;

subplot_tight(1,2,1,spacing)
if exist('ts_LUarea_past0_bl','var')
    tmp = ts_LUarea_past0_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
    hold on
    tmp = ts_LUarea_crop0_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'--k','LineWidth',lineWidth)
end
tmp = ts_LUarea_past_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
if exist('ts_LUarea_past0_bl','var')
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'-','LineWidth',lineWidth,'Color',0.75*ones(1,3))
else
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
end
hold on
tmp = ts_LUarea_crop_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
if exist('ts_LUarea_crop0_bl','var')
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'--','LineWidth',lineWidth,'Color',0.75*ones(1,3))
else
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'--k','LineWidth',lineWidth)
end
set(gca,'ColorOrderIndex',1) ;
plot(yearList_future,conv_fact*tmp_ts_LUarea_past_yr,'-','LineWidth',lineWidth)
set(gca,'ColorOrderIndex',1) ;
plot(yearList_future,conv_fact*tmp_ts_LUarea_crop_yr,'--','LineWidth',lineWidth)

hold off
% legend([strcat(stdLegend,', pasture') strcat(stdLegend,', cropland')], ...
%        'Location','NorthEastOutside') ;
legend('Pasture','Cropland', ...
       'Location','SouthEast') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel(units)
ht = title(['Agricultural area' title_suffix]) ;
letterlabel_align0('A',ht,do_caps) ;

subplot_tight(1,2,2,spacing) ;
tmp = ts_LUarea_ntrl_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,conv_fact*movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,conv_fact*tmp_ts_LUarea_ntrl_yr,'LineWidth',lineWidth)
hold off
legend(stdLegend, ...
       'Location','SouthWest') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel(units)
ht = title(['Natural area' title_suffix]) ;
letterlabel_align0('B',ht,do_caps) ;

if do_save
    export_fig([outDir_ts 'landUse' file_suffix '.pdf'])
    close
end


%% Plot timeseries: N fertilizer, irrigation

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 18 ;
ignYrs = 0 ;
Nsmth = 1 ;
spacing = [0.15 0.1] ;   % vert, horiz
blYears = yearList_baseline ;
% blYears = 1960:yearList_baseline(end) ;
perArea = false ;
thisPos = [1 350 1440 455] ;
%%%%%%%%%%%%%%%%%%%

if perArea
    this_ts_nflux_fert_bl = ts_nflux_fert_bl ...
                         ./ (ts_LUarea_crop_bl * 1e-4) ; % m2 to ha
    this_ts_nflux_fert_yr = ts_nflux_fert_yr ...
                         ./ (ts_LUarea_crop_yr * 1e-4) ; % m2 to ha
%     this_ts_nflux_fert_bl = -(ts_nflux_fert_CerealsC3_bl+ts_nflux_fert_CerealsC4_bl+ts_nflux_fert_Miscanthus_bl+ts_nflux_fert_Oilcrops_bl+ts_nflux_fert_Pulses_bl+ts_nflux_fert_Rice_bl+ts_nflux_fert_StarchyRoots_bl) ...
%                          ./ (ts_LUarea_crop_bl * 1e-4) ; % m2 to ha
%     this_ts_nflux_fert_yr = -(ts_nflux_fert_CerealsC3_yr+ts_nflux_fert_CerealsC4_yr+ts_nflux_fert_Miscanthus_yr+ts_nflux_fert_Oilcrops_yr+ts_nflux_fert_Pulses_yr+ts_nflux_fert_Rice_yr+ts_nflux_fert_StarchyRoots_yr) ...
%                          ./ (ts_LUarea_crop_yr * 1e-4) ; % m2 to ha           
    units_nflux_fert = 'kg/ha' ;
else
    this_ts_nflux_fert_bl = ts_nflux_fert_bl * 1e-9 ; % kgN to TgN
    this_ts_nflux_fert_yr = ts_nflux_fert_yr * 1e-9 ; % kgN to TgN
    units_nflux_fert = 'TgN' ;
end
if perArea
    this_ts_irrig_bl = (ts_irrig_bl * 1e3) ... % m3 to L
                         ./ (ts_LUarea_crop_bl * 1e-4) ; % m2 to ha
    this_ts_irrig_yr = (ts_irrig_yr * 1e3) ... % m3 to L
                         ./ (ts_LUarea_crop_yr * 1e-4) ; % m2 to ha
    units_irrig = 'L/ha' ;
else
    this_ts_irrig_bl = ts_irrig_bl*cf_m3_to_km3*1e-3 ; % m3 to 1000km3
    this_ts_irrig_yr = ts_irrig_yr*cf_m3_to_km3*1e-3 ; % m3 to 1000km3
    units_irrig = '1000 km^3' ;
end


[tmp_this_ts_nflux_fert_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, this_ts_nflux_fert_bl, this_ts_nflux_fert_yr, ignYrs, yearList_future) ;

figure('Position',thisPos,'Color','w') ;
subplot_tight(1,2,1,spacing)
tmp = this_ts_nflux_fert_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    -tmp, [], -tmp_this_ts_nflux_fert_yr, ...
    1, Nsmth, ...
    units_nflux_fert, stdLegend, 'N fertilization', title_suffix, lineWidth, fontSize, ...
    skip3rdColor)

clear tmp_*

%%%%%%%%%%%%%%%%%%%
% Plot timeseries: Irrigation
% Options %%%%%%%%%
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_this_ts_irrig_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, this_ts_irrig_bl, this_ts_irrig_yr, ignYrs, yearList_future) ;

subplot_tight(1,2,2,spacing)
tmp = this_ts_irrig_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    tmp, [], tmp_this_ts_irrig_yr, ...
    1, Nsmth, ...
    units_irrig, stdLegend, 'Irrigation', title_suffix, lineWidth, fontSize, ...
    skip3rdColor)

clear tmp_*

if perArea
    file_suffix = [file_suffix '_perArea'] ;
end

if do_save
    export_fig([outDir_ts 'mgmt_inputs' file_suffix '.pdf'])
    close
end


%% Plot timeseries: N fertilizer on each crop

% Options %%%%%%%%%
thisVar = 'nflux_fert' ;
title_prefix = 'N applied' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = 1e-3*1e-6 ;   % kg to Mt
units = 'Mt N' ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
thisLegend = stdLegend ;

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames);for c = 1:length(cmds); eval(cmds{c}); end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, {}, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
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
if exist('crop0_area_YXBH','var')
    thisVar = 'nflux0_fert' ;
    cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    file_suffix = CFT_timeseries_plot(...
        cell_bl, cell_yr, {}, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
        theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
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
end


%% Plot timeseries: N fertilizer PER HECTARE on each crop

% % Options %%%%%%%%%
% thisVar = 'nflux_fert' ;
% title_prefix = 'N applied per ha' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% ignYrs = 0 ;
% Nsmth = 1 ;
% conv_fact = 1e4 ;   % kg/m2 to kg/ha (kg/m2 * m2/ha)
% units = 'kg ha^{-1}' ;
% %%%%%%%%%%%%%%%%%%%
% 
% theseCFTnames = CFTnames ;
% thisLegend = stdLegend ;
% 
% clear cell_bl cell_yr
% cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_bl_nflux_fert = cell_bl ;
% cell_yr_nflux_fert = cell_yr ;
% clear cell_bl cell_yr
% cmds = get_cell_forPlot(whos, 'croparea', 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_bl_croparea = cell_bl ;
% clear cell_bl
% cmds = get_cell_forPlot(whos, 'croparea0', 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_bl_croparea0 = cell_bl ;
% cmds = get_cell_forPlot(whos, 'croparea', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_yr_croparea = cell_yr ;
% clear cell_bl cell_yr
% % Get kg/m2
% for i = 1:length(cell_bl_nflux_fert)
%     cell_bl{i} = cell_bl_nflux_fert{i} ./ cell_bl_croparea{i} ;
% %     cell_bl{i} = cell_bl_nflux_fert{i} ./ cell_bl_croparea0{i} ;
% %     cell_bl{i} = cell_bl_nflux_fert{i} ./ cell_bl_croparea{i} .* (cell_bl_croparea0{i} ./ cell_bl_croparea{i});
%     cell_yr{i} = cell_yr_nflux_fert{i} ./ cell_yr_croparea{i} ;
% end
% 
% file_suffix = CFT_timeseries_plot(...
%     cell_bl, cell_yr, {}, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'ignYrs',ignYrs, ...
%     'Nsmth',Nsmth, ...
%     'conv_fact',conv_fact) ;
% if do_save
%     export_fig([outDir_ts 'nfertPerHa' file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: Production of each crop

% Options %%%%%%%%%
thisVar = 'cropprod' ;
title_prefix = 'Production' ;
plum_area_adjustment = 1 - unhCropFrac ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = cf_kg2Mt ;
units = 'Mt DM' ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
thisLegend = stdLegend_plusFAO ;

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end

% Adjust for "technology" increase
cell_bl = adj_yield_tech(cell_bl, yearList_baseline) ;
cell_yr = adj_yield_tech(cell_yr, yearList_future) ;

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact, ...
    'plum_area_adjustment', plum_area_adjustment) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries: Production of each crop, EXPECTED

% % Options %%%%%%%%%
% % plum_area_adjustment = 1 ;
% % plum_area_adjustment = 1-0.28 ;
% plum_area_adjustment = 1 - unhCropFrac ;
% lpjg_area_adjustment = 1 ;
% % lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;
% thisVar = 'cropprodExp' ;
% title_prefix = 'Production (exp.)' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% ignYrs = 0 ;
% Nsmth = 1 ;
% conv_fact = cf_kg2Mt ;
% units = 'Mt DM' ;
% %%%%%%%%%%%%%%%%%%%
% 
% theseCFTnames = CFTnames ;
% thisLegend = stdLegend_plusFAO ;
% 
% clear cell_bl cell_yr
% cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% 
% file_suffix = CFT_timeseries_plot(...
%     cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'ignYrs',ignYrs, ...
%     'Nsmth',Nsmth, ...
%     'conv_fact',conv_fact, ...
%     'plum_area_adjustment', plum_area_adjustment, ...
%     'lpjg_area_adjustment', lpjg_area_adjustment) ;
% 
% if do_save
%     export_fig([outDir_ts thisVar file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: Production of each crop, ACTUAL/EXPECTED

% % Options %%%%%%%%%
% % title_prefix = 'Prod. diff: (Act.-Exp./Exp.)' ;
% title_prefix = 'Prod. diff' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% Nsmth = 1 ;
% units = '%' ;
% %%%%%%%%%%%%%%%%%%%
% 
% theseCFTnames = CFTnames ;
% thisLegend = stdLegend(2:end) ;
% 
% clear cell_yr cell_yr_act cell_yr_exp
% cmds = get_cell_forPlot(whos, 'cropprod', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_yr_act = cell_yr ; clear cell_yr
% cell_yr_act = adj_yield_tech(cell_yr_act, yearList_future) ;
% cmds = get_cell_forPlot(whos, 'cropprodExp', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_yr_exp = cell_yr ; clear cell_yr
% CFT_ActVsExp_plot(...
%     cell_yr_act, cell_yr_exp, yearList_future, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'Nsmth',Nsmth) ;
% 
% if do_save
%     export_fig([outDir_ts 'cropProdDiffFromExp.pdf'])
%     close
% end



%% Plot timeseries: Area of each crop

% Options %%%%%%%%%
plum_area_adjustment = 1 ;
% plum_area_adjustment = 1-0.28 ;
% plum_area_adjustment = 1 - unhCropFrac ;
lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;
thisVar = 'croparea' ;
title_prefix = 'Area' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = 1e-6*1e-6 ;
units = 'Million km^2' ;
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
    CFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact, ...
    'plum_area_adjustment', plum_area_adjustment, ...
    'lpjg_area_adjustment', lpjg_area_adjustment) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end



%% Plot timeseries: Yield of each crop

% Options %%%%%%%%%
do_adjYieldTech = true ;
thisVar = 'yield' ;
title_prefix = 'Yield' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = cf_kgPm2_to_tonsPha ;
units = 't ha^{-1}' ;
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

% Adjust for "technology" increase (do not include FAO!)
if do_adjYieldTech
    cell_bl = adj_yield_tech(cell_bl, yearList_baseline) ;
    cell_yr = adj_yield_tech(cell_yr, yearList_future) ;
    title_prefix = [title_prefix ' (techAdj)'] ;
end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;

if do_adjYieldTech
    file_suffix = [file_suffix '_techAdj'] ;
end

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries: Yield of each crop (EXPECTED)

% % Options %%%%%%%%%
% thisVar = 'yieldExp' ;
% title_prefix = 'Yield (PLUM-exp)' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% ignYrs = 0 ;
% Nsmth = 1 ;
% conv_fact = cf_kgPm2_to_tonsPha ;
% units = 't ha^{-1}' ;
% %%%%%%%%%%%%%%%%%%%
% 
% theseCFTnames = CFTnames ;
% if include_fao
%     thisLegend = stdLegend_plusFAO ;
% else
%     thisLegend = stdLegend ;
% end
% 
% clear cell_bl cell_yr
% cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% if include_fao
%     cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% else
%     cell_fao = {} ;
%     fao.tmp_fao_yearList = [] ;
% end
% 
% file_suffix = CFT_timeseries_plot(...
%     cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'ignYrs',ignYrs, ...
%     'Nsmth',Nsmth, ...
%     'conv_fact',conv_fact) ;
% 
% if do_save
%     export_fig([outDir_ts thisVar file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: Yield of each crop, ACTUAL/EXPECTED

% % Options %%%%%%%%%
% do_adjYieldTech = false ;
% title_prefix = 'Yield diff' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% Nsmth = 1 ;
% units = '%' ;
% %%%%%%%%%%%%%%%%%%%
% 
% theseCFTnames = CFTnames ;
% thisLegend = stdLegend(2:end) ;
% 
% clear cell_yr cell_yr_act cell_yr_exp
% cmds = get_cell_forPlot(whos, 'yield', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_yr_act = cell_yr ; clear cell_yr
% cmds = get_cell_forPlot(whos, 'yieldExp', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_yr_exp = cell_yr ; clear cell_yr
% 
% % Adjust for technology change, if doing so
% file_suffix = '' ;
% if do_adjYieldTech
%     cell_yr_act = adj_yield_tech(cell_yr_act, yearList_future) ;
%     title_prefix = [title_prefix ' (techAdj)'] ;
%     file_suffix = [file_suffix '_techAdj'] ;
% end
% 
% CFT_ActVsExp_plot(...
%     cell_yr_act, cell_yr_exp, yearList_future, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'Nsmth',Nsmth) ;
% 
% if do_save
%     export_fig([outDir_ts 'yieldDiffFromExp' file_suffix '.pdf'])
%     close
% end



%% Plot timeseries: total area (FOR TROUBLESHOOTING)

% error('Make this work with SI units!')
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

% % ts_nloss_bl = ts_nflux_flux_bl + ts_nflux_harvest_bl + ts_nflux_leach_bl + ts_nflux_LUch_bl ;
% % ts_nloss_yr = ts_nflux_flux_yr + ts_nflux_harvest_yr + ts_nflux_leach_yr + ts_nflux_LUch_yr ;
% ts_nloss_bl = ts_nflux_flux_bl + ts_nflux_leach_bl ;
% ts_nloss_yr = ts_nflux_flux_yr + ts_nflux_leach_yr ;
% 
[tmp_ts_nflux_flux_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_nflux_flux_bl, ts_nflux_flux_yr, ignYrs, yearList_future) ;
[tmp_ts_nflux_leach_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_nflux_leach_bl, ts_nflux_leach_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_nflux_flux_bl, [], tmp_ts_nflux_flux_yr, ...
    cf_kg2Tg, Nsmth, ...
    'TgN', stdLegend, 'N loss: Gaseous', '', lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('A',ht,do_caps) ;

subplot_tight(1,2,2,spacing)
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_nflux_leach_bl, [], tmp_ts_nflux_leach_yr, ...
    cf_kg2Tg, Nsmth, ...
    'TgN', stdLegend, 'N loss: Dissolved', '', lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('B',ht,do_caps) ;

clear tmp_*

if do_save
    export_fig([outDir_ts 'Nloss' file_suffix '.pdf'])
    close
end


%% Bar graph: N loss

% years_endh = 2000:2009 ;
% years_begf = 2011:2020 ;
% years_endf = 2090:2099 ;
% 
% % Name, code, conversion factor, formatSpec mean, formatSpec SEM
% rowInfo = {'N loss (TgN)', 'nloss', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'N loss: Gaseous (TgN)', 'nflux_flux', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'N loss: Dissolved (TgN)', 'nflux_leach', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            } ;


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
if skip3rdColor
    plot(yearList_future,tmp_ts_albedo1_yr(:,1),'LineWidth',lineWidth)
    h = plot(yearList_future,tmp_ts_albedo1_yr(:,2),'LineWidth',lineWidth) ;
    set(h,'Color',[255 159 56]/255) ;
    set(gca,'ColorOrderIndex',4) ;
    plot(yearList_future,tmp_ts_albedo1_yr(:,3),'LineWidth',lineWidth)
else
    plot(yearList_future,tmp_ts_albedo1_yr,'LineWidth',lineWidth)
end



ax = gca;
ax.ColorOrderIndex = 1;
plot(yearList_baseline,movmean(ts_albedo7_bl,Nsmth),'--k','LineWidth',lineWidth)
% plot(yearList_future,tmp_ts_albedo7_yr,'--','LineWidth',lineWidth)
if skip3rdColor
    plot(yearList_future,tmp_ts_albedo7_yr(:,1),'LineWidth',lineWidth)
    h = plot(yearList_future,tmp_ts_albedo7_yr(:,2),'LineWidth',lineWidth) ;
    set(h,'Color',[255 159 56]/255) ;
    set(gca,'ColorOrderIndex',4) ;
    plot(yearList_future,tmp_ts_albedo7_yr(:,3),'LineWidth',lineWidth)
else
    plot(yearList_future,tmp_ts_albedo7_yr,'LineWidth',lineWidth)
end
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
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_cpool_VegC_bl, [], tmp_ts_cpool_VegC_yr, ...
    cf_kg2Pg, Nsmth, ...
    'PgC', stdLegend, 'Vegetation C', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('A',ht,do_caps) ;

subplot_tight(1,2,2,spacing)
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_cpool_Total_bl, [], tmp_ts_cpool_Total_yr, ...
    cf_kg2Pg, Nsmth, ...
    'PgC', stdLegend, 'Total C', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('B',ht,do_caps) ;

clear tmp_*

if do_save
    export_fig([outDir_ts 'cpools' file_suffix '.pdf'])
    close
end


%% Plot timeseries: BVOCs

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_amon_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_amon_bl, ts_amon_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing) ;
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_amon_bl, [], tmp_ts_amon_yr, ...
    cf_kg2Pg, Nsmth, ...
    'TgC', stdLegend, 'Monoterpene emissions', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('A',ht,do_caps) ;
clear tmp_*

[tmp_ts_aiso_yr, title_suffix, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_aiso_bl, ts_aiso_yr, ignYrs, yearList_future) ;

subplot_tight(1,2,2,spacing) ;
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_aiso_bl, [], tmp_ts_aiso_yr, ...
    cf_kg2Pg, Nsmth, ...
    'TgC', stdLegend, 'Isoprene emissions', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('B',ht,do_caps) ;
set(gca,'FontSize',fontSize)

clear tmp_*

if do_save
    export_fig([outDir_ts 'bvoc' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Evapotranspiration and runoff

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 5 ;
figurePosition = [1 376 1440 429] ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_aevapaaet_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_aevapaaet_bl, ts_aevapaaet_yr, ignYrs, yearList_future) ;

figure('Position',figurePosition,'Color','w') ;
subplot_tight(1,2,1,spacing) ;
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_aevapaaet_bl, [], tmp_ts_aevapaaet_yr, ...
    cf_m3_to_km3*1e-3, Nsmth, ...
    '1000 km^3', stdLegend, 'Evapotranspiration', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('A',ht,do_caps) ;
clear tmp_*

[tmp_ts_tot_runoff_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_tot_runoff_bl, ts_tot_runoff_yr, ignYrs, yearList_future) ;
subplot_tight(1,2,2,spacing) ;
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_tot_runoff_bl, [], tmp_ts_tot_runoff_yr, ...
    cf_m3_to_km3*1e-3, Nsmth, ...
    '1000 km^3', stdLegend, 'Runoff', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('B',ht,do_caps) ;
clear tmp_*

if do_save
    export_fig([outDir_ts 'evapotranspiration_runoff' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Calories

% Options %%%%%%%%%
do_adjYieldTech = true ;

plum_area_adjustment = 1 ;
% plum_area_adjustment = 1-0.28 ;

lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;

lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 1 ;
figurePosition = [1 376 1440 429] ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_kcal_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_kcal_bl, ts_kcal_yr, ignYrs, yearList_future) ;

% Adjust for technology change, if doing so (do not include FAO!)
tmp_ts_kcal_bl = ts_kcal_bl ;
if do_adjYieldTech
    tmp_ts_kcal_bl = adj_yield_tech(tmp_ts_kcal_bl,yearList_baseline) ;
    tmp_ts_kcal_yr = adj_yield_tech(tmp_ts_kcal_yr,yearList_future) ;
    file_suffix = [file_suffix '_techAdj'] ;
    title_suffix = [title_suffix ' (techAdj)'] ;
end

% Adjust for unhandled crops, if doing so (DO include FAO)
tmp_ts_kcal_fao = ts_kcal_fao ;
if ~isequal(lpjg_area_adjustment,1)
    tmp_ts_kcal_bl = tmp_ts_kcal_bl .* lpjg_area_adjustment ;
    [~,IA] = intersect(yearList_baseline,fao.tmp_fao_yearList) ;
    tmp_ts_kcal_fao = tmp_ts_kcal_fao .* lpjg_area_adjustment(IA) ;
    title_suffix = [title_suffix ' (blAdj)'] ;
    file_suffix = [file_suffix '_blAdj'] ;
end
if plum_area_adjustment ~= 1
    tmp_ts_kcal_yr = tmp_ts_kcal_yr * plum_area_adjustment ;
    title_suffix = [title_suffix ' (PLUM\times' num2str(plum_area_adjustment) ')'] ;
    file_suffix = [file_suffix '_plumAdj' num2str(plum_area_adjustment)] ;
end

figure('Position',figurePosition,'Color','w') ;

plot_timeseries(...
    yearList_baseline, fao.tmp_fao_yearList, yearList_future, ...
    tmp_ts_kcal_bl, tmp_ts_kcal_fao, tmp_ts_kcal_yr, ...
    cf_kcal2Ecal, Nsmth, ...
    'Ecal', stdLegend_plusFAO, 'Caloric production', ...
    title_suffix, lineWidth, fontSize, skip3rdColor) ;
clear tmp_*

if do_save
    export_fig([outDir_ts 'calories' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Calories2

% Options %%%%%%%%%
do_adjYieldTech = true ;

plum_area_adjustment = 1 ;
% plum_area_adjustment = 1-0.28 ;

lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;

lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 1 ;
figurePosition = [1 376 1440 429] ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_kcal2_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_kcal2_bl, ts_kcal2_yr, ignYrs, yearList_future) ;

% Adjust for technology change, if doing so (do not include FAO!)
tmp_ts_kcal2_bl = ts_kcal2_bl ;
if do_adjYieldTech
    tmp_ts_kcal2_bl = adj_yield_tech(tmp_ts_kcal2_bl,yearList_baseline) ;
    tmp_ts_kcal2_yr = adj_yield_tech(tmp_ts_kcal2_yr,yearList_future) ;
    file_suffix = [file_suffix '_techAdj'] ;
    title_suffix = [title_suffix ' (techAdj)'] ;
end

% Adjust for unhandled crops, if doing so (DO include FAO)
tmp_ts_kcal2_fao = ts_kcal2_fao ;
if ~isequal(lpjg_area_adjustment,1)
    tmp_ts_kcal2_bl = tmp_ts_kcal2_bl .* lpjg_area_adjustment ;
    [~,IA] = intersect(yearList_baseline,fao.tmp_fao_yearList) ;
    tmp_ts_kcal2_fao = tmp_ts_kcal2_fao .* lpjg_area_adjustment(IA) ;
    title_suffix = [title_suffix ' (blAdj)'] ;
    file_suffix = [file_suffix '_blAdj'] ;
end
if plum_area_adjustment ~= 1
    tmp_ts_kcal2_yr = tmp_ts_kcal2_yr * plum_area_adjustment ;
    title_suffix = [title_suffix ' (PLUM\times' num2str(plum_area_adjustment) ')'] ;
    file_suffix = [file_suffix '_plumAdj' num2str(plum_area_adjustment)] ;
end

figure('Position',figurePosition,'Color','w') ;

plot_timeseries(...
    yearList_baseline, fao.tmp_fao_yearList, yearList_future, ...
    tmp_ts_kcal2_bl, tmp_ts_kcal2_fao, tmp_ts_kcal2_yr, ...
    cf_kcal2Ecal, Nsmth, ...
    'Ecal', stdLegend_plusFAO, 'Caloric production', ...
    title_suffix, lineWidth, fontSize, skip3rdColor) ;
clear tmp_*

if do_save
    export_fig([outDir_ts 'calories2' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Calories, EXPECTED

% % Options %%%%%%%%%
% 
% do_adjYieldTech = true ;
% 
% plum_area_adjustment = 1 ;
% % plum_area_adjustment = 1-0.28 ;
% 
% % lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;
% 
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% ignYrs = 0 ;
% Nsmth = 1 ;
% figurePosition = [1 376 1440 429] ;
% %%%%%%%%%%%%%%%%%%%
% 
% [tmp_ts_kcalExp_yr, title_suffix, file_suffix] = ...
%     rebase_future2baseline(rebase, Nsmth, ts_kcal_bl, ts_kcalExp_yr, ignYrs, yearList_future) ;
% 
% % Adjust for technology change, if doing so (do not include FAO!)
% tmp_ts_kcal_bl = ts_kcal_bl ;
% if do_adjYieldTech
%     tmp_ts_kcal_bl = adj_yield_tech(tmp_ts_kcal_bl,yearList_baseline) ;
%     tmp_ts_kcalExp_yr = adj_yield_tech(tmp_ts_kcalExp_yr,yearList_future) ;
%     file_suffix = [file_suffix '_techAdj'] ;
%     title_suffix = [title_suffix ' (techAdj)'] ;
% end
% 
% % Adjust for unhandled crops, if doing so (DO include FAO)
% tmp_ts_kcal_fao = ts_kcal_fao ;
% if lpjg_area_adjustment ~= 1
%     tmp_ts_kcal_bl = tmp_ts_kcal_bl .* lpjg_area_adjustment ;
%     [~,IA] = intersect(yearList_baseline,fao.tmp_fao_yearList) ;
%     tmp_ts_kcal_fao = tmp_ts_kcal_fao .* lpjg_area_adjustment(IA) ;
%     title_suffix = [title_suffix ' (blAdj)'] ;
%     file_suffix = [file_suffix '_blAdj'] ;
% end
% if plum_area_adjustment ~= 1
%     tmp_ts_kcalExp_yr = tmp_ts_kcalExp_yr * plum_area_adjustment ;
%     title_suffix = [title_suffix ' (PLUM\times' num2str(plum_area_adjustment) ')'] ;
%     file_suffix = [file_suffix '_plumAdj' num2str(plum_area_adjustment)] ;
% end
% 
% figure('Position',figurePosition,'Color','w') ;
% plot(yearList_baseline,cf_kcal2Ecal*movmean(tmp_ts_kcal_bl,Nsmth),'-k','LineWidth',lineWidth)
% hold on
% plot(fao.tmp_fao_yearList,cf_kcal2Ecal*tmp_ts_kcal_fao,'--k','LineWidth',lineWidth)
% plot(yearList_future,cf_kcal2Ecal*tmp_ts_kcalExp_yr,'LineWidth',lineWidth)
% hold off
% legend(stdLegend_plusFAO, ...
%        'Location','NorthWest') ;
% set(gca,'FontSize',fontSize)
% xlabel('Year')
% ylabel('Ecal')
% ht = title(['Caloric production (PLUM-exp)' title_suffix]) ;
% clear tmp_*
% 
% if do_save
%     export_fig([outDir_ts 'caloriesExp' file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: Calories, ACTUAL/EXPECTED

% % Options %%%%%%%%%
% do_adjYieldTech = false ;
% lineWidth = 2 ;
% fontSize = 14 ;
% Nsmth = 1 ;
% figurePosition = [1 376 1440 429] ;
% %%%%%%%%%%%%%%%%%%%
% 
% file_suffix = '' ;
% title_suffix = '' ;
% if do_adjYieldTech
%     file_suffix = [file_suffix '_techAdj'] ;
%     title_suffix = [title_suffix ' (techAdj)'] ;
% end
% 
% figure('Position',figurePosition,'Color','w') ;
% tmp_ts_kcal_yr = ts_kcal_yr ;
% if do_adjYieldTech
%     tmp_ts_kcal_yr = adj_yield_tech(tmp_ts_kcal_yr,yearList_future) ;
% end
% plot(yearList_future,movmean((tmp_ts_kcal_yr - ts_kcalExp_yr)./ts_kcalExp_yr*100,Nsmth,1),'-','LineWidth',lineWidth)
% hold on
% plot(get(gca,'XLim'),[0 0],'--k')
% hold off
% legend(stdLegend(2:end),'Location','SouthWest') ;
% set(gca,'FontSize',fontSize)
% xlabel('Year')
% ylabel('Ecal')
% title(['Calories diff.' title_suffix]) ;
% clear tmp_*
% 
% if do_save
%     export_fig([outDir_ts 'caloriesDiffFromExp' file_suffix '.pdf'])
%     close
% end



%% Map end-of-run mean cropfracs

% % Options %%%%%%%%%
% title_prefix = 'Frac. crop' ;
% whichCFTs = 'lpjg' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% y2include = 65:350 ;
% %%%%%%%%%%%%%%%%%%%
% 
% error('Make this work as PLUMnative!')
% 
% Ny = 2 ;
% Nmaps = length(runList) + 1 ;
% if strcmp(whichCFTs,'lpjg')
%     Ncrops = Ncrops_lpjg ;
%     theseCrops = CFTnames ;
%     Nx = 3 ;
%     theseMaps = maps_cropfracs ;
%     figure_position = figurePos ;
% elseif strcmp(whichCFTs,'plum')
%     Nx = 4 ;
%     theseMaps = maps_cropfracs_plum7 ;
%     figure_position = figurePos ;
% else
%     error(['whichCFTs (' whichCFTs ') not recognized!']) ;
% end
% 
% 
% for c = 1:Ncrops
%     figure('Position',figure_position,'Color','w') ;
%     thisCrop = theseCrops{c} ;
%     for p = 1:Nmaps
%         subplot_tight(Ny,Nx,p,spacing)
%         if p==1
%             tmp = theseMaps.maps_YXvyB(:,:,strcmp(theseMaps.varNames,thisCrop),:) ;
%             tmp(maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'CROPLAND'),:)==0) = 0 ;
%             tmp = mean(tmp,4) ;
%             pcolor(tmp(y2include,:)) ;
%             thisTitle = 'Baseline' ;
%         else
%             pcolor(mean(theseMaps.maps_YXvyr(y2include,:,strcmp(theseMaps.varNames,thisCrop),:,p-1),4)) ;
%             thisTitle = runList{p-1} ;
%         end
%         shading flat ; axis equal tight off
%         caxis([0 1]) ; colorbar('SouthOutside')
%         title([thisCrop ' (' thisTitle ')'])
%         set(gca,'FontSize',fontSize) ;
%     end
%     
%     if do_save
%         export_fig([outDir_maps 'maps_cropfracs_' thisCrop '.png'],['-r' num2str(pngres)])
%         close
%     end
% end
% 
% clear theseMaps Ncrops Nmaps Nx Ny
% 
% 
% %% Map changes in BD hotspot area: CI
% 
% % Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edgecolor = 0.6*ones(3,1) ;
% latlim = [-60,80];
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% textX = 25 ;
% textY_1 = 50 ;
% textY_2 = 20 ;
% cbarOrient = 'SouthOutside' ;
% lineWidth = 0.25 ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% error('Make this work with SI units!')
% 
% if strcmp(version('-release'),'2014b')
%     
%     % Biodiversity hotspots
%     hotspot_shp = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/hotspots_clipByGridlist.shp' ;
%     hotspot_YX = dlmread('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/hotspots_raster.txt',...
%         ' ', 6, 0) ;
%     hotspot_YX(hotspot_YX==-9999) = NaN ;
%     hotspot_YX(:,721) = [] ;
%     hotspot_YX = flipud(hotspot_YX) ;
%     hotspot_YX(nanmask) = NaN ;
%     hotspot_YX(~nanmask & isnan(hotspot_YX)) = 0 ;
%     hotspot_area_YX = hotspot_YX.*land_area_YX ;
%     hotspot_area_YXB = hotspot_area_YX .* maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'NATURAL'),end) ;
%     % hotspot_area_YXB = hotspot_area_YX ;
%     hotspot_area_YXr = repmat(hotspot_area_YX,[1 1 Nruns]) .* squeeze(maps_LU.maps_YXvyr(:,:,strcmp(maps_LU.varNames,'NATURAL'),end,:)) ;
%     
%     hotspot_diff_YXr = hotspot_area_YXr - repmat(hotspot_area_YXB,[1 1 Nruns]) ;
%     
%     map_hotspot_diffs(...
%         hotspot_area_YXB, hotspot_diff_YXr, hotspot_YX, hotspot_shp, ...
%         spacing, latlim, edgecolor, cbarOrient, fontSize, ...
%         textX, textY_1, textY_2, ssp_plot_index, lineWidth, ...
%         yearList_baseline, yearList_future, runList)
%     
%     if do_save
%         export_fig([outDir_maps 'areaDiff_BDhotspots_CI.png'],['-r' num2str(pngres)])
%         close
%     end
%     
% else
%     warning('Skipping hotspots')
% end


%% Map changes in BD hotspot area: glob200

% % Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edgecolor = 0.6*ones(3,1) ;
% latlim = [-60,80];
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% textX = 25 ;
% textY_1 = 50 ;
% textY_2 = 20 ;
% cbarOrient = 'SouthOutside' ;
% lineWidth = 0.25 ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% error('Make this work with SI units!')
% 
% if strcmp(version('-release'),'2014b')
% 
%     % Biodiversity hotspots
%     hotspot_YX = imread('/Users/sam/Geodata/global200ecoregions/g200_terr_raster0.5deg.tif') ;
%     hotspot_shp = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work/g200_terr_from0.5raster.shp' ;
%     hotspot_YX = flipud(hotspot_YX) ;
%     hotspot_area_YX = hotspot_YX.*land_area_YX ;
%     hotspot_area_YXB = hotspot_area_YX .* maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'NATURAL'),end) ;
%     % hotspot_area_YXB = hotspot_area_YX ;
%     hotspot_area_YXr = repmat(hotspot_area_YX,[1 1 Nruns]) .* squeeze(maps_LU.maps_YXvyr(:,:,strcmp(maps_LU.varNames,'NATURAL'),end,:)) ;
%     
%     hotspot_diff_YXr = hotspot_area_YXr - repmat(hotspot_area_YXB,[1 1 Nruns]) ;
%     
%     map_hotspot_diffs(...
%         hotspot_area_YXB, hotspot_diff_YXr, hotspot_YX, hotspot_shp, ...
%         spacing, latlim, edgecolor, cbarOrient, fontSize, ...
%         textX, textY_1, textY_2, ssp_plot_index, lineWidth, ...
%         yearList_baseline, yearList_future, runList)
%     
%     if do_save
%         export_fig([outDir_maps 'areaDiff_BDhotspots_glob200.png'],['-r' num2str(pngres)])
%         close
%     end
% 
% else
%     warning('Skipping hotspots')
% end


%%

disp('All done!')



