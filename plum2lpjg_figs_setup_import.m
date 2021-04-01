%% Setup

cd '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/MATLAB_work'

% Define raster reference object and missing value
R = georasterref('RasterSize', [360 720], ...
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

addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work/')
addpath(genpath('~/Documents/Dropbox/Dissertation/MATLAB work'))

continents_shp = '/Users/sam/Geodata/General/continents_from_countries/continents_from_countries.shp' ;

trimFirstFuture = 0 ;
if strcmp(thisVer,'20180424agmip7')
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
%     runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
    runList = {'SSP1-45','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
%     runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
    runList = {'SSP1-45','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','constLU','constClimCO2'} ;
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
    runList = {'SSP3-60','constLU','constClimCO2'} ;
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
    runList = {'SSP4-60','constLU','constClimCO2'} ;
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
    runList = {'SSP5-85','constLU','constClimCO2'} ;
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
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
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
    runList = {'SSP1-45','constLU','constClimCO2'} ;
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
    runList = {'SSP3-60','constLU','constClimCO2'} ;
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
    runList = {'SSP4-60','constLU','constClimCO2'} ;
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
    runList = {'SSP5-85','constLU','constClimCO2'} ;
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
    runList = {'SSP1-45co2','SSP3-60co2','SSP4-60co2','SSP5-85co2'} ;
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
elseif strcmp(thisVer,'harm3')
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_SSP1_RCP45/output-2019-02-27-103914';
        'LPJGPLUM_2011-2100_harm3_SSP3_RCP60/output-2019-02-27-093027';
        'LPJGPLUM_2011-2100_harm3_SSP4_RCP60/output-2019-02-27-093259';
        'LPJGPLUM_2011-2100_harm3_SSP5_RCP85/output-2019-02-27-104120';
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP4.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm3_constClim')
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_SSP1_RCP45co2_constClim/output-2019-03-07-164547' ;
        'LPJGPLUM_2011-2100_harm3_SSP3_RCP60co2_constClim/output-2019-03-07-164547' ;
        'LPJGPLUM_2011-2100_harm3_SSP4_RCP60co2_constClim/output-2019-03-07-170847' ;
        'LPJGPLUM_2011-2100_harm3_SSP5_RCP85co2_constClim/output-2019-03-07-170908' ;
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP4.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm3_constCO2')
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_SSP1_RCP45_constCO2/output-2019-03-07-164547' ;
        'LPJGPLUM_2011-2100_harm3_SSP3_RCP60_constCO2/output-2019-03-07-170846' ;
        'LPJGPLUM_2011-2100_harm3_SSP4_RCP60_constCO2/output-2019-03-07-170858' ;
        'LPJGPLUM_2011-2100_harm3_SSP5_RCP85_constCO2/output-2019-03-07-170910' ;
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP4.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm3_constLU')
    runList = {'RCP4.5','RCP6.0','RCP8.5'} ;
    runColNames = {'RCP45','RCP60','RCP85'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_constLU_RCP45/output-2019-03-07-164549' ;
        'LPJGPLUM_2011-2100_harm3_constLU_RCP60/output-2019-03-07-164546' ;
        'LPJGPLUM_2011-2100_harm3_constLU_RCP85/output-2019-03-07-164546' ;
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        } ; warning('constLU run, but will display regular PLUM outputs! (from SSP3 for RCP60)')
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = true ;
elseif strcmp(thisVer,'harm3_onlyCO2')
    runList = {'RCP4.5','RCP6.0','RCP8.5'} ;
    runColNames = {'RCP45','RCP60','RCP85'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_constLU_RCP45CO2_constClim/output-2019-06-19-120106';
        'LPJGPLUM_2011-2100_harm3_constLU_RCP60CO2_constClim/output-2019-06-19-120946';
        'LPJGPLUM_2011-2100_harm3_constLU_RCP85CO2_constClim/output-2019-06-19-122314';
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        } ; warning('constLU run, but will display regular PLUM outputs! (from SSP3 for RCP60)')
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = true ;
elseif strcmp(thisVer,'harm3_onlyClim')
    runList = {'RCP4.5','RCP6.0','RCP8.5'} ;
    runColNames = {'RCP45','RCP60','RCP85'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_constLU_RCP45clim_constCO2/output-2019-06-19-122601';
        'LPJGPLUM_2011-2100_harm3_constLU_RCP60clim_constCO2/output-2019-06-19-123448';
        'LPJGPLUM_2011-2100_harm3_constLU_RCP85clim_constCO2/output-2019-06-19-124811';
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        } ; warning('constLU run, but will display regular PLUM outputs! (from SSP3 for RCP60)')
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = true ;
elseif strcmp(thisVer,'harm3_constClimCO2')
    runList = {'SSP1-45','SSP3-60','SSP4-60','SSP5-85'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_SSP1_constCO2_constClim/output-2019-03-12-091426' ;
        'LPJGPLUM_2011-2100_harm3_SSP3_constCO2_constClim/output-2019-03-12-080759' ;
        'LPJGPLUM_2011-2100_harm3_SSP4_constCO2_constClim/output-2019-03-12-083131' ;
        'LPJGPLUM_2011-2100_harm3_SSP5_constCO2_constClim/output-2019-03-12-094118' ;
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP4.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        } ;
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm3_S1R4.5_attr')
    runList = {'SSP1-45','constLU','constClimCO2','constClim'} ;
    runColNames = {'Full','constLU','constClimCO2','constClim'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_SSP1_RCP45/output-2019-02-27-103914' ;
        'LPJGPLUM_2011-2100_harm3_constLU_RCP45/output-2019-03-07-164549' ;
        'LPJGPLUM_2011-2100_harm3_SSP1_constCO2_constClim/output-2019-03-12-091426' ;
        'LPJGPLUM_2011-2100_harm3_SSP1_RCP45co2_constClim/output-2019-03-07-164547' ;
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1';
        } ; warning('constLU run, but will display regular PLUM outputs!')
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm3_S3R6.0_attr')
    runList = {'SSP3-60','constLU','constClimCO2','constClim'} ;
    runColNames = {'Full','constLU','constClimCO2','constClim'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_SSP3_RCP60/output-2019-02-27-093027' ;
        'LPJGPLUM_2011-2100_harm3_constLU_RCP60/output-2019-03-07-164546' ;
        'LPJGPLUM_2011-2100_harm3_SSP3_constCO2_constClim/output-2019-03-12-080759' ;
        'LPJGPLUM_2011-2100_harm3_SSP3_RCP60co2_constClim/output-2019-03-07-164547' ;
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v12.s1';
        } ; warning('constLU run, but will display regular PLUM outputs!')
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm3_S4R6.0_attr')
    runList = {'SSP4-60','constLU','constClimCO2','constClim'} ;
    runColNames = {'Full','constLU','constClimCO2','constClim'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_SSP4_RCP60/output-2019-02-27-093259' ;
        'LPJGPLUM_2011-2100_harm3_constLU_RCP60/output-2019-03-07-164546' ;
        'LPJGPLUM_2011-2100_harm3_SSP4_constCO2_constClim/output-2019-03-12-083131' ;
        'LPJGPLUM_2011-2100_harm3_SSP4_RCP60co2_constClim/output-2019-03-07-170847' ;
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP4.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP4.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP4.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP4.v12.s1';
        } ; warning('constLU run, but will display regular PLUM outputs!')
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
    yearList_baseline = 1850:2010 ;
    skip3rdColor = false ;
elseif strcmp(thisVer,'harm3_S5R8.5_attr')
%     runList = {'SSP5-85','constLU','constClimCO2'} ;
%     runColNames = {'Full','constLU','constClimCO2'} ;
%     runList = {'SSP5-85','constLU','constClimCO2','constClim'} ;
%     runColNames = {'Full','constLU','constClimCO2','constClim'} ;
    runList = {'SSP5-85','constLU','constClimCO2','constCO2'} ;
    runColNames = {'Full','constLU','constClimCO2','constCO2'} ;
    runDirs = {
        'LPJGPLUM_2011-2100_harm3_SSP5_RCP85/output-2019-02-27-104120' ;
        'LPJGPLUM_2011-2100_harm3_constLU_RCP85/output-2019-03-07-164546' ;
        'LPJGPLUM_2011-2100_harm3_SSP5_constCO2_constClim/output-2019-03-12-094118' ;
        'LPJGPLUM_2011-2100_harm3_SSP5_RCP85_constCO2/output-2019-03-07-170910' ;
        } ;
    runDirs_plum = {
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1';
        } ; warning('constLU run, but will display regular PLUM outputs!')
    yearList_future = 2011:2100 ;
    baselineDir = 'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;
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
clear firstdec_tmp
lastdec_tmp = load([baselineDir 'last_decade.mat']) ;
bl_map_fields = fieldnames(lastdec_tmp) ;
for f = 1:length(bl_map_fields)
    thisField = bl_map_fields{f} ;
    eval(['maps_' thisField ' = renameStructField(lastdec_tmp.' thisField ',''maps_YXvy'',''maps_YXvyB'') ;']) ;
end
clear lastdec_tmp
last30yrs_tmp = load([baselineDir 'last_30yrs.mat']) ;
bl_map_fields = fieldnames(last30yrs_tmp) ;
for f = 1:length(bl_map_fields)
    thisField = bl_map_fields{f} ;
    eval(['maps_' thisField ' = renameStructField(last30yrs_tmp.' thisField ',''maps_YXvs'',''maps_YXvsB'') ;']) ;
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
    firstdec_tmp = load([runDirs{r} 'first_decade.mat']) ;
    lastdec_tmp = load([runDirs{r} 'last_decade.mat']) ;
    last30yrs_tmp = load([runDirs{r} 'last_30yrs.mat']) ;
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
            if ~(isfield(firstdec_tmp,thisVar_out) || isfield(lastdec_tmp,thisVar_out) || isfield(last30yrs_tmp,thisVar_out))
                warning([thisVar_out ' not found in future run ' num2str(r) '. Removing.']) ;
                clear(thisVar_in)
                continue
            end
            
            if contains(thisVar_in,'_d1')
                eval([thisVar_in '.maps_YXvyr = nan([size(firstdec_tmp.' thisVar_out '.maps_YXvy) Nruns],''single'') ;']) ;
            elseif contains(thisVar_in,'_d9')
                eval([thisVar_in '.maps_YXvyr = nan([size(lastdec_tmp.' thisVar_out '.maps_YXvy) Nruns],''single'') ;']) ;
            elseif contains(thisVar_in,'_last30')
                eval([thisVar_in '.maps_YXvsr = nan([size(last30yrs_tmp.' thisVar_out '.maps_YXvs) Nruns],''single'') ;']) ;
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
        elseif contains(thisVar_in,'last30')
            if ~isfield(last30yrs_tmp,thisVar_out)
                warning([thisVar_in ' not found in future run ' num2str(r) '. Removing.']) ;
                vars_maps_bl_toRemove(v) = true ;
                eval(['clear ' thisVar_in]) ;
                continue
            end
            eval(['isequal_varNames = isequal(' thisVar_in '.varNames, last30yrs_tmp.' thisVar_out '.varNames) ; ']) ;
            if ~isequal_varNames
                eval(['isequal_varNames_afterSort = isequal(sort(' thisVar_in '.varNames), sort(last30yrs_tmp.' thisVar_out '.varNames)) ; ']) ;
                if ~isequal_varNames_afterSort
                    error(['~isequal_varNames (' thisVar_in '). Not fixable.'])
                end
                warning(['~isequal_varNames (' thisVar_in '). Fixing.'])
                eval(['[~,~,IB] = intersect(' thisVar_in '.varNames,last30yrs_tmp.' thisVar_out '.varNames,''stable'') ;']) ;
                eval([thisVar_in '.maps_YXvsr(:,:,:,:,r) = last30yrs_tmp.' thisVar_out '.maps_YXvs(:,:,IB,:) ;']) ;
            else
                eval([thisVar_in '.maps_YXvsr(:,:,:,:,r) = last30yrs_tmp.' thisVar_out '.maps_YXvs ;']) ;
            end
            eval(['isequal_statList = isequal(' thisVar_in '.statList, last30yrs_tmp.' thisVar_out '.statList) ; ']) ;
            if ~isequal_statList
                eval(['isequal_varNames_afterSort = isequal(sort(' thisVar_in '.statList), sort(last30yrs_tmp.' thisVar_out '.statList)) ; ']) ;
                if ~isequal_varNames_afterSort
                    error(['~isequal_statList (' thisVar_in '). Not fixable.'])
                end
                warning(['~isequal_statList (' thisVar_in '). Fixing.'])
                eval(['[~,~,IB] = intersect(' thisVar_in '.statList,last30yrs_tmp.' thisVar_out '.statList,''stable'') ;']) ;
                eval([thisVar_in '.maps_YXvsr(:,:,:,:,r) = last30yrs_tmp.' thisVar_out '.maps_YXvs(:,:,:,IB) ;']) ;
            else
                eval([thisVar_in '.maps_YXvsr(:,:,:,:,r) = last30yrs_tmp.' thisVar_out '.maps_YXvs ;']) ;
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
    if isfield(last30yrs_tmp,'expyield_last30')
        if r == 1
            maps_yieldExp_last30 = renameStructField(last30yrs_tmp.expyield_last30,'maps_YXvs','maps_YXvsr') ;
        else
            maps_yieldExp_last30.maps_YXvsr(:,:,:,:,r) = last30yrs_tmp.expyield_last30.maps_YXvs ;
        end
    end

    clear firstdec_tmp
    clear lastdec_tmp
    clear last30yrs_tmp

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
land_area_unmasked_YX = land_area_YX ;
land_area_unmasked_weights_YX = land_area_unmasked_YX ./ nansum(nansum(land_area_unmasked_YX)) ;
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


%% Perform secondary calculations

disp('Performing secondary calculations...')

% Adjust for technological improvement, if doing so
if do_adjYieldTech
    
    % Maps
    maps_yield_d1 = adj_yield_tech(maps_yield_d1, yearList_baseline(end-9:end), yearList_future(1:10)) ;
    maps_yield_d9 = adj_yield_tech(maps_yield_d9, yearList_baseline(end-9:end), yearList_future(end-9:end)) ;
    
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
    kcal_per_g = get_kcalDensity(thisCrop) ;
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
        kcal_per_g = get_kcalDensity(thisCrop) ;
        kcal_per_kg = 1e3 * kcal_per_g ;
        eval(['ts_kcalExp2_' thisSuffix ' = ts_kcalExp2_' thisSuffix ' + kcal_per_kg * (' tmp_name{i} '.* ts_croparea_' thisCrop '_yr) ;']) ;
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
    
    % Change in 5th and 95th percentile of peak monthly runoff
    % After Asadieh & Krakauer (2017)
    maps_pk_runoff_d1.maps_p05_YXB = prctile(maps_pk_runoff_d1.maps_YXvyB,5,4) ;
    maps_pk_runoff_d1.maps_p95_YXB = prctile(maps_pk_runoff_d1.maps_YXvyB,95,4) ;
    maps_pk_runoff_d9.maps_p05_YXB = prctile(maps_pk_runoff_d9.maps_YXvyB,5,4) ;
    maps_pk_runoff_d9.maps_p95_YXB = prctile(maps_pk_runoff_d9.maps_YXvyB,95,4) ;
    maps_pk_runoff_d1.maps_p05_YXr = squeeze(prctile(maps_pk_runoff_d1.maps_YXvyr(:,:,:,:,:),5,4)) ;
    maps_pk_runoff_d1.maps_p95_YXr = squeeze(prctile(maps_pk_runoff_d1.maps_YXvyr(:,:,:,:,:),95,4)) ;
    maps_pk_runoff_d9.maps_p05_YXr = squeeze(prctile(maps_pk_runoff_d9.maps_YXvyr(:,:,:,:,:),5,4)) ;
    maps_pk_runoff_d9.maps_p95_YXr = squeeze(prctile(maps_pk_runoff_d9.maps_YXvyr(:,:,:,:,:),95,4)) ;
    below_thresh_YX = mean(maps_awater_d1.maps_YXvyB(:,:,strcmp(maps_awater_d1.varNames,'Runoff'),:),4)/365 < 0.01 ;
    maps_pk_runoff_d1.maps_p05_YXB(below_thresh_YX) = NaN ;
    maps_pk_runoff_d1.maps_p95_YXB(below_thresh_YX) = NaN ;
    maps_pk_runoff_d9.maps_p05_YXB(below_thresh_YX) = NaN ;
    maps_pk_runoff_d9.maps_p95_YXB(below_thresh_YX) = NaN ;
    maps_pk_runoff_d1.maps_p05_YXr(repmat(below_thresh_YX, [1 1 Nruns])) = NaN ;
    maps_pk_runoff_d1.maps_p95_YXr(repmat(below_thresh_YX, [1 1 Nruns])) = NaN ;
    maps_pk_runoff_d9.maps_p05_YXr(repmat(below_thresh_YX, [1 1 Nruns])) = NaN ;
    maps_pk_runoff_d9.maps_p95_YXr(repmat(below_thresh_YX, [1 1 Nruns])) = NaN ;

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
bare_area_YXmean = mean(cat(3,bare_area_YXBH,bare_area_YXr),3) ;
bare_frac_YXmean = bare_area_YXmean ./ gcel_area_YX ;
vegd_frac_YXmean = 1 - bare_frac_YXmean ;
vegd_area_YXmean = vegd_frac_YXmean .* gcel_area_YX ;

disp('Done performing secondary calculations.')


%% Import biodiversity hotspots

disp('Importing biodiversity hotspots...')
hotspot_YX = flipud(imread('/Users/sam/Geodata/BiodiversityHotspotsRevisited_ConservationInternational_2004/data/hotspots_revisited_2004.outerlimit.tif')) ;
hotspot_YX(nanmask) = NaN ;
hotspot_YX = 1==hotspot_YX ;
hotspot_area_YX = hotspot_YX.*gcel_area_YX ;
hotspot_shp = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/hotspots_clipByGridlist.shp' ;

disp('Importing Congolian swamp and lowland forests...')
ecoid_YX = flipud(imread('/Users/sam/Geodata/General/WWF terrestrial ecosystems/wwf_terr_ecos_UnpackClip.halfDeg.ECO_ID.tif')) ;
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
hotspotCSLF_area_YX = (hotspot_YX | cslf_YX).*gcel_area_YX ;
cslf_shp = '/Users/sam/Geodata/General/WWF terrestrial ecosystems/wwf_terr_ecos_UnpackClip.halfDeg.CSLF.shp' ;


%% Import food production units and basins

disp('Importing FPUs...')

fpu_YX = flipud(dlmread('/Users/Shared/PLUM/food_production_units/FPU.asc',' ',6,0)) ;
fpu_YX(fpu_YX==-9999) = NaN ;

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

% Get basin numbers
basin_list = unique(basins_YX(~isnan(basins_YX))) ;
Nbasins = length(basin_list) ;


%% Import population

disp('Importing population...')

pop = readtable('/Users/Shared/PLUM/other_plum_data/SSP/ssp.csv') ;
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
    kcal_per_g = get_kcalDensity(thisCrop) ;
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


disp('Done.')

