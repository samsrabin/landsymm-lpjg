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