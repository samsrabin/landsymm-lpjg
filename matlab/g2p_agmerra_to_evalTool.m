topDir = '/Volumes/Reacher/G2P/outputs_LPJG/remap12_2016' ;
dirList = { ...
    'calibration_Ks3/output-2021-05-03-225549' ;
    'calibration_Ks3_Oilcrops_ownDates/output-2021-05-03-230307' ;
    'calibration_Ks3_Oilcrops_ownDates_TeSo/output-2021-05-03-225929' ;
    'calibration_Ks3_Oilcrops_TeSWfert/output-2021-05-04-172927' ;
    'calibration_Ks3_Oilcrops_ownDates_TeSWfert/output-2021-05-04-172648' ;
    'calibration_Ks3_Oilcrops_ownDates_TeSWfert_TeSo/output-2021-05-04-172908' ;
    } ;

% Setup
addpath(genpath('/Users/Shared/GGCMI/inputs/phase3/ISIMIP3/_MATLAB_ISIMIP3'))

% For evaluation tool
runInfo.phase = 'EV' ;
runInfo.climate_forcing = 'AgMERRA' ;
runInfo.sens_scenario = 'default' ;
runInfo.climate_scenario = 'obsclim' ;
runInfo.yearList_out = 1980:2010 ;
runInfo.soc_scenario = 'socscen' ;
do_checks = false ;
runInfo.gdd_thresh = -1 ;
runInfo = get_more_runInfo(runInfo) ;

for d = 1:length(dirList)
    file_yield = sprintf('%s/%s/yield_plantyear_st.out', topDir, dirList{d}) ;
    dir_out = sprintf('%s/%s/forEvalTool', topDir, dirList{d}) ;
    process_from_yieldfile(file_yield, runInfo, dir_out)
end
disp('All done!')