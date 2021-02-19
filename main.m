%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Master file for crop calibration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Information about this calibration run

% % 2018-02-16
% version_name = 'remap2_PLUM6xtra' ;   % calib_ver = 16
% % filename_guess_yield = '/Volumes/WDMPP_Storage/PLUM/trunk_runs_external/calib.remap2.PLUM6xtra.1901-2005/output-2018-02-16-131139/yield.out.gz' ;
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap2.PLUM6xtra.1901-2005/output-2018-02-16-131139/yield.out.gz' ;
% filename_guess_landuse = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
% filename_guess_cropfrac = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.txt' ;
% calib_ver = 16 ;

% % Ani testing
% version_name = 'ani_test' ;
% filename_guess_yield = '/Users/Shared/ani_testing/yield.out' ;
% filename_guess_landuse = '/Users/Shared/ani_testing/lu_1850_2016_hyde32_maxlncr_brazil_0p25_856.txt' ;
% filename_guess_cropfrac = '/Users/Shared/ani_testing/crop_fractions_MIRCA_aggregated_4CNcrops_fromSam_0p25degRes.txt' ;
% calib_ver = 17 ;

% % 2018-10-07
% version_name = 'remap5_PLUM6xtra' ;   % calib_ver = 16
% % filename_guess_yield = '/Volumes/WDMPP_Storage/PLUM/trunk_runs_external/calib.remap2.PLUM6xtra.1901-2005/output-2018-02-16-131139/yield.out.gz' ;
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap5.PLUM6xtra.1901-2005/output-2018-10-05-024345/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5/LU.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt.gz' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5/cropfracs.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt.gz' ;
% calib_ver = 16 ;

% % 2018-12-07
% version_name = 'remap5_PLUM6xtra_gcArea' ;   % calib_ver = 18
% % filename_guess_yield = '/Volumes/WDMPP_Storage/PLUM/trunk_runs_external/calib.remap2.PLUM6xtra.1901-2005/output-2018-02-16-131139/yield.out.gz' ;
% filename_guess_yield = '/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/calib.remap5.PLUM6xtra.1901-2005/output-2018-10-05-024345/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5/LU.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt.gz' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5/cropfracs.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt.gz' ;
% calib_ver = 18 ;

% % 2019-02-07
% version_name = 'remap5d' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap5d.PLUM6xtra.1901-2005/output-2019-02-07-042215/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5d/LU.remapv5d.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5d/cropfracs.remapv5d.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types

% % Jianyong 2019-01-28 with N-fixing
% version_name = 'jianyong_20190128_Nfix' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/jianyong_calib.20190128/yield_Nfix.out.gz' ;
% filename_guess_landuse = '/Users/Shared/jianyong_calib.20190128/lu_1850_2017_hyde321_maxlncr_halfdeg.txt.gz' ;
% filename_guess_cropfrac = '/Users/Shared/jianyong_calib.20190128/crop_fractions_fromMIRCA_PLUM_1yr_midpoint_Stijn.txt' ;
% calib_ver = 21 ;   % The version of mapping FAO to PLUM crop types

% % Jianyong 2019-01-28 without N-fixing
% version_name = 'jianyong_20190128_noNfix' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/jianyong_calib.20190128/yield_No_Nfix.out.gz' ;
% filename_guess_landuse = '/Users/Shared/jianyong_calib.20190128/lu_1850_2017_hyde321_maxlncr_halfdeg.txt.gz' ;
% filename_guess_cropfrac = '/Users/Shared/jianyong_calib.20190128/crop_fractions_fromMIRCA_PLUM_1yr_midpoint_Stijn.txt' ;
% calib_ver = 21 ;   % The version of mapping FAO to PLUM crop types

% % 2019-02-18
% version_name = 'remap5e' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap5e.1901-2005/output-2019-02-18-122738/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5e/LU.remapv5e.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5e/cropfracs.remapv5e.txt' ;
% calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types

% % 2019-02-19
% version_name = 'remap8a' ;   % calib_ver = 19
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap8a.1901-2005/output-2019-02-18-112224/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v8a/LU.remapv8a.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v8a/cropfracs.remapv8a.txt' ;
% calib_ver = 19 ;   % The version of mapping FAO to PLUM crop types

% % 2019-02-20
% version_name = 'remap8b' ;   % calib_ver = 20
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap8b.1901-2005/output-2019-02-18-102028/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v8b/LU.remapv8b.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v8b/cropfracs.remapv8b.txt' ;
% calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types

% 2021-02-18
% Use Oilcrops as a proxy for FruitVeg and Sugar
version_name = 'remap8b_oilProxy' ;
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap5e.1901-2005/output-2019-02-18-122738/yield.out.gz' ;
filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap5d.PLUM6xtra.1901-2005/output-2019-02-07-042215/yield.out.gz' ;
filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v8b/LU.remapv8b.txt' ;
filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v8b/cropfracs.remapv8b.txt' ;
calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types
oilcrops_proxy_fruitveg_sugar = true ;


%% Other options and setup

if ~exist('oilcrops_proxy_fruitveg_sugar', 'var')
    oilcrops_proxy_fruitveg_sugar = false ;
end
fake_ggcmi_sugar = false ;
fake_ggcmi_oilcrops = false ;

% Country map file
if strcmp(version_name, 'ani_test')
    filename_countriesMap = 'ne_10m_admin_0_countries_ssrIDs_0.25deg.tif' ;
else
    filename_countriesMap = 'country_boundaries62892.noNeg99.extrapd.asc' ;
end

% Years for calibration
year1 = 1995 ;
yearN = 2005 ;

% Paths to calibration code and data
dir_code = '/Users/Shared/PLUM/crop_calib_code/' ;
dir_data = '/Users/Shared/PLUM/crop_calib_data/' ;

% Path to figure output directory
dir_outfigs = '/Users/Shared/PLUM/crop_calibration_for_PLUM/figures/' ;

% Add code files to path (just for this session)
addpath(genpath(dir_code))

% Do slope-only regression
regression_type = 'slope-only' ;

% Set up figure and diary files
out_file = [version_name '_v' num2str(calib_ver)] ;
out_diary = [dir_outfigs out_file '.txt'] ;
out_figure = [dir_outfigs out_file '.pdf'] ;


%% Do it

crop_calibration


%% Save figure

export_fig(out_figure,'-r300')


%% Do Miscanthus calibration kludge


