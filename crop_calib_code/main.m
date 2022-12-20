%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Master file for crop calibration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Information about this calibration run

% % 2018-02-16
% verName_calib = 'remap2_PLUM6xtra' ;   % calib_ver = 16
% % filename_guess_yield = '/Volumes/WDMPP_Storage/PLUM/trunk_runs_external/calib.remap2.PLUM6xtra.1901-2005/output-2018-02-16-131139/yield.out.gz' ;
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap2.PLUM6xtra.1901-2005/output-2018-02-16-131139/yield.out.gz' ;
% filename_guess_landuse = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
% filename_guess_cropfrac = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.txt' ;
% calib_ver = 16 ;

% % Ani testing
% verName_calib = 'ani_test' ;
% filename_guess_yield = '/Users/Shared/ani_testing/yield.out' ;
% filename_guess_landuse = '/Users/Shared/ani_testing/lu_1850_2016_hyde32_maxlncr_brazil_0p25_856.txt' ;
% filename_guess_cropfrac = '/Users/Shared/ani_testing/crop_fractions_MIRCA_aggregated_4CNcrops_fromSam_0p25degRes.txt' ;
% calib_ver = 17 ;

% % 2018-10-07
% verName_calib = 'remap5_PLUM6xtra' ;   % calib_ver = 16
% % filename_guess_yield = '/Volumes/WDMPP_Storage/PLUM/trunk_runs_external/calib.remap2.PLUM6xtra.1901-2005/output-2018-02-16-131139/yield.out.gz' ;
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap5.PLUM6xtra.1901-2005/output-2018-10-05-024345/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5/LU.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt.gz' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5/cropfracs.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt.gz' ;
% calib_ver = 16 ;

% % 2018-12-07
% verName_calib = 'remap5_PLUM6xtra_gcArea' ;   % calib_ver = 18
% % filename_guess_yield = '/Volumes/WDMPP_Storage/PLUM/trunk_runs_external/calib.remap2.PLUM6xtra.1901-2005/output-2018-02-16-131139/yield.out.gz' ;
% filename_guess_yield = '/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/calib.remap5.PLUM6xtra.1901-2005/output-2018-10-05-024345/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5/LU.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt.gz' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5/cropfracs.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt.gz' ;
% calib_ver = 18 ;

% % 2019-02-07
% verName_calib = 'remap5d' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap5d.PLUM6xtra.1901-2005/output-2019-02-07-042215/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5d/LU.remapv5d.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5d/cropfracs.remapv5d.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types

% % Jianyong 2019-01-28 with N-fixing
% verName_calib = 'jianyong_20190128_Nfix' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/jianyong_calib.20190128/yield_Nfix.out.gz' ;
% filename_guess_landuse = '/Users/Shared/jianyong_calib.20190128/lu_1850_2017_hyde321_maxlncr_halfdeg.txt.gz' ;
% filename_guess_cropfrac = '/Users/Shared/jianyong_calib.20190128/crop_fractions_fromMIRCA_PLUM_1yr_midpoint_Stijn.txt' ;
% calib_ver = 21 ;   % The version of mapping FAO to PLUM crop types

% % Jianyong 2019-01-28 without N-fixing
% verName_calib = 'jianyong_20190128_noNfix' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/jianyong_calib.20190128/yield_No_Nfix.out.gz' ;
% filename_guess_landuse = '/Users/Shared/jianyong_calib.20190128/lu_1850_2017_hyde321_maxlncr_halfdeg.txt.gz' ;
% filename_guess_cropfrac = '/Users/Shared/jianyong_calib.20190128/crop_fractions_fromMIRCA_PLUM_1yr_midpoint_Stijn.txt' ;
% calib_ver = 21 ;   % The version of mapping FAO to PLUM crop types

% % 2019-02-18
% verName_calib = 'remap5e' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap5e.1901-2005/output-2019-02-18-122738/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5e/LU.remapv5e.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5e/cropfracs.remapv5e.txt' ;
% calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types

% % 2019-02-19
% verName_calib = 'remap8a' ;   % calib_ver = 19
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap8a.1901-2005/output-2019-02-18-112224/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v8a/LU.remapv8a.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v8a/cropfracs.remapv8a.txt' ;
% calib_ver = 19 ;   % The version of mapping FAO to PLUM crop types

% % 2019-02-20
% verName_calib = 'remap8b' ;   % calib_ver = 20
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap8b.1901-2005/output-2019-02-18-102028/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v8b/LU.remapv8b.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v8b/cropfracs.remapv8b.txt' ;
% calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types

% % 2021-02-18
% % Use Oilcrops as a proxy for FruitVeg and Sugar
% verName_calib = 'remap8b_oilProxy' ;
% % filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap5e.1901-2005/output-2019-02-18-122738/yield.out.gz' ;
% filename_guess_yield = '/Users/Shared/PLUM/trunk_runs/calib.remap5d.PLUM6xtra.1901-2005/output-2019-02-07-042215/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v8b/LU.remapv8b.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v8b/cropfracs.remapv8b.txt' ;
% calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types
% oilcrops_proxy_fruitveg_sugar = true ;

% % 2021-03-31
% verName_calib = 'remap8c' ;   % calib_ver = 20
% filename_guess_yield = '/Users/Shared/PLUM/ssp13/calib.remap8c.1901-2005/output-2021-03-26-215749/yield.out.gz' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v8c/LU.remapv8c.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v8c/cropfracs.remapv8c.txt' ;
% calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types

% 2022-10-30
verName_calib = 'remap10_sai-landsymm' ;   % calib_ver = 18
filename_guess_yield = '/Users/Shared/SAI-LandSyMM/output/remap10/calibration/output-2022-10-31-024406/yield_st.out.gz' ;
filename_guess_landuse = '/Users/Shared/SAI-LandSyMM/input/LU/remaps_v10_f09_g17/LU.remapv10_f09_g17.17249.txt' ;
filename_guess_cropfrac = '/Users/Shared/SAI-LandSyMM/input/LU/remaps_v10_f09_g17/cropfracs.remapv10_f09_g17.17249.txt' ;
calib_ver = 23 ;   % The version of mapping FAO to PLUM crop types
remapVer = '10' ;
drop_northpole = true ;
drop_southpole = true ;
lons_centered_on_180 = true ;
in_prec = 5 ; % Actually 6, but rounding errors (?)
filename_countriesMap = 'country_boundaries_f09_g17.asc' ;


%% Other options and setup

if ~exist('drop_northpole', 'var')
    drop_northpole = false ;
end
if ~exist('drop_southpole', 'var')
    drop_southpole = false ;
end
if ~exist('lons_centered_on_180', 'var')
    lons_centered_on_180 = false ;
end
if ~exist('in_prec', 'var')
    in_prec = 2 ;
end

if ~exist('oilcrops_proxy_fruitveg_sugar', 'var')
    oilcrops_proxy_fruitveg_sugar = false ;
end

% Country map file
if ~exist('filename_countriesMap', 'var')
    if strcmp(verName_calib, 'ani_test')
        filename_countriesMap = 'ne_10m_admin_0_countries_ssrIDs_0.25deg.tif' ;
    else
        filename_countriesMap = 'country_boundaries62892.noNeg99.extrapd.asc' ;
    end
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
out_file = [verName_calib '_v' num2str(calib_ver)] ;
out_diary = [dir_outfigs out_file '.txt'] ;
out_figure = [dir_outfigs out_file '.pdf'] ;


%% Do it

crop_calibration


%% Save figure

export_fig(out_figure,'-r300')


%% Do Miscanthus calibration kludge


