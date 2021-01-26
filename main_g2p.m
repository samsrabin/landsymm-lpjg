%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Master file for GCMI2PLUM crop calibration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Information about this calibration run

% version_name = 'LPJ-GUESS_ggcmi_remap5e' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_calib/LPJ-GUESS.out' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5e/LU.remapv5e.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5e/cropfracs.remapv5e.txt' ;
% filename_countriesMap = 'country_boundaries62892.noNeg99.extrapd.asc' ;
% calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types

% version_name = 'LPJmL_ggcmi_remap5e' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_calib/LPJmL.out' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5e/LU.remapv5e.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5e/cropfracs.remapv5e.txt' ;
% filename_countriesMap = 'country_boundaries62892.noNeg99.extrapd.asc' ;
% calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types

version_name = 'pDSSAT_ggcmi_remap5e' ;   % calib_ver = 18
filename_guess_yield = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_calib/pDSSAT.out' ;
filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5e/LU.remapv5e.txt' ;
filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5e/cropfracs.remapv5e.txt' ;
filename_countriesMap = 'country_boundaries62892.noNeg99.extrapd.asc' ;
calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types



%% Other options and setup

% Years for calibration
year1 = 1980 ;
yearN = 2010 ;

% Paths to calibration code and data
dir_code = '/Users/Shared/PLUM/crop_calib_code/' ;
dir_data = '/Users/Shared/PLUM/crop_calib_data/' ;

% Path to figure output directory
dir_outfigs = '/Volumes/Reacher/G2P/outputs_calibration/' ;

% Add code files to path (just for this session)
addpath(genpath(dir_code))

regression_type = 'slope-only' ;


%% Do it

crop_calibration


%% Save figure

out_file = [dir_outfigs version_name '_v' num2str(calib_ver) '.pdf'] ;
export_fig(out_file,'-r300')


%% Do Miscanthus calibration kludge


