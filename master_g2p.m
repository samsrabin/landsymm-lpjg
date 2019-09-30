%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Master file for GCMI2PLUM crop calibration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Information about this calibration run

% 2019-09-30
version_name = 'remap5e' ;   % calib_ver = 18
filename_emu_yield = '/Users/Shared/GGCMI2PLUM/emulator/Sam/outputs_calib/LPJ-GUESS.out' ;
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
dir_outfigs = '/Users/Shared/GGCMI2PLUM/emulator/Sam/outputs_calib/calib_figs/' ;

% Add code files to path (just for this session)
addpath(genpath(dir_code))


%% Do it

crop_calibration


%% Save figure

out_file = [dir_outfigs version_name '_v' num2str(calib_ver) '.pdf'] ;
export_fig(out_file,'-r300')


%% Notes

% calib_ver = 1 ;    % As Peter Alexander described
% calib_ver = 2 ;    % Looking at each crop separately. NOT ACTUALLY USEFUL FOR CALIBRATION
% calib_ver = 3 ;    % As 1, but with oil+cakes equivalent
% calib_ver = 4 ;    % As 1, but using oil crops as in Commodity Balance "Oilcrops"
% calib_ver = 5 ;    % As 4, but using Commodity Balance for Production
% calib_ver = 6 ;    % As 1, but with pre-processing as 5
% calib_ver = 7 ;    % As 6, but with oil+cakes equivalent
% calib_ver = 8 ;    % As 6, but with updated FAOSTAT data (2017-02-17)
% calib_ver = 9 ;    % As 5, but with Peter's old FAOSTAT data (2015-12)
% calib_ver = 10 ;   % As 9, but without rye
% calib_ver = 11 ;   % As 9, but with Stijn's new crops
% calib_ver = 12 ;   % As 10, but with PLUM6
% calib_ver = 13 ;   % As 9, but with only exact crop matches (for pure6 runs)
% calib_ver = 14 ;   % As 13, but without Soybeans or Pulses (for pure4 runs)
% calib_ver = 15 ;   % ~9 for FAO data, plus Rye, and with PLUM6 types (PLUM6xtra)
% calib_ver = 16 ;   % As 15, but correctly with Rye instead of Oats in Products


%% Do Miscanthus calibration kludge


