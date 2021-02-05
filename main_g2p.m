%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Master file for GCMI2PLUM crop calibration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Information about this calibration run

% model_name = 'EPIC-TAMU' ;
% model_name = 'GEPIC' ;
model_name = 'pDSSAT' ;

% remapVer = '5e' ; calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types
% remapVer = '8b' ; calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types
% remapVer = '9_g2p' ; calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types
remapVer = '10_g2p' ; calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types

ctry_excluded_area_thresh = 0.1 ; % The fraction of a country's excluded
% area of a given crop (due to no simulated yield) above which the country
% will be excluded from that crop's calibration.
ctrymapVer = 1 ;


%% Other options and setup

% Get LU files
if contains(remapVer, '_g2p')
    LUdir = sprintf('/Volumes/Reacher/G2P/inputs/LU/remaps_v%s', ...
        remapVer) ;
else
    LUdir = sprintf('/Users/Shared/PLUM/input/remaps_v%s', ...
        remapVer) ;
end
filename_guess_landuse = sprintf( ...
    '%s/LU.remapv%s.txt', ...
    LUdir, remapVer);
filename_guess_cropfrac = sprintf( ...
    '%s/cropfracs.remapv%s.txt', ...
    LUdir, remapVer);
filename_guess_sugars = sprintf( ...
    '%s/sugar.mat', ...
    LUdir);

% Get countries map
if ctrymapVer == 1
    filename_countriesMap = 'country_boundaries62892.noNeg99.extrapd.asc' ;
else
    error('ctrymapVer %d not recognized', ctrymapVer)
end

% Get version name
version_name = sprintf('ggcmi_%s_%s%d_ctry%d_caet%g', ...
    model_name, remapVer, calib_ver, ctrymapVer, ctry_excluded_area_thresh) ;

% Path to emulated baseline outputs for whatever N inputs Christoph used
dirname_emuBL_yields = sprintf( ...
    '/Volumes/Reacher/GGCMI/AgMIP.output/CMIP_emulated/yields/CMIP6/A1/ssp126/%s', ...
    model_name) ;
if ~exist(dirname_emuBL_yields, 'dir')
    error('Phase 2 yield directory not found (%s)', dirname_emuBL_yields)
end
adaptation = 1 ;

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
addpath(genpath('/Users/sam/Documents/git_repos/g2p_emulation/matlab/emu2plum'))

regression_type = 'slope-only' ;

% Set up figure and diary files
out_file = [version_name '_v' num2str(calib_ver)] ;
out_diary = [dir_outfigs out_file '.txt'] ;
out_figure = [dir_outfigs out_file '.pdf'] ;


%% Do it

crop_calibration


%% Save figure

export_fig(out_figure,'-r300')
close

disp('All done!')


%% Do Miscanthus calibration kludge


