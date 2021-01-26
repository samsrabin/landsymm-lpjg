%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Master file for GCMI2PLUM crop calibration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Information about this calibration run

model_name = 'pDSSAT' ;

remapVer = '20210126' ;
calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types

filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5e/LU.remapv5e.txt' ;
filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5e/cropfracs.remapv5e.txt' ;
filename_countriesMap = 'country_boundaries62892.noNeg99.extrapd.asc' ;


%% Other options and setup

% Get version name
version_name = sprintf('%s_ggcmi_remap%s', model_name, remapVer) ;

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


%% Do it

crop_calibration


%% Save figure

out_file = [dir_outfigs version_name '_v' num2str(calib_ver) '.pdf'] ;
export_fig(out_file,'-r300')


%% Do Miscanthus calibration kludge


