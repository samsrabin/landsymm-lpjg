%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Master file for GCMI2PLUM crop calibration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Information about this calibration run

% model_name = 'EPIC-TAMU' ;
% model_name = 'GEPIC' ;
% model_name = 'pDSSAT' ;
model_name = 'LPJ-GUESS-sim' ;
% lpjg_run_topDir = '/Volumes/Reacher/G2P/outputs_LPJG/remap12_2016/calibration' ;
% lpjg_run_topDir = '/Volumes/Reacher/G2P/outputs_LPJG/remap12_2016/calibration_Ks1' ;
% lpjg_run_topDir = '/Volumes/Reacher/G2P/outputs_LPJG/remap12_2016/calibration_Ks2' ;
% lpjg_run_topDir = '/Volumes/Reacher/G2P/outputs_LPJG/remap12_2016/calibration_Ks3' ;
% lpjg_run_topDir = '/Volumes/Reacher/G2P/outputs_LPJG/remap12_2016/calibration_Ks3_Oilcrops_ownDates' ;
lpjg_run_topDir = '/Volumes/Reacher/G2P/outputs_LPJG/remap12_2016/calibration_Ks3_Oilcrops_ownDates_TeSo' ;

% remapVer = '5e' ; calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types
% remapVer = '8b' ; calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types
% remapVer = '9_g2p' ; calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types
% remapVer = '10_g2p' ; calib_ver = 20 ;   % The version of mapping FAO to PLUM crop types
% remapVer = '11_g2p' ; calib_ver = 23 ;   % The version of mapping FAO to PLUM crop types
% remapVer = '13_g2p' ; calib_ver = 23 ;   % The version of mapping FAO to PLUM crop types
remapVer = '12_g2p' ; calib_ver = 23 ;   % The version of mapping FAO to PLUM crop types

ctry_excluded_area_thresh = 0.1 ; % The fraction of a country's excluded
% area of a given crop (due to no simulated yield) above which the country
% will be excluded from that crop's calibration.
ctrymapVer = 1 ;

% Compare individual years? Or take mean over period?
indiv_years = false ;

% Only doing this for g2p runs!
do_remove_area_dueto_NaNsim = true ;


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
if any(strcmp(remapVer, {'12_g2p', '13_g2p'}))
    % Stripped-down crop list used for running. Complete crop list needed
    % for calibration.
    tmp = sprintf('v%s', strrep(remapVer, '_g2p', '')) ;
    filename_guess_cropfrac = strrep(filename_guess_cropfrac, ...
        tmp, 'v11') ;
    clear tmp
end
% filename_guess_sugars = sprintf( ...
%     '%s/sugar.mat', ...
%     LUdir);

% Get yield file, if LPJ-GUESS simulation. Otherwise, get path to emulated 
% baseline outputs for whatever N inputs Christoph used.
clear filename_guess_yield
specified_topdir = false ;
if strcmp(model_name, 'LPJ-GUESS-sim')
    specified_topdir = exist('lpjg_run_topDir', 'var') ;
    if ~specified_topdir
        lpjg_run_topDir = sprintf('/Volumes/Reacher/G2P/outputs_LPJG/remap%s/calibration', ...
            strrep(remapVer, '_g2p', '')) ;
    end
    if ~exist(lpjg_run_topDir, 'dir')
        error('lpjg_run_topDir not found: %s', lpjg_run_topDir)
    end
    lpjg_run_outDirs = dir(sprintf('%s/output-*', lpjg_run_topDir)) ;
    if isempty(lpjg_run_outDirs)
        error('No output directories found in %s', lpjg_run_topDir)
    end
    filename_guess_yield = sprintf('%s/%s/yield_plantyear_st.out', ...
        lpjg_run_topDir, lpjg_run_outDirs(end).name) ;
else
    dirname_emuBL_yields = sprintf( ...
        '/Volumes/Reacher/GGCMI/AgMIP.output/CMIP_emulated/yields/CMIP6/A1/ssp126/%s', ...
        model_name) ;
    if ~exist(dirname_emuBL_yields, 'dir')
        error('Phase 2 yield directory not found (%s)', dirname_emuBL_yields)
    end
    adaptation = 1 ;
end

oilcrops_proxy_fruitveg_sugar = false ;

% Get countries map
if ctrymapVer == 1
    filename_countriesMap = fullfile(landsymm_lpjg_path(), 'data', 'geodata', 'country_boundaries', 'country_boundaries62892.noNeg99.extrapd.asc') ;
else
    error('ctrymapVer %d not recognized', ctrymapVer)
end

% Get version name
verName_calib = sprintf('ggcmi_%s_%s%d_ctry%d_caet%g', ...
    model_name, remapVer, calib_ver, ctrymapVer, ctry_excluded_area_thresh) ;

% Years for calibration
if indiv_years
    year1 = 1995 ;
    yearN = 2005 ;
else
    year1 = 1980 ;
    yearN = 2010 ;
end

% Paths to calibration code and data
dir_code = '/Users/Shared/PLUM/crop_calib_code/' ;
dir_data = '/Users/Shared/PLUM/crop_calib_data/' ;

% Path to figure output directory
if specified_topdir
    dir_outfigs = [lpjg_run_topDir '/outputs/'] ;
else
    dir_outfigs = '/Volumes/Reacher/G2P/outputs_calibration/' ;
end
if ~exist(dir_outfigs, 'dir')
    mkdir(dir_outfigs)
end

% Add code files to path (just for this session)
addpath(genpath(landsymm_lpjg_path()))

regression_type = 'slope-only' ;

% Set up figure and diary files
out_file = [verName_calib '_v' num2str(calib_ver)] ;
out_diary = [dir_outfigs out_file '.txt'] ;
out_csv = [dir_outfigs out_file '.csv'] ;
out_figure = [dir_outfigs out_file '.pdf'] ;


%% Do it

crop_calibration


%% Save figure and table

fprintf('Saving outputs to %s...\n', dir_outfigs)

export_fig(out_figure,'-r300')
close

T_out = table(shiftdim(listCrops_fa2o), shiftdim(round(calib_factors_u,3)), ...
    'VariableNames', {'Crop', 'calibration_factor'}) ;
writetable(T_out, out_csv)

disp('All done!')


%% Do Miscanthus calibration kludge


