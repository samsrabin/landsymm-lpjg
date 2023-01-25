%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Master file for crop calibration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSR 2022-12-28
% - Given outputs from an LPJ-GUESS run, a way to map LPJ-GUESS/PLUM crop types to FAO
%   crop types, and some other info, this will generate a figure with the "calibration
%   factors" that PLUM will multiply the LPJ-GUESS outputs by.
% - Calculates Miscanthus calibration factors, but doesn't use current LPJ-GUESS outputs
%   if calib_ver < 24.

addpath(genpath(landsymm_lpjg_path()))

% get_calibration_factors_options.m must be somewhere on your path.
% There, specify the following variables:
%     year1, yearN: The first and last years, respectively, to be used for calibration.
%                   Applies to both LPJ-GUESS outputs and FAO data.
%                   E.g.: year1 = 1995; yearN = 2005;
%     dir_data: The path to the directory containing various general files needed by the
%               calibration scripts. 
%               E.g.: dir_data = '/Users/Shared/PLUM/crop_calib_data/' ;
%     dir_outfigs: The path to the directory where you want the output figure to go.
%                  E.g.: dir_outfigs = '/Users/Shared/PLUM/crop_calibration_for_PLUM/figures/' ;
%     verName_calib: Mostly just used for naming files output by calibration scripts.
%                    Can also be used to set version-specific processing options.
%                    E.g.: verName_calib = 'recal_2022-12' ;
%     filename_guess_yield: The LPJ-GUESS output yields you'll be calibrating with. Always 
%                           use yield_st.out if it's present.
%                           E.g.: filename_guess_yield = '/Users/Shared/PLUM/recalibration_2022-12/remap10_yp/calibration/output-2022-12-20-150542/yield_st.out.gz' ;
%     filename_guess_landuse: The land use file to be used for area-weighted averaging in
%                             calibration. This is usually the same as what was used in
%                             the actual runs, but doesn't need to be. It just needs to
%                             have columns for all the crops in the LPJ-GUESS output.
%                             E.g.: filename_guess_landuse = '/Users/Shared/LandSyMM/inputs/LU/remaps_v10_old_62892_gL/LU.remapv10_old_62892_gL.txt' ;
%     filename_guess_cropfrac: As filename_guess_landuse, but for crop fractions instead
%                              of general land use fractions.
%                              E.g.: filename_guess_cropfrac = '/Users/Shared/LandSyMM/inputs/LU/remaps_v10_old_62892_gL/cropfracs.remapv10_old_62892_gL.txt' ;
%     calib_ver: The version of mapping FAO to PLUM crop types. Also sets some other 
%                settings in crop calibration scripts, including FAO files  to use. 
%                Mappings and files are defined in get_FAO_info.m.
%                E.g.: calib_ver = 24 ;
%     remapVer: Used for setting the list of LPJ-GUESS output crops being calibrated.
%               See top of script_import_lpj_yields_noCCy.m for where this is used.
%               E.g.: remapVer = '10' ;

get_calibration_factors_options


%% Setup

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

% Do slope-only regression. There's another option but it should NOT be used.
regression_type = 'slope-only' ;

% SSR 2023-01-25: Setting this to true can cause pretty big differences. Previous behavior
% (i.e., before merging in ggcmi2plum branch with commit 2b6d60e) was false.
do_remove_area_dueto_NaNsim = false ;

% Set up figure and diary files
out_file = [verName_calib '_v' num2str(calib_ver)] ;
out_diary = [dir_outfigs out_file '.txt'] ;
out_figure = [dir_outfigs out_file '.pdf'] ;


%% Do it

crop_calibration


%% Save figure

fprintf('Saving figure to %s...\n', out_figure)
export_fig(out_figure,'-r300')
close
disp('Done.')

