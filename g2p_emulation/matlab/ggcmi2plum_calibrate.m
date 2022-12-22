%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Master file for GCMI2PLUM crop calibration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simple version: Get correct production at global level


%% Information about this calibration run

% version_name = 'LPJ-GUESS_emu_remap5e_simple' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_calib/LPJ-GUESS.out' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5e/LU.remapv5e.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5e/cropfracs.remapv5e.txt' ;
% calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types

version_name = 'LPJmL_emu_remap5e_simple' ;   % calib_ver = 18
filename_guess_yield = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_calib/LPJmL.out' ;
filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5e/LU.remapv5e.txt' ;
filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5e/cropfracs.remapv5e.txt' ;
calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types

% version_name = 'pDSSAT_emu_remap5e_simple' ;   % calib_ver = 18
% filename_guess_yield = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_calib/pDSSAT.out' ;
% filename_guess_landuse = '/Users/Shared/PLUM/input/remaps_v5e/LU.remapv5e.txt' ;
% filename_guess_cropfrac = '/Users/Shared/PLUM/input/remaps_v5e/cropfracs.remapv5e.txt' ;
% calib_ver = 18 ;   % The version of mapping FAO to PLUM crop types



%% Other options and setup

% Years for calibration
year1 = 1980 ;
yearN = 2010 ;

% Paths to calibration code and data
dir_code = '/Users/Shared/PLUM/crop_calib_code/' ;
dir_data = '/Users/Shared/PLUM/crop_calib_data/' ;

% Path to figures/tables output directories
dir_outfigs = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_calib/calib_figs/' ;
dir_outtables = '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/outputs_figs/calibration/calibration_tables' ;

% Add code files to path (just for this session)
addpath(genpath(landsymm_lpjg_path()))

is_ggcmi = true ;
need_countries = false ;
script_setup_cropCalibration

if exist(filename_guess_yield, 'file')
    tmp = dir(filename_guess_yield) ;
else
    tmp = dir([filename_guess_yield '.gz']) ;
end
thisModel_lastModDate_str = datestr(tmp.datenum,'yyyymmddHHMMSS') ;


%% Pre-process simulated yields

script_import_lpj_yields_noCCy

emu_area_c = squeeze(mean(nansum(nansum(croparea_lpj_YXcy_comb, 2), 1), 4)) ;
emu_prod_c = squeeze(mean(nansum(nansum(total_lpj_YXcy_comb, 2), 1), 4)) ;


%% Read FAO data

yearList = year1:yearN ;

[~, fao_area, fao_prod, twofiles] = import_FAO_data(...
    calib_ver, year1, yearN, ...
    need_countries) ;
if ~twofiles
    error('Write code to do this with combined FAO file... If that''s really what you want to do.')
end
%%
[listCrops_fa2o, ...
    listCrops_fa2i, listCrops_fa2i_area, listCrops_fa2i_prod, ...
    FAO_to_FAO_key, fao2fao_key_area, fao2fao_key_prod] = ...
    get_FAO_info(calib_ver, twofiles) ;

fao_area_cy = ...
    FAO_to_Ccy2_oneOnly(fao_area, fao2fao_key_area, ...
    listCrops_fa2i_area, listCrops_fa2o, yearList, ...
    'area', true, '') ;
fao_area_c = mean(fao_area_cy, 2) ;

fao_prod_cy = ...
    FAO_to_Ccy2_oneOnly(fao_prod, fao2fao_key_prod, ...
    listCrops_fa2i_prod, listCrops_fa2o, yearList, ...
    'prod', true, 'Production') ;
fao_prod_c = mean(fao_prod_cy, 2) ;

listCrops_fa2o
cf = fao_prod_c ./ emu_prod_c

%% Get output crop lists
listCrops_lpj_comb_all = listCrops_fa2o ;
listCrops_lpj_comb_all{find(strcmp(listCrops_fa2o,'Wheat'))} = 'CerealsC3' ;
listCrops_lpj_comb_all{find(strcmp(listCrops_fa2o,'Maize'))} = 'CerealsC4' ;
listCrops_lpj_comb_all{find(strcmp(listCrops_fa2o,'Starchy roots'))} = 'StarchyRoots' ;
listCrops_corresp_ggcmi = cell(size(listCrops_lpj_comb_all)) ;
for c = 1:length(listCrops_lpj_comb_all)
    thisCrop_calib = listCrops_lpj_comb_all{c} ;
    switch thisCrop_calib
        case 'CerealsC3'    ; listCrops_corresp_ggcmi{c} = 'max_wheat' ;
        case 'Oilcrops'     ; listCrops_corresp_ggcmi{c} = 'soy' ;
        case 'Pulses'       ; listCrops_corresp_ggcmi{c} = 'soy' ;
        case 'CerealsC4'    ; listCrops_corresp_ggcmi{c} = 'maize' ;
        case 'Rice'         ; listCrops_corresp_ggcmi{c} = 'rice' ;
        case 'StarchyRoots' ; listCrops_corresp_ggcmi{c} = 'spring_wheat' ;
        otherwise ; error([thisCrop_calib ' not recognized!'])
    end
end


%% Get super-simple calibration factors


% Create table
table_out = table(shiftdim(listCrops_lpj_comb_all), ...
                  shiftdim(listCrops_corresp_ggcmi), ...
                  cf) ;
table_out.Properties.VariableNames = {'PLUM_crop','GGCMI_crop','calib_factor'} ;


%% Save results table

disp('Exporting results table...')

% Get filename, appending last modified date/time of most recently modified
% file in this model's directory 
% filename_outtable = [dir_outtables thisModel '_' thisModel_lastModDate_str '_cv' num2str(calib_ver) '.csv'] ;
filename_outtable = sprintf('%s/%s_cv%d_%s.csv', dir_outtables, version_name, calib_ver, thisModel_lastModDate_str) ;

% Export
writetable(table_out, filename_outtable) ;

disp('Done.')

