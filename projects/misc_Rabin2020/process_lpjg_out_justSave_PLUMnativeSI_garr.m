%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Work with LPJ-GUESS outputs for impacts paper: %%%
%%% PLUM-style native %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% All TRUE except really unnecessary ones
% do_save.LU              = true ;
% do_save.crops           = true ;
% do_save.yield           = true ;
% do_save.yield_exp       = false ;
% do_save.yield_map       = true ;
% do_save.yield_exp_map   = false ;
% do_save.irrig           = true ;
% do_save.water           = true ;
% do_save.carbon          = true ;
% do_save.mrunoff         = true ;
% do_save.albedo          = true ;
% do_save.bvocs           = true ;
% do_save.Nflux           = true ;
% do_save.Nfert           = true ;
% do_save.fpc             = false ;

%%% All FALSE except selected
do_save.LU              = false ;
do_save.crops           = false ;
do_save.yield           = true ;
do_save.yield_exp       = true ;
do_save.yield_map       = true ;
do_save.yield_exp_map   = true ;
do_save.irrig           = false ;
do_save.water           = false ;
do_save.carbon          = false ;
do_save.mrunoff         = false ;
do_save.albedo          = false ;
do_save.bvocs           = false ;
do_save.Nflux           = false ;
do_save.Nfert           = false ;
do_save.fpc             = false ;

% Determine which system you're on and set up.
thisSystem = get_system_name() ;
if ~any(strcmp(thisSystem, {'ssr_mac', 'ssr_keal', 'ssr_uc2'}))
    error('thisSystem not recognized: %s (is thisSystem needed anymore?)', thisSystem)
end
addpath(genpath(landsymm_lpjg_path()))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define directories to process %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Rabin et al. (2020)
% inDir_list = {...
% %     'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851'
%     'LPJGPLUM_2011-2100_harm3_SSP1_RCP45/output-2019-02-27-103914';
% %     'LPJGPLUM_2011-2100_harm3_SSP3_RCP60/output-2019-02-27-093027';
% %     'LPJGPLUM_2011-2100_harm3_SSP4_RCP60/output-2019-02-27-093259';
% %     'LPJGPLUM_2011-2100_harm3_SSP5_RCP85/output-2019-02-27-104120';
% %     'LPJGPLUM_2011-2100_harm3_constLU_RCP85/output-2019-03-07-164546';
%     } ;
% gridlist_file = 'PLUMout_gridlist.txt';
% calib_name = 'remap5e_v18' ;

%%% ssp13
if strcmp(thisSystem, 'ssr_mac')
    tmp = '/Users/Shared/PLUM/ssp13/' ;
elseif strcmp(thisSystem, 'ssr_uc2')
    tmp = '/home/kit/imk-ifu/lr8247/PLUM/outputs/ssp13/' ;
else
    error('thisSystem not recognized: %s', thisSystem)
end
inDir_list = { ...
   'LPJGPLUM_1850-2010_remap8c' ;
   'LPJGPLUM_2011-2100_ssp13-1' ;
   'LPJGPLUM_2011-2100_ssp13-2' ;
   'LPJGPLUM_2011-2100_ssp13-3' ;
   'LPJGPLUM_2011-2100_ssp13-4' ;
   'LPJGPLUM_2011-2100_ssp13-5' ;
   } ;
inDir_list = strcat(tmp, inDir_list) ;
gridlist_file = 'gridlist_57790.txt' ;
calib_name = 'remap8b_v20.exclOutliers_10iqr' ;

% Finish up
calib_file = sprintf('%s/%s.csv', paper02_repo_path, calib_name) ;
if ~exist(calib_file, 'file')
    error('calib_file does not exist: %s', calib_file)
end


%% Loop through inDir_list

loop_through_indir_list_garr
