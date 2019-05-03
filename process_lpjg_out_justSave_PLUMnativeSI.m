%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Work with LPJ-GUESS outputs for impacts paper: %%%
%%% PLUM-style native %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calib_file = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/MATLAB_work/remap5e_v18.csv' ;

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
do_save.yield           = false ;
do_save.yield_exp       = false ;
do_save.yield_map       = false ;
do_save.yield_exp_map   = false ;
do_save.irrig           = false ;
do_save.water           = true ;
do_save.carbon          = false ;
do_save.mrunoff         = true ;
do_save.albedo          = false ;
do_save.bvocs           = false ;
do_save.Nflux           = true ;
do_save.Nfert           = false ;
do_save.fpc             = false ;

inDir_list = {...
%     'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851'
    'LPJGPLUM_2011-2100_harm3_SSP1_RCP45/output-2019-02-27-103914';
%     'LPJGPLUM_2011-2100_harm3_SSP3_RCP60/output-2019-02-27-093027';
%     'LPJGPLUM_2011-2100_harm3_SSP4_RCP60/output-2019-02-27-093259';
%     'LPJGPLUM_2011-2100_harm3_SSP5_RCP85/output-2019-02-27-104120';
%     'LPJGPLUM_2011-2100_harm3_constLU_RCP85/output-2019-03-07-164546';
    } ;


%% Loop through inDir_list

loop_through_indir_list
