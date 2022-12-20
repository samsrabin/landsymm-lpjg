function process_lpjg_out_justSave_PLUMnativeSI_SCCfn( ...
    inDir_list, calib_file)

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
do_save.yield           = false ;
do_save.yield_exp       = false ;
do_save.yield_map       = false ;
do_save.yield_exp_map   = false ;
do_save.irrig           = false ;
do_save.water           = false ;
do_save.carbon          = false ;
do_save.mrunoff         = false ;
do_save.albedo          = true ;
do_save.bvocs           = false ;
do_save.Nflux           = false ;
do_save.Nfert           = false ;
do_save.fpc             = false ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop through inDir_list %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loop_through_indir_list

end
