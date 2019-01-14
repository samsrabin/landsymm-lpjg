%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create constant-2010 LU and cropfrac files %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import originals

file_lu_in = '/Users/Shared/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
file_nf_in = '/Users/Shared/PLUM/input/Nfert/LUH2/nfert_1700_2015_luh2_aggregate_sum2x2_midpoint_rescaled_v20.txt' ;

lu_in = lpjgu_matlab_readTable(file_lu_in) ;
nf_in = lpjgu_matlab_readTable(file_nf_in) ;


%% Make output datasets (only 2010)

lu_out = lu_in ;
lu_out(lu_out.Year~=2010,:) = [] ;
lu_out(:,strcmp(lu_out.Properties.VariableNames,'Year')) = [] ;

nf_out = nf_in ;
nf_out(nf_out.Year~=2010,:) = [] ;
nf_out(:,strcmp(nf_out.Properties.VariableNames,'Year')) = [] ;


%% Save output datasets

%%% Options %%%%%%%%%%%%%
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('Saving LU...')
% file_lu_out = '/Users/Shared/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.2010only.txt' ;
% lpjgu_matlab_saveTable(lu_out.Properties.VariableNames, lu_out, file_lu_out,...
%     'outPrec', outPrec, ...
%     'outWidth', outWidth, ...
%     'delimiter', delimiter, ...
%     'overwrite', overwrite, ...
%     'fancy', fancy, ...
%     'progress_step_pct', 20) ;

disp('Saving nfert...')
file_nf_out = '/Users/Shared/PLUM/input/Nfert/LUH2/nfert_2010_luh2_aggregate_sum2x2_midpoint_rescaled_v20.txt' ;
lpjgu_matlab_saveTable(nf_out.Properties.VariableNames, nf_out, file_nf_out,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20) ;

