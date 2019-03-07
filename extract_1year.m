%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract one year of LU inputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisYear = 2010 ;

% file_lu_in = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/LU.remapv6.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% file_cf_in = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/cropfracs.remapv6.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% % file_nf_in = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/nfert.remapv6.20180214.ecFertIrr0.setaside0103.m4.txt' ;

% file_lu_in = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/LU.remapv6.20180214.ecFertIrr0.setaside0103.m4.someOfEachCrop.txt' ;
% file_cf_in = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/cropfracs.remapv6.20180214.ecFertIrr0.setaside0103.m4.someOfEachCrop.txt' ;
% % % file_nf_in = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/nfert.remapv6.20180214.ecFertIrr0.setaside0103.m4.txt' ;

% file_lu_in = '/Users/Shared/PLUM/input/remaps_v6p3/LU.remapv6p3.txt' ;
% file_cf_in = '/Users/Shared/PLUM/input/remaps_v6p3/cropfracs.remapv6p3.txt' ;
% file_nf_in = '/Users/Shared/PLUM/input/remaps_v6p3/nfert.remapv6p3.txt' ;

% file_lu_in = '/Users/Shared/PLUM/input/remaps_v6p3/LU.remapv6p3.someOfEachCrop.txt' ;
% file_cf_in = '/Users/Shared/PLUM/input/remaps_v6p3/cropfracs.remapv6p3.someOfEachCrop.txt' ;
% file_nf_in = '/Users/Shared/PLUM/input/remaps_v6p3/nfert.remapv6p3.txt' ;

% file_lu_in = '/Users/Shared/PLUM/input/remaps_v6p3/LU.remapv6p3.someOfEachCrop.txt' ;
% file_cf_in = '/Users/Shared/PLUM/input/remaps_v6p3/cropfracs.remapv6p3.someOfEachCrop.txt' ;
% file_nf_in = '/Users/Shared/PLUM/input/remaps_v6p3/nfert.remapv6p3.txt' ;

file_lu_in = '/Users/Shared/PLUM/input/remaps_v6p7/LU.remapv6p7.someOfEachCrop.txt' ;
file_cf_in = '/Users/Shared/PLUM/input/remaps_v6p7/cropfracs.remapv6p7.someOfEachCrop.txt' ;
file_nf_in = '/Users/Shared/PLUM/input/remaps_v6p7/nfert.remapv6p7.txt' ;



%% Import

lu = lpjgu_matlab_readTable(file_lu_in) ;
cf = lpjgu_matlab_readTable(file_cf_in) ;
nf = lpjgu_matlab_readTable(file_nf_in) ;
disp('Done')


%% Restrict

lu_out = lu(lu.Year==thisYear,~strcmp(lu.Properties.VariableNames,'Year')) ;
cf_out = cf(cf.Year==thisYear,~strcmp(cf.Properties.VariableNames,'Year')) ;
nf_out = nf(nf.Year==thisYear,~strcmp(nf.Properties.VariableNames,'Year')) ;


%% Save

%%% Options %%%%%%%%%%%%%
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = false ;
fancy = false ;
%%%%%%%%%%%%%%%%%%%%%%%%%

lu_out_array = table2array(lu_out) ;
cf_out_array = table2array(cf_out) ;
nf_out_array = table2array(nf_out) ;

file_lu_out = strrep(file_lu_in,'.txt',sprintf('.%d.txt',thisYear)) ;
file_cf_out = strrep(file_cf_in,'.txt',sprintf('.%d.txt',thisYear)) ;
file_nf_out = strrep(file_nf_in,'.txt',sprintf('.%d.txt',thisYear)) ;

disp('Saving LU...')
lpjgu_matlab_saveTable(lu_out.Properties.VariableNames, lu_out_array, file_lu_out,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20) ;

disp('Saving cropfracs...')
lpjgu_matlab_saveTable(cf_out.Properties.VariableNames, cf_out_array, file_cf_out,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20) ;

disp('Saving nfert...')
lpjgu_matlab_saveTable(nf_out.Properties.VariableNames, nf_out_array, file_nf_out,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20) ;

