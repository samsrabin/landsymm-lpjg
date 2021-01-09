%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert himelo0 file to himelo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file_in = '/Users/Shared/PLUM/input/remaps_v6/cropfracs.remapv6.20180214.ecFertIrr0.setaside0103.m4.himelo0.txt' ;
% out_file = '/Users/Shared/PLUM/input/remaps_v6/cropfracs.remapv6.20180214.ecFertIrr0.setaside0103.m4.himelo1yr.txt' ;

file_in = '/Users/Shared/PLUM/input/remaps_v6p6/cropfracs.remapv6p6.himelo0.txt.gz' ;
out_file = '/Users/Shared/PLUM/input/remaps_v6p6/cropfracs.remapv6p6.himelo1yr.txt.gz' ;


%% Import

cf_in = lpjgu_matlab_readTable(file_in) ;
in_header_cell = cf_in.Properties.VariableNames ;


%% Convert

% Extract first year
cf_out = table2array(cf_in(cf_in.Year==min(cf_in.Year),:)) ;
out_header_cell = in_header_cell ;
cf_out(:,strcmp(in_header_cell,'Year')) = [] ;
out_header_cell(strcmp(in_header_cell,'Year')) = [] ;

% Classify columns
lastCharIsNum = @(x) strcmp(x(end),'0') ;
is_lonlat = contains(out_header_cell,{'Lon','Lat'}) ;
is_himelo = cellfun(lastCharIsNum,out_header_cell) ;
is_regcrop = ~is_himelo & ~is_lonlat ;

% Redistribute crop fractions
Nhimelo = length(find(is_himelo)) ;
cf_out(:,is_regcrop) = 0 ;
cf_out(:,is_himelo) = 1/Nhimelo ;


%% Save

%%% Options %%%%%%%%%%%%%
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
%%%%%%%%%%%%%%%%%%%%%%%%%

lpjgu_matlab_saveTable(out_header_cell, cf_out, out_file,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'justZeroCols', find(is_regcrop)) ;