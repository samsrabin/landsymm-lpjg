%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert LU file to just CPB %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

remapVer = '11_g2p' ;

file_in = sprintf('/Volumes/Reacher/G2P/inputs/LU/remaps_v%s/LU.remapv%s.txt', ...
    remapVer, remapVer) ;


%% Import

disp('Importing...')

lu_in = lpjgu_matlab_readTable(file_in) ;
in_header_cell = lu_in.Properties.VariableNames ;

disp('Done.')


%% Convert

disp('Converting...')

% Extract last year
lu_out = table2array(lu_in(lu_in.Year==max(lu_in.Year),:)) ;
lu_out(:,strcmp(in_header_cell,'Year')) = [] ;
in_header_cell(strcmp(in_header_cell,'Year')) = [] ;

% Redistribute
lu_out(:,contains(in_header_cell,{'CROPLAND','PASTURE'})) = repmat(sum(lu_out(:,contains(in_header_cell,{'CROPLAND','PASTURE','NATURAL'})),2) / 2, [1 2]) ;
lu_out(:,strcmp(in_header_cell,'NATURAL')) = 0 ;

% Ensure sum to 1
tmp = lu_out(:,contains(in_header_cell,{'CROPLAND','PASTURE','NATURAL','BARREN'})) ;
thisSum = sum(tmp,2) ;
lu_out(:,contains(in_header_cell,{'CROPLAND','PASTURE','NATURAL','BARREN'})) = tmp ./ repmat(thisSum,[1 4]) ;

disp('Done.')


%% Save

%%% Options %%%%%%%%%%%%%
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Saving...')

out_file = strrep(file_in, '.txt', '.justCPB_1yr.txt') ;

lpjgu_matlab_saveTable(in_header_cell, lu_out, out_file,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'justZeroCols', find(strcmp(in_header_cell,'NATURAL'))) ;
    
disp('Done.')
    
    
    
    
    
% Warning: Giving some from PASTURE to CROPLAND (1132013 cells).
% Warning: Giving some from NATURAL to CROPLAND (2768641 cells).