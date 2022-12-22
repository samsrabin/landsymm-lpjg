function dbl_out = remove_nul(cell_in)

cell_in(strcmp(cell_in,'nul')) = {'NaN'} ;
dbl_out = str2double(cell_in) ;


end