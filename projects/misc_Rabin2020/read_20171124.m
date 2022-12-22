function varNames = read_20171124(...
    in_file, table_name, varNames, ...
    lons_in, lats_in, years_in)

% Set up & parse input arguments
p = inputParser ;
addRequired(p,'in_file',@isstr) ;
addRequired(p,'table_name',@isstr) ;
addRequired(p,'varNames',@isstruct) ;
parse(p,in_file, table_name, varNames) ;

% Read table
disp('Reading file into table...')
in_table = lpjgu_matlab_readTable(...
    in_file, 'do_save_mat',true,'verbose',false) ;

% Check lons, lats, and years against original
lons = in_table.Lon(in_table.Year==min(in_table.Year)) ;
if ~isequal(lons_in,lons)
    error('Longitudes don''t match!')
end
lats = in_table.Lat(in_table.Year==min(in_table.Year)) ;
if ~isequal(lats_in,lats)
    error('Latitudes don''t match!')
end
years = unique(in_table.Year) ;
if ~isequal(years_in,years)
    error('Years don''t match!')
end

% Save variable names
eval(['varNames.' table_name ' = in_table.Properties.VariableNames ;']) ;

disp('Done.')

end

