function yearList = e2p_get_incl_years(filename)
    
% % Read time axis directly
% varid = find(strcmp(varnames, 'time')) ;
% if length(varid) ~= 1
%     error('Expected to find 1 time variable in %s; found %d', ...
%         filename, length(varid))
% end
% attnames = {finfo.Variables(varid).Attributes.Name} ;
% attid = find(strcmp(attnames, 'units')) ;

% Just parse the filename
year_cellstr = regexp(filename, '\d\d\d\d_\d\d\d\d', 'match') ;
if length(year_cellstr) ~= 1
    error('Expected to find 1 possible year-range string in %s; found %d', ...
        filename, length(year_cellstr))
end
year_strs = strsplit(year_cellstr{1}, '_') ;
year_range = str2double(year_strs) ;
yearList = seq(year_range) ;

end