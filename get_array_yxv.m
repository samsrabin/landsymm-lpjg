function [array_yxv, varNames, lonlats_in] = get_array_yxv(thisFile, lonlats_target, yearList)

Nyears = length(yearList) ;

table_in = lpjgu_matlab_readTable(thisFile) ;
lonlats_in = unique(table2array(table_in(:,1:2)), ...
    'rows', 'stable') ;
if ~isempty(lonlats_target) && ~isequal(lonlats_target, lonlats_in)
    if ~isequal(sortrows(lonlats_target), sortrows(lonlats_in))
        error('Gridlist mismatch (add sort code to fix)')
    else
        error('Gridlist mismatch')
    end
end
Ncells = length(lonlats_in) ;

varNames = setdiff(table_in.Properties.VariableNames, ...
    {'Lon','Lat','Year'}, ...
    'stable') ;
Nvars = length(varNames) ;

has_years = any(strcmp(table_in.Properties.VariableNames,'Year')) ;
if ~has_years
    array_xv = table2array(table_in(:,3:end)) ;
    array_yxv = repmat(permute(array_xv,[3 1 2]), [Nyears 1 1]) ;
else
    yearList_in = unique(table_in.Year) ;
    if ~(min(yearList_in)<=min(yearList) && max(yearList_in)>=max(yearList))
        error('Mismatch between desired and imported yearLists!')
    end
    isok_year = table_in.Year>=min(yearList) & table_in.Year<=max(yearList) ;
    if any(~isok_year)
        array_yxv = lpjgu_matlab_table2array(table_in(isok_year,4:end), [Nyears Ncells Nvars]) ;
    else
        array_yxv = lpjgu_matlab_table2array(table_in(:,4:end), [Nyears Ncells Nvars]) ;
    end
end

end