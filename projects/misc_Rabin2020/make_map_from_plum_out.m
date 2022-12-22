function out_YXyr = make_map_from_plum_out( ...
    in_struct, thisCommodity, thisVar, thisYearList, ...
    countrygroups_YX, countrygroup_list)

% Check inputs
if isempty(thisVar) && isfield(in_struct,'varNames')
    error('You must provide thisVar!')
elseif ~isfield(in_struct,'varNames') && ~isempty(thisVar)
    warning('Ignoring thisVar (%s)', thisVar)
    thisVar = '' ;
end

% Restrict to thisVar
if ~isempty(thisVar)
    data_umyr = squeeze(in_struct.data_umvyr(:,:,strcmp(in_struct.varNames,thisVar),:,:)) ;
else
    data_umyr = in_struct.data_umyr ;
end
Nruns = size(data_umyr,4) ;

% Restrict to thisCommodity, yearList
data_uyr = squeeze(data_umyr( ...
    :, ...
    strcmp(in_struct.commodities,thisCommodity), ...
    in_struct.years>=min(thisYearList) & in_struct.years<=max(thisYearList), ...
    :)) ;

% Make maps
Nyears = length(thisYearList) ;
out_YXyr = nan([size(countrygroups_YX) Nyears Nruns]) ;
for c = 1:length(countrygroup_list)
    isThisCountrygroup_YX = countrygroups_YX==c ;
    isThisCountrygroup_YXyr = repmat(isThisCountrygroup_YX,[1 1 Nyears Nruns]) ;
    theData_YX = permute(squeeze(data_uyr(c,:,:)), [3 4 1 2]) ;
    theData_YXyr = repmat(theData_YX, [size(countrygroups_YX) 1 1]) ;
    try
        out_YXyr(isThisCountrygroup_YXyr) = theData_YXyr(isThisCountrygroup_YXyr) ;
    catch ME
        keyboard
    end
    
end


end