function [varNames, cropList, varNames_basei, cropList_basei, Nlist, maxN] = ...
    e2p_get_names(varNames_bl, varNames_fu, ...
    getbasename, getbasenamei, getN)

% Make sure baseline and future variable names match
if ~isequal(varNames_bl, varNames_fu)
    error('Mismatch between variable names in baseline and future.')
end
varNames = varNames_bl ;

% Make sure any non-experimental names (e.g., CerealsC3 as opposed to
% CerealsC3060) have been stripped
is_nonexp = ~cellfun(@isempty,regexp(varNames, '.*\d\d+')) ;
if any(~is_nonexp)
    error('Non-experimental variables need to already have been stripped')
end

% Get list of crops (not distinguishing rainfed vs. irrigated)
cropList = unique(cellfun(getbasename, varNames, 'UniformOutput', false)) ;
Ncrops = length(cropList) ;

% Get list of rainfed and irrigated crops
varNames_basei = cellfun(getbasenamei, varNames, 'UniformOutput', false) ;
cropList_basei = unique(varNames_basei) ;

% Make sure that each crop has irrigated and rainfed versions
if length(cropList_basei) ~= Ncrops*2
    error('Mismatch of crop lists here')
end

% Get N levels
Nlist = unique(cellfun(getN, varNames, 'UniformOutput', false)) ;
maxN = max(str2double(Nlist)) ;


end