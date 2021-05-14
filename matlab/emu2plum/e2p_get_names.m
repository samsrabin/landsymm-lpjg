function [varNames, cropList, varNames_basei, cropList_basei, Nlist, maxN] = ...
    e2p_get_names(varNames_A, varNames_B, ...
    get_unneeded)

% Make sure A and B variable lists match (if A is included)
if ~isempty(varNames_A) && ~isempty(varNames_B)
    if ~isequal(varNames_A, varNames_B)
        error('Mismatch between variable lists A and B.')
    end
    varNames = varNames_B ;
elseif ~isempty(varNames_A)
    varNames = varNames_A ;
elseif ~isempty(varNames_B)
    varNames = varNames_B ;
else
    error('You must supply at least one of varNames_A or varNames_B')
end


% Make sure any non-experimental names (e.g., CerealsC3 as opposed to
% CerealsC3060) have been stripped
if any(get_unneeded(varNames))
    error('Non-experimental variables need to already have been stripped')
end

% Get list of crops (not distinguishing rainfed vs. irrigated)
cropList = unique(getbasename(varNames)) ;
Ncrops = length(cropList) ;

% Get list of rainfed and irrigated crops
varNames_basei = getbasenamei(varNames) ;
cropList_basei = unique(varNames_basei) ;

% Make sure that each crop has irrigated and rainfed versions
if length(cropList_basei) ~= Ncrops*2
    error('Mismatch of crop lists here')
end

% Get N levels
Nlist = unique(getN_char(varNames)) ;
maxN = max(getN_num(varNames)) ;


end