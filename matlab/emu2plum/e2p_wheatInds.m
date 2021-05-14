function [i_thisWW, i_thisSW, i_thisMW, varNames] = ...
    e2p_wheatInds(thisWW, varNames)

thisSW = strrep(strrep(thisWW, 'winter', 'spring'), 'C3w', 'C3s') ;
thisMW = strrep(strrep(thisWW, 'winter', 'max'), 'C3w', 'C3') ;

if strcmp(thisSW, thisWW)
    error('Error finding spring-wheat counterpart for %s', thisWW)
elseif strcmp(thisMW, thisWW)
    error('Error finding max-wheat counterpart for %s', thisWW)
end

i_thisWW = find(strcmp(varNames,thisWW));
if length(i_thisWW) ~= 1
    error('Error finding i_thisWW (%d found)', length(i_thisWW))
end

i_thisSW = find(strcmp(varNames,thisSW));
if length(i_thisSW) ~= 1
    error('Error finding i_thisSW (%d found)', length(i_thisSW))
end

if ~contains(thisMW, varNames)
    varNames{end+1} = thisMW ;
end

i_thisMW = find(strcmp(varNames,thisMW));
if length(i_thisMW) ~= 1
    error('Error finding i_thisMW (%d found)', length(i_thisMW))
end


end