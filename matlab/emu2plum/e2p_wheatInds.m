function [i_theseABetc, i_thisM, varNames] = ...
    e2p_wheatInds(thisA, varNames, combineCrops_row)

combineCrops_dest = combineCrops_row{1} ;
combineCrops_sources = combineCrops_row{2} ;
Nsources = length(combineCrops_sources) ;

% Get list of source variables
theseABetc = cell(1, Nsources) ;
theseABetc{1} = thisA ;
for s = 2:Nsources
    thisSource = combineCrops_sources{s} ;
    theseABetc{s} = strrep(thisA, getbasename(thisA), thisSource) ;
end
thisM = strrep(thisA, getbasename(thisA), combineCrops_dest) ;

if length(unique(theseABetc)) ~= length(theseABetc)
    error('Error defining counterpart source variable name(s) for %s', thisA)
end

missing_sources = setdiff(theseABetc, varNames) ;
if ~isempty(missing_sources)
    error('%s not found in varNames', missing_sources{1})
end

[C, ~, i_theseABetc] = intersect(theseABetc, varNames, 'stable') ;
if ~isequal(shiftdim(C), shiftdim(theseABetc))
    error('~isequal(shiftdim(C), shiftdim(theseABetc))')
end

if ~contains(thisM, varNames)
    varNames{end+1} = thisM ;
end

i_thisM = find(strcmp(varNames,thisM));
if length(i_thisM) ~= 1
    error('Error finding %s (%d found)', thisM, length(i_thisM))
end


end