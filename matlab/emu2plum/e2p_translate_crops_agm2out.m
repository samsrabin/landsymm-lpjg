function I_out = e2p_translate_crops_agm2out(...
    varNames_agm, varNames_out)

I_out = nan(size(varNames_out)) ;

% Assign equivalents
for c = 1:length(varNames_out)
    thisVar_out = varNames_out{c} ;
    thisCrop_out = getbasename(thisVar_out) ;
    iN = strrep(thisVar_out, thisCrop_out, '') ;
    switch thisCrop_out
        case {'CerealsC3'}
            thisCrop_agm = 'max_wheat' ;
        case {'CerealsC4', 'Miscanthus', 'Sugarcane'}
            thisCrop_agm = 'maize' ;
        case {'Oilcrops','Pulses'}
            thisCrop_agm = 'soy' ;
        case {'StarchyRoots', 'FruitAndVeg', 'Sugarbeet'}
            thisCrop_agm = 'spring_wheat' ;
        case {'Sugar'}
            thisCrop_agm = 'sugar' ;
        case {'Rice'}
            thisCrop_agm = 'rice' ;
        otherwise
            error('GGCMI equivalent of %s not specified', thisCrop_out)
    end
    thisVar_agm = [thisCrop_agm iN] ;
    I = find(strcmp(varNames_agm, thisVar_agm)) ;
    if length(I) ~= 1
        error('Error finding %s in varNames_agm: expected 1, found %d', ...
            thisVar_agm, length(I))
    end
    I_out(c) = I ;
end

% Check
if any(isnan(I_out))
    error('any(isnan(I_out))')
end


end