function [I_agm, I_out] = e2p_translate_crops_agm2out(...
    varNames_agm, varNames_out)

Nvals_agm = getN_num(varNames_agm) ;
has_Nvals = ~isempty(Nvals_agm) ;

basenameIs_agm = getbasenamei(varNames_agm) ;

I_agm = nan(size(varNames_out)) ;
I_out = 1:length(varNames_out) ;

% Assign equivalents
for c = 1:length(varNames_out)
    thisVar_out = varNames_out{c} ;
    thisCrop_out = getbasename(thisVar_out) ;
    thisCropI_out = getbasenamei(thisVar_out) ;
    
    if any(strcmp(basenameIs_agm, thisCropI_out))
        thisVar_agm = thisVar_out ;
        thisCropI_agm = thisCropI_out ;
    else
    
        iN = strrep(thisVar_out, thisCrop_out, '') ;
        if has_Nvals
            thisN = getN_num(thisVar_out) ;
            if ~any(Nvals_agm == thisN)
                continue
            end
        end
        switch thisCrop_out
            case {'CerealsC3'}
                thisCrop_agm = 'max_wheat' ;
            case {'CerealsC4', 'Miscanthus', 'Sugarcane'}
                thisCrop_agm = 'maize' ;
            case {'Oilcrops','Pulses','OilNfix'}
                thisCrop_agm = 'soy' ;
            case {'CerealsC3s', 'StarchyRoots', 'FruitAndVeg', 'Sugarbeet', 'OilOther'}
                thisCrop_agm = 'spring_wheat' ;
            case {'CerealsC3w'}
                thisCrop_agm = 'winter_wheat' ;
            case {'Sugar'}
                thisCrop_agm = 'sugar' ;
            case {'Rice'}
                thisCrop_agm = 'rice' ;
            otherwise
                error('GGCMI equivalent of %s not specified', thisCrop_out)
        end
        thisVar_agm = [thisCrop_agm iN] ;
        thisCropI_agm = [thisCrop_agm strrep(thisCropI_out, thisCrop_out, '')] ;
    end
    
    ismatch = strcmp(basenameIs_agm, thisCropI_agm) ;
    if has_Nvals
        ismatch = ismatch & Nvals_agm==thisN ;
    end
    I = find(ismatch) ;
    if length(I) ~= 1
        error('Error finding %s in varNames_agm: expected 1, found %d', ...
            thisVar_agm, length(I))
    end
    I_agm(c) = I ;
end

% Strip skipped
I_out(isnan(I_agm)) = [] ;
I_agm(isnan(I_agm)) = [] ;

% Check
if length(I_out) ~= length(I_agm)
    error('Length mismatch in e2p_translate_crops_agm2out()')
end


end