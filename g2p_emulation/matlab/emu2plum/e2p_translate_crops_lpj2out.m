function [I_lpj, I_out] = e2p_translate_crops_lpj2out(...
    varNames_lpj, varNames_out)

Nvals_lpj = getN_num(varNames_lpj) ;
basenameIs_lpj = getbasenamei(varNames_lpj) ;

I_lpj = nan(size(varNames_out)) ;
I_out = 1:length(varNames_out) ;

% Assign equivalents
for c = 1:length(varNames_out)
    thisVar_out = varNames_out{c} ;
    thisCrop_out = getbasename(thisVar_out) ;
    thisCropI_out = getbasenamei(thisVar_out) ;
    iN = strrep(thisVar_out, thisCrop_out, '') ;
    thisN = getN_num(thisVar_out) ;
    if ~any(Nvals_lpj == thisN)
        continue
    end
    switch thisCrop_out
        case {'CerealsC3'}
            thisCrop_lpj = 'CerealsC3' ;
        case {'CerealsC4', 'Miscanthus', 'Sugarcane'}
            thisCrop_lpj = 'CerealsC4' ;
        case {'Oilcrops','Pulses','OilNfix'}
            thisCrop_lpj = 'Oilcrops' ;
        case {'CerealsC3s', 'StarchyRoots', 'FruitAndVeg', 'Sugarbeet', 'OilOther'}
            thisCrop_lpj = 'CerealsC3s' ;
        case {'CerealsC3w'}
            thisCrop_lpj = 'CerealsC3w' ;
        case {'Rice'}
            thisCrop_lpj = 'Rice' ;
        otherwise
            error('GGCMI equivalent of %s not specified', thisCrop_out)
    end
    thisVar_lpj = [thisCrop_lpj iN] ;
    thisCropI_lpj = [thisCrop_lpj strrep(thisCropI_out, thisCrop_out, '')] ;
    I = find(strcmp(basenameIs_lpj, thisCropI_lpj) & Nvals_lpj==thisN) ;
    if length(I) ~= 1
        error('Error finding %s in varNames_lpj: expected 1, found %d', ...
            thisVar_lpj, length(I))
    end
    I_lpj(c) = I ;
end

% Strip skipped
I_out(isnan(I_lpj)) = [] ;
I_lpj(isnan(I_lpj)) = [] ;

% Check
if length(I_out) ~= length(I_lpj)
    error('Length mismatch in e2p_translate_crops_lpj2out()')
end


end