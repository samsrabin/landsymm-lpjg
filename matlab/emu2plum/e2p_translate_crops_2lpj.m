function [cropList_in_asLpj, used_lpjCrops] = e2p_translate_crops_2lpj(...
    cropList_in, cropList_lpj, colname, varargin)

verbose = true ;
if ~isempty(varargin)
    verbose = varargin{1} ;
    if length(varargin)>1
        error('At most 1 optional argument (verbose) accepted in e2p_translate_crops_2lpj')
    end
end

cropList_in_asLpj = cell(size(cropList_in)) ;

% Assign equivalents
for c = 1:length(cropList_in)
    thisCrop_in = cropList_in{c} ;
    if contains(thisCrop_in, cropList_lpj)
        cropList_in_asLpj{c} = thisCrop_in ;
    elseif contains(lower(thisCrop_in), cropList_lpj)
        cropList_in_asLpj{c} = lower(thisCrop_in) ;
    elseif contains(thisCrop_in, {'Sugarcane'})
        cropList_in_asLpj{c} = 'CerealsC4' ;
    elseif contains(thisCrop_in,{'OilNfix', 'Pulses'})
        cropList_in_asLpj{c} = 'Oilcrops' ;
    elseif contains(thisCrop_in,{'StarchyRoots', 'OilOther', 'Sugarbeet'})
        cropList_in_asLpj{c} = 'CerealsC3s' ;
    else
        error('LPJ-GUESS equivalent of %s not specified', thisCrop_in)
    end
end

% Check
used_lpjCrops = false(size(cropList_lpj)) ;
for c = 1:length(cropList_in_asLpj)
    thisCrop_in = cropList_in_asLpj{c} ;
    if ~contains(thisCrop_in,cropList_lpj)
        error('%s assigned as translation, but does not exist in cropList_lpj', thisCrop_in)
    end
    used_lpjCrops(strcmp(cropList_lpj, thisCrop_in)) = true ;
end

% Display results
if verbose && ~isequal(cropList_in, cropList_in_asLpj)
    disp(' ')
    disp(table(cropList_in', cropList_in_asLpj', 'VariableNames', {colname, 'LPJGequiv'}))
end


end