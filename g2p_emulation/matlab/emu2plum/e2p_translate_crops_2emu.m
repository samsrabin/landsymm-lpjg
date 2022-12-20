function [cropList_in_asEmu, used_emuCrops] = e2p_translate_crops_2emu(...
    cropList_in, cropList_emu, colname, varargin)

verbose = true ;
if ~isempty(varargin)
    verbose = varargin{1} ;
    if length(varargin)>1
        error('At most 1 optional argument (verbose) accepted in e2p_translate_crops_2emu')
    end
end

cropList_in_asEmu = cell(size(cropList_in)) ;

% Assign equivalents
for c = 1:length(cropList_in)
    thisCrop_in = cropList_in{c} ;
    if contains(lower(thisCrop_in),cropList_emu)
        cropList_in_asEmu{c} = lower(thisCrop_in) ;
    elseif strcmp(thisCrop_in,'CerealsC3')
        cropList_in_asEmu{c} = 'max_wheat' ;
    elseif contains(thisCrop_in, {'CerealsC4', 'Sugarcane'})
        cropList_in_asEmu{c} = 'maize' ;
    elseif contains(thisCrop_in,{'Oilcrops', 'Pulses', 'OilNfix'})
        if contains('soy',cropList_emu)
            cropList_in_asEmu{c} = 'soy' ;
        else
            cropList_in_asEmu{c} = 'spring_wheat' ;
        end
    elseif contains(thisCrop_in, ...
            {'CerealsC3s', 'StarchyRoots', 'Rice', 'OilOther', ...
            'Sugarbeet', 'FruitAndVeg'})
        cropList_in_asEmu{c} = 'spring_wheat' ;
    elseif contains(thisCrop_in, 'CerealsC3w')
        cropList_in_asEmu{c} = 'winter_wheat' ;
    else
        error('GGCMI equivalent of %s not specified', thisCrop_in)
    end
end

% Check
used_emuCrops = false(size(cropList_emu)) ;
for c = 1:length(cropList_in_asEmu)
    thisCrop_in = cropList_in_asEmu{c} ;
    if ~contains(thisCrop_in,cropList_emu)
        error('%s assigned as translation, but does not exist in cropList_emu', thisCrop_in)
    end
    used_emuCrops(strcmp(cropList_emu, thisCrop_in)) = true ;
end

% Display results
if verbose && ~isequal(cropList_in, cropList_in_asEmu)
    disp(' ')
    disp(table(cropList_in', cropList_in_asEmu', 'VariableNames', {colname, 'GCCMIequiv'}))
end


end