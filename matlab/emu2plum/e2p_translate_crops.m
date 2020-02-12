function [cropList_lpj_asEmu, used_emuCrops] = e2p_translate_crops(...
    cropList_lpj, cropList_emu)

cropList_lpj_asEmu = cell(size(cropList_lpj)) ;

% Assign equivalents
for c = 1:length(cropList_lpj)
    thisCrop_lpj = cropList_lpj{c} ;
    if contains(lower(thisCrop_lpj),cropList_emu)
        cropList_lpj_asEmu{c} = lower(thisCrop_lpj) ;
    elseif strcmp(thisCrop_lpj,'CerealsC3')
        cropList_lpj_asEmu{c} = 'max_wheat' ;
    elseif strcmp(thisCrop_lpj,'CerealsC4')
        cropList_lpj_asEmu{c} = 'maize' ;
    elseif contains(thisCrop_lpj,{'Oilcrops', 'Pulses'})
        if contains('soy',cropList_emu)
            cropList_lpj_asEmu{c} = 'soy' ;
        else
            cropList_lpj_asEmu{c} = 'spring_wheat' ;
        end
    elseif strcmp(thisCrop_lpj,'StarchyRoots')
        cropList_lpj_asEmu{c} = 'spring_wheat' ;
    else
        error('GGCMI equivalent of %s not specified', thisCrop_lpj)
    end
end

% Check
used_emuCrops = false(size(cropList_emu)) ;
for c = 1:length(cropList_lpj_asEmu)
    thisCrop_lpj = cropList_lpj_asEmu{c} ;
    if ~contains(thisCrop_lpj,cropList_emu)
        error('%s assigned as translation, but does not exist in cropList_emu', thisCrop_lpj)
    end
    used_emuCrops(strcmp(cropList_emu, thisCrop_lpj)) = true ;
end

% Display results
disp(' ')
disp(table(cropList_lpj', cropList_lpj_asEmu', 'VariableNames', {'LPJ', 'GCCMIequiv'}))


end