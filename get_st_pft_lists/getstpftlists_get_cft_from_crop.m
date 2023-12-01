function cft = getstpftlists_get_cft_from_crop(crop)

switch crop
    case {'CerealsC3', 'CerealsC3w'}
        cft = 'TeWW' ;
    case {'CerealsC3s', 'StarchyRoots', 'FruitAndVeg', 'Sugar', ...
            'Sugarbeet', 'OilOther', 'ExtraCrop'}
        cft = 'TeSW' ;
    case {'CerealsC4', 'Sugarcane'}
        cft = 'TeCo' ;
    case {'Rice'}
        cft = 'TrRi' ;
    case {'Oilcrops', 'Pulses', 'OilNfix'}
        cft = 'TeSo' ;
    otherwise
        error('Crop %s not recognized in get_cft_from_crop()', crop)
end


end