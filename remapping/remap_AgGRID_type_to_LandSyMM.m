function aggrid_type = remap_AgGRID_type_to_LandSyMM(landsymm_type)

switch landsymm_type
    case {'CerealsC3', 'StarchyRoots', 'FruitAndVeg', ...
            'Sugarbeet', 'OilOther', 'ExtraCrop'}
        aggrid_type = 'wheat' ;
    case {'CerealsC4', 'Miscanthus', 'Sugarcane'}
        aggrid_type = 'maize' ;
    case {'Oilcrops', 'Pulses', 'OilNfix'}
        aggrid_type = 'soybean' ;
    case {'Rice'}
        aggrid_type = 'rice' ;
    case {'Sugar'}
        aggrid_type = 'combined_sugars' ;
    otherwise
        error('What crop from AgGRID Nfert inputs should I use for %s?', landsymm_type)
end

end