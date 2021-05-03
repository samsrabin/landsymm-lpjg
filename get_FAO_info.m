function [listCrops_fa2o, ...
    listCrops_fa2i, listCrops_fa2i1, listCrops_fa2i2, ...
    FAO_to_FAO_key, FAO_to_FAO_key1, FAO_to_FAO_key2] = ...
    get_FAO_info ( calib_ver, twofiles )

listCrops_fa2i = {} ;
listCrops_fa2i1 = {} ;
listCrops_fa2i2 = {} ;
FAO_to_FAO_key = {} ;
FAO_to_FAO_key1 = {} ;
FAO_to_FAO_key2 = {} ;

% Get FAO_to_FAO_key
if calib_ver==1 || calib_ver==6 || calib_ver==8
    listCrops_fa2o = {'Wheat','Maize','Rice','Oilcrops','Pulses','Starchy roots'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    FAO_to_FAO_key{getFi('Wheat')}          = {'Wheat','Barley','Oats'} ;
    FAO_to_FAO_key{getFi('Maize')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key{getFi('Rice')}           = {'Rice, paddy'} ;
    FAO_to_FAO_key{getFi('Oilcrops')}       = {'Oilcrops Primary'} ;
    FAO_to_FAO_key{getFi('Pulses')}         = {'Pulses,Total'} ;
    FAO_to_FAO_key{getFi('Starchy roots')}  = {'Roots and Tubers,Total'} ;
elseif calib_ver==2
    listCrops_fa2o = {'Wheat','Barley','Oats','Maize','Millet','Sorghum','Rice','Oilcrops','Pulses','Starchy roots'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    FAO_to_FAO_key{getFi('Wheat')}          = {'Wheat'} ;
    FAO_to_FAO_key{getFi('Barley')}         = {'Barley'} ;
    FAO_to_FAO_key{getFi('Oats')}           = {'Oats'} ;
    FAO_to_FAO_key{getFi('Maize')}          = {'Maize'} ;
    FAO_to_FAO_key{getFi('Millet')}         = {'Millet'} ;
    FAO_to_FAO_key{getFi('Sorghum')}        = {'Sorghum'} ;
    FAO_to_FAO_key{getFi('Rice')}           = {'Rice, paddy'} ;
    FAO_to_FAO_key{getFi('Oilcrops')}       = {'Oilcrops Primary'} ;
    FAO_to_FAO_key{getFi('Pulses')}         = {'Pulses,Total'} ;
    FAO_to_FAO_key{getFi('Starchy roots')}  = {'Roots and Tubers,Total'} ;
elseif calib_ver==3 || calib_ver==7
    listCrops_fa2o = {'Wheat','Maize','Rice','Oilcrops','Pulses','Starchy roots'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    FAO_to_FAO_key{getFi('Wheat')}          = {'Wheat','Barley','Oats'} ;
    FAO_to_FAO_key{getFi('Maize')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key{getFi('Rice')}           = {'Rice, paddy'} ;
    FAO_to_FAO_key{getFi('Oilcrops')}       = {'Oilcrops Primary','Oilcakes Equivalent'} ;
    FAO_to_FAO_key{getFi('Pulses')}         = {'Pulses,Total'} ;
    FAO_to_FAO_key{getFi('Starchy roots')}  = {'Roots and Tubers,Total'} ;
elseif calib_ver==4
    listCrops_fa2o = {'Wheat','Maize','Rice','Oilcrops','Pulses','Starchy roots'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    FAO_to_FAO_key{getFi('Wheat')}          = {'Wheat','Barley','Oats'} ;
    FAO_to_FAO_key{getFi('Maize')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key{getFi('Rice')}           = {'Rice, paddy'} ;
    FAO_to_FAO_key{getFi('Oilcrops')}       = {'Coconuts','Cottonseed',...
        'Groundnuts, with shell','Karite nuts (sheanuts)',...
        'Castor oil seed','Tung nuts','Jojoba seed','Safflower seed',...
        'Poppy seed','Melonseed','Tallowtree seed','Kapok fruit',...
        'Kapokseed in shell','Kapokseed shelled','Linseed','Hempseed',...
        'Oilseeds nes','Flour, oilseeds','Olives','Palm kernels',...
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed'} ;
    FAO_to_FAO_key{getFi('Pulses')}         = {'Pulses,Total'} ;
    FAO_to_FAO_key{getFi('Starchy roots')}  = {'Roots and Tubers,Total'} ;
elseif calib_ver==5
%     error('Do this')
    listCrops_fa2o = {'Wheat','Maize','Rice','Oilcrops','Pulses','Starchy roots'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key1{getFi('Wheat')}          = {'Wheat','Barley','Oats'} ;
    FAO_to_FAO_key1{getFi('Maize')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key1{getFi('Rice')}           = {'Rice, paddy'} ;
    FAO_to_FAO_key1{getFi('Oilcrops')}       = {'Coconuts',...
        ...%'Cottonseed',...
        'Seed cotton',... % This includes cotton grown for fiber...
        'Groundnuts, with shell','Karite nuts (sheanuts)',...
        'Castor oil seed','Tung nuts','Jojoba seed','Safflower seed',...
        'Poppy seed','Melonseed','Tallowtree seed','Kapok fruit',...
        'Kapokseed in shell','Kapokseed shelled','Linseed','Hempseed',...
        'Oilseeds nes','Flour, oilseeds','Olives',...
        ...%'Palm kernels',...
        'Oil, palm fruit',... % "Palm kernels" not included in Area
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed'} ;
    FAO_to_FAO_key1{getFi('Pulses')}         = {'Pulses,Total'} ;
    FAO_to_FAO_key1{getFi('Starchy roots')}  = {'Roots and Tubers,Total'} ;
    % Get item names for COMMODITY BALANCE: CROPS PRIMARY EQUIVALENT dataset
    FAO_to_FAO_key2{getFi('Wheat')}          = {'Wheat and products','Barley and products','Oats'} ;
    FAO_to_FAO_key2{getFi('Maize')}          = {'Maize and products','Millet and products','Sorghum and products'} ;
    FAO_to_FAO_key2{getFi('Rice')}           = {'Rice (Paddy Equivalent)'} ;
    FAO_to_FAO_key2{getFi('Oilcrops')}       = {'Oilcrops'} ;
    FAO_to_FAO_key2{getFi('Pulses')}         = {'Pulses'} ;
    FAO_to_FAO_key2{getFi('Starchy roots')}  = {'Starchy Roots'} ;
elseif calib_ver==9
%     error('Do this')
    listCrops_fa2o = {'Wheat','Maize','Rice','Oilcrops','Pulses','Starchy roots'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key1{getFi('Wheat')}          = {'Wheat','Barley','Oats'} ;
    FAO_to_FAO_key1{getFi('Maize')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key1{getFi('Rice')}           = {'Rice paddy'} ;
    FAO_to_FAO_key1{getFi('Oilcrops')}       = {'Coconuts',...
        ...%'Cottonseed',...
        'Seed cotton',... % This includes cotton grown for fiber...
        'Groundnuts with shell','Karite nuts (sheanuts)',...
        'Castor oil seed','Tung nuts','Jojoba seed','Safflower seed',...
        'Poppy seed','Melonseed','Tallowtree seed','Kapok fruit',...
        ...%'Kapokseed in shell','Kapokseed shelled',... % Not included in Area
        'Linseed','Hempseed',...
        'Oilseeds nes',...
        ...%'Flour, oilseeds', % Not included in Area
        'Olives',...
        ...%'Palm kernels',... % "Palm kernels" not included in Area
        'Oil palm fruit',...
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed'} ;
    FAO_to_FAO_key1{getFi('Pulses')}         = {'Pulses Total'} ;
    FAO_to_FAO_key1{getFi('Starchy roots')}  = {'Roots and Tubers Total'} ;
    % Get item names for COMMODITY BALANCE: CROPS PRIMARY EQUIVALENT dataset
    FAO_to_FAO_key2{getFi('Wheat')}          = {'Wheat and products','Barley and products','Oats'} ;
    FAO_to_FAO_key2{getFi('Maize')}          = {'Maize and products','Millet and products','Sorghum and products'} ;
    FAO_to_FAO_key2{getFi('Rice')}           = {'Rice (Paddy Equivalent)'} ;
    FAO_to_FAO_key2{getFi('Oilcrops')}       = {'Oilcrops'} ;
    FAO_to_FAO_key2{getFi('Pulses')}         = {'Pulses'} ;
    FAO_to_FAO_key2{getFi('Starchy roots')}  = {'Starchy Roots'} ;
elseif calib_ver==10 || calib_ver==12
%     error('Do this')
    listCrops_fa2o = {'Wheat','Maize','Rice','Oilcrops','Pulses','Starchy roots'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key1{getFi('Wheat')}          = {'Wheat','Barley'} ;
    FAO_to_FAO_key1{getFi('Maize')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key1{getFi('Rice')}           = {'Rice paddy'} ;
    FAO_to_FAO_key1{getFi('Oilcrops')}       = {'Coconuts',...
        ...%'Cottonseed',...
        'Seed cotton',... % This includes cotton grown for fiber...
        'Groundnuts with shell','Karite nuts (sheanuts)',...
        'Castor oil seed','Tung nuts','Jojoba seed','Safflower seed',...
        'Poppy seed','Melonseed','Tallowtree seed','Kapok fruit',...
        ...%'Kapokseed in shell','Kapokseed shelled',... % Not included in Area
        'Linseed','Hempseed',...
        'Oilseeds nes',...
        ...%'Flour, oilseeds', % Not included in Area
        'Olives',...
        ...%'Palm kernels',... % "Palm kernels" not included in Area
        'Oil palm fruit',...
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed'} ;
    FAO_to_FAO_key1{getFi('Pulses')}         = {'Pulses Total'} ;
    FAO_to_FAO_key1{getFi('Starchy roots')}  = {'Roots and Tubers Total'} ;
    % Get item names for COMMODITY BALANCE: CROPS PRIMARY EQUIVALENT dataset
    FAO_to_FAO_key2{getFi('Wheat')}          = {'Wheat and products','Barley and products','Oats'} ;
    FAO_to_FAO_key2{getFi('Maize')}          = {'Maize and products','Millet and products','Sorghum and products'} ;
    FAO_to_FAO_key2{getFi('Rice')}           = {'Rice (Paddy Equivalent)'} ;
    FAO_to_FAO_key2{getFi('Oilcrops')}       = {'Oilcrops'} ;
    FAO_to_FAO_key2{getFi('Pulses')}         = {'Pulses'} ;
    FAO_to_FAO_key2{getFi('Starchy roots')}  = {'Starchy Roots'} ;
elseif calib_ver==11
%     error('Do this')
    listCrops_fa2o = {'Faba bean','Sorghum','Soybean'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    FAO_to_FAO_key{getFi('Faba bean')}         = {'Pulses,Total'} ;
    FAO_to_FAO_key{getFi('Sorghum')}           = {'Sorghum'} ;
    FAO_to_FAO_key{getFi('Soybean')}           = {'Soybeans'} ;
elseif calib_ver==13
%     error('Do this')
    listCrops_fa2o = {'Wheat','Maize','Sorghum','Pulses','Soybeans','Rice'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key1{getFi('Wheat')}          = {'Wheat'} ;
    FAO_to_FAO_key1{getFi('Maize')}          = {'Maize'} ;
    FAO_to_FAO_key1{getFi('Sorghum')}        = {'Sorghum'} ;
    FAO_to_FAO_key1{getFi('Rice')}           = {'Rice paddy'} ;
    FAO_to_FAO_key1{getFi('Soybeans')}       = {'Soybeans'} ;
    FAO_to_FAO_key1{getFi('Pulses')}         = {'Pulses Total'} ;
    % Get item names for COMMODITY BALANCE: CROPS PRIMARY EQUIVALENT dataset
    FAO_to_FAO_key2{getFi('Wheat')}          = {'Wheat and products'} ;
    FAO_to_FAO_key2{getFi('Maize')}          = {'Maize and products'} ;
    FAO_to_FAO_key2{getFi('Sorghum')}        = {'Sorghum and products'} ;
    FAO_to_FAO_key2{getFi('Rice')}           = {'Rice (Paddy Equivalent)'} ;
    FAO_to_FAO_key2{getFi('Soybeans')}       = {'Soyabeans'} ;
    FAO_to_FAO_key2{getFi('Pulses')}         = {'Pulses'} ;
elseif calib_ver==14
%     error('Do this')
    listCrops_fa2o = {'Wheat','Maize','Sorghum','Rice'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key1{getFi('Wheat')}          = {'Wheat'} ;
    FAO_to_FAO_key1{getFi('Maize')}          = {'Maize'} ;
    FAO_to_FAO_key1{getFi('Sorghum')}        = {'Sorghum'} ;
    FAO_to_FAO_key1{getFi('Rice')}           = {'Rice paddy'} ;
    % Get item names for COMMODITY BALANCE: CROPS PRIMARY EQUIVALENT dataset
    FAO_to_FAO_key2{getFi('Wheat')}          = {'Wheat and products'} ;
    FAO_to_FAO_key2{getFi('Maize')}          = {'Maize and products'} ;
    FAO_to_FAO_key2{getFi('Sorghum')}        = {'Sorghum and products'} ;
    FAO_to_FAO_key2{getFi('Rice')}           = {'Rice (Paddy Equivalent)'} ;
elseif calib_ver==15
%     error('Do this')
    listCrops_fa2o = {'Wheat','Maize','Rice','Oilcrops','Pulses','Starchy roots'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key1{getFi('Wheat')}          = {'Wheat','Barley','Rye'} ;
    FAO_to_FAO_key1{getFi('Maize')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key1{getFi('Rice')}           = {'Rice paddy'} ;
    FAO_to_FAO_key1{getFi('Oilcrops')}       = {'Coconuts',...
        ...%'Cottonseed',...
        'Seed cotton',... % This includes cotton grown for fiber...
        'Groundnuts with shell','Karite nuts (sheanuts)',...
        'Castor oil seed','Tung nuts','Jojoba seed','Safflower seed',...
        'Poppy seed','Melonseed','Tallowtree seed','Kapok fruit',...
        ...%'Kapokseed in shell','Kapokseed shelled',... % Not included in Area
        'Linseed','Hempseed',...
        'Oilseeds nes',...
        ...%'Flour, oilseeds', % Not included in Area
        'Olives',...
        ...%'Palm kernels',... % "Palm kernels" not included in Area
        'Oil palm fruit',...
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed'} ;
    FAO_to_FAO_key1{getFi('Pulses')}         = {'Pulses Total'} ;
    FAO_to_FAO_key1{getFi('Starchy roots')}  = {'Roots and Tubers Total'} ;
    % Get item names for COMMODITY BALANCE: CROPS PRIMARY EQUIVALENT dataset
    FAO_to_FAO_key2{getFi('Wheat')}          = {'Wheat and products','Barley and products','Oats'} ;
    FAO_to_FAO_key2{getFi('Maize')}          = {'Maize and products','Millet and products','Sorghum and products'} ;
    FAO_to_FAO_key2{getFi('Rice')}           = {'Rice (Paddy Equivalent)'} ;
    FAO_to_FAO_key2{getFi('Oilcrops')}       = {'Oilcrops'} ;
    FAO_to_FAO_key2{getFi('Pulses')}         = {'Pulses'} ;
    FAO_to_FAO_key2{getFi('Starchy roots')}  = {'Starchy Roots'} ;
elseif calib_ver==16 || calib_ver==18
%     error('Do this')
    listCrops_fa2o = {'Wheat','Maize','Rice','Oilcrops','Pulses','Starchy roots'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key1{getFi('Wheat')}          = {'Wheat','Barley','Rye'} ;
    FAO_to_FAO_key1{getFi('Maize')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key1{getFi('Rice')}           = {'Rice, paddy'} ;
    FAO_to_FAO_key1{getFi('Oilcrops')}       = {'Coconuts',...
        ...%'Cottonseed',...
        'Seed cotton',... % This includes cotton grown for fiber...
        'Groundnuts with shell','Karite nuts (sheanuts)',...
        'Castor oil seed','Tung nuts','Jojoba seed','Safflower seed',...
        'Poppy seed','Melonseed','Tallowtree seed','Kapok fruit',...
        ...%'Kapokseed in shell','Kapokseed shelled',... % Not included in Area
        'Linseed','Hempseed',...
        'Oilseeds nes',...
        ...%'Flour, oilseeds', % Not included in Area
        'Olives',...
        ...%'Palm kernels',... % "Palm kernels" not included in Area
        'Oil palm fruit',...
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed'} ;
    FAO_to_FAO_key1{getFi('Pulses')}         = {'Pulses Total'} ;
    FAO_to_FAO_key1{getFi('Starchy roots')}  = {'Roots and Tubers Total'} ;
    % Get item names for COMMODITY BALANCE: CROPS PRIMARY EQUIVALENT dataset
    FAO_to_FAO_key2{getFi('Wheat')}          = {'Wheat and products','Barley and products','Rye and products'} ;
    FAO_to_FAO_key2{getFi('Maize')}          = {'Maize and products','Millet and products','Sorghum and products'} ;
    FAO_to_FAO_key2{getFi('Rice')}           = {'Rice (Paddy Equivalent)'} ;
    FAO_to_FAO_key2{getFi('Oilcrops')}       = {'Oilcrops'} ;
    FAO_to_FAO_key2{getFi('Pulses')}         = {'Pulses'} ;
    FAO_to_FAO_key2{getFi('Starchy roots')}  = {'Starchy Roots'} ;
elseif calib_ver==17
%     error('Do this')
    listCrops_fa2o = {'CerealsC3','CerealsC4','Rice','OtherC3'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key{getFi('CerealsC3')}          = {'Wheat','Barley','Rye'} ;
    FAO_to_FAO_key{getFi('CerealsC4')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key{getFi('Rice')}           = {'Rice, paddy'} ;
    FAO_to_FAO_key{getFi('OtherC3')}       = {'Coconuts',...
        ...%'Cottonseed',...
        'Seed cotton',... % This includes cotton grown for fiber...
        'Groundnuts, with shell','Karite nuts (sheanuts)',...
        'Castor oil seed','Tung nuts','Jojoba seed','Safflower seed',...
        'Poppy seed','Melonseed','Tallowtree seed','Kapok fruit',...
        ...%'Kapokseed in shell','Kapokseed shelled',... % Not included in Area
        'Linseed','Hempseed',...
        'Oilseeds nes',...
        ...%'Flour, oilseeds', % Not included in Area
        'Olives',...
        ...%'Palm kernels',... % "Palm kernels" not included in Area
        'Oil, palm fruit',...
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed', ...
        'Pulses,Total', ...
        'Roots and Tubers,Total'} ;
elseif calib_ver==19
%     error('Do this')
    listCrops_fa2o = {'Wheat','Maize','Rice','Oilcrops','Pulses', ...
        'Starchy roots','Sugar','DateCitGrape'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key1{getFi('Wheat')}          = {'Wheat','Barley','Rye'} ;
    FAO_to_FAO_key1{getFi('Maize')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key1{getFi('Rice')}           = {'Rice, paddy'} ;
    FAO_to_FAO_key1{getFi('Oilcrops')}       = {'Coconuts',...
        ...%'Cottonseed',...
        'Seed cotton',... % This includes cotton grown for fiber...
        'Groundnuts with shell','Karite nuts (sheanuts)',...
        'Castor oil seed','Tung nuts','Jojoba seed','Safflower seed',...
        'Poppy seed','Melonseed','Tallowtree seed','Kapok fruit',...
        ...%'Kapokseed in shell','Kapokseed shelled',... % Not included in Area
        'Linseed','Hempseed',...
        'Oilseeds nes',...
        ...%'Flour, oilseeds', % Not included in Area
        'Olives',...
        ...%'Palm kernels',... % "Palm kernels" not included in Area
        'Oil palm fruit',...
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed'} ;
    FAO_to_FAO_key1{getFi('Pulses')}         = {'Pulses Total'} ;
    FAO_to_FAO_key1{getFi('Starchy roots')}  = {'Roots and Tubers Total'} ;
    FAO_to_FAO_key1{getFi('Sugar')}          = {'Sugar beet', 'Sugar cane'} ;
        % Sugar (production) could also include 'Sugar crops nes'
    FAO_to_FAO_key1{getFi('DateCitGrape')}    = {'Dates', 'Citrus Fruit Total', 'Grapes'} ;
    % Get item names for COMMODITY BALANCE: CROPS PRIMARY EQUIVALENT dataset
    FAO_to_FAO_key2{getFi('Wheat')}          = {'Wheat and products','Barley and products','Rye and products'} ;
    FAO_to_FAO_key2{getFi('Maize')}          = {'Maize and products','Millet and products','Sorghum and products'} ;
    FAO_to_FAO_key2{getFi('Rice')}           = {'Rice (Paddy Equivalent)'} ;
    FAO_to_FAO_key2{getFi('Oilcrops')}       = {'Oilcrops'} ;
    FAO_to_FAO_key2{getFi('Pulses')}         = {'Pulses'} ;
    FAO_to_FAO_key2{getFi('Starchy roots')}  = {'Starchy Roots'} ;
    FAO_to_FAO_key2{getFi('Sugar')}          = {'Sugar beet', 'Sugar cane'} ;
        % Sugar (commodity balance) could also include:
        %    Sugar & Sweeteners
        %    Sugar (Raw Equivalent)
        %    Sugar Crops
        %    Sugar Raw Equivalent
        %    Sugar Refined Equiv
        %    Sugar beet
        %    Sugar cane
        %    Sugar non-centrifugal
%     FAO_to_FAO_key2{getFi('DateCitGrape')}    = {'Dates', ...
%         'Citrus Other','Grapefruit and products','Lemons Limes and products','Oranges Mandarines', ...
%         'Grapes and products (excl wine)', 'Wine'} ;
    FAO_to_FAO_key2{getFi('DateCitGrape')}    = {'Dates', ...
        'Citrus Other','Grapefruit and products','Lemons Limes and products','Oranges Mandarines', ...
        'Grapes and products (excl wine)'} ;
elseif calib_ver==20
%     error('Do this')
    listCrops_fa2o = {'Wheat','Maize','Rice','Oilcrops','Pulses', ...
        'Starchy roots','Sugar','FruitAndVeg'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key1{getFi('Wheat')}          = {'Wheat','Barley','Rye'} ;
    FAO_to_FAO_key1{getFi('Maize')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key1{getFi('Rice')}           = {'Rice, paddy'} ;
    FAO_to_FAO_key1{getFi('Oilcrops')}       = {'Coconuts',...
        ...%'Cottonseed',...
        'Seed cotton',... % This includes cotton grown for fiber...
        'Groundnuts with shell','Karite nuts (sheanuts)',...
        'Castor oil seed','Tung nuts','Jojoba seed','Safflower seed',...
        'Poppy seed','Melonseed','Tallowtree seed','Kapok fruit',...
        ...%'Kapokseed in shell','Kapokseed shelled',... % Not included in Area
        'Linseed','Hempseed',...
        'Oilseeds nes',...
        ...%'Flour, oilseeds', % Not included in Area
        'Olives',...
        ...%'Palm kernels',... % "Palm kernels" not included in Area
        'Oil palm fruit',...
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed'} ;
    FAO_to_FAO_key1{getFi('Pulses')}         = {'Pulses Total'} ;
    FAO_to_FAO_key1{getFi('Starchy roots')}  = {'Roots and Tubers Total'} ;
    FAO_to_FAO_key1{getFi('Sugar')}          = {'Sugar beet', 'Sugar cane'} ;
        % Sugar (production) could also include 'Sugar crops nes'
    FAO_to_FAO_key1{getFi('FruitAndVeg')}    = {...
        'Apples', 'Apricots', 'Avocados', 'Bananas', 'Berries nes', 'Blueberries', ...
        'Carobs', 'Cashewapple', 'Cherries', 'Cherries, sour', 'Cranberries', ...
        'Currants', 'Dates', 'Figs', 'Fruit, citrus nes', 'Fruit, fresh nes', ...
        'Fruit, pome nes', 'Fruit, stone nes', 'Fruit, tropical fresh nes', ...
        'Gooseberries', 'Grapefruit (inc. pomelos)', 'Grapes', 'Kiwi fruit', ...
        'Lemons and limes', 'Mangoes, mangosteens, guavas', 'Melons, other (inc.cantaloupes)', ...
        'Oranges', 'Papayas', 'Peaches and nectarines', 'Pears', 'Persimmons', ...
        'Pineapples', 'Plantains and others', 'Plums and sloes', 'Quinces', ...
        'Raspberries', 'Strawberries', 'Tangerines, mandarins, clementines, satsumas', 'Watermelons' ...
        'Artichokes', 'Asparagus', 'Beans, green', 'Cabbages and other brassicas', ...
        'Carrots and turnips', 'Cassava leaves', 'Cauliflowers and broccoli', ...
        'Chillies and peppers, green', 'Cucumbers and gherkins', 'Eggplants (aubergines)', ...
        'Garlic', 'Leeks, other alliaceous vegetables', 'Lettuce and chicory', 'Maize, green', ...
        'Mushrooms and truffles', 'Okra', 'Onions, dry', 'Onions, shallots, green', ...
        'Peas, green', 'Pumpkins, squash and gourds', 'Spinach', 'String beans', 'Tomatoes', ...
        'Vegetables, fresh nes', 'Vegetables, leguminous nes' ...
        } ;
    % Get item names for COMMODITY BALANCE: CROPS PRIMARY EQUIVALENT dataset
    FAO_to_FAO_key2{getFi('Wheat')}          = {'Wheat and products','Barley and products','Rye and products'} ;
    FAO_to_FAO_key2{getFi('Maize')}          = {'Maize and products','Millet and products','Sorghum and products'} ;
    FAO_to_FAO_key2{getFi('Rice')}           = {'Rice (Paddy Equivalent)'} ;
    FAO_to_FAO_key2{getFi('Oilcrops')}       = {'Oilcrops'} ;
    FAO_to_FAO_key2{getFi('Pulses')}         = {'Pulses'} ;
    FAO_to_FAO_key2{getFi('Starchy roots')}  = {'Starchy Roots'} ;
    FAO_to_FAO_key2{getFi('Sugar')}          = {'Sugar beet', 'Sugar cane'} ;
        % Sugar (commodity balance) could also include:
        %    Sugar & Sweeteners
        %    Sugar (Raw Equivalent)
        %    Sugar Crops
        %    Sugar Raw Equivalent
        %    Sugar Refined Equiv
        %    Sugar beet
        %    Sugar cane
        %    Sugar non-centrifugal
    FAO_to_FAO_key2{getFi('FruitAndVeg')}    = {'Fruits - Excluding Wine','Vegetables'} ;
elseif calib_ver==21
    listCrops_fa2o = {'Faba bean','Sorghum','Soybean'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    FAO_to_FAO_key{getFi('Faba bean')}         = {'Pulses Total'} ;
%     FAO_to_FAO_key{getFi('Faba bean')}         = {'Broad beans horse beans dry','Vegetables leguminous nes'} ;
    FAO_to_FAO_key{getFi('Sorghum')}           = {'Sorghum'} ;
    FAO_to_FAO_key{getFi('Soybean')}           = {'Soybeans'} ;
elseif calib_ver==22
    % GGCMI phase 3 / ISIMIP3
    listCrops_fa2o = {'bea', 'mai', 'mil', 'ric', 'sor', 'soy', 'whe'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    FAO_to_FAO_key{getFi('bea')}         = {'Beans, dry'} ;
%     FAO_to_FAO_key{getFi('bea')}         = {'Beans, dry', 'Beans, green'} ;
    FAO_to_FAO_key{getFi('mai')}           = {'Maize'} ;
%     FAO_to_FAO_key{getFi('mai')}           = {'Maize', 'Maize, green'} ;
    FAO_to_FAO_key{getFi('mil')}           = {'Millet'} ;
    FAO_to_FAO_key{getFi('ric')}           = {'Rice, paddy'} ;
%     FAO_to_FAO_key{getFi('ric')}           = {'Rice, paddy (rice milled equivalent'} ;
    FAO_to_FAO_key{getFi('sor')}           = {'Sorghum'} ;
    FAO_to_FAO_key{getFi('soy')}           = {'Soybeans'} ;
    FAO_to_FAO_key{getFi('whe')}           = {'Wheat'} ;
elseif calib_ver==23
%     error('Do this')
    listCrops_fa2o = {'CerealsC3','CerealsC4','Rice','OilNfix','OilOther','Pulses', ...
        'StarchyRoots','Sugarbeet','Sugarcane','FruitAndVeg'} ;
    getFi = @(x)find(strcmp(listCrops_fa2o,x)) ;
    % Get item names for PRODUCTION: CROPS dataset
    FAO_to_FAO_key1{getFi('CerealsC3')}          = {'Wheat','Barley','Rye'} ;
    FAO_to_FAO_key1{getFi('CerealsC4')}          = {'Maize','Millet','Sorghum'} ;
    FAO_to_FAO_key1{getFi('Rice')}           = {'Rice, paddy'} ;
    allOilcrops_key1 = { ...
        'Coconuts',...
        ...%'Cottonseed',...
        'Seed cotton',... % This includes cotton grown for fiber...
        'Groundnuts with shell','Karite nuts (sheanuts)',...
        'Castor oil seed','Tung nuts','Jojoba seed','Safflower seed',...
        'Poppy seed','Melonseed','Tallowtree seed','Kapok fruit',...
        ...%'Kapokseed in shell','Kapokseed shelled',... % Not included in Area
        'Linseed','Hempseed',...
        'Oilseeds nes',...
        ...%'Flour, oilseeds', % Not included in Area
        'Olives',...
        ...%'Palm kernels',... % "Palm kernels" not included in Area
        'Oil palm fruit',...
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed'} ;
    FAO_to_FAO_key1{getFi('OilNfix')}        = { ...
        'Groundnuts with shell', 'Soybeans'} ;
    FAO_to_FAO_key1{getFi('OilOther')}        = setdiff( ...
        allOilcrops_key1, FAO_to_FAO_key1{getFi('OilNfix')}) ;
        
    FAO_to_FAO_key1{getFi('Pulses')}         = {'Pulses Total'} ;
    FAO_to_FAO_key1{getFi('StarchyRoots')}  = {'Roots and Tubers Total'} ;
    FAO_to_FAO_key1{getFi('Sugarbeet')}          = {'Sugar beet'} ;
    FAO_to_FAO_key1{getFi('Sugarcane')}          = {'Sugar cane'} ;
        % Sugar (production) could also include 'Sugar crops nes'
    FAO_to_FAO_key1{getFi('FruitAndVeg')}    = {...
        'Apples', 'Apricots', 'Avocados', 'Bananas', 'Berries nes', 'Blueberries', ...
        'Carobs', 'Cashewapple', 'Cherries', 'Cherries, sour', 'Cranberries', ...
        'Currants', 'Dates', 'Figs', 'Fruit, citrus nes', 'Fruit, fresh nes', ...
        'Fruit, pome nes', 'Fruit, stone nes', 'Fruit, tropical fresh nes', ...
        'Gooseberries', 'Grapefruit (inc. pomelos)', 'Grapes', 'Kiwi fruit', ...
        'Lemons and limes', 'Mangoes, mangosteens, guavas', 'Melons, other (inc.cantaloupes)', ...
        'Oranges', 'Papayas', 'Peaches and nectarines', 'Pears', 'Persimmons', ...
        'Pineapples', 'Plantains and others', 'Plums and sloes', 'Quinces', ...
        'Raspberries', 'Strawberries', 'Tangerines, mandarins, clementines, satsumas', 'Watermelons' ...
        'Artichokes', 'Asparagus', 'Beans, green', 'Cabbages and other brassicas', ...
        'Carrots and turnips', 'Cassava leaves', 'Cauliflowers and broccoli', ...
        'Chillies and peppers, green', 'Cucumbers and gherkins', 'Eggplants (aubergines)', ...
        'Garlic', 'Leeks, other alliaceous vegetables', 'Lettuce and chicory', 'Maize, green', ...
        'Mushrooms and truffles', 'Okra', 'Onions, dry', 'Onions, shallots, green', ...
        'Peas, green', 'Pumpkins, squash and gourds', 'Spinach', 'String beans', 'Tomatoes', ...
        'Vegetables, fresh nes', 'Vegetables, leguminous nes' ...
        } ;
    % Get item names for COMMODITY BALANCE: CROPS PRIMARY EQUIVALENT dataset
    FAO_to_FAO_key2{getFi('CerealsC3')}          = {'Wheat and products','Barley and products','Rye and products'} ;
    FAO_to_FAO_key2{getFi('CerealsC4')}          = {'Maize and products','Millet and products','Sorghum and products'} ;
    FAO_to_FAO_key2{getFi('Rice')}           = {'Rice (Paddy Equivalent)'} ;
    allOilcrops_key2 = { ...
        'Coconuts - Incl Copra', 'Cottonseed', 'Groundnuts (Shelled Eq)', ...
        'Oilcrops, Other', 'Olives (including preserved)', 'Palm kernels', ...
        'Rape and Mustardseed', 'Sesame seed', 'Soyabeans', 'Sunflower seed'} ;
    FAO_to_FAO_key2{getFi('OilNfix')}        = { ...
        'Groundnuts (Shelled Eq)', 'Soyabeans'} ;
    FAO_to_FAO_key2{getFi('OilOther')}       = setdiff( ...
        allOilcrops_key2, FAO_to_FAO_key2{getFi('OilNfix')}) ;
    FAO_to_FAO_key2{getFi('Pulses')}         = {'Pulses'} ;
    FAO_to_FAO_key2{getFi('StarchyRoots')}  = {'Starchy Roots'} ;
    FAO_to_FAO_key2{getFi('Sugarbeet')}          = {'Sugar beet'} ;
    FAO_to_FAO_key2{getFi('Sugarcane')}          = {'Sugar cane'} ;
        % Sugar (commodity balance) could also include:
        %    Sugar & Sweeteners
        %    Sugar (Raw Equivalent)
        %    Sugar Crops
        %    Sugar Raw Equivalent
        %    Sugar Refined Equiv
        %    Sugar beet
        %    Sugar cane
        %    Sugar non-centrifugal
    FAO_to_FAO_key2{getFi('FruitAndVeg')}    = {'Fruits - Excluding Wine','Vegetables'} ;
else
    error(['calib_ver not recognized: ' num2str(calib_ver)])
end

% Get list of all crops needed from FAO data, 2017-04-10
if ~twofiles
    listCrops_fa2i = {} ;
    for c1 = 1:length(listCrops_fa2o)
        thisCell = FAO_to_FAO_key{c1} ;
        if ~iscell(thisCell)
            listCrops_fa2i{end+1} = thisCell ;
        else
            for c2 = 1:length(thisCell)
                listCrops_fa2i{end+1} = thisCell{c2} ;
            end ; clear c2 thisCrop
        end
        clear thisCell
    end ; clear c1
else
    listCrops_fa2i1 = {} ;
    for c1 = 1:length(listCrops_fa2o)
        thisCell = FAO_to_FAO_key1{c1} ;
        if ~iscell(thisCell)
            listCrops_fa2i1{end+1} = thisCell ;
        else
            for c2 = 1:length(thisCell)
                listCrops_fa2i1{end+1} = thisCell{c2} ;
            end ; clear c2 thisCrop
        end
        clear thisCell
    end ; clear c1
    listCrops_fa2i2 = {} ;
    for c1 = 1:length(listCrops_fa2o)
        thisCell = FAO_to_FAO_key2{c1} ;
        if ~iscell(thisCell)
            listCrops_fa2i2{end+1} = thisCell ;
        else
            for c2 = 1:length(thisCell)
                listCrops_fa2i2{end+1} = thisCell{c2} ;
            end ; clear c2 thisCrop
        end
        clear thisCell
    end ; clear c1
end


end