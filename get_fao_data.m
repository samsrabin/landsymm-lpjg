function [total_fa2_Ccy, croparea_fa2_Ccy, yield_fa2_Ccy, ...
    fao_itemNames, ...
    listCrops_fa2o, Ncrops_fa2o, ...
    listCountries_map_present_all, ...
    is_tropical, is_xtratrop, ...
    calib_ver_used, twofiles] ...
    = get_fao_data(year1,yearN,calib_ver,...
    Ncountries, listCountries_map_present, countries_YX, countries_key)

% UNITS
%    Area harvested: ha
%    Production:     metric tons
%    Yield:          Hg/ha

check_country_names = true ;

if calib_ver<9
    error('Currently get_fao_data() is only tested for calib_ver>=9.')
end

[strip_fao_nans, fix_cotedivoire, combine_sudans,...
    combine_subChinas, ~, combine_serbmont] ...
    = get_FAOread_options(calib_ver) ;

[fao, fao1, fao2, twofiles, ...
    listCountries_map_present_all, is_tropical, is_xtratrop] = import_FAO_data(...
    calib_ver, year1, yearN, ...
    countries_YX, countries_key, listCountries_map_present, ...
    strip_fao_nans, fix_cotedivoire, combine_sudans, ...
    combine_subChinas, combine_serbmont) ;

x=1;

% Process FAO data into country-crop-year arrays
%%%verbose = false ;
% % % % % Define crop map for FAO-->PLUM (based on MIRCA mapping but with names
% % % % % changed to match those used by FAO)
%%%getPi = @(x) find(strcmp(listCrops_lpj_comb,x)) ;
% % % % FAO_to_PLUM_key{getPi('TeWW')}  = {'Wheat','Barley','Rye'} ;
% % % % FAO_to_PLUM_key{getPi('TeSW')}  = {'Pulses,Total','Potatoes','Sugar beet','Cassava','Sunflower seed','Soybeans','Groundnuts, with shell','Rapeseed'} ;
% % % % % FAO_to_PLUM_key{getPi('TeCo')}  = {'Maize','Millet','Sorghum'} ;
% % % % FAO_to_PLUM_key{getPi('TeCo')}  = {'Maize'} ;
% % % % % FAO_to_PLUM_key{getPi('TrRi')}  = {'Rice, paddy'} ;
% % % % FAO_to_PLUM_key{getPi('TrRi')}  = {'Rice, paddy'} ;
% % % % % [croparea_fao_Ccy,total_fao_Ccy,yield_fao_Ccy] = ...
% % % % %             FAO_to_Ccy(fao,FAO_to_PLUM_key,listCountries_map_present,...
% % % % %             listCrops_fao,Ncrops_lpj_comb,listYears_fao,verbose) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get raw FAO data for PLUM-relevant crops %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Crops from FAO
% listCrops_fa2i = {'Wheat','Barley','Oats','Maize','Millet','Sorghum','Rice, paddy','Oilcrops Primary','Pulses,Total','Roots and Tubers,Total'} ;

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
elseif calib_ver==16
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
    FAO_to_FAO_key{getFi('Rice')}           = {'Rice paddy'} ;
    FAO_to_FAO_key{getFi('OtherC3')}       = {'Coconuts',...
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
        'Rapeseed','Mustard seed','Sesame seed','Soybeans','Sunflower seed', ...
        'Pulses Total', ...
        'Roots and Tubers Total'} ;
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
%%%    Ncrops_fa2i = length(listCrops_fa2i) ;
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
%%%    Ncrops_fa2i1 = length(listCrops_fa2i1) ;
%%%    Ncrops_fa2i2 = length(listCrops_fa2i2) ;
end

% Finish up
Ncrops_fa2o = length(listCrops_fa2o) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get item names
if ~twofiles
    warning('Add code to get itemNames!')
    fao_itemNames = [] ;
else
    fao_itemNames.a = unique(fao1.ItemName(strcmp(fao1.ElementName,'Area harvested'))) ;
    fao_itemNames.p = unique(fao2.ItemName(strcmp(fao2.ElementName,'Production'))) ;
end

% Get _Ccy arrays
ignoreInfYield = true ;
ignoreNoData = true ;
% if ~exist('yield_fa2_Ccy','var') || calib_ver~=calib_ver_used
    verbose = false ;
    if ~twofiles
        disp('Getting FA2 _Ccy arrays...')
        [total_fa2_Ccy,croparea_fa2_Ccy,yield_fa2_Ccy] = ...
            FAO_to_Ccy2(fao,FAO_to_FAO_key,listCountries_map_present,...
            listCrops_fa2i,listCrops_fa2o,year1:yearN,verbose,...
            ignoreInfYield,ignoreNoData) ;
    else
        disp('Getting FA2 _Ccy array: Area...')
        croparea_fa2_Ccy = ...
            FAO_to_Ccy2_areaOnly(fao1,FAO_to_FAO_key1,listCountries_map_present_all,...
            listCrops_fa2i1,listCrops_fa2o,year1:yearN,verbose,...
            ignoreNoData) ;
        disp('Done.')
        disp(' ')
        disp('Getting FA2 _Ccy array: Production...')
        total_fa2_Ccy = ...
            FAO_to_Ccy2_totalOnly(fao2,FAO_to_FAO_key2,listCountries_map_present_all,...
            listCrops_fa2i2,listCrops_fa2o,year1:yearN,verbose,...
            ignoreInfYield,ignoreNoData) ;
        % Reconcile FAO NaNs
        either_fa2_nan = (isnan(total_fa2_Ccy) | isnan(croparea_fa2_Ccy)) ;
        total_fa2_Ccy(either_fa2_nan) = NaN ;
        croparea_fa2_Ccy(either_fa2_nan) = NaN ;
        % There shouldn't be any FAO production if no FAO data
        if any(total_fa2_Ccy>0 & croparea_fa2_Ccy==0)
            error('At least one cell has production without crop area.')
        end
        % Calculate yield (tDM/ha)
        yield_fa2_Ccy = total_fa2_Ccy ./ croparea_fa2_Ccy ;
        % Sanity check
        if any(isinf(yield_fa2_Ccy(:)))
            if ~ignoreInfYield
                error('At least one member of yield_fa2_Ccy is Inf!')
            else
                warning('At least one member of yield_fa2_Ccy is Inf! IGNORING')
                total_fa2_Ccy(isinf(yield_fa2_Ccy)) = NaN ;
                croparea_fa2_Ccy(isinf(yield_fa2_Ccy)) = NaN ;
                yield_fa2_Ccy(isinf(yield_fa2_Ccy)) = NaN ;
            end
        end
    end
    disp('Done.')
    calib_ver_used = calib_ver ;
% end
% 

end