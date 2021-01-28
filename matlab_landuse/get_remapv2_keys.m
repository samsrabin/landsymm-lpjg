function [THISlist_cropsCombined_out, THISin2out_keyCombined_frac, ...
    THISlist_ignore_frac, ...
    out_names, out_keys, out_ignores, out_ignore_types] ...
    = get_remapv2_keys(in_version)

out_names = {'20180105b', '20180206', '20180210', '20180212', ...
    '20180214', '20190216', ...
    '20180301ani', ...
    'WithFruitVegSugar_a', 'WithFruitVegSugar_b', ...
    'WithFruitVeg_sepSugar', ...
    'jianyong01', 'jianyong01b', ...
    'ani01', '20200928'} ;
if ~any(strcmp(out_names,in_version))
    error(['thisVer ' in_version ' not recognized!'])
end

% 20180105b
thisOne = '20180105b' ;
list_cropsCombined_out = {'TeWW','TeSW','TeCo','TrRi'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('TeWW')}  =  {'Wheat','Barley'} ;
in2out_keyCombined_frac{getOci('TeSW')}  =  {'Pulses','Potatoes','Sugarbeet','Cassava','Sunflower','Soybeans','GroundnutsPeanuts','RapeseedCanola'} ;
in2out_keyCombined_frac{getOci('TeCo')}  =  {'Maize','Millet','Sorghum'} ;
in2out_keyCombined_frac{getOci('TrRi')}  =  {'Rice'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Sugarcane';'Oilpalm';'Citrus';'Datepalm';'GrapesVine';
    'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials';
    'FodderGrasses';'Rye'} ;

% 20180206
thisOne = '20180206' ;
list_cropsCombined_out = {'Wheat','Maize','Sorghum','Rice','Soybeans','Pulses'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('Wheat')}  = {'Wheat'} ;
in2out_keyCombined_frac{getOci('Maize')}  = {'Maize'} ;
in2out_keyCombined_frac{getOci('Sorghum')}  = {'Sorghum'} ;
in2out_keyCombined_frac{getOci('Rice')}  = {'Rice'} ;
in2out_keyCombined_frac{getOci('Soybeans')}  = {'Soybeans'} ;
in2out_keyCombined_frac{getOci('Pulses')}  = {'Pulses'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Sugarcane';'Oilpalm';'Citrus';'Datepalm';'GrapesVine';
    'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials';
    'FodderGrasses';'Rye';'Barley';'Millet';
    'Sunflower';'GroundnutsPeanuts';'RapeseedCanola';
    'Potatoes';'Sugarbeet';'Cassava'} ;

% 20180210
thisOne = '20180210' ;
list_cropsCombined_out = {'CerealsC3','CerealsC4','Rice','Oilcrops','Pulses','StarchyRoots'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('CerealsC3')}   = {'Wheat','Barley'} ;
in2out_keyCombined_frac{getOci('CerealsC4')}   = {'Maize','Millet','Sorghum'} ;
in2out_keyCombined_frac{getOci('Rice')}        = {'Rice'} ;
in2out_keyCombined_frac{getOci('Oilcrops')}    = {'Sunflower','Soybeans','GroundnutsPeanuts','RapeseedCanola'} ;
in2out_keyCombined_frac{getOci('Pulses')}      = {'Pulses'} ;
in2out_keyCombined_frac{getOci('StarchyRoots')}= {'Potatoes','Sugarbeet','Cassava'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Sugarcane';'Oilpalm';'Citrus';'Datepalm';'GrapesVine';
    'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials';
    'FodderGrasses';'Rye'} ;

% 20180212
thisOne = '20180212' ;
list_cropsCombined_out = {'Wheat','Maize','Sorghum','Rice'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('Wheat')}  = {'Wheat'} ;
in2out_keyCombined_frac{getOci('Maize')}  = {'Maize'} ;
in2out_keyCombined_frac{getOci('Sorghum')}  = {'Sorghum'} ;
in2out_keyCombined_frac{getOci('Rice')}  = {'Rice'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Sugarcane';'Oilpalm';'Citrus';'Datepalm';'GrapesVine';
    'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials';
    'FodderGrasses';'Rye';'Barley';'Millet';
    'Sunflower';'GroundnutsPeanuts';'RapeseedCanola';
    'Potatoes';'Sugarbeet';'Cassava';
    'Soybeans';'Pulses'} ;

% 20180214
thisOne = '20180214' ;
list_cropsCombined_out = {'CerealsC3','CerealsC4','Rice','Oilcrops','Pulses','StarchyRoots'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('CerealsC3')}   = {'Wheat','Barley','Rye'} ;
in2out_keyCombined_frac{getOci('CerealsC4')}   = {'Maize','Millet','Sorghum'} ;
in2out_keyCombined_frac{getOci('Rice')}        = {'Rice'} ;
in2out_keyCombined_frac{getOci('Oilcrops')}    = {'Sunflower','Soybeans','GroundnutsPeanuts','RapeseedCanola','Oilpalm'} ;
in2out_keyCombined_frac{getOci('Pulses')}      = {'Pulses'} ;
in2out_keyCombined_frac{getOci('StarchyRoots')}= {'Potatoes','Sugarbeet','Cassava'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Sugarcane';'Citrus';'Datepalm';'GrapesVine';
    'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials';
    'FodderGrasses'} ;

% 20190216
%%% As 20180214, but with Sugarbeet ignored instead of in StarchyRoots
thisOne = '20190216' ;
list_cropsCombined_out = {'CerealsC3','CerealsC4','Rice','Oilcrops','Pulses','StarchyRoots'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('CerealsC3')}   = {'Wheat','Barley','Rye'} ;
in2out_keyCombined_frac{getOci('CerealsC4')}   = {'Maize','Millet','Sorghum'} ;
in2out_keyCombined_frac{getOci('Rice')}        = {'Rice'} ;
in2out_keyCombined_frac{getOci('Oilcrops')}    = {'Sunflower','Soybeans','GroundnutsPeanuts','RapeseedCanola','Oilpalm'} ;
in2out_keyCombined_frac{getOci('Pulses')}      = {'Pulses'} ;
in2out_keyCombined_frac{getOci('StarchyRoots')}= {'Potatoes','Cassava'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Sugarcane';'Citrus';'Datepalm';'GrapesVine';
    'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials';
    'FodderGrasses';'Sugarbeet'} ;

% 20180301ani
thisOne = '20180301ani' ;
list_cropsCombined_out = {'TeWW','TeSW','TeCo','TrRi'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('TeWW')} = {'Wheat','Barley','Rye'} ;
in2out_keyCombined_frac{getOci('TeCo')} = {'Maize','Millet','Sorghum'} ;
in2out_keyCombined_frac{getOci('TrRi')} = {'Rice'} ;
in2out_keyCombined_frac{getOci('TeSW')} = {'Sunflower','Soybeans','GroundnutsPeanuts','RapeseedCanola','Oilpalm','Pulses','Potatoes','Sugarbeet','Cassava'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Sugarcane';'Citrus';'Datepalm';'GrapesVine';
    'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials';
    'FodderGrasses'} ;

% WithFruitVegSugar_a
% NOTE: Sugarbeet moved from Starchy Roots to Sugar. Should have never been
% included in Starchy Roots, probably, especially considering that it is
% its own crop in the FAO data (i.e., not lumped in with Roots And Tubers
% Total, which was the only item in the calibration for v16/v18
% previously).
thisOne = 'WithFruitVegSugar_a' ;
list_cropsCombined_out = {'CerealsC3','CerealsC4','Rice','Oilcrops','Pulses','StarchyRoots','DateCitGrape','Sugar'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('CerealsC3')}   = {'Wheat','Barley','Rye'} ;
in2out_keyCombined_frac{getOci('CerealsC4')}   = {'Maize','Millet','Sorghum'} ;
in2out_keyCombined_frac{getOci('Rice')}        = {'Rice'} ;
in2out_keyCombined_frac{getOci('Oilcrops')}    = {'Sunflower','Soybeans','GroundnutsPeanuts','RapeseedCanola','Oilpalm'} ;
in2out_keyCombined_frac{getOci('Pulses')}      = {'Pulses'} ;
in2out_keyCombined_frac{getOci('StarchyRoots')}= {'Potatoes','Cassava'} ;
in2out_keyCombined_frac{getOci('DateCitGrape')} = {'Datepalm','Citrus','GrapesVine'} ;
in2out_keyCombined_frac{getOci('Sugar')}       = {'Sugarbeet','Sugarcane'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Cotton';'Cocoa';'Coffee';'FodderGrasses';
     'OtherAnnuals';'OtherPerennials'} ;

% WithFruitVegSugar_b
% As WithFruitVegSugar_a, but with Other Annual and Other Perennial included in
% FruitAndVeg.
% Other annual:
%   Almonds; Apples; Apricots; Areca Nuts (Betel); Avocados; Bananas; 
%   Blueberries; Brazil Nuts; Carobs; Cashew Nuts; Cashewapple; Cherries; 
%   Chestnuts; Cinnamon (Canella); Cloves; Coconuts; Coir; Cranberries; 
%   Currants; Figs; Fruit Fresh, other; Fruit Tropical Fresh, other; 
%   Ginger; Gooseberries; Hazelnuts (Filberts); Hops; Kapok Fiber; 
%   Kapokseed in Shell; Karite Nuts (Sheanuts); Kiwi Fruit; Kolanuts; 
%   Mangoes; Mate; Natural Gums; Natural Rubber; 
%   Nutmeg, Mace and Cardamons; Nuts, other; Olives; Papayas; 
%   Peaches and Nectarines; Pears; Pepper; Persimmons; Pimento; Pineapples;
%   Pistachios; Plantains; Plums; Quinces; Raspberries; Sour Cherries; 
%   Stone Fruit other; Strawberries; Tea; Tung Nuts; Vanilla; Walnuts
% Other perennial:
%   Abaca (Manila Hemp); Agave Fibers, other; Anise, Badian and Fennel; 
%   Artichokes; Asparagus; Beans, Green; Beets for Fodder; Berries, other; 
%   Broad Beans, Green; Buckwheat; Cabbage for Fodder; Cabbages; 
%   Canary Seed; Cantaloupes & other Melons; Carrots; Carrots for Fodder; 
%   Castor Beans; Cauliflower; Cereals, other; Chicory Roots; 
%   Chillies & Peppers, Green; Cucumbers and Gherkins; Eggplants; 
%   Fiber Crops, other; Flax Fiber and Tow; Fonio; Forage Products, other; 
%   Garlic; Green Corn (Maize); Green Oilseeds for Fodder; 
%   Hemp Fiber and Tow; Hempseed; Jute; Jute-Like Fibers; Legumes, other; 
%   Lettuce; Linseed; Melonseed; Mixed Grain; Mushrooms; Mustard Seed; 
%   Oats; Oilseeds, other; Okra; Onions & Shallots, Green; Onions, Dry; 
%   Peas, Green; Peppermint; Poppy Seed; Pumpkins, Squash, Gourds; 
%   Pyrethrum, Dried Flowers; Quinoa; Ramie; Roots and Tubers, other; 
%   Safflower Seed; Sesame Seed; Sisal; Spices, other; Spinach; 
%   String Beans; Sugar Crops, other; Swedes for Fodder; Sweet Potatoes; 
%   Taro; Tobacco Leaves; Tomatoes; Triticale; Turnips for Fodder; 
%   Vegetables & Roots for Fodder; Vegetables Fresh, other; Watermelons; 
%   Yams; Yautia
thisOne = 'WithFruitVegSugar_b' ;
list_cropsCombined_out = {'CerealsC3','CerealsC4','Rice','Oilcrops','Pulses','StarchyRoots','FruitAndVeg','Sugar'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('CerealsC3')}   = {'Wheat','Barley','Rye'} ;
in2out_keyCombined_frac{getOci('CerealsC4')}   = {'Maize','Millet','Sorghum'} ;
in2out_keyCombined_frac{getOci('Rice')}        = {'Rice'} ;
in2out_keyCombined_frac{getOci('Oilcrops')}    = {'Sunflower','Soybeans','GroundnutsPeanuts','RapeseedCanola','Oilpalm'} ;
in2out_keyCombined_frac{getOci('Pulses')}      = {'Pulses'} ;
in2out_keyCombined_frac{getOci('StarchyRoots')}= {'Potatoes','Cassava'} ;
in2out_keyCombined_frac{getOci('FruitAndVeg')} = {'Datepalm','Citrus','GrapesVine','OtherAnnuals','OtherPerennials'} ;
in2out_keyCombined_frac{getOci('Sugar')}       = {'Sugarbeet','Sugarcane'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Cotton';'Cocoa';'Coffee';'FodderGrasses'} ;

% WithFruitVeg_sepSugar
% As WithFruitVegSugar_b, but with Sugar separated into Sugarbeet and
% Sugarcane
thisOne = 'WithFruitVeg_sepSugar' ;
list_cropsCombined_out = {'CerealsC3','CerealsC4','Rice','Oilcrops', ...
    'Pulses','StarchyRoots','FruitAndVeg','Sugarbeet','Sugarcane'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('CerealsC3')}   = {'Wheat','Barley','Rye'} ;
in2out_keyCombined_frac{getOci('CerealsC4')}   = {'Maize','Millet','Sorghum'} ;
in2out_keyCombined_frac{getOci('Rice')}        = {'Rice'} ;
in2out_keyCombined_frac{getOci('Oilcrops')}    = {'Sunflower','Soybeans','GroundnutsPeanuts','RapeseedCanola','Oilpalm'} ;
in2out_keyCombined_frac{getOci('Pulses')}      = {'Pulses'} ;
in2out_keyCombined_frac{getOci('StarchyRoots')}= {'Potatoes','Cassava'} ;
in2out_keyCombined_frac{getOci('FruitAndVeg')} = {'Datepalm','Citrus','GrapesVine','OtherAnnuals','OtherPerennials'} ;
in2out_keyCombined_frac{getOci('Sugarbeet')}       = {'Sugarbeet'} ;
in2out_keyCombined_frac{getOci('Sugarcane')}       = {'Sugarcane'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Cotton';'Cocoa';'Coffee';'FodderGrasses'} ;

% Jianyong v01
thisOne = 'jianyong01' ;
list_cropsCombined_out = {'Wheat','Maize','Sorghum','Rice','Soybean','FabaBean'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('Wheat')}  = {'Wheat'} ;
in2out_keyCombined_frac{getOci('Maize')}  = {'Maize'} ;
in2out_keyCombined_frac{getOci('Sorghum')}  = {'Sorghum'} ;
in2out_keyCombined_frac{getOci('Rice')}  = {'Rice'} ;
in2out_keyCombined_frac{getOci('Soybean')}  = {'Soybeans'} ;
in2out_keyCombined_frac{getOci('FabaBean')}  = {'Pulses'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Sugarcane';'Oilpalm';'Citrus';'Datepalm';'GrapesVine';
    'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials';
    'FodderGrasses';'Rye';'Barley';'Millet';
    'Sunflower';'GroundnutsPeanuts';'RapeseedCanola';
    'Potatoes';'Sugarbeet';'Cassava'} ;

% Jianyong v01b
%%% Same mapping, just different MIRCA names
thisOne = 'jianyong01b' ;
list_cropsCombined_out = {'Wheat','Maize','Sorghum','Rice','Soybean','FabaBean'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('Wheat')}  = {'Wheat'} ;
in2out_keyCombined_frac{getOci('Maize')}  = {'Maize'} ;
in2out_keyCombined_frac{getOci('Sorghum')}  = {'Sorghum'} ;
in2out_keyCombined_frac{getOci('Rice')}  = {'Rice'} ;
in2out_keyCombined_frac{getOci('Soybean')}  = {'Soybean'} ;
in2out_keyCombined_frac{getOci('FabaBean')}  = {'Pulses'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Sugar cane';'Oil palm';'Citrus';'Date palm';'Grapes / Vine';
    'Cotton';'Cocoa';'Coffee';'Others annual';'Others perennial';
    'Fodder grasses';'Rye';'Barley';'Millet';
    'Sunflower';'Groundnuts / Peanuts';'Rape seed / Canola';
    'Potatoes';'Sugar beet';'Cassava'} ;

% Ani v01
thisOne = 'ani01' ;
list_cropsCombined_out = {'TeWW','TeCo','TrRi','TeSW'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('TeWW')}  = {'Wheat', 'Rye', 'Barley'} ;
in2out_keyCombined_frac{getOci('TeCo')}  = {'Maize', 'Millet', 'Sorghum', 'Sugar cane'} ;
in2out_keyCombined_frac{getOci('TrRi')}  = {'Rice'} ;
in2out_keyCombined_frac{getOci('TeSW')}  = {'Soybean','Pulses', ...
    'Oil palm','Citrus','Date palm','Grapes / Vine', ...
    'Cotton','Cocoa','Coffee','Others annual','Others perennial', ...
    'Fodder grasses', ...
    'Sunflower','Groundnuts / Peanuts','Rape seed / Canola', ...
    'Potatoes','Sugar beet','Cassava'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = {} ;

% 20200928
% Initial version for GGCMI phase 3

% in2out_keyCombined_frac{getOci('CerealsC3')}   = {'Wheat','Barley','Rye'} ;
% in2out_keyCombined_frac{getOci('CerealsC4')}   = {'Maize','Millet','Sorghum'} ;
% in2out_keyCombined_frac{getOci('Rice')}        = {'Rice'} ;
% in2out_keyCombined_frac{getOci('Oilcrops')}    = {'Sunflower','Soybeans','GroundnutsPeanuts','RapeseedCanola','Oilpalm'} ;
% in2out_keyCombined_frac{getOci('Pulses')}      = {'Pulses'} ;
% in2out_keyCombined_frac{getOci('StarchyRoots')}= {'Potatoes','Cassava'} ;
% out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
% list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
% clear in2out_keyCombined_frac
% out_ignores{strcmp(out_names,thisOne)} = ...
%     {'Sugarcane';'Citrus';'Datepalm';'GrapesVine';
%     'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials';
%     'FodderGrasses';'Sugarbeet'} ;

thisOne = '20200928' ;
list_cropsCombined_out = {'whe','mai','ric','soy','bea','sor','mil'} ;
getOci = @(x) find(strcmp(list_cropsCombined_out,x)) ;
in2out_keyCombined_frac{getOci('whe')}  = {'Wheat'} ;
in2out_keyCombined_frac{getOci('mai')}  = {'Maize'} ;
in2out_keyCombined_frac{getOci('ric')}  = {'Rice'} ;
in2out_keyCombined_frac{getOci('soy')}  = {'Soybeans'} ;
in2out_keyCombined_frac{getOci('bea')}  = {'Pulses'} ;
in2out_keyCombined_frac{getOci('sor')}  = {'Sorghum'} ;
in2out_keyCombined_frac{getOci('mil')}  = {'Millet'} ;
out_keys{strcmp(out_names,thisOne)} = in2out_keyCombined_frac ;
list_cropsCombined_out_ALL{strcmp(out_names,thisOne)} = list_cropsCombined_out ;
clear in2out_keyCombined_frac
out_ignores{strcmp(out_names,thisOne)} = ...
    {'Sugarcane';'Citrus';'Datepalm';'GrapesVine';
    'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials';
    'FodderGrasses';'Sugarbeet';'Barley';'Rye';
    'Sunflower';'GroundnutsPeanuts';'RapeseedCanola';'Oilpalm';
    'Potatoes';'Cassava'} ;

% Get ignore types
out_ignore_types = zeros(size(out_ignores)) ;
for i = 1:length(out_ignores)
    if out_ignore_types(i)==0
        out_ignore_types(i) = i ;
        thisOne = out_ignores{i} ;
        for j = (i+1):length(out_ignores)
            if isequal(thisOne,out_ignores{j})
                out_ignore_types(j) = i ;
            end
        end
    end
end

% Get maps for this version
THISlist_cropsCombined_out = list_cropsCombined_out_ALL{strcmp(out_names,in_version)} ;
THISin2out_keyCombined_frac = out_keys{strcmp(out_names,in_version)} ;
THISlist_ignore_frac = out_ignores{strcmp(out_names,in_version)} ;

end