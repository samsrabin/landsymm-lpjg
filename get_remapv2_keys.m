function [THISlist_cropsCombined_out, THISin2out_keyCombined_frac, ...
    THISlist_ignore_frac, ...
    out_names, out_keys, out_ignores, out_ignore_types] ...
    = get_remapv2_keys(in_version)

out_names = {'20180105b', '20180206', '20180210', '20180212', '20180214', '20180301ani'} ;
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