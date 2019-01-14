function kcal_per_g = get_kcalDensity(thisCrop)

% Values taken from Franck et al. (2011, Ecological Modelling)
% Same method as Krause et al. (2017)

switch thisCrop
    case 'CerealsC3';     kcal_per_g = 3.34 ; % Wheat
    case 'CerealsC4';     kcal_per_g = mean([3.56 3.40]) ; % Maize, millet
    case 'Oilcrops';      kcal_per_g = mean([3.08 3.35 4.14 4.94]) ; % Sunflower, soybean, groundnut, rapeseed
    case 'Pulses';        kcal_per_g = 3.41 ; % Pulses
    case 'Rice';          kcal_per_g = 2.80 ; % Rice
    case 'StarchyRoots';  kcal_per_g = mean([0.70 1.09]) ; % Sugar beet, cassava
    case {'CC3G','CC4G'}; kcal_per_g = 0 ;
    otherwise ; error(['thisCrop (' thisCrop ') not recognized for calorie density!']) ;
end

end