function kcal_per_g = get_kcalDensity2(thisCrop)

% Values taken from Franck et al. (2011, Ecological Modelling)
% Same method as Krause et al. (2017)
% THIS VERSION ACTUALLY ADJUSTS DRY TO WET MATTER
% kcal_per_g = kcal_per_g_wet/dry_matter_frac

switch thisCrop
    case 'CerealsC3';     kcal_per_g = 3.34/0.88 ; % Wheat
    case 'CerealsC4';     kcal_per_g = mean([3.56/0.88 3.40/0.88]) ; % Maize, millet
    case 'Oilcrops';      kcal_per_g = mean([3.08/0.93 3.35/0.90 4.14/0.94 4.94/0.92]) ; % Sunflower, soybean, groundnut, rapeseed
    case 'Pulses';        kcal_per_g = 3.41/0.90 ; % Pulses
    case 'Rice';          kcal_per_g = 2.80/0.87 ; % Rice
    case 'StarchyRoots';  kcal_per_g = mean([0.70/0.24 1.09/0.35]) ; % Sugar beet, cassava
    case {'CC3G','CC4G'}; kcal_per_g = 0 ;
    otherwise ; error(['thisCrop (' thisCrop ') not recognized for calorie density!']) ;
end

end