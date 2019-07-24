function kcal_per_g = get_kcalDensity2(thisCrop_orig)

% Values taken from Franck et al. (2011, Ecological Modelling)
% Same method as Krause et al. (2017)
% THIS VERSION ACTUALLY ADJUSTS DRY TO WET MATTER
% kcal_per_g = kcal_per_g_wet/dry_matter_frac

thisCrop = lower(thisCrop_orig) ;
if strcmp(thisCrop(end-2:end),'irr')
    thisCrop = thisCrop(1:end-3) ;
end

switch thisCrop
    case {'cerealsc3', 'teww', 'tesw'}
        kcal_per_g = 3.34/0.88 ; % Wheat
    case {'cerealsc4', 'teco'}
        kcal_per_g = mean([3.56/0.88 3.40/0.88]) ; % Maize, millet
    case 'cereals';       kcal_per_g = mean([3.34/0.88 3.56/0.88 3.40/0.88]) ; % Wheat, maize, millet
    case 'oilcrops';      kcal_per_g = mean([3.08/0.93 3.35/0.90 4.14/0.94 4.94/0.92]) ; % Sunflower, soybean, groundnut, rapeseed
    case 'pulses';        kcal_per_g = 3.41/0.90 ; % Pulses
    case {'rice', 'trri'}
        kcal_per_g = 2.80/0.87 ; % Rice
    case 'starchyroots';  kcal_per_g = mean([0.70/0.24 1.09/0.35]) ; % Sugar beet, cassava
    case {'cc3g','cc4g'}; kcal_per_g = 0 ;
    otherwise
        error(['thisCrop (' thisCrop_orig ') not recognized for calorie density!']) ;
end

end