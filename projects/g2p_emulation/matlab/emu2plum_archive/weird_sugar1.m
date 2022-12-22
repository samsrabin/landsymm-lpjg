% thisCrop = 'maize' ;
thisCrop = 'spring_wheat' ;
irrFert = '0010' ;

map_bl_YX = lpjgu_vector2map( ...
    data_bl_emu0.garr_xv(:,strcmp(data_bl_emu0.varNames, [thisCrop irrFert(2:end)])), ...
    [360 720], data_bl_emu0.list2map) ;
map_fu_YX = lpjgu_vector2map( ...
    max( ...
    data_fu_emu2.garr_xvt(:,strcmp(data_fu_emu2.varNames, [thisCrop irrFert]),:), ...
    [], 3), ...
    [360 720], data_bl_emu0.list2map) ;

map_YX = map_fu_YX ./ map_bl_YX ;
incl_x = ~isnan(map_YX) & ~isinf(map_YX) ;
data_x = map_YX(incl_x) ;
[B,I] = rmoutliers(data_x) ;
if ~isempty(I)
    fprintf('Removing %d outliers (min. %0.1f)\n', ...
        length(I), min(data_x(I))) ;
    data_x(I) = NaN ;
    map_YX(incl_x) = data_x ;
end

shademap()
axis equal tight
title(sprintf('max(fu2) / bl0: %s%s', strrep(thisCrop, '_', '\_'), irrFert))

%%

% thisCrop0 = 'maize'; thisCrop3 = 'CerealsC4' ;
% thisCrop0 = 'spring_wheat'; thisCrop3 = 'FruitAndVeg' ;
% thisCrop0 = 'spring_wheat'; thisCrop3 = 'Sugar' ;
thisCrop0 = 'maize'; thisCrop3 = 'Sugar' ;
irrFert = '0010' ;

map_bl_YX = lpjgu_vector2map( ...
    data_bl_emu0.garr_xv(:,strcmp(data_bl_emu0.varNames, [thisCrop0 irrFert(2:end)])), ...
    [360 720], data_bl_emu0.list2map) ;
map_fu_YX = lpjgu_vector2map( ...
    max( ...
    data_fu_emu3.garr_xvt(:,strcmp(data_fu_emu3.varNames, [thisCrop3 irrFert]),:), ...
    [], 3), ...
    [360 720], data_bl_emu0.list2map) ;

shademap(map_fu_YX ./ map_bl_YX)
axis equal tight
title(strrep(sprintf('max(fu3) / bl0: (%s/%s)%s', thisCrop3, thisCrop0, irrFert), ...
    '_', '\_'))