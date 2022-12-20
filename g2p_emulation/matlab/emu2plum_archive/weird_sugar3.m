% map_mai_YX = get_map(data_fu_emu2, ['maize' irrFert]) ;
% map_swh_YX = get_map(data_fu_emu2, ['spring_wheat' irrFert]) ;
map_mai_YX = get_map(data_fu_emu3, ['CerealsC4' irrFert]) ;
map_swh_YX = get_map(data_fu_emu3, ['FruitAndVeg' irrFert]) ;

map_sgc_YX = map_mai_YX * cf_emu0.calibration_factor( ...
    strcmp(cf_emu0.Crop, 'Sugarcane')) ;
map_sgb_YX = map_swh_YX * cf_emu0.calibration_factor( ...
    strcmp(cf_emu0.Crop, 'Sugarbeet')) ;
map_sug_YX = max(map_sgc_YX, map_sgb_YX);

% chose_sgc_YX = double(map_sgc_YX > map_sgb_YX | isnan(map_sgb_YX)) ;
% chose_sgc_YX(isnan(map_sug_YX)) = NaN ;
% shademap(chose_sgc_YX)

tmp_sgc = map_sgc_YX ;
tmp_sgc(isnan(tmp_sgc)) = 0 ;
tmp_sgb = map_sgb_YX ;
tmp_sgb(isnan(tmp_sgb)) = 0 ;
sgc_minus_sgb_YX = tmp_sgc - tmp_sgb ;
sgc_minus_sgb_YX(isnan(map_sgc_YX) & isnan(map_sgb_YX)) = NaN ;

caxis_sug = [0 max(max(map_sug_YX))] ;

spacing = [0.025 0.025] ;
fontSize = 14 ;


figure('Color', 'w', 'Position', figurePos)

subplot_tight(2, 2, 1, spacing)
pcolor(map_mai_YX) ; shading flat ; axis equal tight off
colorbar
title('maize'); set(gca, 'FontSize', fontSize)

subplot_tight(2, 2, 2, spacing)
pcolor(map_swh_YX) ; shading flat ; axis equal tight off
colorbar
title('spring wheat'); set(gca, 'FontSize', fontSize)

subplot_tight(2, 2, 3, spacing)
pcolor(map_sgc_YX) ; shading flat ; axis equal tight off
colorbar
title('sugarcane'); set(gca, 'FontSize', fontSize)

subplot_tight(2, 2, 4, spacing)
pcolor(map_sgb_YX) ; shading flat ; axis equal tight off
colorbar
title('sugarbeet'); set(gca, 'FontSize', fontSize)

figure('Color', 'w', 'Position', figurePos)

subplot_tight(2, 2, 1, spacing)
pcolor(map_sug_YX) ; shading flat ; axis equal tight off
caxis(caxis_sug); colorbar
title('sugar'); set(gca, 'FontSize', fontSize)

subplot_tight(2, 2, 2, spacing)
pcolor(sgc_minus_sgb_YX) ; shading flat ; axis equal tight off
colorbar
title('cane minus beet'); set(gca, 'FontSize', fontSize)

subplot_tight(2, 2, 3, spacing)
pcolor(map_sgc_YX) ; shading flat ; axis equal tight off
caxis(caxis_sug); colorbar
title('sugarcane'); set(gca, 'FontSize', fontSize)

subplot_tight(2, 2, 4, spacing)
pcolor(map_sgb_YX) ; shading flat ; axis equal tight off
caxis(caxis_sug); colorbar
title('sugarbeet'); set(gca, 'FontSize', fontSize)



function map_YX = get_map(data_in, thisVar)

I = find(strcmp(data_in.varNames, thisVar)) ;
if length(I) ~= 1
    error('Expected to find 1 match of %s in data_in.varNames; found %d', ...
        thisVar, length(I))
end

data_x = max(data_in.garr_xvt(:,I), ...
    [], 3) ;

% Convert kg/m2 to tons/ha
data_x = 10 * data_x ;

map_YX = lpjgu_vector2map( ...
    data_x, ...
    [360 720], data_in.list2map) ;


end