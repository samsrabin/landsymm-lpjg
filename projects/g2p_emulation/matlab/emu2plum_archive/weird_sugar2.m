thisVar = 'spring_wheat010' ;
% thisVar = 'maize010' ;

doit(data_bl_agm, thisVar)
title(strrep(...
    sprintf('data_bl_agm %s', thisVar), ...
    '_', '\_'))
doit(data_bl_emu0, thisVar)
title(strrep(...
    sprintf('data_bl_emu0 %s', thisVar), ...
    '_', '\_'))
%%
tmp = data_bl_emu0 ;
tmp.garr_xv = max(deltas_emu_xvt, [], 3) ;
doit(tmp, thisVar)
title(strrep(...
    sprintf('max deltas %s', thisVar), ...
    '_', '\_'))

function doit(data_in, thisVar)

data_x = data_in.garr_xv(:,strcmp(data_in.varNames, thisVar)) ;
map_YX = lpjgu_vector2map( ...
    data_x, ...
    [360 720], data_in.list2map) ;
figure('Color', 'w', 'Position', figurePos)
pcolor(map_YX)
% pcolor(log10(map_YX))
% pcolor(map_YX > 5)
caxis([0 max(data_x(data_x<5))]) ;
shading flat ; axis equal tight off
colorbar
set(gca, 'FontSize', 14)

end