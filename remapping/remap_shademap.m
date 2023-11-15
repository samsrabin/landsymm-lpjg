function remap_shademap(map_YX, fig_title, fig_outfile)

figure(...
    'Position', figurePos, ...
    'Color', 'w' ...
    )
pcolor(map_YX) ;
shading flat ;
axis equal tight off ;
caxis([0 1])
colorbar

set(gca, ...
    'FontSize', 14 ...
    )

colormap(gca, 'cool')

fig_title = replace(fig_title, '_', '\_') ;
title(fig_title)

out_res = 150 ;
out_res_str = sprintf('-r%d', out_res) ;
export_fig(fig_outfile, out_res_str)
close

end