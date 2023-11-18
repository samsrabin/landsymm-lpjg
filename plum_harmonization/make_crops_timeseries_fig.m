function make_crops_timeseries_fig(ts_base_cy, ts_orig_cyr, ts_harm_cyr, ...
    LPJGcrops, legend_ts, yearList_baselineLU, yearList_orig, units, ...
    titleWord, fileWord, out_dir, timeseries_legend_loc)

Ncrops_lpjg = length(LPJGcrops) ;
ny = ceil(Ncrops_lpjg/2) ;

if Ncrops_lpjg <= 8
    spacing = [0.05 0.1] ;
    thisPos = [1 41 1440 764] ;
    fontSize = 14 ;
else
    spacing = [0.05 0.1] ;
    thisPos = [1 41 1440 (764/4)*ny] ;
    fontSize = 14 ;
end

figure('Color', 'w', 'Resize', 'off', 'Position', thisPos)

for v = 1:Ncrops_lpjg
    subplot_tight(ny,2,v,spacing) ;
    if ~isempty(ts_base_cy)
        plot(yearList_baselineLU,ts_base_cy(v,:),'-k','LineWidth',2) ;
    end
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList_orig,squeeze(ts_orig_cyr(v,:,:)),'--','LineWidth',1)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList_orig,squeeze(ts_harm_cyr(v,:,:)),'-','LineWidth',1)
    hold off
    title(sprintf('%s: %s', titleWord, LPJGcrops{v})) ;
    set(gca,'FontSize',fontSize)
    ylabel(units)
    if ~isempty(legend_ts)
        legend(legend_ts, 'Location', timeseries_legend_loc)
    end
end

% Save
warning off export_fig:exportgraphics
export_fig(sprintf('%stimeSeries_%s.pdf', out_dir, fileWord))
close

end