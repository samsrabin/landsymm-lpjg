function map_hotspot_diffs(...
    hotspot_area_YXB, hotspot_diff_YXr, hotspot_YX, hotspot_shp, ...
    spacing, latlim, edgecolor, cbarOrient, fontSize, ...
    textX, textY_1, textY_2, ssp_plot_index, lineWidth, ...
    yearList_baseline, yearList_future, runList, ...
    conv_fact_map, conv_fact_total, units_map, units_total, do_caps)

bground = nan(size(hotspot_diff_YXr(:,:,1))) ;
bground(~isnan(mean(hotspot_diff_YXr,3))) = 1 ;

lat = -89.75:0.5:89.75+0.5 ;
lon = -179.75:0.5:179.75+0.5 ;
lonlim = [-180,180];

Nruns = length(runList) ;

figure('Color','w','Position',figurePos) ;
h1 = subplot_tight(2,3,1,spacing) ;
hotspot_BL_2map = conv_fact_map*hotspot_area_YXB ;
hotspot_BL_2map(hotspot_YX==0) = NaN ;
map_with_SHPoverlay_v2(hotspot_BL_2map,[], ...
                        'bground',bground, ...
                        'lonlim', lonlim, ...
                        'latlim', latlim, ...
                        'thisColormap', 'parula', ...
                        'fontSize', fontSize, ...
                        'edgeColor', edgecolor, ...
                        'lineWidth', lineWidth, ...
                        'cbarOrient', cbarOrient, ...
                        'units_map', units_map) ;
total_bl = conv_fact_total*nansum(nansum(hotspot_area_YXB)) ;
t0 = text(textX,textY_1,sprintf('%0.1f %s', total_bl, units_total),'FontSize',fontSize+2) ;
set(t0,'Position',[textX/720 textY_1/360 0],'Units','normalized')
set(t0,'Position',[textX/720 textY_1/360 0],'Units','normalized')
ht = title(['BD hotspot area, ' num2str(yearList_baseline(end)) ' (km^2)']) ;
letterlabel_align0(char(1 + 64),ht,do_caps) ;
colorlim = 0 ;

caxis_lims = max(max(max(conv_fact_map*abs(hotspot_diff_YXr)))) ;
caxis_lims = [-caxis_lims caxis_lims] ;

for r = 1:Nruns
    subplot_tight(2,3,ssp_plot_index(r),spacing) ;
    map_with_SHPoverlay_v2(conv_fact_map*hotspot_diff_YXr(:,:,r),hotspot_shp, ...
                'bground',bground,...
                'lonlim', lonlim, ...
                'latlim', latlim, ...
                'thisColormap', 'rdbu_ssr', ...
                'fontSize', fontSize, ...
                'edgeColor', edgecolor, ...
                'lineWidth', lineWidth, ...
                'cbarOrient', cbarOrient, ...
                'flip', true, ...
                'caxis_lims', caxis_lims, ...
                'units_map', units_map) ;
%     title(['\Delta BD hotspot area, ' num2str(yearList_future(end)) ': ' runList{r} ' (km^2)'])
    ht = title(runList{r}) ;
    letterlabel_align0(char(ssp_plot_index(r) + 64),ht,do_caps) ;
    total_yr = conv_fact_total*nansum(nansum(hotspot_diff_YXr(:,:,r))) ;
    if total_yr>0
        t1 = text(textX,textY_1,sprintf('+%0.1f %s', total_yr, units_total),'FontSize',fontSize+2) ;
    else
        t1 = text(textX,textY_1,sprintf('%0.1f %s', total_yr, units_total),'FontSize',fontSize+2) ;
    end
    set(t1,'Position',[textX/720 textY_1/360 0],'Units','normalized')
    set(t1,'Position',[textX/720 textY_1/360 0],'Units','normalized')
    pctDiff = round(100*total_yr./total_bl,1) ;
    pctDiff_str = num2str(pctDiff) ;
    if pctDiff>0
        pctDiff_str = ['+' pctDiff_str] ;
    end
    t2 = text(textX,textY_2,['(' pctDiff_str ' %)'],'FontSize',fontSize+2) ;
    set(t2,'Position',[textX/720 textY_2/360 0],'Units','normalized')
    set(t2,'Position',[textX/720 textY_2/360 0],'Units','normalized')
    % Set up for changing color axis limits
    gcas{r} = gca ;
    colorlim = max(colorlim, max(abs(caxis))) ;
end
for r = 1:Nruns
    caxis(gcas{r},[-colorlim colorlim])
end

% colormap(h1,'parula')


end
