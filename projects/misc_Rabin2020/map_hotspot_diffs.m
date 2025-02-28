function map_hotspot_diffs(...
    hotspot_area_YXB, hotspot_diff_YXr, hotspot_YX, hotspot_shp, cslf_shp, ...
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

multiple_baselines = ndims(hotspot_area_YXB)==3 ;
if multiple_baselines
    % Every run has a different baseline, so don't bother plotting (a).
    % Do, however, calculate baseline totals
    hotspot_area_YXrB = hotspot_area_YXB ;
    hotspot_area_YXB = nan(size(hotspot_YX)) ;
    total_bl_r = nan(Nruns,1) ;
    for r = 1:Nruns
        total_bl_r(r) = conv_fact_total*nansum(nansum(hotspot_area_YXrB(:,:,r))) ;
    end
elseif ndims(hotspot_area_YXB)~=2
    error('ndims(hotspot_area_YXB) = %d', ndims(hotspot_area_YXB))
end

figure('Color','w','Position',[1   211   720   594]) ;
% h1 = subplot_tight(2,3,1,spacing) ;
% h1 = subplot_tight(3,1,1,spacing) ;
h1 = subplot_tight(3,1,1,spacing.*[1 6]) ;
hotspot_BL_2map = conv_fact_map*hotspot_area_YXB ;
hotspot_BL_2map(hotspot_YX==0) = NaN ;
[~,hcb_bl] = map_with_SHPoverlay_v2(hotspot_BL_2map,hotspot_shp,...[], ...
                        'bground',bground, ...
                        'lonlim', lonlim, ...
                        'latlim', latlim, ...
                        'thisColormap', 'parula', ...
                        'fontSize', fontSize, ...
                        'edgeColor', edgecolor, ...
                        'lineWidth', lineWidth, ...
                        'cbarOrient', cbarOrient, ...
                        'units_map', units_map, ...
                        'shapedata2', cslf_shp) ;
if ~multiple_baselines
    total_bl = conv_fact_total*nansum(nansum(hotspot_area_YXB)) ;
%     t0 = text(textX,textY_1,sprintf('%0.1f %s', total_bl, units_total),'FontSize',fontSize+2) ;
    t0 = text(textX,textY_1,sprintf('%0.1f %s', total_bl, units_total),'FontSize',fontSize) ;
    set(t0,'Position',[textX/720 textY_1/360+0.2 0],'Units','normalized')
    set(t0,'Position',[textX/720 textY_1/360+0.2 0],'Units','normalized')
    ht = title(['BD hotspot area, ' num2str(yearList_baseline(end)) ' (km^2)']) ;
    letterlabel_align0(char(1 + 64),ht,do_caps) ;
end
colorlim = 0 ;

caxis_lims = max(max(max(conv_fact_map*abs(hotspot_diff_YXr)))) ;
caxis_lims = [-caxis_lims caxis_lims] ;

for r = 1:Nruns
%     subplot_tight(2,3,ssp_plot_index(r),spacing) ;
    subplot_tight(3,2,r+2,spacing) ;
    shiftup = -0.05 + (ceil(r/2)-1)*0.1 ;
    [h, cb_lims] = map_with_SHPoverlay_v2(conv_fact_map*hotspot_diff_YXr(:,:,r),hotspot_shp, ...
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
                'units_map', units_map, ...
                'shapedata2', cslf_shp, ...
                'shiftup', shiftup) ;
%     title(['\Delta BD hotspot area, ' num2str(yearList_future(end)) ': ' runList{r} ' (km^2)'])
    ht = title(runList{r}) ;
%     letterlabel_align0(char(ssp_plot_index(r) + 64),ht,do_caps) ;
    letterlabel_align0(char(r+1 + 64),ht,do_caps) ;
    total_yr = conv_fact_total*nansum(nansum(hotspot_diff_YXr(:,:,r))) ;
    if total_yr>0
        t1 = text(textX,textY_1,sprintf('+%0.1f %s', total_yr, units_total),'FontSize',fontSize) ;
    else
        t1 = text(textX,textY_1,sprintf('%0.1f %s', total_yr, units_total),'FontSize',fontSize) ;
    end
    set(t1,'Position',[textX/720 textY_1/360+0.2 0],'Units','normalized')
    set(t1,'Position',[textX/720 textY_1/360+0.2 0],'Units','normalized')
    if multiple_baselines
        pctDiff = round(100*total_yr./total_bl_r(r),1) ;
    else
        pctDiff = round(100*total_yr./total_bl,1) ;
    end
    pctDiff_str = num2str(pctDiff) ;
    if pctDiff>0
        pctDiff_str = ['+' pctDiff_str] ;
    end
%     t2 = text(textX,textY_2,['(' pctDiff_str ' %)'],'FontSize',fontSize+2) ;
    t2 = text(textX,textY_2,['(' pctDiff_str ' %)'],'FontSize',fontSize) ;
    set(t2,'Position',[textX/720 textY_2/360+0.13 0],'Units','normalized')
    set(t2,'Position',[textX/720 textY_2/360+0.13 0],'Units','normalized')
    
    % Set up for changing color axis limits
    gcas{r} = gca ;
    colorlim = max(colorlim, max(abs(caxis))) ;
end
for r = 1:Nruns
    caxis(gcas{r},[-colorlim colorlim])
end

thisPos = get(h,'Position') ;
hcb = colorbar(h,'Location','SouthOutside') ;
set(hcb,'XLim',cb_lims) ;
set(h,'Position',thisPos)
hcb.Position(2) = 0.075 ;
ylabel(hcb, units_map)
hcb.FontSize = fontSize ;
hcb.TickDirection = 'out' ;
x=1 ;
thisPos = hcb_bl.Position ;
thisPos(2) = 0.15 ;
hcb.Position = thisPos ;


% colormap(h1,'parula')


end
