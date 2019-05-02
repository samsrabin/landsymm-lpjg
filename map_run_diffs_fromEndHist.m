function map_run_diffs_fromEndHist( ...
    maps_d9, title_text, sumvars, ...
    do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, land_area_YX, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim, prctile_clim)

% Get missing info
Nruns = length(runList) ;
if isempty(land_area_YX)
    land_area_YX = ones(size(maps_d9.maps_YXvyr,1),size(maps_d9.maps_YXvyr,2)) ;
end

figure('Color','w','Position',thisPos) ;

[~,IA] = intersect(maps_d9.varNames,sumvars) ;
if isempty(IA)
    error('sumvars not found in maps_d9.varNames!')
end
endh_YXy = squeeze(sum(maps_d9.maps_YXvyB(:,:,IA,:),3)) ;
endh_YXmean = mean(endh_YXy,3) ;

clim_max = 0 ;
hs = [] ;
hcbs = [] ;
for r = 1:Nruns
    hs(r) = subplot_tight(ny,nx,r,spacing) ;
    
    % Get data
    endf_YXy = squeeze(sum(maps_d9.maps_YXvyr(:,:,IA,:,r),3)) ;
    endf_YXmean = mean(endf_YXy,3) ;
    if do_pct
        map = (endf_YXmean - endh_YXmean) ./ endh_YXmean * 100;
        units_map = '%' ;
    else
        map = (endf_YXmean - endh_YXmean) .* conv_fact_map ;
    end
    
    % Make plot
    pcolor(map(69:end,:)) ; shading flat ; axis equal tight off
    if do_pct
        if isempty(pct_clim)
            this_clim_max = min(this_clim_max,100) ;
        else
            this_clim_max = pct_clim ;
        end
    elseif ~isempty(prctile_clim)
        this_clim_max = prctile(abs(map(~isnan(map))), prctile_clim) ;
    else
        this_clim_max = max(abs(caxis)) ;
    end
    caxis(this_clim_max*[-1 1])
    clim_max = max(clim_max,this_clim_max) ;
    hcbs(r) = colorbar(colorBarLoc) ;
    hcb = hcbs(r) ;
    colormap(gca,brighten(brewermap(64,'RdBu_ssr'),-0.3)) ;
    
    % Add labels
    set(get(hcb,'XLabel'),'String',['2090s minus 2000s (' units_map ')'])
    ht = title(['Diff. in ' title_text ', 2001-2010 to 2091-2100 (' runList{r} ')']) ;
    set(gca,'FontSize',fontSize)
    letterlabel_align0(char(r + 64),ht,do_caps) ;
    
    % Add data
    if ~isempty(conv_fact_total)
        mean_endh = conv_fact_total * nansum(nansum(endh_YXmean .* land_area_YX)) ;
        sd_endh = std(nansum(nansum(endh_YXy .* repmat(land_area_YX,[1 1 10]),1),2),...
            0,3) * conv_fact_total ;
        mean_endf = conv_fact_total * nansum(nansum(endf_YXmean .* land_area_YX)) ;
        sd_endf = std(nansum(nansum(endf_YXy .* repmat(land_area_YX,[1 1 10]),1),2),...
            0,3) * conv_fact_total ;
%         text(textX,textY_1,['2000s: ' num2str(round(mean_endh)) '±' num2str(round(sd_endh)) ' ' units_total],'FontSize',fontSize-2) ;
%         text(textX,textY_2,['2090s: ' num2str(round(mean_endf)) '±' num2str(round(sd_endf)) ' ' units_total],'FontSize',fontSize-2) ;
        if contains(title_text, 'albedo')
            text(textX,textY_1, ...
                sprintf('2000s: %0.3f±%0.3f %s', mean_endh, sd_endh, units_total), ...
                'FontSize',fontSize-2) ;
            text(textX,textY_2, ...
                sprintf('2090s: %0.3f±%0.3f %s', mean_endf, sd_endf, units_total), ...
                'FontSize',fontSize-2) ;
        else
            text(textX,textY_1, ...
                sprintf('2000s: %d±%d %s', round(mean_endh), round(sd_endh), units_total), ...
                'FontSize',fontSize-2) ;
            text(textX,textY_2, ...
                sprintf('2090s: %d±%d %s', round(mean_endf), round(sd_endf), units_total), ...
                'FontSize',fontSize-2) ;
        end
    end
end
if equalize_cbars && (~do_pct || (do_pct && isempty(pct_clim)))
    for r = 1:Nruns
%         caxis(hs(r),[-this_clim_max this_clim_max]) ;
        caxis(hs(r),[-clim_max clim_max]) ;
        if ~do_pct && ~isempty(prctile_clim)
            hcb_ticks = get(hcbs(r),'Ticks') ;
            hcb_limits = get(hcbs(r),'Limits') ;
            hcb_ticklabels = get(hcbs(r),'TickLabels') ;
            new_ticks = [hcb_limits(1) hcb_ticks hcb_limits(2)] ;
            new_ticklabels = [{'<'};hcb_ticklabels;{'>'}] ;
            set(hcbs(r), ...
                'Ticks', new_ticks, ...
                'TickLabels',  new_ticklabels)
        end
    end
end

end