function map_run_diffs_fromEndHist( ...
    maps_d9, title_text, sumvars, ...
    do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, land_area_YX, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim)

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
for r = 1:Nruns
    hs(r) = subplot_tight(ny,nx,r,spacing) ;
    
    % Get data
    endf_YXy = squeeze(sum(maps_d9.maps_YXvyr(:,:,IA,:,r),3)) ;
    endf_YXmean = mean(endf_YXy,3) ;
    if do_pct
        map = (endf_YXmean - endh_YXmean) ./ endh_YXmean * 100;
        units_map = '%' ;
    else
        map = (endf_YXmean - endh_YXmean) * conv_fact_map ;
    end
    
    % Make plot
    pcolor(map(69:end,:)) ; shading flat ; axis equal tight off
    this_clim_max = max(abs(caxis)) ;
    if do_pct
        if isempty(pct_clim)
            this_clim_max = min(this_clim_max,100) ;
        else
            caxis(pct_clim*[-1 1])
        end
    else
        caxis([-this_clim_max this_clim_max]) ;
    end
    clim_max = max(clim_max,this_clim_max) ;
    hcb = colorbar(colorBarLoc) ;
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
        caxis(hs(r),[-this_clim_max this_clim_max]) ;
    end
end


end