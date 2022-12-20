function map_run_diffs( ...
    maps_d1, maps_d9, title_text, sumvars, ...
    do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, land_area_YX, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim)

% Get missing info
Nruns = length(runList) ;
if isempty(land_area_YX)
    land_area_YX = ones(size(maps_d1.maps_YXvyr,1),size(maps_d1.maps_YXvyr,2)) ;
end

figure('Color','w','Position',thisPos) ;

clim_max = 0 ;
hs = [] ;
for r = 1:Nruns
    hs(r) = subplot_tight(ny,nx,r,spacing) ;
    
    % Get data
    [~,IA] = intersect(maps_d1.varNames,sumvars) ;
    if isempty(IA)
        error('sumvars not found in maps_d1.varNames!')
    end
    begf_YXy = squeeze(sum(maps_d1.maps_YXvyr(:,:,IA,:,r),3)) ;
    begf_YXmean = mean(begf_YXy,3) ;
    [~,IA] = intersect(maps_d9.varNames,sumvars) ;
    if isempty(IA)
        error('sumvars not found in maps_d9.varNames!')
    end
    endf_YXy = squeeze(sum(maps_d9.maps_YXvyr(:,:,IA,:,r),3)) ;
    endf_YXmean = mean(endf_YXy,3) ;
    if do_pct
        map = (endf_YXmean - begf_YXmean) ./ begf_YXmean * 100;
        units_map = '%' ;
    else
        map = (endf_YXmean - begf_YXmean) * conv_fact_map ;
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
    set(get(hcb,'XLabel'),'String',['2090s minus 2010s (' units_map ')'])
    ht = title(['Diff. in ' title_text ', 2011-2020 to 2091-2100 (' runList{r} ')']) ;
    set(gca,'FontSize',fontSize)
    letterlabel_align0(char(r + 64),ht,do_caps) ;
    
    % Add data
    if ~isempty(conv_fact_total)
        mean_begf = conv_fact_total * nansum(nansum(begf_YXmean .* land_area_YX)) ;
        sd_begf = std(nansum(nansum(begf_YXy .* repmat(land_area_YX,[1 1 10]),1),2),...
            0,3) * conv_fact_total ;
        mean_endf = conv_fact_total * nansum(nansum(endf_YXmean .* land_area_YX)) ;
        sd_endf = std(nansum(nansum(endf_YXy .* repmat(land_area_YX,[1 1 10]),1),2),...
            0,3) * conv_fact_total ;
        text(textX,textY_1,['2010s: ' num2str(round(mean_begf)) '±' num2str(round(sd_begf)) ' ' units_total],'FontSize',fontSize-2) ;
        text(textX,textY_2,['2090s: ' num2str(round(mean_endf)) '±' num2str(round(sd_endf)) ' ' units_total],'FontSize',fontSize-2) ;
    end
end
if equalize_cbars && (~do_pct || (do_pct && isempty(pct_clim)))
    for r = 1:Nruns
        caxis(hs(r),[-this_clim_max this_clim_max]) ;
    end
end


end