function maps_YXr = map_run_diffs_fromEndHist_oneCol( ...
    data_d9, title_text, sumvars, ...
    fontSize, spacing, textX, textY_1, textY_2, ...
    thisPos, colorBarLoc, runList, do_caps, land_area_YX, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, varargin)

already_maps = true ;
if ~isempty(varargin)
    if length(varargin)==2
        map_size = varargin{1} ;
        list2map = varargin{2} ;
        already_maps = false ;
    else
        error('map_run_diffs_fromEndHist_oneCol() accepts 0 or 2 optional arguments: map_size and list2map')
    end
end

% Get missing info
Nruns = length(runList) ;
if isempty(land_area_YX)
    if already_maps
        land_area_YX = ones(size(data_d9.maps_YXvyr,1),size(data_d9.maps_YXvyr,2)) ;
    else
        land_area_YX = ones(map_size) ;
    end
end

figure('Color','w','Position',thisPos) ;

[~,IA] = intersect(data_d9.varNames,sumvars) ;
if isempty(IA)
    error('sumvars not found in data_d9.varNames!')
end
if already_maps
    endh_YXy = squeeze(sum(data_d9.maps_YXvyB(:,:,IA,:),3)) ;
    endh_YXmean = mean(endh_YXy,3) ;
else
    endh_xy = squeeze(sum(data_d9.garr_xvyB(:,IA,:),2)) ;
    endh_xmean = mean(endh_xy,2) ;
    endh_YXmean = lpjgu_vector2map(endh_xmean, map_size, list2map) ;
end

clim_max = 0 ;
hs = [] ;
hcbs = [] ;
if already_maps
    maps_YXr = nan([size(land_area_YX) Nruns]) ;
else
    maps_YXr = nan([map_size Nruns]) ;
end
for r = 1:Nruns
    hs(r) = subplot_tight(Nruns,1,r,spacing) ;
    
    % Get data
    if already_maps
        endf_YXy = squeeze(sum(data_d9.maps_YXvyr(:,:,IA,:,r),3)) ;
        endf_YXmean = mean(endf_YXy,3) ;
    else
        endf_xy = squeeze(sum(data_d9.garr_xvyr(:,IA,:,r),2)) ;
        endf_xmean = mean(endf_xy,2) ;
        endf_YXmean = lpjgu_vector2map(endf_xmean, map_size, list2map) ;
    end
    map = (endf_YXmean - endh_YXmean) .* conv_fact_map ;
    maps_YXr(:,:,r) = map ;
    
    % Make plot
    pcolor(map(69:end,:)) ; shading flat ; axis equal tight off
    if ~isempty(prctile_clim)
        this_clim_max = prctile(abs(map(~isnan(map))), prctile_clim) ;
    else
        this_clim_max = max(abs(caxis)) ;
    end
    caxis(this_clim_max*[-1 1])
    clim_max = max(clim_max,this_clim_max) ;
    colormap(gca,brighten(brewermap(64,'RdBu_ssr'),-0.3)) ;
        
    % Add labels
    if r==1
%         hsgt = sgtitle(sprintf('Diff. in %s, 2000s to 2090s', title_text)) ;
%         set(hsgt, 'FontSize', fontSize+4, 'FontWeight', 'bold')
        htmp = text(0, 0, sprintf('Diff. in %s, 2000s to 2090s', title_text), ...
            'FontSize', fontSize+4, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center') ;
        htmp.Units = 'normalized' ;
        htmp.Position = [0.5 1.25 0] ;
    end
    ht = title(runList{r}) ;
    set(gca,'FontSize',fontSize)
    letterlabel_align0(char(r + 64),ht,do_caps) ;
    
    % Add data
    if ~isempty(conv_fact_total)
        mean_endh = conv_fact_total * nansum(nansum(endh_YXmean .* land_area_YX)) ;
        mean_endf = conv_fact_total * nansum(nansum(endf_YXmean .* land_area_YX)) ;
        if already_maps
            sd_endh = std(nansum(nansum(endh_YXy .* repmat(land_area_YX,[1 1 10]),1),2),...
                0,3) * conv_fact_total ;
            sd_endf = std(nansum(nansum(endf_YXy .* repmat(land_area_YX,[1 1 10]),1),2),...
                0,3) * conv_fact_total ;
        else
            sd_endh = std(nansum(endh_xy .* repmat(land_area_YX(list2map),[1 10]),1),...
                0,2) * conv_fact_total ;
            sd_endf = std(nansum(endf_xy .* repmat(land_area_YX(list2map),[1 10]),1),...
                0,2) * conv_fact_total ;
        end
        data_fontSize = fontSize ;
        if contains(title_text, 'albedo')
            text(textX,textY_1, ...
                sprintf('2000s: %0.3f�%0.3f %s', mean_endh, sd_endh, units_total), ...
                'FontSize',data_fontSize) ;
            text(textX,textY_2, ...
                sprintf('2090s: %0.3f�%0.3f %s', mean_endf, sd_endf, units_total), ...
                'FontSize',data_fontSize) ;
        else
            text(textX,textY_1, ...
                sprintf('2000s: %d�%d %s', round(mean_endh), round(sd_endh), units_total), ...
                'FontSize',data_fontSize) ;
            text(textX,textY_2, ...
                sprintf('2090s: %d�%d %s', round(mean_endf), round(sd_endf), units_total), ...
                'FontSize',data_fontSize) ;
        end
    end
end

    
% Equalize color axes
for r = 1:Nruns
    caxis(hs(r),[-clim_max clim_max]) ;
end

% Move subplots
for r = 1:Nruns
    shiftit = -0.01 + (r-1)*0.0125 ;
    set(hs(r),'Position', get(hs(r),'Position') + [0 shiftit 0 0])
end

% Add big colorbar at bottom
cb_pos = get(hs(end),'Position') ;
cb_pos(2) = 0.05 ; % Y position in figure
cb_pos(4) = 0.015 ; % Height
hcb = colorbar('Location','SouthOutside','Position',cb_pos);
set(get(hcb,'XLabel'),'String',['2090s minus 2000s (' units_map ')'])
set(hcb, 'FontSize', fontSize)
if ~isempty(prctile_clim)
    hcb_ticks = get(hcb,'Ticks') ;
    hcb_limits = get(hcb,'Limits') ;
    hcb_ticklabels = get(hcb,'TickLabels') ;
    new_ticks = [hcb_limits(1) hcb_ticks hcb_limits(2)] ;
    new_ticklabels = [{'<'};hcb_ticklabels;{'>'}] ;
    set(hcb, ...
        'Ticks', new_ticks, ...
        'TickLabels',  new_ticklabels)
end





end