function maps_YXr = map_run_diffs_fromEndHist_oneCol_v2( ...
    data_d9, title_text, sumvars, ...
    fontSize, spacing, textX, textY_1, textY_2, ...
    thisPos, runList, do_caps, land_area_YX, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    bins_lowBnds, this_colormap_name, ...
    map_size, list2map, precision_total)

if ~isint(precision_total)
    error('precision_total must be int')
end

% Get missing info
Nruns = length(runList) ;
if isempty(land_area_YX)
    land_area_YX = ones(map_size) ;
end

figure('Color','w','Position',thisPos) ;

[~,IA] = intersect(data_d9.varNames,sumvars) ;
if isempty(IA)
    error('sumvars not found in data_d9.varNames!')
end
endh_xy = squeeze(sum(data_d9.garr_xvyB(:,IA,:),2)) ;
endh_xmean = mean(endh_xy,2) ;
endh_YXmean = lpjgu_vector2map(endh_xmean, map_size, list2map) ;

clim_max = 0 ;
hs = [] ;
hcbs = [] ;
maps_YXr = nan([map_size Nruns]) ;
for r = 1:Nruns
    hs(r) = subplot_tight(Nruns,1,r,spacing) ;
    
    % Get data
    endf_xy = squeeze(sum(data_d9.garr_xvyr(:,IA,:,r),2)) ;
    endf_xmean = mean(endf_xy,2) ;
    endf_YXmean = lpjgu_vector2map(endf_xmean, map_size, list2map) ;
    map = (endf_YXmean - endh_YXmean) .* conv_fact_map ;
    maps_YXr(:,:,r) = map ;
    
    % Make plot
%     Nbins = 64 ;
%     map_YX_bin = map ;
    Nbins = length(bins_lowBnds) ;
    map_YX_bin = nan(size(map)) ;
    for b = 1:Nbins
        map_YX_bin(map >= bins_lowBnds(b)) = b+0.1 ;
    end
        
    % Plot map
    R = georefcells([-90 90],[-180 180],[360 720]) ;
    hwm = worldmap(map_YX_bin, R) ;
    setm(hwm,'MapLatLimit',[-55 90])
    setm(hwm,'mapprojection','Robinson', 'frame','off','grid','on')
    set(findall(hwm,'Tag','PLabel'),'visible','off')
    set(findall(hwm,'Tag','MLabel'),'visible','off')
    geoshow(map_YX_bin, R, ...
        'DisplayType', 'surface')
    colormap(gca,brighten(brewermap(Nbins,this_colormap_name),-0.3)) ;
    if Nbins<64
        new_caxis = [1 Nbins+1] ;
    else
        new_caxis = [-1 1]*max(abs(caxis)) ;
    end
    caxis(new_caxis)
    
    % Plot continent outlines
    hold on
    geoshow('landareas.shp', 'FaceAlpha', 0)
    hold off
    
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
        sd_endh = std(nansum(endh_xy .* repmat(land_area_YX(list2map),[1 10]),1),...
            0,2) * conv_fact_total ;
        sd_endf = std(nansum(endf_xy .* repmat(land_area_YX(list2map),[1 10]),1),...
            0,2) * conv_fact_total ;
        data_fontSize = fontSize ;
%         if contains(title_text, 'albedo')
%             text(textX,textY_1, ...
%                 sprintf('2000s: %0.3f�%0.3f %s', mean_endh, sd_endh, units_total), ...
%                 'FontSize',data_fontSize, 'Units', 'normalized') ;
%             text(textX,textY_2, ...
%                 sprintf('2090s: %0.3f�%0.3f %s', mean_endf, sd_endf, units_total), ...
%                 'FontSize',data_fontSize, 'Units', 'normalized') ;
%         else
%             text(textX,textY_1, ...
%                 sprintf('2000s: %d�%d %s', round(mean_endh), round(sd_endh), units_total), ...
%                 'FontSize',data_fontSize, 'Units', 'normalized') ;
%             text(textX,textY_2, ...
%                 sprintf('2090s: %d�%d %s', round(mean_endf), round(sd_endf), units_total), ...
%                 'FontSize',data_fontSize, 'Units', 'normalized') ;
%         end
        the_diff = mean_endf - mean_endh ;
        pct_diff = 100*(the_diff)/mean_endh ;
        if the_diff >= 0
            plussign = '+' ;
        else
            plussign = '' ;
        end
        if precision_total == 0
            the_diff = round(the_diff) ;
            data_char = '%d' ;
        elseif precision_total>0
            data_char = sprintf('%%0.%df', precision_total) ;
        else
            error('Problem with precision_total')
        end
        the_text = sprintf( ...
            ['%s' data_char ' %s\n(%s%0.1f%%)'], ...
            plussign, the_diff, units_total, plussign, pct_diff) ;
        ht = text(textX,textY_1, ...
            the_text, ...
            'FontSize',data_fontSize, 'Units', 'normalized') ;
    end
    pause(0.1)
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
if Nbins<64
    hcb.Ticks = 1:(Nbins+1) ;
    bins_lowBnds_str = strrep(cellstr(num2str(bins_lowBnds')), ' ', '') ;
    bins_lowBnds_str = [bins_lowBnds_str ; {num2str(-bins_lowBnds(1))}] ;
    hcb.TickLabels = bins_lowBnds_str ;
end
if any(any(any(maps_YXr>max(caxis))))
    hcb.TickLabels{1} = ['\leq' hcb.TickLabels{1}] ;
    hcb.TickLabels{end} = ['\geq' hcb.TickLabels{end}] ;
end




end