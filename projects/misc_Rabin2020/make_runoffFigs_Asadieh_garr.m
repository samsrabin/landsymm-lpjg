function pctDiff_xr = make_runoffFigs_Asadieh_garr( ...
    garr_mon_runoff_last30, garr_awater_last30, runList, ...
    droughtOrFlood, land_area_weights_x, ...
    do_norm, spacing, norm_ticks, fontSize, fontSize_text, Ystart, ...
    basins_x, do_aggregate_basins, map_size, list2map)

if ~strcmp(droughtOrFlood, 'drought') && ~strcmp(droughtOrFlood, 'flood')
    error('droughtOrFlood (%s) must be either ''drought'' or ''flood''', droughtOrFlood)
end

Nruns = length(runList) ;
BoverA_to_dq = @(BoverA) (BoverA-1)./(BoverA+1) ;
dq_to_BoverA = @(dq) -(dq+1)./(dq-1) ;
this_colormap = brighten(brewermap(64,'RdBu_ssr'),-0.3) ;
Ncolors = size(this_colormap,1) ;

if do_aggregate_basins
    [garr_mon_runoff_last30, garr_awater_last30] = ...
        aggregate_to_basins(...
        garr_mon_runoff_last30, garr_awater_last30, ...
        basins_x) ;
end

% After Asadieh & Krakauer (2017): Ignore cells with runoff < 0.01mm/day in
% last 10 years of baseline
ind_var = strcmp(garr_awater_last30.varNames,'Runoff') ;
ind_stat = strcmp(garr_awater_last30.statList,'last30_mean') ;
below_thresh_annMean_x = garr_awater_last30.garr_xvsB(:,ind_var,ind_stat)/365 < 0.01 ; 

if strcmp(droughtOrFlood, 'drought')
    ind_var = strcmp(garr_awater_last30.varNames,'Runoff') ;
    ind_stat = strcmp(garr_awater_last30.statList,'last30_prctile05') ;
    pXX_q20c_x = garr_awater_last30.garr_xvsB(:,ind_var,ind_stat) ;
    pXX_q21c_xr = squeeze(garr_awater_last30.garr_xvsr(:,ind_var,ind_stat,:)) ;
    below_thresh_x = below_thresh_annMean_x ;   % | pXX_q20c_x / 365 == 0 ; % Not in Asadieh & Krakauer, but otherwise you get Infs...
    pXX_q20c_x(below_thresh_x) = NaN ;
    pXX_q20c_YX = lpjgu_vector2map(pXX_q20c_x, map_size, list2map) ;
    pXX_q21c_xr(repmat(below_thresh_x, [1 Nruns])) = NaN ;
    pXX_q20c_xr = repmat(pXX_q20c_x,[1 Nruns]) ;
    pXX_dq_xr = (pXX_q21c_xr - pXX_q20c_xr) ./ (pXX_q21c_xr + pXX_q20c_xr) ;
    thisTitle = '10-year lowest annual runoff' ;
else
    ind_stat = strcmp(garr_mon_runoff_last30.statList,'last30_prctile95') ;
    pXX_q20c_x = garr_mon_runoff_last30.garr_xvsB(:,:,ind_stat) ;
    pXX_q21c_xr = squeeze(garr_mon_runoff_last30.garr_xvsr(:,:,ind_stat,:)) ;
    below_thresh_x = below_thresh_annMean_x ;     % | pXX_q20c_x * 12/365 < 0.01 ; % Not in Asadieh & Krakauer, but otherwise you get Infs...
    pXX_q20c_x(below_thresh_x) = NaN ;
    pXX_q20c_YX = lpjgu_vector2map(pXX_q20c_x, map_size, list2map) ;
    pXX_q21c_xr(repmat(below_thresh_x, [1 Nruns])) = NaN ;
    pXX_q20c_xr = repmat(pXX_q20c_x,[1 Nruns]) ;
    pXX_dq_xr = (pXX_q21c_xr - pXX_q20c_xr) ./ (pXX_q21c_xr + pXX_q20c_xr) ;
    thisTitle = '10-year highest monthly runoff' ;
end

% Asadieh & Krakauer (2017) do not re-weight (their Table 1 sums to 75.9%
% of land area)
land_area_weights_x(below_thresh_x) = NaN ;
land_area_weights_YX = lpjgu_vector2map(land_area_weights_x, map_size, list2map) ;
fprintf('%s: %0.1f%% of land area included\n', ...
    droughtOrFlood, 100*nansum(nansum(land_area_weights_YX)))
% land_area_weights_YX = land_area_weights_YX ./ nansum(nansum(land_area_weights_YX)) ; % Re-weight

pctDiff_xr = nan(size(pXX_q21c_xr)) ;
figure('Color','w','Position',figurePos) ;
for r = 1:Nruns
    subplot_tight(2,2,r,spacing)
    
    % Get percent difference
    pctDiff_x = 100*(pXX_q21c_xr(:,r) - pXX_q20c_x) ./ pXX_q20c_x ;
    pctDiff_x(pXX_q20c_x==0 & pXX_q21c_xr(:,r)==0) = 0 ;
    pctDiff_xr(:,r) = pctDiff_x ;
    
    % Make map
    if do_norm
        pXX_dq_YX = lpjgu_vector2map(pXX_dq_xr(:,r), map_size, list2map) ;
        pcolor(pXX_dq_YX(Ystart:end,:)) ;
        shading flat
        axis equal tight off
        caxis([-1 1]) ;
        colormap(gca,this_colormap)
        set(gca,'FontSize',fontSize)
        hcb = colorbar ;
        if ~isempty(norm_ticks)
            cb_ticks_new = sort([-BoverA_to_dq(norm_ticks(~isinf(norm_ticks))) 0 BoverA_to_dq(norm_ticks(~isinf(norm_ticks)))]) ;
            cb_ticklabels_new_num = [flip(norm_ticks,2) 1 norm_ticks] ;
            cb_ticklabels_new = compose('%d', cb_ticklabels_new_num) ;
            yep = ~isint(cb_ticklabels_new_num) & ~isinf(cb_ticklabels_new_num) ;
            cb_ticklabels_new(yep) = compose('%0.1f', cb_ticklabels_new_num(yep)) ;
            if any(isinf(norm_ticks))
                cb_ticks_new = [-1 cb_ticks_new 1] ;
            end
            cb_ticklabels_new = strcat('\times', cb_ticklabels_new) ;
            cb_ticklabels_new(cb_ticks_new<0) = strcat(cb_ticklabels_new(cb_ticks_new<0), '^{-1}') ;
            hcb.Ticks = cb_ticks_new ;
            clear cb_ticks_new yep
        else
            cb_ticks = hcb.Ticks ;
            cb_ticks_new = (1 + abs(cb_ticks)) ./ (1 - abs(cb_ticks)) ;
            cb_ticklabels_new = compose('%0.1f', cb_ticks_new) ;
            cb_ticklabels_new = strcat('\times', cb_ticklabels_new) ;
            cb_ticklabels_new(cb_ticks<0) = strcat(cb_ticklabels_new(cb_ticks<0),'^{-1}') ;
        end
        hcb.TickLabels = cb_ticklabels_new ;
        clear cb_ticklabels_new
    else
        pctDiff_YX = lpjgu_vector2map(pctDiff_x, map_size, list2map) ;
        pcolor(pctDiff_YX(Ystart:end,:)) ;
        shading flat
        axis equal tight off
        colormap(gca,this_colormap)
        set(gca,'FontSize',fontSize)
        hcb = colorbar ;
        caxis([-100 100]) ;
        cb_ticklabels = hcb.TickLabels ;
        cb_ticklabels = strcat(cb_ticklabels, '%') ;
        cb_ticklabels(hcb.Ticks>0) = strcat('+', cb_ticklabels(hcb.Ticks>0)) ;
        cb_ticklabels{end} = [char(8805) cb_ticklabels{end}] ;
        hcb.TickLabels = cb_ticklabels ;
    end
    title(runList{r})
    
    
    % Add stats
    pct_inc = 100*nansum(land_area_weights_x(pctDiff_x>0)) ;
    pct_dec = 100*nansum(land_area_weights_x(pctDiff_x<0)) ;
    median_inc = median(pctDiff_x(pctDiff_x>0)) ;
    median_dec = median(pctDiff_x(pctDiff_x<0)) ;
    pXX_dq_x = pXX_dq_xr(:,r) ;
    mean_inc = 100*(dq_to_BoverA(mean(pXX_dq_x(pctDiff_x>0)))-1) ;
    mean_dec = 100*(dq_to_BoverA(mean(pXX_dq_x(pctDiff_x<0)))-1) ;
    
    xlims = get(gca,'XLim') ;
    ylims = get(gca,'YLim') ;
    xrange = xlims(2) - xlims(1) ;
    yrange = ylims(2) - ylims(1) ;
    x = 0.05*xrange ;
    text(x, 0.35*yrange, ...
        sprintf('%0.0f%% increasing\n(median +%0.0f%%,\nmean +%0.0f%%)', pct_inc, median_inc, mean_inc), ...
        'FontSize', fontSize_text, 'Color', this_colormap(round(0.85*Ncolors),:)) ;
    text(x, 0.1*yrange, ...
        sprintf('%0.0f%% decreasing\n(median %0.0f%%,\nmean %0.0f%%)', pct_dec, median_dec, mean_dec), ...
        'FontSize', fontSize_text, 'Color', this_colormap(round(0.15*Ncolors),:)) ;
end
sgtitle(thisTitle,'FontSize',fontSize*1.5,'FontWeight','bold')


end


function [garr_mon_runoff_last30, garr_awater_last30] = ...
    aggregate_to_basins(...
    garr_mon_runoff_last30, garr_awater_last30, ...
    basins_x)

basin_list = unique(basins_x) ;
Nbasins = length(basin_list) ;

fprintf('Aggregating basins... ')
pct_done = 0 ;
dpf = 0.1 ;

for b = 1:Nbasins
    
%     fprintf('Basin %d of %d\n', b, Nbasins)
    if b>1 && rem(b/Nbasins,dpf) < rem((b-1)/Nbasins,dpf)
        pct_done = pct_done + 100*dpf ;
        fprintf('%d%%...', pct_done)
    end
    
    isThisBasin = basins_x==basin_list(b) ;
    
    % garr_mon_runoff_last30.garr_xvsB
    size_vsB = size(garr_mon_runoff_last30.garr_xvsB) ;
    size_vsB = size_vsB(2:3) ;
    for v = 1:size_vsB(1)
        for s = 1:size_vsB(2)
            garr_mon_runoff_last30.garr_xvsB(:,v,s) = ...
                    average_this_basin(isThisBasin, ...
                    garr_mon_runoff_last30.garr_xvsB(:,v,s)) ;
        end ; clear s
    end ; clear v
    
    % garr_mon_runoff_last30.garr_xvsr
    size_vsr = size(garr_mon_runoff_last30.garr_xvsr) ;
    size_vsr = size_vsr(2:4) ;
    for v = 1:size_vsr(1)
        for s = 1:size_vsr(2)
            for r = 1:size_vsr(3)
                garr_mon_runoff_last30.garr_xvsr(:,v,s,r) = ...
                    average_this_basin(isThisBasin, ...
                    garr_mon_runoff_last30.garr_xvsr(:,v,s,r)) ;
            end ; clear r
        end ; clear s
    end ; clear v
    
    % garr_awater_last30.garr_xvsB
    size_vsB = size(garr_awater_last30.garr_xvsB) ;
    size_vsB = size_vsB(2:3) ;
    for v = 1:size_vsB(1)
        for s = 1:size_vsB(2)
            garr_awater_last30.garr_xvsB(:,v,s) = ...
                    average_this_basin(isThisBasin, ...
                    garr_awater_last30.garr_xvsB(:,v,s)) ;
        end ; clear s
    end ; clear v
    
    % garr_awater_last30.garr_xvsr
    size_vsr = size(garr_awater_last30.garr_xvsr) ;
    size_vsr = size_vsr(2:4) ;
    for v = 1:size_vsr(1)
        for s = 1:size_vsr(2)
            for r = 1:size_vsr(3)
                garr_awater_last30.garr_xvsr(:,v,s,r) = ...
                    average_this_basin(isThisBasin, ...
                    garr_awater_last30.garr_xvsr(:,v,s,r)) ;
            end ; clear r
        end ; clear s
    end ; clear v
    
end
fprintf('Done.\n')

end


function tmp_YX = average_this_basin(isThisBasin, tmp_YX)

tmp = tmp_YX(isThisBasin) ;
NinThisBasin = length(tmp(~isnan(tmp))) ;
if NinThisBasin == 0
    return
end
tmp_YX(isThisBasin) = nansum(tmp) / NinThisBasin ;


end

