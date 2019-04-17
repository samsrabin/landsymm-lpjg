function make_runoffFigs_Asadieh( ...
    maps_mon_runoff_d9, maps_awater_d9, runList, ...
    droughtOrFlood, ...
    do_norm, spacing, norm_ticks, fontSize, fontSize_text, Ystart)

if ~strcmp(droughtOrFlood, 'drought') && ~strcmp(droughtOrFlood, 'flood')
    error('droughtOrFlood (%s) must be either ''drought'' or ''flood''', droughtOrFlood)
end

size_YX = [size(maps_mon_runoff_d9.maps_YXvyB,1) size(maps_mon_runoff_d9.maps_YXvyB,2)] ;
Nruns = length(runList) ;
BoverA_to_dq = @(x) (x-1)./(x+1) ;
this_colormap = brighten(brewermap(64,'RdBu_ssr'),-0.3) ;
Ncolors = size(this_colormap,1) ;

% After Asadieh & Krakauer (2017): Ignore cells with runoff < 0.01mm/day in
% last 10 years of baseline
ii = strcmp(maps_awater_d9.varNames,'Runoff') ;
below_thresh_annMean_YX = mean(maps_awater_d9.maps_YXvyB(:,:,ii,:),4)/365 < 0.01 ; 

if strcmp(droughtOrFlood, 'drought')
    pXX_q20c_YX = prctile(sum(maps_mon_runoff_d9.maps_YXvyB,3),5,4) ;
    pXX_q21c_YXr = squeeze(prctile(sum(maps_mon_runoff_d9.maps_YXvyr,3),5,4)) ;
    below_thresh_YX = below_thresh_annMean_YX ;   % | pXX_q20c_YX / 365 == 0 ; % Not in Asadieh & Krakauer, but otherwise you get Infs...
    pXX_q20c_YX(below_thresh_YX) = NaN ;
    pXX_q21c_YXr(repmat(below_thresh_YX, [1 1 Nruns])) = NaN ;
    pXX_q20c_YXr = repmat(pXX_q20c_YX,[1 1 Nruns]) ;
    pXX_dq_YXr = (pXX_q21c_YXr - pXX_q20c_YXr) ./ (pXX_q21c_YXr + pXX_q20c_YXr) ;
    thisTitle = '10-year lowest annual runoff' ;
else
    pXX_q20c_YX = prctile(reshape(maps_mon_runoff_d9.maps_YXvyB, [size_YX 120]),95,3) ;
    pXX_q21c_YXr = squeeze(prctile(reshape(maps_mon_runoff_d9.maps_YXvyr,[size_YX 120 Nruns]),95,3)) ;
    below_thresh_YX = below_thresh_annMean_YX ;     % | pXX_q20c_YX * 12/365 < 0.01 ; % Not in Asadieh & Krakauer, but otherwise you get Infs...
    pXX_q20c_YX(below_thresh_YX) = NaN ;
    pXX_q21c_YXr(repmat(below_thresh_YX, [1 1 Nruns])) = NaN ;
    pXX_q20c_YXr = repmat(pXX_q20c_YX,[1 1 Nruns]) ;
    pXX_dq_YXr = (pXX_q21c_YXr - pXX_q20c_YXr) ./ (pXX_q21c_YXr + pXX_q20c_YXr) ;
    thisTitle = '10-year highest monthly runoff' ;
end

figure('Color','w','Position',figurePos) ;
for r = 1:Nruns
    subplot_tight(2,2,r,spacing)
    
    % Get percent difference
    pctDiff_YX = 100*(pXX_q21c_YXr(:,:,r) - pXX_q20c_YX) ./ pXX_q20c_YX ;
    pctDiff_YX(pXX_q20c_YX==0 & pXX_q21c_YXr(:,:,r)==0) = 0 ;
    
    % Make map
    if do_norm
        pcolor(pXX_dq_YXr(Ystart:end,:,r)) ;
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
    pct_inc = 100*length(find(pctDiff_YX>0)) ./ length(find(~isnan(pctDiff_YX))) ;
    median_inc = median(pctDiff_YX(pctDiff_YX>0)) ;
    pct_dec = 100*length(find(pctDiff_YX<0)) ./ length(find(~isnan(pctDiff_YX))) ;
    median_dec = median(pctDiff_YX(pctDiff_YX<0)) ;
    xlims = get(gca,'XLim') ;
    ylims = get(gca,'YLim') ;
    xrange = xlims(2) - xlims(1) ;
    yrange = ylims(2) - ylims(1) ;
    x = 0.05*xrange ;
    text(x, 0.3*yrange, ...
        sprintf('%0.0f%% increasing\n(median +%0.0f%%)', pct_inc, median_inc), ...
        'FontSize', fontSize_text, 'Color', this_colormap(round(0.85*Ncolors),:)) ;
    text(x, 0.1*yrange, ...
        sprintf('%0.0f%% decreasing\n(median %0.0f%%)', pct_dec, median_dec), ...
        'FontSize', fontSize_text, 'Color', this_colormap(round(0.15*Ncolors),:)) ;
end
sgtitle(thisTitle,'FontSize',fontSize*1.5,'FontWeight','bold')




end