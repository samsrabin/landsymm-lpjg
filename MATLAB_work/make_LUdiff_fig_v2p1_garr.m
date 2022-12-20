function [h, caxis_max] = make_LUdiff_fig_v2p1_garr(...
    maps_xB, maps_xr, ...
    y1, yN, LUname, runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, i1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, ...
    do_caps, thisStat, as_pct_change, map_size, list2map)

total_bl = nansum(nansum(maps_xB)) ;

% Get maps
maps_xBr = repmat(maps_xB,[1 Nruns]) ;
this_diff_xr = maps_xr - maps_xBr ;
if as_pct_change
    this_diff_xr = 100*(this_diff_xr ./ maps_xBr) ;
    this_diff_xr(maps_xBr==0) = NaN ;
else
    this_diff_xr = this_diff_xr*conv_fact_map ;
end
maps_xB = maps_xB*conv_fact_map ;
clear maps_xBr

% Set outliers to median + 2*IQR
iqr_mult = 2 ;
for r = 1:Nruns
    tmp_x = this_diff_xr(:,r) ;
    diff_prctile_25 = prctile(tmp_x, 25) ;
    diff_prctile_50 = prctile(tmp_x, 50) ;
    diff_prctile_75 = prctile(tmp_x, 75) ;
    thisWidth = iqr_mult*(diff_prctile_75 - diff_prctile_25) ;
    tmp_x(tmp_x>diff_prctile_50+thisWidth) = diff_prctile_50+thisWidth ;
    tmp_x(tmp_x<diff_prctile_50-thisWidth) = diff_prctile_50-thisWidth ;
    this_diff_xr(:,r) = tmp_x ;
end

figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if only1bl && r>1
        maps_xB = [] ;
    end
    runName2 = runList{r} ;
    if only1bl
        i2 = ssp_plot_index(r) ;
        MAPS_XB = maps_xB ;
        runName1 = '' ;
    else
        i1 = (r-1)*2+1  ;
        i2 = (r-1)*2+2  ;
        MAPS_XB = maps_xB(:,r) ;
        total_bl = nansum(MAPS_XB/conv_fact_map) ;
        runName1 = runName2 ;
    end
    if isempty(MAPS_XB)
        MAPS_YXB = [] ;
    else
        MAPS_YXB = lpjgu_vector2map(MAPS_XB, map_size, list2map) ;
    end
    this_diff_YX = lpjgu_vector2map(this_diff_xr(:,r), map_size, list2map) ;
    [h, caxis_max] = actually_make_fig(...
        MAPS_YXB, this_diff_YX, total_bl, ...
        y1, yN, LUname, runName1, runName2, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, i1, i2, colorBarLoc, ...
        conv_fact_map, conv_fact_total, units_map, units_total, ...
        do_caps, only1bl, thisStat, as_pct_change) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
%     caxis(gcas(r),[-colorlim colorlim])
    if as_pct_change && colorlim>100
%         caxis(gcas(r),[-colorlim colorlim])
%         caxis(gcas(r),100*[-1 1])
        percentage_colormap(gcas, 'rdbu_ssr', colorlim)
    else
        caxis(gcas(r),[-colorlim colorlim])
    end
end

end


function [h, caxis_max, ht] = actually_make_fig(...
    maps_YXB, this_diff_YX, total_bl, ...
    y1, yN, LUname, runName1, runName2, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, i1, i2, colorBarLoc, ...
    conv_fact_map, conv_fact_total, units_map, units_total, ...
    do_caps, only1bl, thisStat, as_pct_change)

% Baseline
if ~isempty(maps_YXB)
    subplot_tight(ny,nx,i1,spacing) ;
    pcolor(maps_YXB(69:end,:)) ; shading flat ; axis equal tight off
    colorbar(colorBarLoc)
    set(gca,'XTick',[],'YTick',[])
    set(gca,'FontSize',fontSize)
    if only1bl
        thisTitle = sprintf('%s %s, %d (%s)', LUname, thisStat, y1, units_map) ;
    else
        thisTitle = sprintf('%s %s, %d: %s (%s)', LUname, thisStat, y1, runName1, units_map) ;
    end
    ht = title(thisTitle) ;
    if isempty(total_bl)
        total_bl = nansum(nansum(maps_YXB)) ;
    end
    text(textX,textY_1,[num2str(round(total_bl*conv_fact_total,1)) ' ' units_total],'FontSize',fontSize-2) ;
    letterlabel_align0(char(i1 + 64),ht,do_caps) ;
end

% Future
subplot_tight(ny,nx,i2,spacing) ;
pcolor(this_diff_YX(69:end,:)) ; shading flat ; axis equal tight off
colorbar(colorBarLoc)
colormap(gca,brighten(brewermap(64,'RdBu_ssr'),-0.3))
set(gca,'XTick',[],'YTick',[])
set(gca,'FontSize',fontSize)
if as_pct_change
    ht = title(sprintf('%s %s %s, %d: %s (%%)', '\Delta', LUname, thisStat, yN, runName2)) ;
else
    ht = title(sprintf('%s %s %s, %d: %s (%s)', '\Delta', LUname, thisStat, yN, runName2, units_map)) ;
end

% Total
if ~isempty(conv_fact_total)
    total_yr = nansum(nansum(this_diff_YX/conv_fact_map)) ;
    pctDiff = round(100*total_yr./total_bl,1) ;
    pctDiff_str = num2str(pctDiff) ;
    if pctDiff>0
        pctDiff_str = ['+' pctDiff_str] ;
    end
    text(textX,textY_1,[num2str(round(total_yr*conv_fact_total,1)) ' ' units_total],'FontSize',fontSize-2) ;
    text(textX,textY_2,['(' pctDiff_str ' %)'],'FontSize',fontSize-2) ;
    if i2>3 && only1bl
        i2 = i2 - 1 ;
    end
    letterlabel_align0(char(i2 + 64),ht,do_caps) ;
end

% Set up for changing color axis limits
h = gca ;
caxis_max = max(abs(caxis)) ;

end