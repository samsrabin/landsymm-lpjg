function [h, caxis_max] = make_LUdiff_fig_v2(...
    this_area_YXB, this_diff_YXr, ...
    y1, yN, LUname, runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, i1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, ...
    do_caps)

total_bl = nansum(nansum(this_area_YXB)) ;

this_area_YXB = this_area_YXB*conv_fact_map ;
this_diff_YXr = this_diff_YXr*conv_fact_map ;

figure('Color','w','Position',thisPos) ;
gcas = {} ;
colorlim = 0 ;
for r = 1:Nruns
    if only1bl && r>1
        this_area_YXB = [] ;
    end
    runName2 = runList{r} ;
    if only1bl
        i2 = ssp_plot_index(r) ;
        THIS_area_YXB = this_area_YXB ;
        runName1 = '' ;
    else
        i1 = (r-1)*2+1  ;
        i2 = (r-1)*2+2  ;
        THIS_area_YXB = this_area_YXB(:,:,r) ;
        total_bl = nansum(nansum(THIS_area_YXB/conv_fact_map)) ;
        runName1 = runName2 ;
    end
    [h, caxis_max] = actually_make_fig(...
        THIS_area_YXB, this_diff_YXr(:,:,r), total_bl, ...
        y1, yN, LUname, runName1, runName2, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, i1, i2, colorBarLoc, ...
        conv_fact_map, conv_fact_total, units_map, units_total, ...
        do_caps, only1bl) ;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    gcas = [gcas h] ;
    colorlim = max(colorlim,caxis_max) ;
end
for r = 1:Nruns
    caxis(gcas(r),[-colorlim colorlim])
end

end


function [h, caxis_max, ht] = actually_make_fig(...
    this_area_YXB, this_diff_YX, total_bl, ...
    y1, yN, LUname, runName1, runName2, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, i1, i2, colorBarLoc, ...
    conv_fact_map, conv_fact_total, units_map, units_total, ...
    do_caps, only1bl)

% Baseline
if ~isempty(this_area_YXB)
    subplot_tight(ny,nx,i1,spacing) ;
    pcolor(this_area_YXB(69:end,:)) ; shading flat ; axis equal tight off
    colorbar(colorBarLoc)
    set(gca,'XTick',[],'YTick',[])
    set(gca,'FontSize',fontSize)
    if only1bl
        ht = title([LUname ' area, ' num2str(y1) ' (' units_map ')']) ;
    else
        ht = title([LUname ' area, ' num2str(y1) ': ' runName1 ' (' units_map ')']) ;
    end
    if isempty(total_bl)
        total_bl = nansum(nansum(this_area_YXB)) ;
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
ht = title(['\Delta ' LUname ' area, ' num2str(yN) ': ' runName2 ' (' units_map ')']) ;
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

% Set up for changing color axis limits
h = gca ;
caxis_max = max(abs(caxis)) ;

end