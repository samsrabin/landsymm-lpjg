function [h, caxis_max] = make_LUdiff_fig(...
    this_area_YXB, this_diff_YX, total_bl, ...
    y1, yN, LUname, runName1, runName2, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, i1, i2, colorBarLoc)

% Baseline
if ~isempty(this_area_YXB)
    subplot_tight(ny,nx,i1,spacing) ;
    pcolor(this_area_YXB(69:end,:)) ; shading flat ; axis equal tight off
    colorbar(colorBarLoc)
    set(gca,'XTick',[],'YTick',[])
    set(gca,'FontSize',fontSize)
    title([LUname ' area, ' num2str(y1) ': ' runName1 ' (km^2)'])
    if isempty(total_bl)
        total_bl = 1e-6*nansum(nansum(this_area_YXB)) ;
    end
    text(textX,textY_1,[num2str(round(total_bl,1)) ' Mkm^2'],'FontSize',fontSize-2) ;
end

% Future
subplot_tight(ny,nx,i2,spacing) ;
pcolor(this_diff_YX(69:end,:)) ; shading flat ; axis equal tight off
colorbar(colorBarLoc)
colormap(gca,brighten(brewermap(64,'RdBu_ssr'),-0.3))
set(gca,'XTick',[],'YTick',[])
set(gca,'FontSize',fontSize)
title(['\Delta ' LUname ' area, ' num2str(yN) ': ' runName2 ' (km^2)'])
total_yr = 1e-6*nansum(nansum(this_diff_YX)) ;
pctDiff = round(100*total_yr./total_bl,1) ;
pctDiff_str = num2str(pctDiff) ;
if pctDiff>0
    pctDiff_str = ['+' pctDiff_str] ;
end
text(textX,textY_1,[num2str(round(total_yr,1)) ' Mkm^2'],'FontSize',fontSize-2) ;
text(textX,textY_2,['(' pctDiff_str ' %)'],'FontSize',fontSize-2) ;

% Set up for changing color axis limits
h = gca ;
caxis_max = max(abs(caxis)) ;

end