function [diff_crop_YXr, diff_past_YXr] = make_LUdiff_fig_v4(...
    crop_area_YXB, past_area_YXB, diff_crop_YXr, diff_past_YXr, ...
    y1, yN, LUname, runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, i1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, ...
    do_caps, same_caxis, Nbins)

% this_colormap_name = 'RdBu_ssr' ;
% this_colormap_name = '-PRGn_ssr' ;
this_colormap_name = '-PiYG_ssr' ;

total_crop_bl = nansum(nansum(crop_area_YXB)) ;
total_past_bl = nansum(nansum(past_area_YXB)) ;

crop_area_YXB = crop_area_YXB*conv_fact_map ;
past_area_YXB = past_area_YXB*conv_fact_map ;
diff_crop_YXr = diff_crop_YXr*conv_fact_map ;
diff_past_YXr = diff_past_YXr*conv_fact_map ;

figure('Color','w','Position',thisPos) ;
gcas_crop = {} ;
gcas_past = {} ;
colorlim = 0 ;
for r = 1:Nruns
%     if only1bl && r>1
%         crop_area_YXB = [] ;
%         past_area_YXB = [] ;
%     end
    runName2 = runList{r} ;
    i1 = (r-1)*2 + 1 ;
    i2 = i1 + 1 ;
    if only1bl
%         i2 = ssp_plot_index(r) ;
        THIScrop_area_YXB = crop_area_YXB ;
        THISpast_area_YXB = past_area_YXB ;
        runName1 = '' ;
    else
%         i1 = (r-1)*2+1  ;
%         i2 = (r-1)*2+2  ;
        THIScrop_area_YXB = crop_area_YXB(:,:,r) ;
        THISpast_area_YXB = past_area_YXB(:,:,r) ;
        total_crop_bl = nansum(nansum(THIScrop_area_YXB/conv_fact_map)) ;
        total_past_bl = nansum(nansum(THISpast_area_YXB/conv_fact_map)) ;
        runName1 = runName2 ;
    end
    [hc, hp, caxis_max_crop, caxis_max_past] = actually_make_fig(...
        THIScrop_area_YXB, THISpast_area_YXB, diff_crop_YXr(:,:,r), diff_past_YXr(:,:,r), ...
        total_crop_bl, total_past_bl, ...
        y1, yN, LUname, runName1, runName2, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, i1, i2, colorBarLoc, ...
        conv_fact_map, conv_fact_total, units_map, units_total, ...
        do_caps, only1bl, this_colormap_name, Nbins) ;
    caxis(hc,[-max(abs(caxis(hc))) max(abs(caxis(hc)))])
    caxis(hp,[-max(abs(caxis(hp))) max(abs(caxis(hp)))])
    gcas_crop = [gcas_crop hc] ;
    gcas_past = [gcas_past hp] ;
    colorlim_crop = max(colorlim,caxis_max_crop) ;
    colorlim_past = max(colorlim,caxis_max_past) ;
end
if same_caxis
    colorlim_crop = max(colorlim_crop, colorlim_past) ;
    manual_caxis = colorlim_crop>2900 && colorlim_crop<3000 ;
    if manual_caxis
        colorlim_crop = 2750 ;
    end
    colorlim_past = colorlim_crop ;
end
for r = 1:Nruns
    hc = gcas_crop(r) ;
    hp = gcas_past(r) ;
    caxis(hc,[-colorlim_crop colorlim_crop])
    caxis(hp,[-colorlim_past colorlim_past])
    Yshift = 0.02*(r-3) ;
    hc.Position(2) = hc.Position(2) + Yshift ;
    hp.Position(2) = hp.Position(2) + Yshift ;
end

hcb = add_big_colorbar(gcas_crop(Nruns), fontSize, units_map, manual_caxis) ;
if same_caxis
    hcb.Position(1) = 0.5 - hcb.Position(3)/2 ;
else
    add_big_colorbar(gcas_past(Nruns), fontSize, units_map, false) ;
end

end


function cmap = get_colormap(colormap_name, Nbins)

flip_colormap = strcmp(colormap_name(1),'-') ;
if flip_colormap
    colormap_name = colormap_name(2:end) ;
end
cmap = brighten(brewermap(Nbins,colormap_name),-0.3) ;
if flip_colormap
    cmap = flipud(cmap) ;
end

end


function [hc, hp, caxis_max_crop, caxis_max_past, ht] = actually_make_fig(...
    crop_area_YXB, past_area_YXB, diff_crop_YX, diff_past_YX, ...
    total_crop_bl, total_past_bl, ...
    y1, yN, LUname, runName1, runName2, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, i1, i2, colorBarLoc, ...
    conv_fact_map, conv_fact_total, units_map, units_total, ...
    do_caps, only1bl, this_colormap_name, Nbins)

en_dash = char(8211) ;
% flip_colormap = false ;
% if strcmp(this_colormap_name(1),'-')
%     flip_colormap = true ;
%     this_colormap_name = this_colormap_name(2:end) ;
% end
% this_colormap = brighten(brewermap(Nbins,this_colormap_name),-0.3) ;
% if flip_colormap
%     this_colormap = flipud(this_colormap) ;
% end
this_colormap = get_colormap(this_colormap_name, Nbins) ;
if ny ~= 4
    error('Code for one big colorbar (per column) only tested with ny==4')
end

% Set up for single colorbar and overarching title for each column
shiftup = 0.01*(i2/2 - 1) ;
if i1 > 2
    shiftup = shiftup * 2 ;
end
colTitle_yPos = 1.07 ;
colTitle_fontSize = fontSize + 4 ;

% Cropland
if ~isempty(crop_area_YXB) && isempty(total_crop_bl)
    total_crop_bl = nansum(nansum(crop_area_YXB)) ;
end
subplot_tight(ny,nx,i1,spacing) ;
if i1 > 1
    set(gca,'Position', get(gca,'Position') + [0 shiftup 0 0])
end
pcolor(diff_crop_YX(69:end,:)) ; shading flat ; axis equal tight off
% colorbar(colorBarLoc)
colormap(gca,this_colormap)
set(gca,'XTick',[],'YTick',[])
set(gca,'FontSize',fontSize)
ht = text(-0.01, 0.4, ...
    runName2, ...
    'Units', 'normalized', ...
    'Rotation', 90, ...
    'HorizontalAlignment', 'center', ...
    'FontSize', fontSize+4, 'FontWeight', 'bold') ;
add_totals_v3( ...
    diff_crop_YX, total_crop_bl, ...
    conv_fact_map, conv_fact_total, units_total, ...
    textX, textY_1, textY_2, fontSize, ...
    i1, do_caps)

% Set up for changing color axis limits
hc = gca ;
caxis_max_crop = max(abs(caxis)) ;

if i1==1
    htmp = text(0, 0, sprintf('%s cropland area, %d%s%d','\Delta',y1,en_dash,yN), ...
        'FontSize', colTitle_fontSize, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center') ;
    htmp.Units = 'normalized' ;
    htmp.Position = [0.5 colTitle_yPos 0] ;
end

% Pasture
if ~isempty(past_area_YXB) && isempty(total_past_bl)
    total_past_bl = nansum(nansum(past_area_YXB)) ;
end
subplot_tight(ny,nx,i2,spacing) ;
if i1 > 1
    set(gca,'Position', get(gca,'Position') + [0 shiftup 0 0])
end
pcolor(diff_past_YX(69:end,:)) ; shading flat ; axis equal tight off
colormap(gca,this_colormap)
set(gca,'XTick',[],'YTick',[])
set(gca,'FontSize',fontSize)
% ht = title(runName2) ;
add_totals_v3( ...
    diff_past_YX, total_past_bl, ...
    conv_fact_map, conv_fact_total, units_total, ...
    textX, textY_1, textY_2, fontSize, ...
    i1+1, do_caps)

% Set up for changing color axis limits
hp = gca ;
caxis_max_past = max(abs(caxis)) ;

if i1==1
    htmp = text(0, 0, sprintf('%s pasture area, %d%s%d','\Delta',y1,en_dash,yN), ...
        'FontSize', colTitle_fontSize, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center') ;
    htmp.Units = 'normalized' ;
    htmp.Position = [0.5 colTitle_yPos 0] ;
end

end


function add_totals_v3( ...
    diff_this_YX, total_this_bl, ...
    conv_fact_map, conv_fact_total, units_total, ...
    textX, textY_1, textY_2, fontSize, ...
    i1, do_caps)

% this_fontSize = fontSize - 2 ;
this_fontSize = fontSize ;

total_yr = nansum(nansum(diff_this_YX/conv_fact_map)) ;
pctDiff = round(100*total_yr./total_this_bl,1) ;
pctDiff_str = num2str(pctDiff) ;
diff_format = '%0.1f %s' ;
if pctDiff>0
    pctDiff_str = ['+' pctDiff_str] ;
    diff_format = ['+' diff_format] ;
end
text(textX,textY_1,sprintf(diff_format, ...
    total_yr*conv_fact_total,units_total),...
    'FontSize',this_fontSize) ;
text(textX,textY_2,['(' pctDiff_str '%)'],'FontSize',this_fontSize) ;
textY_3 = max(textY_1, textY_2) + abs(textY_1 - textY_2) ;

thisChar = char(i1 + 64) ;
if do_caps==1
    thisChar = upper(thisChar) ;
elseif do_caps==-1
    thisChar = lower(thisChar) ;
end
thisLabel = sprintf('(%s)', thisChar) ;
text(textX,textY_3,thisLabel, ...
    'FontSize',this_fontSize, 'FontWeight', 'bold') ;


end


function [hcb, hyl] = add_big_colorbar(h, fontSize, units_map, manual_caxis)

thisPos = get(h,'Position') ;
hcb = colorbar(h,'Location','SouthOutside') ;
set(h,'Position',thisPos)
hcb.Position(2) = 0.075 ;
hyl = ylabel(hcb, sprintf('(%s)',units_map)) ;
hcb.FontSize = fontSize ;

% Mess with ticks
hcb.TickDirection = 'out' ;
if manual_caxis
    if max(hcb.Limits)~=2750
        error('Deal with tick marks here')
    end
    hcb.Ticks = -2750:500:2750 ;
%     keyboard
    for t = 1:length(hcb.Ticks)
        thisTick = hcb.Ticks(t) ;
        if abs(thisTick)>250 && rem(abs(thisTick)-250,1000)~=0
            hcb.TickLabels{t} = '' ;
        end
    end
end


% Move units label
hyl.Units = 'normalized' ;
hyl.Position(1) = 1.05 ;
hyl.Position(2) = -0.2 ;


end