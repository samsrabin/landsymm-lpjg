function make_LUfrac_fig_v5(...
    area_orig_r, area_harm_r, ...
    orig_frac_YXrH, harm_frac_YXrH, ...
    runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, ...
    Nruns, thisPos, units_map, units_total, ...
    do_caps, bins_lowBnds, this_colormap_name, col_titles, ...
    lines_overlay)

figure('Color','w','Position',thisPos) ;
gcas_lu1 = {} ;
gcas_lu2 = {} ;
colorlim = 0 ;
Nbins = length(bins_lowBnds) ;
for r = 1:Nruns
    runName = runList{r} ;
    i1 = (r-1)*2 + 1 ;
    i2 = i1 + 1 ;
    [hc, hp] = actually_make_fig(...
        orig_frac_YXrH(:,:,r), harm_frac_YXrH(:,:,r), ...
        area_orig_r(min(r, length(area_orig_r))), ...
        area_harm_r(min(r, length(area_harm_r))), ...
        runName, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, i1, i2, ...
        units_total, ...
        do_caps, this_colormap_name, bins_lowBnds, col_titles, ...
        lines_overlay) ;
    gcas_lu1 = [gcas_lu1 hc] ;
    gcas_lu2 = [gcas_lu2 hp] ;
end
Xshift = 0.05 ;
for r = 1:Nruns
    hc = gcas_lu1(r) ;
    hp = gcas_lu2(r) ;
    hc.Position(1) = hc.Position(1) + Xshift ;
    Yshift = 0.04*(r-3) ;
    hc.Position(2) = hc.Position(2) + Yshift ;
    hp.Position(2) = hp.Position(2) + Yshift ;
end

hcb = add_big_colorbar(gcas_lu1(Nruns), fontSize, units_map, bins_lowBnds) ;
hcb.Position(1) = 0.5 - hcb.Position(3)/2 ;

end


function cmap = get_colormap(colormap_name, Nbins)

flip_colormap = strcmp(colormap_name(1),'-') ;
if flip_colormap
    colormap_name = colormap_name(2:end) ;
end
try 
    cmap = colormap(colormap_name) ;
catch ME
    cmap = brighten(brewermap(Nbins,colormap_name),-0.3) ;
end
if flip_colormap
    cmap = flipud(cmap) ;
end

end


function [hc, hp, ht] = actually_make_fig(...
    orig_frac_YX, harm_frac_YX, ...
    area_orig, area_harm, ...
    runName, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, i1, i2, ...
    units_total, ...
    do_caps, this_colormap_name, bins_lowBnds, col_titles, ...
    lines_overlay)

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

% lu1
subplot_tight(ny,nx,i1,spacing) ;
if i1 > 1
    set(gca,'Position', get(gca,'Position') + [0 shiftup 0 0])
end
is_lu1 = true ;
ht = plot_map(orig_frac_YX, bins_lowBnds, this_colormap_name, ...
    is_lu1, i1, runName, fontSize, ...
    area_orig, units_total, ...
    textX, textY_1, textY_2, ...
    do_caps, lines_overlay) ;
hc = gca ;

% Add lu1 column title
if i1==1
    htmp = text(0, 0, col_titles{1}, ...
        'FontSize', colTitle_fontSize, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center') ;
    htmp.Units = 'normalized' ;
    htmp.Position = [0.5 colTitle_yPos 0] ;
end

% lu2
subplot_tight(ny,nx,i2,spacing) ;
if i1 > 1
    set(gca,'Position', get(gca,'Position') + [0 shiftup 0 0])
end
is_lu1 = false ;
plot_map(harm_frac_YX, bins_lowBnds, this_colormap_name, ...
    is_lu1, i1+1, runName, fontSize, ...
    area_harm, units_total, ...
    textX, textY_1, textY_2, ...
    do_caps, lines_overlay) ;
hp = gca ;

% Add lu2 column title
if i1==1
    htmp = text(0, 0, col_titles{2}, ...
        'FontSize', colTitle_fontSize, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center') ;
    htmp.Units = 'normalized' ;
    htmp.Position = [0.5 colTitle_yPos 0] ;
end

pause(0.1)

end


function add_totals( ...
    area_in, ...
    units_total, ...
    textX, textY_1, textY_2, fontSize, ...
    i1, do_caps)

% this_fontSize = fontSize - 2 ;
this_fontSize = fontSize ;

% pctDiff = round(100*total_diff./area_bl,1) ;
% pctDiff_str = num2str(pctDiff) ;
diff_format = '%0.1f %s' ;
% if pctDiff>0
%     pctDiff_str = ['+' pctDiff_str] ;
%     diff_format = ['+' diff_format] ;
% end
theText = sprintf(diff_format, area_in, units_total) ;

htext = text(0, 0, theText,...
    'FontSize',this_fontSize) ;
htext.Units = 'normalized' ;
htext.Position = [textX textY_1 0] ;

% htext2 = text(0, 0, ['(' pctDiff_str '%)'], 'FontSize', this_fontSize) ;
% htext2.Units = 'normalized' ;
% htext2.Position = [textX textY_2 0] ;

textY_3 = max(textY_1, textY_2) + abs(textY_1 - textY_2) ;

thisChar = char(i1 + 64) ;
if do_caps==1
    thisChar = upper(thisChar) ;
elseif do_caps==-1
    thisChar = lower(thisChar) ;
end
thisLabel = sprintf('(%s)', thisChar) ;
text(textX,textY_3,thisLabel, ...
    'FontSize',this_fontSize, 'FontWeight', 'bold', ...
    'Units', 'normalized') ;


end


function [hcb, hyl] = add_big_colorbar(h, fontSize, units_map, bins_lowBnds)

thisPos = get(h,'Position') ;
hcb = colorbar(h,'Location','SouthOutside') ;
set(h,'Position',thisPos)
hcb.Position(2) = 0.075 ;
hyl = ylabel(hcb, sprintf('(%s)',units_map)) ;
hcb.FontSize = fontSize ;

% Mess with ticks
hcb.TickDirection = 'out' ;
% Nbins = length(bins_lowBnds) ;
% hcb.Ticks = 1:(Nbins+1) ;
% bins_lowBnds_str = strrep(cellstr(num2str(bins_lowBnds')), ' ', '') ;
% bins_lowBnds_str = [bins_lowBnds_str ; {'100'}] ;
% hcb.TickLabels = bins_lowBnds_str ;

% Move units label
hyl.Units = 'normalized' ;
hyl.Position(1) = 1.09 ;
hyl.Position(2) = -0.5 ;


end


function ht = plot_map(map_YX, bins_lowBnds, this_colormap_name, ...
    is_lu1, i_in, runName, fontSize, ...
    area_in, units_total, ...
    textX, textY_1, textY_2, ...
    do_caps, lines_overlay)

% Create binned version of map
Nbins = length(bins_lowBnds) ;
map_YX_bin = nan(size(map_YX)) ;
for b = 1:Nbins
    map_YX_bin(map_YX >= bins_lowBnds(b)) = b+0.1 ;
end

% Get colormap
this_colormap = get_colormap(this_colormap_name, Nbins) ;

% Plot map
R = georefcells([-90 90],[-180 180],[360 720]) ;
hwm = worldmap(map_YX_bin, R) ;
setm(hwm,'MapLatLimit',[-55 90])
setm(hwm,'mapprojection','Robinson', 'frame','off','grid','on')
set(findall(hwm,'Tag','PLabel'),'visible','off')
set(findall(hwm,'Tag','MLabel'),'visible','off')
geoshow(map_YX_bin, R, ...
    'DisplayType', 'surface')
colormap(gca, this_colormap)
% caxis([1 Nbins+1])
caxis([0 100])

% Plot continent outlines
if ~isempty(lines_overlay)
    hold on
    geoshow(lines_overlay, 'FaceAlpha', 0)
    hold off
end

% Add text
set(gca,'FontSize',fontSize)
if is_lu1
    ht = text( ...
        ...-0.01, ...
        0.02, ...
        0.4, ...
        runName, ...
        'Units', 'normalized', ...
        'Rotation', 90, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', fontSize+4, 'FontWeight', 'bold') ;
end
add_totals(area_in, units_total, ...
    textX, textY_1, textY_2, fontSize, ...
    i_in, do_caps)


end