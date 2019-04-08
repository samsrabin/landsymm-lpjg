Nvars = size(rowInfo,1) ;
mean_endh_v = nan(Nvars,1) ;
mean_endf_vr = nan(Nvars,Nruns) ;
errb_endh_v = nan(Nvars,1) ;
errb_endf_vr = nan(Nvars,Nruns) ;
any_notFirstDecade = false ;
for c = 1:Nvars
    
    % Get values
    thisVar = rowInfo{c,2} ;
    thisConv = rowInfo{c,3} ;
    if strcmp(thisVar, 'hotspot_area')
        % End-hist
        hotspot_area_YXy = repmat(hotspot_area_YX, [1 1 1 length(years_endh)]) ;
        area_endh_y = thisConv * squeeze(nansum(nansum( ...
            hotspot_area_YXy .* maps_LU_d9.maps_YXvyB(:,:,strcmp(maps_LU_d9.varNames,'NATURAL'),:), ...
            1), 2)) ;
        mean_endh_v(c) = mean(area_endh_y) ;
        % End-fut
        hotspot_area_YXy = repmat(hotspot_area_YX, [1 1 1 length(years_endf)]) ;
        area_endf_yr = thisConv * squeeze(nansum(nansum( ...
            hotspot_area_YXy .* maps_LU_d9.maps_YXvyr(:,:,strcmp(maps_LU_d9.varNames,'NATURAL'),:,:), ...
            1), 2)) ;
        mean_endf_vr(c,:) = mean(area_endf_yr,1) ;
        % Error bars
        if strcmp(sd_or_sem,'st. dev.')
            errb_endh_v(c) = std(area_endh_y) ;
            errb_endf_vr(c,:) = std(area_endf_yr, 1) ;
        else
            errb_endh_v(c) = sem_ssr(area_endh_y, years_endh) ;
            errb_endf_vr(c,:) = sem_ssr(area_endf_yr, years_endf) ;
        end
        
    elseif contains(thisVar,'demand')
        % Which commodity
        thisCommod = strsplit(thisVar,'.') ;
        thisCommod = thisCommod{2} ;
        i = find(strcmp(commods,thisCommod)) ;
        % Get difference from first year to last decade
        if length(unique(ts_commodDemand_yvr(1,i,:))) > 1
            error('This code assumes all runs have identical 2010 demand!')
        end
        mean_endh_v(c) = thisConv*ts_commodDemand_yvr(1,i,1) ;
        errb_endh_v(c) = 0 ;
        ok_years = yearList_PLUMout>=min(years_endf) & yearList_PLUMout<=max(years_endf) ;
        mean_endf_vr(c,:) = thisConv*mean(ts_commodDemand_yvr(ok_years,i,:),1) ;
        if strcmp(sd_or_sem,'st. dev.')
            errb_endf_vr(c,:) = thisConv*std(ts_commodDemand_yvr(ok_years,i,:), 1) ;
        else
            errb_endf_vr(c,:) = thisConv*sem_ssr(ts_commodDemand_yvr(ok_years,i,:), years_endf) ;
        end
        % Add marker indicating difference is from firstYear_PLUM, not
        % years_endh
        if ~any_notFirstDecade
            firstYear_PLUM = yearList_PLUMout(1) ;
            any_notFirstDecade = true ;
        end
        rowInfo{c,1} = [rowInfo{c,1} '*'] ;
    else
        if contains(thisVar,'+')
            theseVars = strsplit(thisVar,'+') ;
            thisVar = theseVars{1} ;
            ts_thisVar_bl = eval(['ts_' thisVar '_bl']) ;
            ts_thisVar_yr = eval(['ts_' thisVar '_yr']) ;
            for v = 2:length(theseVars)
                thisVar = theseVars{v} ;
                ts_thisVar_bl = ts_thisVar_bl + eval(['ts_' thisVar '_bl']) ;
                ts_thisVar_yr = ts_thisVar_yr + eval(['ts_' thisVar '_yr']) ;
            end
            thisVar = 'thisVar' ;
        end
        mean_endh_v(c) = thisConv*eval(['mean(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
        mean_endf_vr(c,:) = thisConv*eval(['mean(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
        if strcmp(sd_or_sem,'st. dev.')
            errb_endh_v(c) = thisConv*eval(['std(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
            errb_endf_vr(c,:) = thisConv*eval(['std(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:), 1)']) ;
        else
            errb_endh_v(c) = thisConv*eval(['sem_ssr(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)), years_endh)']) ;
            errb_endf_vr(c,:) = thisConv*eval(['sem_ssr(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:), years_endf)']) ;
        end
        clear ts_thisVar_yr % Does nothing if does not exist
    end
end

% Calculate % difference
mean_endh_vr = repmat(mean_endh_v,[1 Nruns]) ;
errb_endh_vr = repmat(errb_endh_v,[1 Nruns]) ;
mean_diff_vr = mean_endf_vr - mean_endh_vr ;
mean_diffPct_vr = 100 * (mean_diff_vr ./ mean_endh_vr) ;
errb_diff_vr = sqrt(errb_endh_vr.^2 + errb_endf_vr.^2) ;
% EQUIVALENT TO BELOW % errb_diffPct_vr = 100 * ((mean_endf_vr+errb_diff_vr-mean_endh_vr)./mean_endh_vr - (mean_endf_vr-mean_endh_vr)./mean_endh_vr) ;
errb_diffPct_vr = 100 * (errb_diff_vr ./ mean_endh_vr) ;

% Make figure
if strcmp(orientation,'v')
    figure('Color','w','Position',figurePos) ;
    subplot_tight(1,1,1,[0.2 0.04])
    bar(mean_diffPct_vr, 'grouped') ;
    set(gca, 'XTickLabel', rowInfo(:,1)) ;
    xtickangle(45) ;
elseif strcmp(orientation,'h')
    figure('Color','w','Position',[1    33   720   772]) ;
    subplot_tight(1,1,1,[0.04 0.2])
    barh(mean_diffPct_vr, 'grouped') ;
    set(gca, 'YTickLabel', rowInfo(:,1)) ;
else
    error('orientation (%s) not recognized', orientation)
end

%%%%%%%%%%%%%%%%%%%%%%
%%% Add error bars %%%
%%%%%%%%%%%%%%%%%%%%%%
% https://www.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab

% Finding the number of groups and the number of bars in each group
ngroups = size(mean_diffPct_vr, 1);
nbars = size(mean_diffPct_vr, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
hold on
if strcmp(orientation,'v')
    lims = get(gca,'YLim') ;
    thisPad_factor = 0.005 ;
else
    lims = get(gca,'XLim') ;
    thisPad_factor = 0.01 ;
end
lims_diff = lims(2) - lims(1) ;
thisPad = thisPad_factor * lims_diff ;
for i = 1:nbars
    % Calculate center of each bar
    if strcmp(orientation,'v')
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars) ;
        he = errorbar(x, mean_diffPct_vr(:,i), errb_diffPct_vr(:,i), 'linestyle', 'none', 'linewidth', 0.5) ;
    else
        y = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars) ;
        he = errorbar(mean_diffPct_vr(:,i), y, errb_diffPct_vr(:,i), 'horizontal', 'linestyle', 'none', 'linewidth', 0.5) ;
    end
    set(he, 'Color', errbar_color) ;
    for g = 1:ngroups
        if strcmp(orientation,'v')
            thisX = he.XData(g) ;
            thisY = he.YData(g) + he.YPositiveDelta(g) ;
            if thisY < 0
                thisY = he.YData(g) - he.YNegativeDelta(g) - thisPad ;
            else
                thisY = thisY + thisPad ;
            end
        else
            thisY = he.YData(g) ;
            thisX = he.XData(g) + he.XPositiveDelta(g) ;
            if thisX < 0
                thisX = he.XData(g) - he.XNegativeDelta(g) - thisPad ;
            else
                thisX = thisX + thisPad ;
            end
        end
        % Mean diff (%)
        thisText1 = sprintf('%.0f%%', mean_diffPct_vr(g,i)) ;
        if mean_diffPct_vr(g,i) > 0
            thisText1 = ['+' thisText1] ;
        end
        % Mean diff (absolute units)
        thisText2 = sprintf([rowInfo{g,4} ' %s'], mean_diff_vr(g,i), rowInfo{g,6}) ;
        if mean_diffPct_vr(g,i) > 0
            thisText2 = ['+' thisText2] ;
        end
        % Combine and add
        thisText = sprintf('%s (%s)', thisText1, thisText2) ;
        ht = text(thisX, thisY, thisText) ;
        if strcmp(orientation,'v')
            ht.Rotation = 90 ;
        else
            ht.FontSize = 8 ;
        end
        if (strcmp(orientation,'v') && thisY < 0) ...
        || (strcmp(orientation,'h') && thisX < 0)
            ht.HorizontalAlignment = 'right' ;
        end
    end
end
hold off


%%%%%%%%%%%%%%%%%%%%%%
%%% Add separators %%%
%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(where2sep)
    hold on
    if strcmp(orientation,'v')
        ylims = get(gca,'YLim') ;
        for i = 1:length(where2sep)
            x = where2sep(i) ;
            line([x x], ylims, 'Color', 0.8*(ones(3,1)), 'LineWidth', 10) ;
        end
    else
        xlims = get(gca,'XLim') ;
        for i = 1:length(where2sep)
            y = where2sep(i) ;
            line(xlims, [y y], 'Color', 0.8*(ones(3,1)), 'LineWidth', 5) ;
        end
    end
    hold off
end

%%%%%%%%%%%%%%
%%% Finish %%%
%%%%%%%%%%%%%%

% Add legend, and labels
if strcmp(orientation,'v')
    legend(runList, 'Location', 'Northwest')
    hxl = xlabel('Indicator') ;
    hyl = ylabel(['Change +/- ' sd_or_sem ' (%)']) ;
else
    legend(runList, 'Location', 'Northeast')
    hyl = ylabel('Indicator') ;
    hxl = xlabel(['Change +/- ' sd_or_sem ' (%)']) ;
end
set(gca, 'FontSize', fontSize) ;
hxl.FontWeight = 'bold' ;
hyl.FontWeight = 'bold' ;

% Reposition axexs
h = gca ;
h.Units = 'normalized' ;
op = get(h, 'OuterPosition') ;
if strcmp(orientation,'v')
    op(2) = 0.02 ;
    op(4) = 1.05 ;
    op(3) = 1.24 ;
    set(h, 'OuterPosition', op) ;
else
    op(1) = 0 ;
    op(3) = 1 ;
    op(2) = -0.17 ;
    op(4) = 1.2678 ;
    set(h, 'OuterPosition', op) ;
    set(h,'YDir','reverse')
end

% Remove surrounding chartjunk
box off