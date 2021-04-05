%% Import things if needed

if any(any(strcmp(rowInfo,'co2'))) && ~exist('ts_co2_yr','var')
    import_co2
end

if any(any(strcmp(rowInfo,'temp'))) && ~exist('ts_temp_yr','var') ...
|| any(any(strcmp(rowInfo,'prec'))) && ~exist('ts_prec_yr','var')
    if any(strcmp({'harm3', 'ssp13'}, thisVer))
        warning('Reading LPJ-GUESS output temps instead of actual CMIP5 data')
        import_lpjgOut_temp_precip
    elseif any(any(strcmp(rowInfo,'prec')))
        error('Make script for precip equivalent to import_temperature.m (or set %s to read LPJ-GUESS output temps instead of actual CMIP5 data)', ...
            thisVer)
    else
        import_temperature
    end
end


%% Start bar graph

disp('Making bar graph...')

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
    if strcmp(thisVar, 'hotspot_area') ...
            || strcmp(thisVar, 'hotspotCSLF_area') ...
            || strcmp(thisVar, 'hotspotCSLF_area_nobioenergy')
        % End-hist
        if strcmp(thisVar, 'hotspot_area')
            hotspot_area_xy = repmat(hotspot_area_x, [1 1 length(years_endh)]) ;
        else
            hotspot_area_xy = repmat(hotspotCSLF_area_x, [1 1 length(years_endh)]) ;
        end
        tmpfrac_ntrl_x1y = garr_LU_d9.garr_xvyB(:,strcmp(garr_LU_d9.varNames,'NATURAL'),:) ;
        if strcmp(thisVar, 'hotspotCSLF_area_nobioenergy')
            tmpfrac_ntrl_x1y = tmpfrac_ntrl_x1y + ...
                garr_LU_d9.garr_xvyB(:,strcmp(garr_LU_d9.varNames,'CROPLAND'),:) ...
                .* garr_cropfracs_d9.garr_xvyB(:,strcmp(garr_cropfracs_d9.varNames,'Miscanthus'),:) ;
        end
        area_endh_y = thisConv * squeeze(nansum( ...
            hotspot_area_xy .* tmpfrac_ntrl_x1y, ...
            1)) ;
        mean_endh_v(c) = mean(area_endh_y) ;
        clear tmpfrac_ntrl_x1y
        % End-fut
        if strcmp(thisVar, 'hotspot_area')
            hotspot_area_xy = repmat(hotspot_area_x, [1 1 length(years_endf)]) ;
        else
            hotspot_area_xy = repmat(hotspotCSLF_area_x, [1 1 length(years_endf)]) ;
        end
        tmpfrac_ntrl_x1yr = garr_LU_d9.garr_xvyr(:,strcmp(garr_LU_d9.varNames,'NATURAL'),:,:) ;
        if strcmp(thisVar, 'hotspotCSLF_area_nobioenergy')
            tmpfrac_ntrl_x1yr = tmpfrac_ntrl_x1yr + ...
                garr_LU_d9.garr_xvyr(:,strcmp(garr_LU_d9.varNames,'CROPLAND'),:,:) ...
                .* garr_cropfracs_d9.garr_xvyr(:,strcmp(garr_cropfracs_d9.varNames,'Miscanthus'),:,:) ;
        end
        area_endf_yr = thisConv * squeeze(nansum( ...
            hotspot_area_xy .* tmpfrac_ntrl_x1yr, ...
            1)) ;
        mean_endf_vr(c,:) = mean(area_endf_yr,1) ;
        clear tmpfrac_ntrl_x1yr
        % Error bars
        if strcmp(sd_or_sem,'st. dev.')
            errb_endh_v(c) = std(area_endh_y) ;
            errb_endf_vr(c,:) = std(area_endf_yr, 1) ;
        elseif strcmp(sd_or_sem,'SEM')
            errb_endh_v(c) = sem_ssr(area_endh_y, years_endh) ;
            errb_endf_vr(c,:) = sem_ssr(area_endf_yr, years_endf) ;
        elseif ~isempty(sd_or_sem)
            error('sd_or_sem not recognized: %s', sd_or_sem) ;
        end
        
    elseif contains(thisVar,'emand')
        % Which commodity
        thisVar_spl = strsplit(thisVar,'.') ;
        thisDemand = thisVar_spl{1} ;
        thisCommod = thisVar_spl{2} ;
        clear thisVar_spl
        i = find(strcmp(commods,thisCommod)) ;
        % Get difference from first year to last decade
        eval(sprintf(...
            'ts_commodDemandTMP_1ir = ts_commod%s_yvr(1,i,:) ;', ...
            thisDemand)) ;
        eval(sprintf(...
            'ts_commodDemandTMP_Fir = ts_commod%s_yvr(ok_years,i,:) ;', ...
            thisDemand)) ;
        if length(unique(ts_commodDemandTMP_1ir)) > 1
            diffPct = 100 * (max(ts_commodDemandTMP_1ir) - min(ts_commodDemandTMP_1ir)) ...
                / min(ts_commodDemandTMP_1ir) ;
            tolerancePct = 0.1 ;
            if diffPct >= tolerancePct
                error('2010 demand for %s differs among runs by %0.2g%% (tolerance %0.2g%%)', ...
                    thisCommod, diffPct, tolerancePct)
            end
            warning('2010 demand for %s differs among runs by %0.2g%%; using mean as baseline', ...
                thisCommod, diffPct) ;
            mean_endh_v(c) = thisConv*mean(ts_commodDemandTMP_1ir,3) ;
        else
            mean_endh_v(c) = thisConv*ts_commodDemandTMP_1ir(:,:,1) ;
        end
        
        clear ts_commodDemandTMP_1ir
        errb_endh_v(c) = 0 ;
        ok_years = yearList_PLUMout>=min(years_endf) & yearList_PLUMout<=max(years_endf) ;
        mean_endf_vr(c,:) = thisConv*mean(ts_commodDemandTMP_Fir,1) ;
        if strcmp(sd_or_sem,'st. dev.')
            errb_endf_vr(c,:) = thisConv*std(ts_commodDemandTMP_Fir, 1) ;
        else
            errb_endf_vr(c,:) = thisConv*sem_ssr(ts_commodDemandTMP_Fir, years_endf) ;
        end
        clear ts_commodDemandTMP_Fir
        % Add marker indicating difference is from firstYear_PLUM, not
        % years_endh
        if ~any_notFirstDecade
            firstYear_PLUM = yearList_PLUMout(1) ;
            any_notFirstDecade = true ;
        end
        rowInfo{c,1} = [rowInfo{c,1} '*'] ;
        
    elseif strcmp(thisVar, 'pop')
        eval(sprintf('thisVar_yr = %s_yr ;', thisVar)) ;
        % Get difference from first year to last decade
        if length(unique(thisVar_yr(1,:))) > 1
            error('This code assumes all runs have identical 2010 %s!', thisVar)
        end
        mean_endh_v(c) = thisConv*thisVar_yr(1,1) ;
        errb_endh_v(c) = 0 ;
        ok_years = yearList_PLUMout>=min(years_endf) & yearList_PLUMout<=max(years_endf) ;
        mean_endf_vr(c,:) = thisConv*mean(thisVar_yr(ok_years,:),1) ;
        if strcmp(sd_or_sem,'st. dev.')
            errb_endf_vr(c,:) = thisConv*std(thisVar_yr(ok_years,:), 1) ;
        else
            errb_endf_vr(c,:) = thisConv*sem_ssr(thisVar_yr(ok_years,:), years_endf) ;
        end
        % Add marker indicating difference is from firstYear_PLUM, not
        % years_endh
        if ~any_notFirstDecade
            firstYear_PLUM = yearList_PLUMout(1) ;
            any_notFirstDecade = true ;
        end
        rowInfo{c,1} = [rowInfo{c,1} '*'] ;
        clear thisVar_yr
        
    elseif any(strcmp({'tot_runoff','aevapaaet','prec'}, thisVar))
        % Use last 30 years of each century
        ts_thisVarTMP_bl = eval(['ts_' thisVar '_bl(yearList_baseline>=1971 & yearList_baseline<=2000) ;']) ;
        ts_thisVarTMP_yr = eval(['ts_' thisVar '_yr(yearList_future>=2071 & yearList_future<=2100,:) ;']) ;
        mean_endh_v(c) = thisConv*mean(ts_thisVarTMP_bl) ;
        mean_endf_vr(c,:) = thisConv*mean(ts_thisVarTMP_yr,1) ;
        clear tmp_endh_y tmp_endf_yr
        if strcmp(sd_or_sem,'st. dev.')
            errb_endh_v(c) = thisConv*std(ts_thisVarTMP_bl) ;
            errb_endf_vr(c,:) = thisConv*std(ts_thisVarTMP_yr, 1) ;
        elseif strcmp(sd_or_sem,'SEM')
            errb_endh_v(c) = thisConv*sem_ssr(ts_thisVarTMP_bl) ;
            errb_endf_vr(c,:) = thisConv*sem_ssr(ts_thisVarTMP_yr, years_endf) ;
        elseif ~isempty(sd_or_sem)
            error('sd_or_sem not recognized: %s', sd_or_sem) ;
        end
        % Add marker indicating difference is from last 3 decades, not
        % years_endh and years_endf
        if ~any_notFirstDecade
            any_notFirstDecade = true ;
        end
        rowInfo{c,1} = [rowInfo{c,1} '^#'] ;
    else
        % Is this per-capita?
        is_percapita = contains(thisVar,'PC') ;
        thisVar = strrep(thisVar, 'PC', '') ;
        
        % Is this a BVOC variable?
        is_bvoc = contains(thisVar,{'aiso','amon'}) ;
        
        % Add multiple variables, if specified
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
        
        % Get means and variation
        try
            ts_thisVarTMP_bl = eval(['ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)) ;']) ;
        catch ME
            keyboard
        end

        ts_thisVarTMP_yr = eval(['ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:) ;']) ;
        if is_percapita
            [~,IA] = intersect(yearList_pop, years_endh) ;
            if length(IA) ~= 1
                error('This code assumes that there is only one year of overlap between population data and years_endh')
            elseif abs(max(pop_yr(IA,:)) - min(pop_yr(IA,:))) > 0
                error('This code assumes that all runs have identical starting population!')
            end
            ts_thisVarTMP_bl = ts_thisVarTMP_bl / pop_yr(IA,1) ;
            clear IA
            [tmp_years,IA] = intersect(yearList_pop, years_endf) ;
            if ~isequal(shiftdim(tmp_years), shiftdim(years_endf))
                error('Problem with overlap of yearList_pop and years_endf')
            end
            ts_thisVarTMP_yr = ts_thisVarTMP_yr ./ pop_yr(IA,:) ;
            clear IA tmp_years
            % Add marker indicating difference is from firstYear_PLUM, not
            % years_endh
            if ~any_notFirstDecade
                firstYear_PLUM = yearList_PLUMout(1) ;
                any_notFirstDecade = true ;
            end
            rowInfo{c,1} = [rowInfo{c,1} '*'] ;
        end
        mean_endh_v(c) = thisConv*mean(ts_thisVarTMP_bl) ;
        mean_endf_vr(c,:) = thisConv*mean(ts_thisVarTMP_yr,1) ;
        if strcmp(thisVar,'temp') && strcmp(rowInfo{c,6},[char(176) 'C'])
            mean_endh_v(c) = mean_endh_v(c) - 273.15 ;
            mean_endf_vr(c,:) = mean_endf_vr(c,:) - 273.15 ;
        end
        clear tmp_endh_y tmp_endf_yr
        if strcmp(sd_or_sem,'st. dev.')
            errb_endh_v(c) = thisConv*std(ts_thisVarTMP_bl) ;
            errb_endf_vr(c,:) = thisConv*std(ts_thisVarTMP_yr, 1) ;
        elseif strcmp(sd_or_sem,'SEM')
            errb_endh_v(c) = thisConv*sem_ssr(ts_thisVarTMP_bl) ;
            errb_endf_vr(c,:) = thisConv*sem_ssr(ts_thisVarTMP_yr, years_endf) ;
        elseif ~isempty(sd_or_sem)
            error('sd_or_sem not recognized: %s', sd_or_sem) ;
        end
        clear ts_thisVar_yr % Does nothing if does not exist
        
        % Do not show bars for BVOC emissions with constant climate (daily
        % temperature range is not detrended)
        if is_bvoc && any(contains(runList,'constClim'))
            mean_endf_vr(c,contains(runList,'constClim')) = NaN ;
            errb_endf_vr(c,contains(runList,'constClim')) = NaN ;
        end
    end
end

% Calculate % difference
mean_endh_vr = repmat(mean_endh_v,[1 Nruns]) ;
mean_diff_vr = mean_endf_vr - mean_endh_vr ;
mean_diffPct_vr = 100 * (mean_diff_vr ./ mean_endh_vr) ;
errb_endh_vr = repmat(errb_endh_v,[1 Nruns]) ;
errb_diff_vr = sqrt(errb_endh_vr.^2 + errb_endf_vr.^2) ;
% EQUIVALENT TO BELOW % errb_diffPct_vr = 100 * ((mean_endf_vr+errb_diff_vr-mean_endh_vr)./mean_endh_vr - (mean_endf_vr-mean_endh_vr)./mean_endh_vr) ;
errb_diffPct_vr = 100 * (errb_diff_vr ./ mean_endh_vr) ;
if isempty(sd_or_sem)
    errb_diff_vr = zeros(size(errb_diff_vr)) ;
    errb_diffPct_vr = zeros(size(errb_diff_vr)) ;
end



%% Make figure

groupnames = rowInfo(:,1) ;
% unitnames = rowInfo(:,6) ;
% groupnames = strcat(groupnames, strcat(strcat('\n(', unitnames), ')')) ;

figure('Color','w','Position',figure_position) ;

if strcmp(orientation,'v')
    subplot_tight(1,1,1,[0.2 0.04])
    bar(mean_diffPct_vr, 'grouped') ;
    set(gca, 'XTickLabel', groupnames) ;
    xtickangle(45) ;
elseif strcmp(orientation,'h')
    subplot_tight(1,1,1,[0.04 0.2])
    barh(mean_diffPct_vr, 'grouped') ;
    set(gca, 'YTickLabel', groupnames) ;
else
    error('orientation (%s) not recognized', orientation)
end


%% Add error bars and/or data labels
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
h_barLabels = {} ;
x = 0 ;
for i = 1:nbars
    % Calculate center of each bar
    if strcmp(orientation,'v')
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars) ;
        he = errorbar(x, mean_diffPct_vr(:,i), errb_diffPct_vr(:,i), 'linestyle', 'none', 'linewidth', 0.5) ;
    else
        y = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars) ;
        he = errorbar(mean_diffPct_vr(:,i), y, errb_diffPct_vr(:,i), 'horizontal', 'linestyle', 'none', 'linewidth', 0.5) ;
    end
    if ~isempty(sd_or_sem)
        set(he, 'Color', errbar_color) ;
    else
        set(he, 'Visible', 'off') ;
    end
    for g = 1:ngroups
        x = x + 1 ;
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
        elseif ngroups*nbars<=40
            ht.FontSize = 10 ;
        else
            ht.FontSize = 8 ;
        end
        if strcmp(orientation,'v')
            error('Write realignment code for vertical orientation')
        end
        if strcmp(orientation,'h') && thisX < 0
            ht.Position(1) = thisPad ;
        end
        h_barLabels{x} = ht ;
    end
end
hold off
clear x


%% Add separators

if ~isempty(where2sep)
    hold on
    if strcmp(orientation,'v')
        error('Deal with this for vertically-oriented bars!')
        ylims = get(gca,'YLim') ;
        for i = 1:length(where2sep)
            x = where2sep(i) ;
            line([x x], ylims, 'Color', 0.8*(ones(3,1)), 'LineWidth', 10) ;
        end
    else
        xlims = get(gca,'XLim') ;
        sepwidth = groupwidth * 0.2 ;
        for i = 1:length(where2sep)
            plot(xlims,where2sep(i)*ones(size(xlims)), '--k')
            ht = text(xlims(1), where2sep(i), sep_labels{i}) ;
            if ngroups*nbars<=40
                sep_label_fontSize = 12 ;
            else
                sep_label_fontSize = 12 ;
            end
            set(ht, ...
                'HorizontalAlignment', 'center', ...
                'FontSize', sep_label_fontSize, ...
                'FontWeight', 'bold', ...
                'FontAngle', 'italic', ...
                'BackgroundColor', 'w')
        end
    end
    hold off
end


%% Finish

% Add legend, and labels
if ~isempty(sd_or_sem)
    axis_label = ['Change ' plusminus ' across-year ' sd_or_sem ' (%)'] ;
else
    axis_label = 'Change (%)' ;
end
if strcmp(orientation,'v')
    legend(runList, 'Location', 'Northwest')
%     hxl = xlabel('Indicator') ;
%     hxl.FontWeight = 'bold' ;
    hyl = ylabel(axis_label) ;
    hyl.FontWeight = 'bold' ;
else
    legend(runList, 'Location', 'Northeast')
%     hyl = ylabel('Indicator') ;
%     hyl.FontWeight = 'bold' ;
    hxl = xlabel(axis_label) ;
    hxl.FontWeight = 'bold' ;
end
set(gca, 'FontSize', fontSize) ;

% Reposition axes
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

% Get rid of extraneous space at top
ylims = get(gca,'YLim') ;
set(gca,'YLim',[0.4 max(ylims)])

% Tweak bar label locations
xlims = get(gca,'XLim') ;
hf = gcf ;
hf.Units = 'pixels' ;
ha = gca ;
axis_position_norm = ha.Position ;
ha.Units = 'pixels' ;
for x = 1:length(h_barLabels)
    thisH = h_barLabels{x} ;
    thisH_units_orig = thisH.Units ;
    thisH.Units = 'pixels' ;
    % Move bar labels that have fallen off the right edge
    if ha.Position(1) + thisH.Extent(1) + thisH.Extent(3) > hf.Position(3)
        tmp = strsplit(thisH.String,' ') ;
        thisH.String = {tmp{1}, strjoin(tmp(2:end))} ;
        clear tmp
        if ha.Position(1) + thisH.Extent(1) + thisH.Extent(3) > hf.Position(3)
            thisH.Units = thisH_units_orig ;
            thisH.String = strjoin(thisH.String) ;
            newPosition = thisH.Position ;
            if contains(thisH.String,'^')
                newPosition(2) = newPosition(2) + 0.13*ngroups/9 ; % Not minus because we've flipped the Y axis
            else
                newPosition(2) = newPosition(2) + 0.15*ngroups/9 ; % Not minus because we've flipped the Y axis
            end
            excess = thisH.Extent(1)+thisH.Extent(3) - xlims(2) ;
            newPosition(1) = newPosition(1) - excess ;
            thisH.Position = newPosition ;
        end
    end
    thisH.Units = thisH_units_orig ;
    % Nudge up bar labels that were pushed down because of a superscript
    if contains(thisH.String,'^')
        newPosition = thisH.Position ;
        newPosition(2) = newPosition(2) - 0.2*thisH.Extent(4) ; % Not plus because we've flipped the Y axis
        thisH.Position = newPosition ;
    end
end
ha.Units = 'normalized' ;
ha.Position = axis_position_norm ;

% Nudge legend into a good position
ylims = get(gca,'YLim') ;
nudge = [0 -(groupwidth+sepwidth*1.5)/(ylims(2)-ylims(1)) 0 0] ;
[does_overlap, hlegend] = test_for_overlap(hf) ;
while does_overlap
    if hlegend.Position(2) + nudge(2) > ylims(1)
        hlegend.Position = hlegend.Position + nudge ;
    else
        warning('Setting legend position to ''Best''')
        hlegend.Location = 'Best' ;
        break
    end
    [does_overlap, hlegend] = test_for_overlap(hf) ;
end


disp('Done.')