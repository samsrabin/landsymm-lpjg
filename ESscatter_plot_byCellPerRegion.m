function ESscatter_plot_byCellPerRegion( ...
    mapstructX, thisVarX, labelX, ...
    mapstructY, thisVarY, labelY, ...
    regions_YX, region_list, runList)

%%%%% Options
fontSize = 14 ;
markerSize = 36*2 ;
markerEdgeAlpha = 0.25 ;

% Check inputs
ESscatter_check_inputs(...
    mapstructX, thisVarX, ...
    mapstructY, thisVarY, ...
    regions_YX) ;

% Get information about inputs
regions_codes = unique(regions_YX(~isnan(regions_YX))) ;
Nregions = length(regions_codes) ;
Nruns = size(mapstructX.maps_YXvyr,5) ;
if Nregions==7
    nx = 4 ;
    ny = 2 ;
    spacing = [0.1 0.05] ;   % [v h]
else
    error('Set up for %d regions', Nregions)
end

% Get mean for variable in question
[mapsX_YXB, mapsX_YXr, mapsY_YXB, mapsY_YXr]= ...
    ESscatter_getMeanOverPeriod(...
        mapstructX, thisVarX, ...
        mapstructY, thisVarY) ;

% Make figures
for r = 1:Nruns
    hs = [] ; % Axis handles
    xmin = NaN ;
    xmax = NaN ;
    ymin = NaN ;
    ymax = NaN ;
    
    figure('Color','w','Position',figurePos) ;
    for c = 1:Nregions
        thisCountry = regions_codes(c) ;
        hs(c) = subplot_tight(ny,nx,c,spacing) ;
        tmp_YX = mapsX_YXr(:,:,r) ; 
        thisDiffX = tmp_YX(regions_YX==thisCountry) - mapsX_YXB(regions_YX==thisCountry) ;
        tmp_YX = mapsY_YXr(:,:,r) ;
        thisDiffY = tmp_YX(regions_YX==thisCountry) - mapsY_YXB(regions_YX==thisCountry) ;
        h = scatter(thisDiffX, thisDiffY) ;
        h.SizeData = markerSize ;
        h.MarkerEdgeAlpha = markerEdgeAlpha ;
        title(sprintf('%s: %s', runList{r}, region_list{c}))
        xlabel(['\Delta ' labelX])
        ylabel(['\Delta ' labelY])
        set(gca,'FontSize',fontSize)
        xmin = min(xmin, min(get(gca,'XLim'))) ;
        ymin = min(ymin, min(get(gca,'YLim'))) ;
        xmax = max(xmax, max(get(gca,'XLim'))) ;
        ymax = max(ymax, max(get(gca,'YLim'))) ;
    end
    
    % Equalize axes
    xlims = [xmin xmax] ;
    ylims = [ymin ymax] ;
    for c = 1:Nregions
        set(hs(c), 'XLim', xlims, 'YLim', ylims)
    end
    
    % Add zero-lines
    if xlims(1)<0 && xlims(2)>0
        for c = 1:Nregions
            hold(hs(c),'on')
            plot(hs(c),[0 0],ylims,'--k')
            hold(hs(c),'off')
        end
    end
    if ylims(1)<0 && ylims(2)>0
        for c = 1:Nregions
            hold(hs(c),'on')
            plot(hs(c),xlims,[0 0],'--k')
            hold(hs(c),'off')
        end
    end
end



end