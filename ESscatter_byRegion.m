function ESscatter_byRegion( ...
    mapsX_YXyr, labelX, ...
    mapsY_YXyr, labelY, ...
    countries_YX, runList)

%%%%% Options
fontSize = 14 ;
markerSize = 15 ;
spacing = [0.1 0.1] ;   % [vertical horizontal]

% Get information about inputs
countries_codes = unique(countries_YX(~isnan(countries_YX))) ;
Ncountries = length(countries_codes) ;
Nruns = length(runList) ;

% Get mean maps across all years
mapsX_YXr = squeeze(mean(mapsX_YXyr,3)) ;
mapsY_YXr = squeeze(mean(mapsY_YXyr,3)) ;

% Get means for each country
resX_cr = nan(Ncountries,Nruns) ;
resY_cr = nan(Ncountries,Nruns) ;
for c = 1:Ncountries
    thisCountry = countries_codes(c) ;
    for r = 1:Nruns
        tmp_YX = mapsX_YXr(:,:,r) ;        
        resX_cr(c,r) = nanmean(tmp_YX(countries_YX==thisCountry)) ;
        tmp_YX = mapsY_YXr(:,:,r) ;        
        resY_cr(c,r) = nanmean(tmp_YX(countries_YX==thisCountry)) ;
    end
end

% Trim countries that are all NaN in both datasets
bad_countries = find(all(isnan(resX_cr),2) & all(isnan(resY_cr),2)) ;
resX_cr(bad_countries,:) = [] ;
resY_cr(bad_countries,:) = [] ;
% Trim countries that are all NaN in only one dataset
bad_countries = find(all(isnan(resX_cr),2)) ;
if ~isempty(bad_countries)
    warning('Some countries were bad in X but not Y. Removing from both.')
    resX_cr(bad_countries,:) = [] ;
    resY_cr(bad_countries,:) = [] ;
end
bad_countries = find(all(isnan(resY_cr),2)) ;
if ~isempty(bad_countries)
    warning('Some countries were bad in Y but not X. Removing from both.')
    resX_cr(bad_countries,:) = [] ;
    resY_cr(bad_countries,:) = [] ;
end
if any(any(isnan(resX_cr)))
    error('NaN remaining in resX_cr')
elseif any(any(isnan(resY_cr)))
    error('NaN remaining in resY_cr')
end

% Make figure
figure('Color','w','Position',figurePos) ;
hs = [] ; % Axis handles
xmin = NaN ;
xmax = NaN ;
ymin = NaN ;
ymax = NaN ;
for r = 1:Nruns
    hs(r) = subplot_tight(2,2,r,spacing) ;
    plot(resX_cr(:,r),resY_cr(:,r),'.', ...
        'MarkerSize',markerSize)
    title(runList{r})
    xlabel(labelX)
    ylabel(labelY)
    set(gca,'FontSize',fontSize)
    xmin = min(xmin, min(get(gca,'XLim'))) ;
    ymin = min(ymin, min(get(gca,'YLim'))) ;
    xmax = max(xmax, max(get(gca,'XLim'))) ;
    ymax = max(ymax, max(get(gca,'YLim'))) ;
end

% Equalize axes
xlims = [xmin xmax] ;
ylims = [ymin ymax] ;
for r = 1:Nruns
    set(hs(r), 'XLim', xlims, 'YLim', ylims)
end

% Add zero-lines
if xlims(1)<0 && xlims(2)>0
    for r = 1:Nruns
        hold(hs(r),'on')
        plot(hs(r),[0 0],ylims,'--k')
        hold(hs(r),'off')
    end
end
if ylims(1)<0 && ylims(2)>0
    for r = 1:Nruns
        hold(hs(r),'on')
        plot(hs(r),xlims,[0 0],'--k')
        hold(hs(r),'off')
    end
end



end