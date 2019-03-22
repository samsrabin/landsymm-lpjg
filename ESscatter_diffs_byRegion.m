function ESscatter_diffs_byRegion( ...
    mapstructX, thisVarX, labelX, ...
    mapstructY, thisVarY, labelY, ...
    countries_YX, runList)

%%%%% Options
fontSize = 14 ;
markerSize = 15 ;
spacing = [0.1 0.1] ;   % [vertical horizontal]

% Check inputs
ESscatter_check_inputs(...
    mapstructX, thisVarX, ...
    mapstructY, thisVarY, ...
    countries_YX) ;

% Get information about inputs
countries_codes = unique(countries_YX(~isnan(countries_YX))) ;
Ncountries = length(countries_codes) ;
Nruns = size(mapstructX.maps_YXvyr,5) ;

% Get mean for variable in question
[mapsX_YXB, mapsX_YXr, mapsY_YXB, mapsY_YXr]= ...
    ESscatter_getMeanOverPeriod(...
        mapstructX, thisVarX, ...
        mapstructY, thisVarY) ;

% Get means for each country
resX_cB = nan(Ncountries,1) ;
resY_cB = nan(Ncountries,1) ;
resX_cr = nan(Ncountries,Nruns) ;
resY_cr = nan(Ncountries,Nruns) ;
for c = 1:Ncountries
    thisCountry = countries_codes(c) ;
    resX_cB(c) = nanmean(mapsX_YXB(countries_YX==thisCountry)) ;
    resY_cB(c) = nanmean(mapsY_YXB(countries_YX==thisCountry)) ;
    for r = 1:Nruns
        tmp_YX = mapsX_YXr(:,:,r) ;        
        resX_cr(c,r) = nanmean(tmp_YX(countries_YX==thisCountry)) ;
        tmp_YX = mapsY_YXr(:,:,r) ;        
        resY_cr(c,r) = nanmean(tmp_YX(countries_YX==thisCountry)) ;
    end
end

% Trim countries that are all NaN
bad_countries = find(isnan(resX_cB) & all(isnan(resX_cr),2) & isnan(resY_cB) & all(isnan(resY_cr),2)) ;
resX_cB(bad_countries) = [] ;
resY_cB(bad_countries) = [] ;
resX_cr(bad_countries,:) = [] ;
resY_cr(bad_countries,:) = [] ;
if any(isnan(resX_cB))
    error('NaN remaining in resX_cB')
elseif any(any(isnan(resX_cr)))
    error('NaN remaining in resX_cr')
elseif any(isnan(resY_cB))
    error('NaN remaining in resY_cB')
elseif any(any(isnan(resY_cr)))
    error('NaN remaining in resY_cr')
end

% Get diffs (future minus historical)
diffX_cr = resX_cr - repmat(resX_cB,[1 Nruns]) ;
diffY_cr = resY_cr - repmat(resY_cB,[1 Nruns]) ;

% Make figure
figure('Color','w','Position',figurePos) ;
hs = [] ; % Axis handles
xmin = NaN ;
xmax = NaN ;
ymin = NaN ;
ymax = NaN ;
for r = 1:Nruns
    hs(r) = subplot_tight(2,2,r,spacing) ;
    plot(diffX_cr(:,r),diffY_cr(:,r),'.', ...
        'MarkerSize',markerSize)
    title(runList{r})
    xlabel(['? ' labelX])
    ylabel(['? ' labelY])
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