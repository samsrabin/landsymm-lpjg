function es_relationship_scatter( ...
    mapstruct1, thisVar1, ...
    mapstruct2, thisVar2, ...
    countries_YX, runList)

% Get information about inputs
countries_list = unique(countries_YX(~isnan(countries_YX))) ;
Ncountries = length(countries_list) ;
Nruns = size(mapstruct1.maps_YXvyr,5) ;

% Get mean for variable in question
maps1_YXB = mean(mapstruct1.maps_YXvyB(:,:,strcmp(mapstruct1.varNames,thisVar1),:),4) ;
maps1_YXr = squeeze(mean(mapstruct1.maps_YXvyr(:,:,strcmp(mapstruct1.varNames,thisVar1),:,:),4)) ;
maps2_YXB = mean(mapstruct2.maps_YXvyB(:,:,strcmp(mapstruct2.varNames,thisVar2),:),4) ;
maps2_YXr = squeeze(mean(mapstruct2.maps_YXvyr(:,:,strcmp(mapstruct2.varNames,thisVar2),:,:),4)) ;

% Get means for each country
res1_cB = nan(Ncountries,1) ;
res2_cB = nan(Ncountries,1) ;
res1_cr = nan(Ncountries,Nruns) ;
res2_cr = nan(Ncountries,Nruns) ;
for c = 1:Ncountries
    thisCountry = countries_list(c) ;
    res1_cB(c) = nanmean(maps1_YXB(countries_YX==thisCountry)) ;
    res2_cB(c) = nanmean(maps2_YXB(countries_YX==thisCountry)) ;
    for r = 1:Nruns
        tmp_YX = maps1_YXr(:,:,r) ;        
        res1_cr(c,r) = nanmean(tmp_YX(countries_YX==thisCountry)) ;
        tmp_YX = maps2_YXr(:,:,r) ;        
        res2_cr(c,r) = nanmean(tmp_YX(countries_YX==thisCountry)) ;
    end
end

% Trim countries that are all NaN
bad_countries = find(isnan(res1_cB) & all(isnan(res1_cr),2) & isnan(res2_cB) & all(isnan(res2_cr),2)) ;
res1_cB(bad_countries) = [] ;
res2_cB(bad_countries) = [] ;
res1_cr(bad_countries,:) = [] ;
res2_cr(bad_countries,:) = [] ;
if any(isnan(res1_cB))
    error('NaN remaining in res1_cB')
elseif any(any(isnan(res1_cr)))
    error('NaN remaining in res1_cr')
elseif any(isnan(res2_cB))
    error('NaN remaining in res2_cB')
elseif any(any(isnan(res2_cr)))
    error('NaN remaining in res2_cr')
end

% Get diffs
diff1_cr = res1_cr - repmat(res1_cB,[1 Nruns]) ;
diff2_cr = res2_cr - repmat(res2_cB,[1 Nruns]) ;

% Make figure
figure('Color','w','Position',figurePos) ;
spacing = 0.05 ;
for r = 1:Nruns
    subplot_tight(2,2,r,spacing)
    plot(diff1_cr(:,r),diff2_cr(:,r),'.') ;
    title(runList{r})
end




end