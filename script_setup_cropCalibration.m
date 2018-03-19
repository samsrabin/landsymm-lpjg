% Import country map and key
if strcmp(filename_countriesMap,'country_boundaries62892.noNeg99.extrapd.asc')
    PLUM_countries = true ;
else
    PLUM_countries = false ;
end
if PLUM_countries
    countries_YX = flipud(dlmread(filename_countriesMap,'',6,0)) ;
    countries_YX(countries_YX<=0) = NaN ;
    countries_key = readtable('country_boundaries_codes4.csv') ;
else
    countries_YX = flipud(imread(filename_countriesMap)) ;
    countries_YX(countries_YX==0) = NaN ; % Water
    countries_key_tmp = readtable('ne_10m_admin_0_countries_ssrIDs.csv') ;
    countries_key = table(countries_key_tmp.ADMIN, countries_key_tmp.ADM_ID_SSR) ;
    countries_key.Properties.VariableNames = {'Country','numCode'} ;
end

% Import land area (km2)
gcel_area_YXqd = transpose(ncread('staticData_quarterdeg.nc','carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread('staticData_quarterdeg.nc','icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear tmp
%%%%% Convert to ha
land_area_YX = land_area_YX * 100 ;