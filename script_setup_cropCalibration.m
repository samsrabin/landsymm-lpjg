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

% Is this a GGCMI run?
is_ggcmi = contains(version_name,'ggcmi') | contains(version_name,'GGCMI') ;

% Import land area (km2)
xres = 360/size(countries_YX,2) ;
yres = 180/size(countries_YX,1) ;
if ~is_ggcmi
    if calib_ver==17 % Put your calib_ver here if you want it to use MCD12C1-derived land area.
        if xres ~= yres
            error('To use this land area map, xres must == yres.')
        end
        landarea_file = ['land_area_km2_fromMCD12C1_' num2str(xres) 'deg.mat'] ;
        load(landarea_file) ;
    elseif calib_ver <= 16
        gcel_area_YXqd = transpose(ncread('staticData_quarterdeg.nc','carea')) ;
        land_frac_YXqd = 1 - flipud(transpose(ncread('staticData_quarterdeg.nc','icwtr'))) ;
        land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
        land_area_YX = aggregate_land_area(land_area_YXqd,xres,yres) ;
    elseif calib_ver == 18
        % Use gridcell area instead of land area
        land_area_YXqd = transpose(ncread('staticData_quarterdeg.nc','carea')) ;
        land_area_YX = aggregate_land_area(land_area_YXqd,xres,yres) ;
    else
        error(['calib_ver ' num2str(calib_ver) ' not recognized in "Import land area"'])
    end
    %%%%% Convert from km2 to ha
    land_area_YX = land_area_YX * 100 ;
end