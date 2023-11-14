% Make sure directory names end with filesep
if ~strcmp(dir_data(end), filesep)
    dir_data = [dir_data filesep] ;
end
if exist('dir_outfigs', 'var') && ~strcmp(dir_outfigs(end), filesep)
    dir_outfigs = [dir_outfigs filesep] ;
end

do_temp_vs_trop = false ;

% Years for calibration
listYears_fao = year1:yearN ;
Nyears_fao = length(listYears_fao) ;

% Add data files to path (just for this session)
addpath(genpath(dir_data))


% Import countries, if needed. Also set up X and Y resolution.
if need_countries
    [strip_fao_nans, fix_cotedivoire, combine_sudans,...
        combine_subChinas, combine_subChinas_map, combine_serbmont] ...
        = get_FAOread_options(calib_ver) ;
    % Import country map and key
    if any(strcmp(filename_countriesMap, ...
            {fullfile(landsymm_lpjg_path(), 'data', 'geodata', 'country_boundaries', 'country_boundaries62892.noNeg99.extrapd.asc'), ...
            fullfile(landsymm_lpjg_path(), 'data', 'geodata', 'country_boundaries', 'country_boundaries_f09_g17.asc')}))
        PLUM_countries = true ;
    else
        PLUM_countries = false ;
    end
    if PLUM_countries
        if contains(filename_countriesMap, 'f09_g17')
            NheaderLines = 7 ;
        else
            NheaderLines = 6 ;
        end
        countries_YX = flipud(dlmread(filename_countriesMap,'',NheaderLines,0)) ;
        countries_YX(countries_YX<=0) = NaN ;
        countries_key = readtable(fullfile(landsymm_lpjg_path(), 'data', 'geodata', 'country_boundaries', 'country_boundaries_codes4.csv') ;
    else
        countries_YX = flipud(imread(filename_countriesMap)) ;
        countries_YX(countries_YX==0) = NaN ; % Water
        countries_key_tmp = readtable('ne_10m_admin_0_countries_ssrIDs.csv') ;
        countries_key = table(countries_key_tmp.ADMIN, countries_key_tmp.ADM_ID_SSR) ;
        countries_key.Properties.VariableNames = {'Country','numCode'} ;
    end
    xres = 360/size(countries_YX,2) ;
    yres = 180/size(countries_YX,1) ;
else
    xres = 0.5 ;
    yres = 0.5 ;
end

% Import land area (km2)
if contains(filename_countriesMap, 'f09_g17')
    landarea_file = '/Users/Shared/CESM_work/CropEvalData_ssr/landuse.timeseries_0.9x1.25_hist_78pfts_CMIP6_simyr1850-2015_c170824.nc4' ;
    land_area_YX = transpose(ncread(landarea_file, 'AREA')) ;
    land_area_YX = land_area_YX(2:end-1,:) ;
elseif calib_ver==17 % Put your calib_ver here if you want it to use MCD12C1-derived land area.
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
elseif calib_ver >= 18 && calib_ver <= 24
    % Use gridcell area instead of land area
    land_area_YXqd = transpose(ncread('staticData_quarterdeg.nc','carea')) ;
    land_area_YX = aggregate_land_area(land_area_YXqd,xres,yres) ;
else
    error(['calib_ver ' num2str(calib_ver) ' not recognized in "Import land area"'])
end
%%%%% Convert from km2 to ha
land_area_YX = land_area_YX * 100 ;


% Is this a GGCMI run?
if ~exist('is_ggcmi', 'var')
    is_ggcmi = contains(verName_calib,'ggcmi') ...
        | contains(verName_calib,'GGCMI') ...
        | contains(verName_calib,'emu') ;
end
if ~exist('indiv_years', 'var')
    indiv_years = true ;
end
