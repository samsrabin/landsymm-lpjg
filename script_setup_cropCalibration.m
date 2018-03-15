% Import country map and key
countries_YX = flipud(dlmread('country_boundaries62892.noNeg99.extrapd.asc','',6,0)) ;
countries_YX(countries_YX<=0) = NaN ;
countries_key = readtable('country_boundaries_codes4.csv') ;

if exist('filename_guess_gridlist','var')
    gl = lpjgu_matlab_readTable_then2map(filename_guess_gridlist) ;
    countries_YX(~gl.mask_YX) = NaN ;
end

% Merge South Sudan with Sudan if not needed
if yearN < 2011 && combine_sudans
    countries_YX(countries_YX==countries_key.numCode(strcmp(countries_key.Country,'South Sudan'))) = countries_key.numCode(strcmp(countries_key.Country,'Sudan')) ;
    countries_key(strcmp(countries_key.Country,'South Sudan'),:) = [] ;
end

% Merge Serbia and Montenegro, if doing so
if combine_serbmont
    countries_YX(countries_YX==countries_key.numCode(strcmp(countries_key.Country,'Montenegro'))) = countries_key.numCode(strcmp(countries_key.Country,'Serbia')) ;
    countries_key(strcmp(countries_key.Country,'Montenegro'),:) = [] ;
    countries_key.Country{strcmp(countries_key.Country,'Serbia')} = 'Serbia and Montenegro' ;
end

% Merge sub-Chinas, if doing so
if combine_subChinas_map
    countries_YX(countries_YX==158) = countries_key.numCode(strcmp(countries_key.Country,'China')) ;
    countries_YX(countries_YX==countries_key.numCode(strcmp(countries_key.Country,...
        'Hong Kong'))) = countries_key.numCode(strcmp(countries_key.Country,'China')) ;
    countries_YX(countries_YX==countries_key.numCode(strcmp(countries_key.Country,...
        'Macau'))) = countries_key.numCode(strcmp(countries_key.Country,'China')) ;
    % Remove from key
    countries_key(strcmp(countries_key.Country,'Hong Kong'),:) = [] ;
    countries_key(strcmp(countries_key.Country,'Macau'),:) = [] ;
end

Ncountries_key = length(countries_key.numCode) ;
countries_map = unique(countries_YX(~isnan(countries_YX))) ;
Ncountries_map = length(countries_map) ;
missing = ones(size(countries_YX)) ;
missing(isnan(countries_YX)) = 0 ;
missing_list = [] ;
for c = 1:Ncountries_map
    thisCountry_fromMap = countries_map(c) ;
    if ~any(countries_key.numCode==thisCountry_fromMap)
        missing(countries_YX==thisCountry_fromMap) = 2 ;
        missing_list = [missing_list thisCountry_fromMap] ;
        warning(['Code ' num2str(thisCountry_fromMap) ' not found in countries_key!'])
    end
end ; clear c
% figure ; pcolor(missing) ; shading flat ; axis equal tight
% disp(missing_list)
if isequal(missing_list,[158 238])
    warning('These are Taiwan and the Falkland Islands. Not sure if FAO data includes them in Argentina/GB/China.')
elseif missing_list == 238
    warning('These are the Falkland Islands. Not sure if FAO data includes them in Argentina/GB.')
else
    warning('Some unexpected countries are missing!')
end
listCountries_map_present = {} ;
for c = 1:Ncountries_key
    thisCountryName_fromKey = countries_key.Country{c} ;
    thisCountryCode_fromKey = countries_key.numCode(c) ;
    if any(countries_YX(~isnan(countries_YX))==thisCountryCode_fromKey)
        listCountries_map_present{end+1} = thisCountryName_fromKey ;
    end
end ; clear c
Ncountries = length(listCountries_map_present) ;

% % % % Classify countries as tropical and/or extratropical
% % % center_lats_YX = repmat((-89.75:0.5:89.75)',[1 720]) ;
% % % is_tropical_YX = (abs(center_lats_YX)<23.5) ;
% % % N_tropical = nan(Ncountries,1) ;
% % % N_xtratrop = nan(Ncountries,1) ;
% % % is_tropical = false(Ncountries,1) ;
% % % is_xtratrop = false(Ncountries,1) ;
% % % getCi = @(x) find(strcmp(countries_key.Country,x)) ;
% % % for C = 1:Ncountries
% % %     thisCountry = listCountries_map_present{C} ;
% % %     thisKey = countries_key.numCode(getCi(thisCountry)) ;
% % %     N_inThisCountry = length(find(countries_YX==thisKey)) ;
% % %     N_tropical(C) = length(find(countries_YX==thisKey & is_tropical_YX)) ;
% % %     N_xtratrop(C) = length(find(countries_YX==thisKey & ~is_tropical_YX)) ;
% % %     if N_tropical(C) > N_xtratrop(C)
% % %         is_tropical(C) = true ;
% % %     elseif N_tropical(C) < N_xtratrop(C)
% % %         is_xtratrop(C) = true ;
% % %     elseif N_tropical(C) == N_xtratrop(C)
% % %         is_tropical(C) = true ;
% % %         is_xtratrop(C) = true ;
% % %     end
% % % end ; clear C
[is_tropical,is_xtratrop] = classify_tropical_countries(...
    listCountries_map_present,countries_YX,countries_key) ;

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