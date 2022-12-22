function [is_tropical,is_xtratrop] = classify_tropical_countries(...
    listCountries_map_present,countries_YX,countries_key)

% Get info
Ncountries = length(listCountries_map_present) ;
lat_size = size(countries_YX,1) ;
lon_size = size(countries_YX,2) ;
lat_res = 180/lat_size ;

% Setup
center_lats_YX = repmat((-90+lat_res/2:lat_res:90-lat_res/2)',[1 lon_size]) ;
is_tropical_YX = (abs(center_lats_YX)<23.5) ;
N_tropical = nan(Ncountries,1) ;
N_xtratrop = nan(Ncountries,1) ;
is_tropical = false(Ncountries,1) ;
is_xtratrop = false(Ncountries,1) ;
getCi = @(x) find(strcmp(countries_key.Country,x)) ;
for C = 1:Ncountries
    thisCountry = listCountries_map_present{C} ;
    thisKey = countries_key.numCode(getCi(thisCountry)) ;
    N_tropical(C) = length(find(countries_YX==thisKey & is_tropical_YX)) ;
    N_xtratrop(C) = length(find(countries_YX==thisKey & ~is_tropical_YX)) ;
    if N_tropical(C) > N_xtratrop(C)
        is_tropical(C) = true ;
    elseif N_tropical(C) < N_xtratrop(C)
        is_xtratrop(C) = true ;
    elseif N_tropical(C) == N_xtratrop(C)
        is_tropical(C) = true ;
        is_xtratrop(C) = true ;
    end
end ; clear C

end