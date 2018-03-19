% Mask countries map to LPJ-GUESS gridlist
gl.mask_YX = ~isnan(mean(mean(croparea_lpj_YXcy_comb,4),3)) ;
if ~isequal(size(gl.mask_YX),size(countries_YX))
    error('~isequal(size(gl.mask_YX),size(countries_YX))')
end
if false % Put calib_ver==XX here if you want to use a mask
    if exist('filename_guess_gridlist','var')
        gl = lpjgu_matlab_readTable_then2map(filename_guess_gridlist) ;
    end
    countries_YX(~gl.mask_YX) = NaN ;
elseif calib_ver<=16
    if exist('filename_guess_gridlist','var')
        warning(['filename_guess_gridlist ignored for calib_ver ' num2str(calib_ver) '.'])
    end
%     countries_YX(~gl.mask_YX) = NaN ;
else
    error(['calib_ver ' num2str(calib_ver) ' not recognized in "Mask countries map to LPJ-GUESS gridlist"'])
end

% Where does LPJ-GUESS have gridcells but countries_YX does not?
inGlNotCtries_YX = gl.mask_YX & isnan(countries_YX) ;
if false % Put calib_ver==XX here if you want to try and fix this
    % First, try to fill with values from "all touched" gdal_rasterize. We do
    % not use this at first because it tends to screw up whenever two
    % countries touch each other. However, where land borders ocean--i.e.,
    % where a lot of cells in this situation are--it should be pretty
    % unambiguously accurate.
    if any(inGlNotCtries_YX(:)) && ~PLUM_countries
        [PATHSTR,NAME,EXT] = fileparts(filename_countriesMap) ;
        filename_countriesMap_at = [PATHSTR NAME '.at' EXT ] ;
        countries_YX_at = flipud(imread(filename_countriesMap_at)) ;
        countries_YX_at(countries_YX_at==0) = NaN ; % Water
        countries_YX(inGlNotCtries_YX) = countries_YX_at(inGlNotCtries_YX) ;
    end
    % Warn if any such cells remain, and mask LPJ-GUESS outputs.
    inGlNotCtries_YX = gl.mask_YX & isnan(countries_YX) ;
    if any(inGlNotCtries_YX(:))
        warning([num2str(length(find(inGlNotCtries_YX))) ' cells in LPJ-GUESS output but not countries map! These will be masked.'])
        inGlNotCtries_YXcy = repmat(inGlNotCtries_YX,[1 1 size(croparea_lpj_YXcy_comb,3) size(croparea_lpj_YXcy_comb,4)]) ;
        croparea_lpj_YXcy_comb(inGlNotCtries_YXcy) = NaN ;
        cropfrac_lpj.maps_YXvy(inGlNotCtries_YXcy) = NaN ;
        cropfrac_lpj_comb.maps_YXvy(inGlNotCtries_YXcy) = NaN ;
        landuse_lpj.maps_YXvy(inGlNotCtries_YXcy) = NaN ;
        total_lpj_YXcy_comb(inGlNotCtries_YXcy) = NaN ;
        yield_lpj.maps_YXvy(inGlNotCtries_YXcy) = NaN ;
        yield_lpj_comb.maps_YXvy(inGlNotCtries_YXcy) = NaN ;
        gl.mask_YX = ~isnan(mean(mean(croparea_lpj_YXcy_comb,4),3));
        if any(any(gl.mask_YX & isnan(countries_YX)))
            error('How did that masking not work???')
        end
    end
elseif calib_ver<=16
    if any(inGlNotCtries_YX(:))
        warning([num2str(length(find(inGlNotCtries_YX))) ' cells in LPJ-GUESS output but not countries map! These will be ignored in calibration. To view: figure;pcolor(inGlNotCtries_YX);shading flat'])
    end
else
    error(['calib_ver ' num2str(calib_ver) ' not recognized in "Where does LPJ-GUESS have gridcells but countries_YX does not"'])
end


% Adjust map
if ~PLUM_countries
    countries_YX(countries_YX==countries_key.numCode(strcmp(countries_key.Country,'Antarctica'))) = NaN ;
    countries_key(strcmp(countries_key.Country,'Antarctica'),:) = [] ;
    countries_YX(countries_YX==countries_key.numCode(strcmp(countries_key.Country,'French Southern and Antarctic Lands'))) = NaN ;
    countries_key(strcmp(countries_key.Country,'French Southern and Antarctic Lands'),:) = [] ;
end

% Adjust names in countries table, if not using PLUM countries
if ~PLUM_countries
    country_name_adj = {...
        'North Korea',        'Democratic People''s Republic of Korea' ;
        'South Korea',        'Republic of Korea' ;
        'Bolivia',            'Bolivia (Plurinational State of)' ;
        'Republic of Congo',  'Congo' ;
        'Ivory Coast',        'Cote d''Ivoire' ;
        'Guinea Bissau',      'Guinea-Bissau' ;
        'Iran',               'Iran (Islamic Republic of)' ;
        'Laos',               'Lao People''s Democratic Republic' ;
        'Palestine',          'Occupied Palestinian Territory' ;
        'Vietnam',            'Viet Nam' ;
        'Republic of Serbia', 'Serbia' ;
        'Syria',              'Syrian Arab Republic' ;
        'Macedonia',          'The former Yugoslav Republic of Macedonia' ;
        'East Timor',         'Timor-Leste' ;
        'Venezuela',          'Venezuela (Bolivarian Republic of)' ;
        'Brunei',             'Brunei Darussalam' ;
        'Moldova',            'Republic of Moldova' ;
        'Russia',             'Russian Federation' ;
        'The Bahamas',        'Bahamas' ; ...
        'Cabo Verde',         'Cabo Verde' ; ...
        } ;
    for c = 1:size(country_name_adj,1)
        countries_key.Country(strcmp(countries_key.Country,country_name_adj{c,1})) ...
            = {country_name_adj{c,2}} ;
    end
    clear country_name_adj
end ; clear c

% Merge names in countries table, if not using PLUM countries
if ~PLUM_countries
    country_name_adj = {...
        'Aland',                            'Finland' ;
        'Baykonur Cosmodrome',              'Kazakhstan' ;
        'Cyprus No Mans Area',              'Cyprus' ;
        'Northern Cyprus',                  'Cyprus' ;
        'Akrotiri Sovereign Base Area',     'Cyprus' ;
        'Somaliland',                       'Somalia' ;
        } ;
    for c = 1:size(country_name_adj,1)
        thisCode_1 = countries_key.numCode(strcmp(countries_key.Country,country_name_adj{c,1})) ;
        thisCode_2 = countries_key.numCode(strcmp(countries_key.Country,country_name_adj{c,2})) ;
        if isempty(thisCode_1)
            error('isempty(thisCode_1)')
        elseif isempty(thisCode_2)
            error('isempty(thisCode_2)')
        end
        countries_YX(countries_YX==thisCode_1) = thisCode_2 ;
        countries_key(strcmp(countries_key.Country,country_name_adj{c,1}),:) = [] ;
    end
    clear country_name_adj
end ; clear c

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
    % Taiwan
    if PLUM_countries
        taiwan_code = 158 ;
    else
        taiwan_code = countries_key.numCode(strcmp(countries_key.Country,'Taiwan')) ;
    end
    countries_YX(countries_YX==taiwan_code) = countries_key.numCode(strcmp(countries_key.Country,'China')) ;
    countries_key(strcmp(countries_key.Country,'Taiwan'),:) = [] ;
    clear taiwan_code
    
    % Hong Kong
    hk_code = countries_key.numCode(strcmp(countries_key.Country,'Hong Kong')) ;
    if ~isempty(hk_code)
        countries_YX(countries_YX==hk_code) = ...
            countries_key.numCode(strcmp(countries_key.Country,'China')) ;
        countries_key(strcmp(countries_key.Country,'Hong Kong'),:) = [] ;
    end
    clear hk_code
    
    % Macau
    macau_code = countries_key.numCode(strcmp(countries_key.Country,'Macau')) ;
    if ~isempty(macau_code)
        countries_YX(countries_YX==macau_code) = ...
            countries_key.numCode(strcmp(countries_key.Country,'China')) ;
        countries_key(strcmp(countries_key.Country,'Macau'),:) = [] ;
    end
    clear macau_code
end

% Check that every mapped code is found in the key
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
if PLUM_countries
    if isequal(missing_list,[158 238])
        warning('These are Taiwan and the Falkland Islands. Not sure if FAO data includes them in Argentina/GB/China.')
    elseif missing_list == 238
        warning('These are the Falkland Islands. Not sure if FAO data includes them in Argentina/GB.')
    else
        error('Some unexpected countries are missing!')
    end
else
    if ~isempty(missing_list)
        error('Some unexpected countries are missing!')
    end
end

% Get list of countries that are present in map
listCountries_map_present = {} ;
for c = 1:Ncountries_key
    thisCountryName_fromKey = countries_key.Country{c} ;
    thisCountryCode_fromKey = countries_key.numCode(c) ;
    if any(countries_YX(~isnan(countries_YX))==thisCountryCode_fromKey)
        listCountries_map_present{end+1} = thisCountryName_fromKey ;
    end
end ; clear c
Ncountries = length(listCountries_map_present) ;

% Classify as tropical or not
[is_tropical,is_xtratrop] = classify_tropical_countries(...
    listCountries_map_present,countries_YX,countries_key) ;






