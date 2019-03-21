%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLUM2LPJG figures: Multiple runs: %%%
%%% PLUM-style native %%%%%%%%%%%%%%%%%%%
%%% Ecosystem service relationships %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


thisVer = 'harm3' ;
% thisVer = 'harm3_constLU' ;
% thisVer = 'harm3_constClim' ;
% thisVer = 'harm3_constCO2' ;
% thisVer = 'harm3_constClimCO2' ;
% thisVer = 'harm3_S1R4.5_attr' ;
% thisVer = 'harm3_S3R6.0_attr' ;
% thisVer = 'harm3_S4R6.0_attr' ;
% thisVer = 'harm3_S5R8.5_attr' ;

do_adjYieldTech = true ; % Apply annual tech. change increase to yields?

unhCropFrac = 0 ; % Set to zero for previous behavior. v10 = 0.177

% ignored_crops = {'CC3G','CC4G'} ;
% ignored_crops = {'CC3G','CC4G','Miscanthus'} ;
ignored_crops = {'CC3G','CC4G','ExtraCrop'} ;

do_save = true ;
rebase = false ;
pngres = 150 ;

do_caps = -1 ;

       
%% Setup and import

run('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/MATLAB_work/plum2lpjg_figs_setup_import.m') ;


%% Import countries

countries_YX = flipud(dlmread('/Users/Shared/PLUM/crop_calib_data/countries/PLUM-specific/country_boundaries62892.noNeg99.extrapd.asc','',6,0)) ;
countries_YX(countries_YX<=0) = NaN ;

% Get country list; remove countries not in map
countries_key = readtable('/Users/Shared/PLUM/crop_calib_data/countries/PLUM-specific/country_boundaries_codes4.csv') ;
country_list_orig = unique(countries_key.Country) ;
countries_notinmap = setdiff(countries_key.numCode, unique(countries_YX(~isnan(countries_YX)))) ;
[~,IA] = intersect(countries_key.numCode,countries_notinmap,'stable') ;
countries_key(IA,:) = [] ;
country_list = unique(countries_key.Country) ;
Ncountries = length(country_list) ;

% Get region list and map
region_list = unique(countries_key.Region) ;
Nregions = length(region_list) ;
regions_YX = nan(size(countries_YX)) ;
for r = 1:Nregions
    thisRegion = region_list{r} ;
    theseCountries = countries_key.numCode(strcmp(countries_key.Region,thisRegion)) ;
    for c = 1:length(theseCountries)
        regions_YX(countries_YX==theseCountries(c)) = r ;
    end
end

% Get country GROUP list; harmonize with country list
countrygroups_key = readtable('/Users/Shared/PLUM/crop_calib_data/countries/PLUM-specific/country_groups.csv') ;
countrygroup_list = unique(countrygroups_key.PlumGroup) ;
countrygroups_key_Country = countrygroups_key.Country ;
for c = 1:length(countrygroups_key_Country)
    thisCountry = countrygroups_key_Country{c} ;
    if ~any(strcmp(country_list,thisCountry))
        if any(strcmp(country_list_orig,thisCountry))
            continue
        end
        warning('%s not found in country_list', thisCountry)
        if strcmp(thisCountry,'Kosovo')
            warning('Kosovo will be assumed part of Serbia.')
            countrygroups_key(strcmp(countrygroups_key.Country,thisCountry),:) = [] ;
            countrygroup_list(strcmp(countrygroup_list,thisCountry)) = [] ;
        else
            error('Aborting!')
        end
    end
end
for c = 1:Ncountries
    thisCountry = country_list{c} ;
    if ~any(strcmp(countrygroups_key.Country, thisCountry))
        warning('%s not found in countrygroups_key.Country', thisCountry)
        if any(strcmp({ ...
                'American Samoa','Guam','Northern Mariana Islands' ...
                }, thisCountry))
            thisSub = 'East Asia & Pacific_other' ;
        elseif strcmp(thisCountry, 'Andorra')
            thisSub = 'France Netherlands & Benlex' ;
        elseif any(strcmp({ ...
                'Cayman Islands', 'Turks and Caicos Islands', 'United States Virgin Islands', ...
                }, thisCountry))
            thisSub = 'Caribbean' ;
        elseif any(strcmp({ ...
                'Occupied Palestinian Territory','Syrian Arab Republic' ...
                }, thisCountry))
            thisSub = 'Middle East other' ;
        elseif any(strcmp({ ...
                'Timor-Leste' ...
                }, thisCountry))
            thisSub = 'Indonesia' ;
        elseif any(strcmp({ ...
                'Western Sahara' ...
                }, thisCountry))
            thisSub = 'North Africa other' ;
        else
            error('Aborting!')
        end
        warning('%s will be grouped with "%s".', thisCountry, thisSub)
        warning('off','MATLAB:table:RowsAddedExistingVars')
        countrygroups_key.Country(end+1) = {thisCountry} ;
        countrygroups_key.PlumGroup(end) = {thisSub} ;
        warning('on','MATLAB:table:RowsAddedExistingVars')
    end
end

% Get country group map
Ngroups = length(countrygroup_list) ;
countrygroups_YX = nan(size(countries_YX)) ;
for g = 1:Ngroups
    thisGroup = countrygroup_list{g} ;
    theseCountries_names = countrygroups_key.Country(strcmp(countrygroups_key.PlumGroup,thisGroup)) ;
    [~,~,IB] = intersect(theseCountries_names, countries_key.Country,'stable') ;
    theseCountries_codes = countries_key.numCode(IB) ;
    if length(theseCountries_codes)>1
        x=1;
    end
    for c = 1:length(theseCountries_codes)
        countrygroups_YX(countries_YX==theseCountries_codes(c)) = g ;
    end
end


%% By-country CARBON vs MONOTERPENES

es_relationship_scatter( ...
    maps_cpool_d9, 'Total', 'Total carbon (kgC m^{-2})', ...
    maps_amon_d9, 'Total', 'Monoterpene emissions (kgC m^{-2})', ...
    countries_YX, runList)


%% By-country CARBON vs ALBEDO

es_relationship_scatter( ...
    maps_cpool_d9, 'Total', 'Total carbon (kgC m^{-2})', ...
    maps_albedo_d9, 'January', 'January albedo', ...
    countries_YX, runList)


%% Per-region by-cell CARBON vs ALBEDO

ESscatter_plot_byCellPerRegion( ...
    maps_cpool_d9, 'Total', 'Total carbon (kgC m^{-2})', ...
    maps_albedo_d9, 'January', 'January albedo', ...
    regions_YX, region_list, runList)












