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

       
%% Setup; import LPJ-GUESS results

do_adjYieldTech = true ; % Apply annual tech. change increase to yields?

unhCropFrac = 0 ; % Set to zero for previous behavior. v10 = 0.177

% ignored_crops = {'CC3G','CC4G'} ;
% ignored_crops = {'CC3G','CC4G','Miscanthus'} ;
ignored_crops = {'CC3G','CC4G','ExtraCrop'} ;

do_save = true ;
rebase = false ;
pngres = 150 ;

do_caps = -1 ;

% Import LPJ-GUESS results
run('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/MATLAB_work/plum2lpjg_figs_setup_import.m') ;


%% Import PLUM outputs

if exist('runDirs_plum', 'var')
    fprintf('Reading PLUM outputs: ')
    countryDemand = struct ; 
    domestic = struct ; 
    for r = 1:Nruns
        thisDir = runDirs_plum{r} ;
        fprintf('%d... ', r)
        table_in_countryDemand = readtable(sprintf('%s/countryDemand.txt', thisDir)) ;
        table_in_countryDemand = sortrows(table_in_countryDemand, ...
            {'Year','Commodity','Country'}) ; % Reverse of dimension order
        
        table_in_domestic = readtable(sprintf('%s/domestic.txt', thisDir)) ;
        table_in_domestic.Import_price = remove_nul(table_in_domestic.Import_price) ;
        table_in_domestic.Export_price = remove_nul(table_in_domestic.Export_price) ;
        table_in_domestic.Net_imports = remove_nul(table_in_domestic.Net_imports) ;
        table_in_domestic.Net_import_cost = remove_nul(table_in_domestic.Net_import_cost) ;
        if r==1
            countryDemand.years = unique(table_in_countryDemand.Year) ;
            countryDemand.countrygroups = unique(table_in_countryDemand.Country) ;
            countryDemand.commodities = unique(table_in_countryDemand.Commodity) ;
            countryDemand.data_umyr = nan( ...
                length(countryDemand.countrygroups), ...
                length(countryDemand.commodities), ...
                length(countryDemand.years), ...
                Nruns) ;
            size_countryDemand_data = size(countryDemand.data_umyr) ;
            domestic.header = table_in_domestic.Properties.VariableNames ;
            domestic.years = unique(table_in_domestic.Year) ;
            domestic.countrygroups = unique(table_in_domestic.Country) ;
            domestic.commodities = unique(table_in_domestic.Crop) ;
            [domestic.varNames, IA] = setdiff(domestic.header, {'Year','Country','Crop'},'stable') ;
            domestic.data_umvyr = nan( ...
                length(domestic.countrygroups), ...
                length(domestic.commodities), ...
                length(domestic.varNames), ...
                length(domestic.years), ...
                Nruns) ;
            size_domestic_data = size(domestic.data_umvyr) ;
        else
            if ~isequal(countryDemand.years, unique(table_in_countryDemand.Year))
                error('countryDemand: YEAR column does not match')
            elseif ~isequal(countryDemand.countrygroups, unique(table_in_countryDemand.Country))
                error('countryDemand: COUNTRYGROUP column does not match')
            elseif ~isequal(countryDemand.commodities, unique(table_in_countryDemand.Commodity))
                error('countryDemand: COMMODITY column does not match')
            end
            if ~isequal(domestic.years, unique(table_in_domestic.Year))
                error('domestic: YEAR column does not match')
            elseif ~isequal(domestic.countrygroups, unique(table_in_domestic.Country))
                error('domestic: COUNTRYGROUP column does not match')
            elseif ~isequal(domestic.commodities, unique(table_in_domestic.Crop))
                error('domestic: COMMODITY column does not match')
            elseif ~isequal(domestic.header, table_in_domestic.Properties.VariableNames)
                error('domestic: HEADER does not match')
            end
        end
        countryDemand.data_umyr(:,:,:,r) = reshape(table_in_countryDemand.Demand,size_countryDemand_data(1:end-1)) ;
        
        table_in_domestic = stack(table_in_domestic, ...
            domestic.varNames, ...
            'NewDataVariableName','Value',...
            'IndexVariableName','Var') ;
        table_in_domestic = sortrows(table_in_domestic, ...
            {'Year','Var','Crop','Country'}) ; % Reverse of dimension order
        domestic.data_umvyr(:,:,:,:,r) = reshape(table_in_domestic.Value,size_domestic_data(1:end-1)) ;
        clear table_in_countryDemand table_in_domestic
    end
    clear IA
    fprintf('Done.\n')
else
    warning('Not importing any PLUM results.')
end


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

ESscatter_diffs_byRegion( ...
    maps_cpool_d9, 'Total', 'Total carbon (kgC m^{-2})', ...
    maps_amon_d9, 'Total', 'Monoterpene emissions (kgC m^{-2})', ...
    countries_YX, runList)


%% By-country CARBON vs ALBEDO

ESscatter_diffs_byRegion( ...
    maps_cpool_d9, 'Total', 'Total carbon (kgC m^{-2})', ...
    maps_albedo_d9, 'January', 'January albedo', ...
    countries_YX, runList)


%% Per-region by-cell CARBON vs ALBEDO

ESscatter_diffs_byCellPerRegion( ...
    maps_cpool_d9, 'Total', 'Total carbon (kgC m^{-2})', ...
    maps_albedo_d9, 'January', 'January albedo', ...
    regions_YX, region_list, runList)



%% By-countrygroup RUMINANT PRODUCTION vs CARBON

tmp_rumprod_2090s_YXyr = make_map_from_plum_out( ...
    domestic, 'ruminants', 'Production', 2091:2100, ...
    countrygroups_YX, countrygroup_list) ;

tmp_carbon_2090s_YXyr = squeeze(maps_cpool_d9.maps_YXvyr(:,:,strcmp(maps_cpool_d9.varNames,'Total'),:,:)) ;

warning('Normalize by area?')
ESscatter_byRegion( ...
    tmp_rumprod_2090s_YXyr, 'Ruminant production (??)', ...
    tmp_carbon_2090s_YXyr, 'Total carbon (kgC m^{-2})', ...
    countries_YX, runList)

clear tmp*

