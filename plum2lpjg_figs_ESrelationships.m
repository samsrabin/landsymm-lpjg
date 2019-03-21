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


%% Countries
countries_YX = flipud(dlmread('/Users/Shared/PLUM/crop_calib_data/countries/PLUM-specific/country_boundaries62892.noNeg99.extrapd.asc','',6,0)) ;
countries_YX(countries_YX<=0) = NaN ;
countries_key = readtable('/Users/Shared/PLUM/crop_calib_data/countries/PLUM-specific/country_boundaries_codes4.csv') ;
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


%% By-country CARBON vs ISOPRENE



es_relationship_scatter( ...
    maps_cpool_d9, 'Total', ...
    maps_amon_d9, 'Total', ...
    countries_YX, runList)













