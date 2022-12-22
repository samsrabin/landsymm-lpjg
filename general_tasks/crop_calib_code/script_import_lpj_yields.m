script_import_lpj_yields_noCCy


%% Make _Ccy datasets
disp('Making _Ccy datasets... ')
verbose = false ;
[croparea_lpj_Ccy,total_lpj_Ccy,yield_lpj_Ccy] = ...
    YXcy_to_Ccy(total_lpj_YXcy_comb,croparea_lpj_YXcy_comb,...
                countries_YX,countries_key,listCountries_map_present,...
                verbose) ;
disp('Done.')