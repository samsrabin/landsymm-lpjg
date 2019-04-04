%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLUM2LPJG figures: Multiple runs: %%%
%%% PLUM-style native %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % thisVer = '20180424agmip7' ;
% % thisVer = '20180424agmip7_asPLUMout2011-2015' ;
% thisVer = 'v4s1_v20180426' ;
% thisVer = 'v4s1_v20180426_asPLUMout2011' ;
% thisVer = 'v4s1_v20180426_asLUH2_2010' ;
% thisVer = 'v6s1_v20180703' ;
% thisVer = 'v10s1_v20180801' ;
% thisVer = 'v10s1_v20180801' ;
% thisVer = 'harm2' ;
% thisVer = 'harm2_constLU' ;
% thisVer = 'harm2.1' ;
% thisVer = 'harm2.1_constClimCO2' ;
% thisVer = 'harm2.1_constLU' ;
% thisVer = 'harm2.1_S1R4.5_attr' ;
% thisVer = 'harm2.1_S3R6.0_attr' ;
% thisVer = 'harm2.1_S4R6.0_attr' ;
% thisVer = 'harm2.1_S5R8.5_attr' ;
% thisVer = 'harm2.2' ;
% thisVer = 'harm2.3' ;
% thisVer = 'harm2.3_constClimCO2' ;
% thisVer = 'harm2.3_constLU' ;
% thisVer = 'harm2.3_S1R4.5_attr' ;
% thisVer = 'harm2.3_S3R6.0_attr' ;
% thisVer = 'harm2.3_S4R6.0_attr' ;
% thisVer = 'harm2.3_S5R8.5_attr' ;
% thisVer = 'harm2.3_constClim' ;
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


%% Save outputs

if false
    disp('Saving outputs...')
    save(['/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/fig_script_outputs/' thisVer '.mat'], ...
        'agri_*', 'bare_*', 'bl_*', 'crop*', 'firstdec_tmp', 'lastdec_tmp', ...
        'land_area_YX', 'gcel_area_YX', 'maps_*', 'nanmask', 'ntrl_*', 'past*', ...
        'ts_*_bl','ts_*_yr')
    disp('Done.')
end


%% Figure-ified table after Krause et al. (2017) Table 2

fontSize = 14 ;

years_endh = 2001:2010 ;
years_endf = 2091:2100 ;

% Name, code, conversion factor, formatSpec mean, formatSpec SEM, units
where2sep = [4.5 9.5 11.5] ;
rowInfo = { ...
           % LU variables
%            'Area: Agriculture', 'LUarea_crop+LUarea_past', 1e6*1e-6, '%0.1f', '%0.1f', 'Mkm^2' ;
           'Area: Non-agri.', 'LUarea_ntrl', 1e-6*1e-6, '%0.1f', '%0.1f', 'Mkm^2' ;
           'Fertilizer', 'nflux_fert', -1e-9, '%0.1f', '%0.1f', 'TgN' ;
           'Irrigation', 'irrig', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f', '1000 km^3' ;
           'Crop prod.', 'kcal', cf_kcalEcal, '%0.1f', '%0.1f', 'Ecal' ;
           % "Higher is better"
           'Veg. C', 'cpool_VegC', cf_kg2Pg, '%d', '%d', 'GtC' ;
%            'Soil/litter C', 'cpool_LitterSoilC', cf_kg2Pg, '%d', '%d', 'GtC' ;
           'Total C', 'cpool_Total', cf_kg2Pg, '%d', '%d', 'GtC' ;
           'Jan. albedo', 'albedo1', 1, '%0.3f', '%0.3f', '' ;
           'Jul. albedo', 'albedo7', 1, '%0.3f', '%0.3f', '' ;
           'Jan. albedo, borfor+tundra', 'albedo1_borfor+albedo1_tundra', 1, '%0.3f', '%0.3f', '' ;
           % "Lower is better"
           'N loss', 'nloss', cf_kg2Tg, '%0.1f', '%0.1f', 'TgN' ;
           'BVOC emis.', 'aiso+amon', cf_kg2Tg, '%0.1f', '%0.1f', 'TgC' ;
           % "Neutral"
           'ET', 'aevapaaet', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f', '1000 km^3' ;
           'Runoff', 'tot_runoff', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f', '1000 km^3' ;
           } ;

Nvars = size(rowInfo,1) ;
mean_endh_v = nan(Nvars,1) ;
mean_endf_vr = nan(Nvars,Nruns) ;
sem_endh_v = nan(Nvars,1) ;
sem_endf_vr = nan(Nvars,Nruns) ;
for c = 1:Nvars
    
    % Get values
    thisVar = rowInfo{c,2} ;
    thisConv = rowInfo{c,3} ;
    if contains(thisVar,'+')
        theseVars = strsplit(thisVar,'+') ;
        thisVar = theseVars{1} ;
        ts_thisVar_bl = eval(['ts_' thisVar '_bl']) ;
        ts_thisVar_yr = eval(['ts_' thisVar '_yr']) ;
        for v = 2:length(theseVars)
            thisVar = theseVars{v} ;
            ts_thisVar_bl = ts_thisVar_bl + eval(['ts_' thisVar '_bl']) ;
            ts_thisVar_yr = ts_thisVar_yr + eval(['ts_' thisVar '_yr']) ;
        end
        thisVar = 'thisVar' ;
    end
    mean_endh_v(c) = thisConv*eval(['mean(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
    sem_endh_v(c) = thisConv*eval(['std(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
    mean_endf_vr(c,:) = thisConv*eval(['mean(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
    sem_endf_vr(c,:) = thisConv*eval(['std(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
    clear ts_thisVar_yr % Does nothing if does not exist
end

% Calculate % difference
mean_endh_vr = repmat(mean_endh_v,[1 Nruns]) ;
sem_endh_vr = repmat(sem_endh_v,[1 Nruns]) ;
mean_diff_vr = mean_endf_vr - mean_endh_vr ;
mean_diffPct_vr = 100 * (mean_diff_vr ./ mean_endh_vr) ;
sem_diff_vr = sqrt(sem_endh_vr.^2 + sem_endf_vr.^2) ;
% EQUIVALENT TO BELOW % sem_diffPct_vr = 100 * ((mean_endf_vr+sem_diff_vr-mean_endh_vr)./mean_endh_vr - (mean_endf_vr-mean_endh_vr)./mean_endh_vr) ;
sem_diffPct_vr = 100 * (sem_diff_vr ./ mean_endh_vr) ;

% Make figure
figure('Color','w','Position',figurePos) ; 
h = bar(mean_diffPct_vr, 'grouped') ;
set(gca, 'XTickLabel', rowInfo(:,1)) ;
xtickangle(45) ;

%%%%%%%%%%%%%%%%%%%%%%
%%% Add error bars %%%
%%%%%%%%%%%%%%%%%%%%%%
% https://www.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab

% Finding the number of groups and the number of bars in each group
ngroups = size(mean_diffPct_vr, 1);
nbars = size(mean_diffPct_vr, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
hold on
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars) ;
    errorbar(x, mean_diffPct_vr(:,i), sem_diffPct_vr(:,i), 'k', 'linestyle', 'none', 'linewidth', 0.5);
end
hold off


%%%%%%%%%%%%%%%%%%%%%%
%%% Add separators %%%
%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(where2sep)
    hold on
    ylims = get(gca,'YLim') ;
    for i = 1:length(where2sep)
%         x = floor(where2sep(i)) + groupwidth/2 ;
        x = where2sep(i) ;
        line([x x], ylims, 'Color', 0.8*(ones(3,1)), 'LineWidth', 10) ;
    end
    hold off
end

%%%%%%%%%%%%%%
%%% Finish %%%
%%%%%%%%%%%%%%

legend(runList)
title('Change in ecosystem service indicators, 2001-2010 to 2091-2100')
xlabel('Indicator')
ylabel('Change (%)')
set(gca, 'FontSize', fontSize) ;


%% Map differences from end of historical to end of future: Isoprene emissions 

% Options %%%%%%%%%
filename_base = [outDir_maps 'bvoc_iso'] ;
title_text = 'isoprene emissions' ;
sumvars = 'Total' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1e-3 ;   % mgC/m2 to gC/m2
units_map = 'gC m^{-2} yr ^{-1}' ;
conv_fact_total = 1e-6*cf_kg2Tg ;   % mgC to TgC
units_total = 'TgC yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_aiso_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: Monoterpene emissions

% Options %%%%%%%%%
filename_base = [outDir_maps 'bvoc_mon'] ;
title_text = 'monoterpene emissions' ;
sumvars = 'Total' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1e-3 ;   % mgC/m2 to gC/m2
units_map = 'gC m^{-2} yr ^{-1}' ;
conv_fact_total = 1e-6*cf_kg2Tg ;   % mgC to TgC
units_total = 'TgC yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_amon_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: Total N loss

% Options %%%%%%%%%
filename_base = [outDir_maps 'nloss'] ;
title_text = 'total N loss' ;
sumvars = {'flux','leach'} ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % kgN/m2
units_map = 'kgN m^{-2} yr ^{-1}' ;
conv_fact_total = cf_kg2Tg ;   % kgN to TgN
units_total = 'TgN yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_nflux_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: Gaseous N loss

% Options %%%%%%%%%
filename_base = [outDir_maps 'nloss_gas'] ;
title_text = 'gaseous N loss' ;
sumvars = 'flux' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % kgN/m2
units_map = 'kgN m^{-2} yr ^{-1}' ;
conv_fact_total = cf_kg2Tg ;   % kgN to TgN
units_total = 'TgN yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_nflux_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: Dissolved N loss

% Options %%%%%%%%%
filename_base = [outDir_maps 'nloss_liq'] ;
title_text = 'dissolved N loss' ;
sumvars = 'leach' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % kgN/m2
units_map = 'kgN m^{-2} yr ^{-1}' ;
conv_fact_total = cf_kg2Tg ;   % kgN to TgN
units_total = 'TgN yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_nflux_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: January albedo

% warning('Skipping January albedo: pass per-year non-bare area as weighting!')

% Options %%%%%%%%%
filename_base = [outDir_maps 'albedo_jan'] ;
title_text = 'January albedo' ;
sumvars = 'January' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % unitless
units_map = 'unitless' ;
conv_fact_total = [] ;   % unitless
units_total = '' ;
this_land_area_map = land_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_albedo_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: July albedo

% warning('Skipping July albedo: pass per-year non-bare area as weighting!')

% Options %%%%%%%%%
filename_base = [outDir_maps 'albedo_jul'] ;
title_text = 'July albedo' ;
sumvars = 'July' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % unitless
units_map = '' ;
conv_fact_total = [] ;   % unitless
units_total = '' ;
this_land_area_map = land_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_albedo_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: Vegetation C

% Options %%%%%%%%%
filename_base = [outDir_maps 'cpool_veg'] ;
title_text = 'vegetation C' ;
sumvars = 'VegC' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1e-6*1e4 ;   % kgC/m2 to tonsC/ha
units_map = 'tons C ha^{-1}' ;
conv_fact_total = cf_kg2Pg ;   % kgC to PgC
units_total = 'PgC' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_cpool_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: Total C

% Options %%%%%%%%%
filename_base = [outDir_maps 'cpool_tot'] ;
title_text = 'total C' ;
sumvars = 'Total' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1e-6*1e4 ;   % kgC/m2 to tonsC/ha
units_map = 'tons C ha^{-1}' ;
conv_fact_total = cf_kg2Pg ;   % kgC to PgC
units_total = 'PgC' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_cpool_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: Evapotranspiration

% Options %%%%%%%%%
filename_base = [outDir_maps 'water_evapotransp'] ;
title_text = 'evapotranspiration' ;
sumvars = {'Evap','Transp'} ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % mm
units_map = 'mm yr^{-1}' ;
conv_fact_total = 1e-3*cf_m3_to_km3*1e-3 ;   % mm to 1000 km3
units_total = '1000 km^3 yr^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_awater_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: Annual runoff

% Options %%%%%%%%%
filename_base = [outDir_maps 'water_runoff'] ;
title_text = 'annual runoff' ;
sumvars = 'Runoff' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % mm
units_map = 'mm yr^{-1}' ;
conv_fact_total = 1e-3*cf_m3_to_km3*1e-3 ;   % mm to 1000 km3
units_total = '1000 km^3 yr^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_awater_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;


%% Map differences from end of historical to end of future: Peak runoff

% Options %%%%%%%%%
filename_base = [outDir_maps 'water_runoff_peak'] ;
title_text = 'peak monthly runoff' ;
sumvars = 'Max' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % mm
units_map = 'mm mon^{-1}' ;
conv_fact_total = [] ;   % Left empty to avoid calculation of total
units_total = '' ;
this_land_area_map = [] ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, maps_pk_runoff_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;

%% stop

stop


%% Map changes in LU area: End-Historical to End-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = figurePos ;
nx = 3 ;
ny = 2 ;
colorBarLoc = 'SouthOutside' ;
conv_fact_map = 1e-6 ;   % m2 to km2
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
units_map = 'km^2' ;
units_total = 'Mkm^2' ;
only1bl = true ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;

% Natural area
make_LUdiff_fig_v2(...
    ntrl_area_YXBH, ntrl_diff_YXrH, ...
    thisY1, thisYN, '"Natural"', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
    close
end

% Cropland area
make_LUdiff_fig_v2(...
    crop_area_YXBH, crop_diff_YXrH, ...
    thisY1, thisYN, 'Cropland', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
    close
end

% "Pasture" area
make_LUdiff_fig_v2(...
    past_area_YXBH, past_diff_YXrH, ...
    thisY1, thisYN, '"Pasture"', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
    close
end

% Agricultural area
make_LUdiff_fig_v2(...
    agri_area_YXBH, agri_diff_YXrH, ...
    thisY1, thisYN, 'Agricultural', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_agri.png'],['-r' num2str(pngres)])
    close
end

% Cropland0 area
if exist('crop0_area_YXBH','var')
    make_LUdiff_fig_v2(...
        crop0_area_YXBH, crop0_diff_YXrH, ...
        thisY1, thisYN, 'Cropland0', runList, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
        Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
    if do_save
        export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop0.png'],['-r' num2str(pngres)])
        close
    end
end

% "Pasture0" area
if exist('past0_area_YXBH','var')
    make_LUdiff_fig_v2(...
        past0_area_YXBH, past0_diff_YXrH, ...
        thisY1, thisYN, '"Pasture0"', runList, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
        Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
    if do_save
        export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past0.png'],['-r' num2str(pngres)])
        close
    end
end


%% Map changes in LU area: End-Historical to Begin-Future

% % Options %%%%%%%%%
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% textX = 25 ;
% textY_1 = 50 ;
% textY_2 = 20 ;
% thisPos = figurePos ;
% nx = 3 ;
% ny = 2 ;
% colorBarLoc = 'SouthOutside' ;
% conv_fact_map = 1e-6 ;   % m2 to km2
% conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
% units_map = 'km^2' ;
% units_total = 'Mkm^2' ;
% only1bl = true ;
% %%%%%%%%%%%%%%%%%%%
% 
% thisY1 = yearList_baseline(end) ;
% thisYN = yearList_future(1) ;
% 
% % Natural area
% make_LUdiff_fig_v2(...
%     ntrl_area_YXBH, ntrl_area_YXBFr - repmat(ntrl_area_YXBH,[1 1 Nruns]), ...
%     thisY1, thisYN, '"Natural"', runList, ...
%     spacing, fontSize, textX, textY_1, textY_2, ...
%     nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%     Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
% if do_save
%     export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
%     close
% end
% 
% % Cropland area
% make_LUdiff_fig_v2(...
%     crop_area_YXBH, crop_area_YXBFr - repmat(crop_area_YXBH,[1 1 Nruns]), ...
%     thisY1, thisYN, 'Cropland', runList, ...
%     spacing, fontSize, textX, textY_1, textY_2, ...
%     nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%     Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
% if do_save
%     export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
%     close
% end
% 
% % "Pasture" area
% make_LUdiff_fig_v2(...
%     past_area_YXBH, past_area_YXBFr - repmat(past_area_YXBH,[1 1 Nruns]), ...
%     thisY1, thisYN, '"Pasture"', runList, ...
%     spacing, fontSize, textX, textY_1, textY_2, ...
%     nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%     Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
% if do_save
%     export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
%     close
% end
% 
% % Agricultural area
% make_LUdiff_fig_v2(...
%     agri_area_YXBH, agri_area_YXBFr - repmat(agri_area_YXBH,[1 1 Nruns]), ...
%     thisY1, thisYN, 'Agricultural', runList, ...
%     spacing, fontSize, textX, textY_1, textY_2, ...
%     nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%     Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
% if do_save
%     export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_agri.png'],['-r' num2str(pngres)])
%     close
% end
% 
% % Cropland0 area
% if exist('crop0_area_YXBH','var')
%     make_LUdiff_fig_v2(...
%         crop0_area_YXBH, crop_area_YXBFr - repmat(crop0_area_YXBH,[1 1 Nruns]), ...
%         thisY1, thisYN, 'Cropland0', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop0.png'],['-r' num2str(pngres)])
%         close
%     end
% end
% 
% % "Pasture0" area
% if exist('past0_area_YXBH','var')
%     make_LUdiff_fig_v2(...
%         past0_area_YXBH, past_area_YXBFr - repmat(past0_area_YXBH,[1 1 Nruns]), ...
%         thisY1, thisYN, '"Pasture0"', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past0.png'],['-r' num2str(pngres)])
%         close
%     end
% end


%% Map changes in LU area: Begin-Future to End-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.05 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = [1 33 935 772] ;
nx = 2 ;
ny = 4 ;
colorBarLoc = 'EastOutside' ;
conv_fact_map = 1e-6 ;   % m2 to km2
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
units_map = 'km^2' ;
units_total = 'Mkm^2' ;
only1bl = false ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_future(1) ;
thisYN = yearList_future(end) ;

% Natural area
make_LUdiff_fig_v2(...
    ntrl_area_YXBFr, ntrl_diff_YXrF, ...
    thisY1, thisYN, '"Natural"', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
    close
end

% Cropland area
make_LUdiff_fig_v2(...
    crop_area_YXBFr, crop_diff_YXrF, ...
    thisY1, thisYN, 'Cropland', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
    close
end

% "Pasture" area
make_LUdiff_fig_v2(...
    past_area_YXBFr, past_diff_YXrF, ...
    thisY1, thisYN, '"Pasture"', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
    close
end

% Agricultural area
make_LUdiff_fig_v2(...
    agri_area_YXBFr, agri_diff_YXrF, ...
    thisY1, thisYN, 'Agricultural', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_agri.png'],['-r' num2str(pngres)])
    close
end


%% Map changes in cropland and pasture area: End-historical to End-Future


%% Map changes in cropland and pasture area: Begin-Future to End-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.05 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = [1 33 935 772] ;
nx = 2 ;
ny = 4 ;
colorBarLoc = 'EastOutside' ;
conv_fact_map = 1e-6 ;   % m2 to km2
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
units_map = 'km^2' ;
units_total = 'Mkm^2' ;
only1bl = false ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_future(1) ;
thisYN = yearList_future(end) ;

make_LUdiff_fig_v3(...
    crop_area_YXBFr, past_area_YXBFr, ...
    crop_diff_YXrF, past_diff_YXrF, ...
    thisY1, thisYN, 'Cropland', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
if do_save
    export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_croppast.png'],['-r' num2str(pngres)])
    close
end



%% Map changes in each crop area: End-Historical to End-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
thisPos = figurePos ;
nx = 3 ;
ny = 2 ;
colorBarLoc = 'SouthOutside' ;
conv_fact_map = 1e-6 ;   % m2 to km2
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
units_map = 'km^2' ;
units_total = 'Mkm^2' ;
only1bl = true ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;
for c = 1:Ncrops
    make_LUdiff_fig_v2(...
        maps_cropareas_YXvBH(:,:,c), squeeze(maps_cropareasDiffs_YXvrH(:,:,c,:)), ...
        thisY1, thisYN, CFTnames_maps{c}, runList, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
        Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
    if do_save
        export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_maps{c} '.png'],['-r' num2str(pngres)])
        close
    end
end


%% Map changes in each crop area: End-Historical to Begin-Future

% % Options %%%%%%%%%
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% textX = 25 ;
% textY_1 = 50 ;
% textY_2 = 20 ;
% thisPos = figurePos ;
% nx = 3 ;
% ny = 2 ;
% colorBarLoc = 'SouthOutside' ;
% conv_fact_map = 1e-6 ;   % m2 to km2
% conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
% units_map = 'km^2' ;
% units_total = 'Mkm^2' ;
% only1bl = true ;
% %%%%%%%%%%%%%%%%%%%
% 
% thisY1 = yearList_baseline(end) ;
% thisYN = yearList_future(1) ;
% 
% for c = 1:Ncrops
%     make_LUdiff_fig_v2(...
%         maps_cropareas_YXvBH(:,:,c), squeeze(maps_cropareas_YXvBFr(:,:,c,:))-maps_cropareas_YXvBH(:,:,c), ...
%         thisY1, thisYN, CFTnames_maps{c}, runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_maps{c} '.png'],['-r' num2str(pngres)])
%         close
%     end
% end


%% Map changes in each crop area: Begin-Future to End-Future

% % Options %%%%%%%%%
% fontSize = 14 ;
% spacing = [0.05 0.05] ;   % [vert, horz]
% textX = 25 ;
% textY_1 = 50 ;
% textY_2 = 20 ;
% thisPos = [1 33 935 772] ;
% nx = 2 ;
% ny = 4 ;
% colorBarLoc = 'EastOutside' ;
% conv_fact_map = 1e-6 ;   % m2 to km2
% conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
% units_map = 'km^2' ;
% units_total = 'Mkm^2' ;
% only1bl = false ;
% %%%%%%%%%%%%%%%%%%%
% 
% thisY1 = yearList_future(1) ;
% thisYN = yearList_future(end) ;
% 
% for c = 1:Ncrops
%     make_LUdiff_fig_v2(...
%         squeeze(maps_cropareas_YXvBFr(:,:,c,:)), squeeze(maps_cropareasDiffs_YXvrF(:,:,c,:)), ...
%         thisY1, thisYN, CFTnames_maps{c}, runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_maps{c} '.png'],['-r' num2str(pngres)])
%         close
%     end
% end


%% Plot timeseries: Land uses

rebase = false ;

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 18 ;
spacing = [0.15 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
blYears = yearList_baseline ;
% blYears = 1960:yearList_baseline(end) ;
conv_fact = 1e-6*1e-6 ;   % m2 to Mkm2
units = 'Million km^2' ;
% thisPos = figurePos ;
thisPos = [1 383 1440 422] ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_LUarea_crop_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_crop_bl, ts_LUarea_crop_yr, ignYrs, yearList_future) ;
[tmp_ts_LUarea_past_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_past_bl, ts_LUarea_past_yr, ignYrs, yearList_future) ;
[tmp_ts_LUarea_ntrl_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_LUarea_ntrl_bl, ts_LUarea_ntrl_yr, ignYrs, yearList_future) ;

figure('Position',thisPos,'Color','w') ;

subplot_tight(1,2,1,spacing)
if exist('ts_LUarea_past0_bl','var')
    tmp = ts_LUarea_past0_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
    hold on
    tmp = ts_LUarea_crop0_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'--k','LineWidth',lineWidth)
end
tmp = ts_LUarea_past_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
if exist('ts_LUarea_past0_bl','var')
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'-','LineWidth',lineWidth,'Color',0.75*ones(1,3))
else
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
end
hold on
tmp = ts_LUarea_crop_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
if exist('ts_LUarea_crop0_bl','var')
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'--','LineWidth',lineWidth,'Color',0.75*ones(1,3))
else
    plot(blYears,conv_fact*movmean(tmp,Nsmth),'--k','LineWidth',lineWidth)
end
set(gca,'ColorOrderIndex',1) ;
plot(yearList_future,conv_fact*tmp_ts_LUarea_past_yr,'-','LineWidth',lineWidth)
set(gca,'ColorOrderIndex',1) ;
plot(yearList_future,conv_fact*tmp_ts_LUarea_crop_yr,'--','LineWidth',lineWidth)

hold off
% legend([strcat(stdLegend,', pasture') strcat(stdLegend,', cropland')], ...
%        'Location','NorthEastOutside') ;
legend('Pasture','Cropland', ...
       'Location','SouthEast') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel(units)
ht = title(['Agricultural area' title_suffix]) ;
letterlabel_align0('A',ht,do_caps) ;

subplot_tight(1,2,2,spacing) ;
tmp = ts_LUarea_ntrl_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot(blYears,conv_fact*movmean(tmp,Nsmth),'-k','LineWidth',lineWidth)
hold on
plot(yearList_future,conv_fact*tmp_ts_LUarea_ntrl_yr,'LineWidth',lineWidth)
hold off
legend(stdLegend, ...
       'Location','SouthWest') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel(units)
ht = title(['Natural area' title_suffix]) ;
letterlabel_align0('B',ht,do_caps) ;
stop
if do_save
    export_fig([outDir_ts 'landUse' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Land uses compared to comparables from Alexander et al. (2017)

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 18 ;
spacing = 0.07 ;
legloc = 'Southeast' ;
%%%%%%%%%%%%%%%%%%%

% Import
table_in = readtable('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/AlexanderEtAl2017_allNonGriddedData.csv') ;

% Restrict to global runs
table_in(~strcmp(table_in.location,'Global'),:) = [] ;

% Only care about cropland and pasture
table_in(~(strcmp(table_in.land_type,'c') | strcmp(table_in.land_type,'p')),:) = [] ;

% Only care about part of the time period
yearList_tmp = 2010:10:2100 ;
Nyears_tmp = length(yearList_tmp) ;
table_in(table_in.year<min(yearList_tmp) | table_in.year>max(yearList_tmp),:) = [] ;

% Make plots
skip_this = false ;
if strcmp(thisVer,'harm3')
    combos = {'SSP1', 'RCP 4.5' ;
              'SSP3', 'RCP 6.0' ;
              'SSP4', 'RCP 6.0' ;
              'SSP5', 'RCP 8.5' ;
             } ;
else
    warning('thisVer (%s) not recognized for comparing land uses to Alexander et al. (2017); skipping.', thisVer) ;
    skip_this = true ;
end
if ~skip_this
    if Nruns ~= 4
        error('Nruns ~= 4')
    end
    for c = 1:length(combos)
        
        % Get matches
        thisEcon = combos{c,1} ;
        thisClim = combos{c,2} ;
        thisTitle_tmp = sprintf('(%s-%s)',thisEcon, thisClim([5 7])) ;
        isMatch = strcmp(table_in.economicScenario,thisEcon) & strcmp(table_in.rcp,thisClim) ;
        column_runs = strcat(strcat(table_in.model(isMatch),':'), table_in.scenarioId(isMatch)) ;
        matching_runs = unique(column_runs) ;
        matching_runs_4legend = strrep(matching_runs,':',': ') ;
        Nmatch = length(matching_runs) ;
        if Nmatch==0
            continue
        end
        
        % Get data
        data_vyc = nan(2,Nyears_tmp,Nmatch) ;
        for i = 1:Nmatch
            thisRun = matching_runs{i} ;
            tmp = strsplit(thisRun, ':') ;
            thisModel = tmp{1} ;
            thisScenID = tmp{2} ;
            isThisRun = strcmp(table_in.model,thisModel) & strcmp(table_in.scenarioId,thisScenID) ;
            data_x = table_in.value_rebased(isThisRun) ;
            years_x = table_in.year(isThisRun) ;
            if isempty(intersect(years_x,yearList_tmp))
                error('No matching years found!')
            end
            lt_x = table_in.land_type(isThisRun) ;
            thisRun_isCrop = strcmp(lt_x,'c') ;
            for y = 1:Nyears_tmp
                thisYear = yearList_tmp(y) ;
                isThisYear = years_x==thisYear ;
                if any(isThisYear)
                    data_vyc(1,y,i) = data_x(isThisYear & thisRun_isCrop) ;
                    data_vyc(2,y,i) = data_x(isThisYear & ~thisRun_isCrop) ;
                end
            end
            
        end
        data_vyc = data_vyc * 1e10 ;   % Mha to m2
        % Rebase to LPJ-GUESS-PLUM 2011 level assuming linear between 2010-2020
        adj = squeeze(ts_LUarea_crop_yr(yearList_future==2011,c) - ...
            (0.9*data_vyc(1,1,:)+0.1*data_vyc(1,2,:))) ;
        adj = permute(adj,[3 2 1]) ;
        data_vyc(1,:,:) = data_vyc(1,:,:) + repmat(adj,[1 length(yearList_tmp) 1]) ;
        adj = squeeze(ts_LUarea_past_yr(yearList_future==2011,c) - ...
            (0.9*data_vyc(2,1,:)+0.1*data_vyc(2,2,:))) ;
        data_vyc(2,:,:) = data_vyc(2,:,:) + repmat(permute(adj,[3 2 1]),[1 length(yearList_tmp) 1]) ;
        data_vyc = data_vyc * 1e-12 ;   % m2 to Mkm2
        data_ycv = permute(data_vyc,[2 3 1]) ;
        
        % Plot
        matching_runs_4legend{end+1} = 'LPJ-GUESS-PLUM' ;
        figure('Position',figurePos,'Color','w') ;
        % Cropland
        subplot_tight(1,2,1,spacing)
        plot(yearList_tmp, data_ycv(:,:,1), 'LineWidth', lineWidth)
        hold on
        plot(yearList_future, ts_LUarea_crop_yr(:,c)*1e-12, 'k', 'LineWidth', lineWidth)
        hold off
        legend(matching_runs_4legend,'Location',legloc)
        set(gca,'FontSize',fontSize) ;
        thisTitle = sprintf('Cropland %s',thisTitle_tmp) ;
        title(thisTitle)
        xlabel('Year')
        ylabel('Area (Mkm^2)')
        % Pasture
        subplot_tight(1,2,2,spacing)
        plot(yearList_tmp, data_ycv(:,:,2), 'LineWidth', lineWidth)
        hold on
        plot(yearList_future, ts_LUarea_past_yr(:,c)*1e-12, 'k', 'LineWidth', lineWidth)
        hold off
        legend(matching_runs_4legend,'Location',legloc)
        set(gca,'FontSize',fontSize) ;
        thisTitle = sprintf('Pasture %s',thisTitle_tmp) ;
        title(thisTitle)
        xlabel('Year')
        ylabel('Area (Mkm^2)')
        
        % Save
        if do_save
            filename_out = sprintf('%s/landUse_compAlexander2017_%s.pdf', ...
                outDir_ts, thisTitle_tmp(2:end-1)) ;
            export_fig(filename_out)
            close
        end
        
    end
    
    clear table_in
end


%% Plot timeseries: N fertilizer, irrigation

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 18 ;
ignYrs = 0 ;
Nsmth = 1 ;
spacing = [0.15 0.1] ;   % vert, horiz
blYears = yearList_baseline ;
% blYears = 1960:yearList_baseline(end) ;
perArea = false ;
thisPos = [1 350 1440 455] ;
%%%%%%%%%%%%%%%%%%%

if perArea
    this_ts_nflux_fert_bl = ts_nflux_fert_bl ...
                         ./ (ts_LUarea_crop_bl * 1e-4) ; % m2 to ha
    this_ts_nflux_fert_yr = ts_nflux_fert_yr ...
                         ./ (ts_LUarea_crop_yr * 1e-4) ; % m2 to ha
%     this_ts_nflux_fert_bl = -(ts_nflux_fert_CerealsC3_bl+ts_nflux_fert_CerealsC4_bl+ts_nflux_fert_Miscanthus_bl+ts_nflux_fert_Oilcrops_bl+ts_nflux_fert_Pulses_bl+ts_nflux_fert_Rice_bl+ts_nflux_fert_StarchyRoots_bl) ...
%                          ./ (ts_LUarea_crop_bl * 1e-4) ; % m2 to ha
%     this_ts_nflux_fert_yr = -(ts_nflux_fert_CerealsC3_yr+ts_nflux_fert_CerealsC4_yr+ts_nflux_fert_Miscanthus_yr+ts_nflux_fert_Oilcrops_yr+ts_nflux_fert_Pulses_yr+ts_nflux_fert_Rice_yr+ts_nflux_fert_StarchyRoots_yr) ...
%                          ./ (ts_LUarea_crop_yr * 1e-4) ; % m2 to ha           
    units_nflux_fert = 'kg/ha' ;
else
    this_ts_nflux_fert_bl = ts_nflux_fert_bl * 1e-9 ; % kgN to TgN
    this_ts_nflux_fert_yr = ts_nflux_fert_yr * 1e-9 ; % kgN to TgN
    units_nflux_fert = 'TgN' ;
end
if perArea
    this_ts_irrig_bl = (ts_irrig_bl * 1e3) ... % m3 to L
                         ./ (ts_LUarea_crop_bl * 1e-4) ; % m2 to ha
    this_ts_irrig_yr = (ts_irrig_yr * 1e3) ... % m3 to L
                         ./ (ts_LUarea_crop_yr * 1e-4) ; % m2 to ha
    units_irrig = 'L/ha' ;
else
    this_ts_irrig_bl = ts_irrig_bl*cf_m3_to_km3*1e-3 ; % m3 to 1000km3
    this_ts_irrig_yr = ts_irrig_yr*cf_m3_to_km3*1e-3 ; % m3 to 1000km3
    units_irrig = '1000 km^3' ;
end


[tmp_this_ts_nflux_fert_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, this_ts_nflux_fert_bl, this_ts_nflux_fert_yr, ignYrs, yearList_future) ;

figure('Position',thisPos,'Color','w') ;
subplot_tight(1,2,1,spacing)
tmp = this_ts_nflux_fert_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    -tmp, [], -tmp_this_ts_nflux_fert_yr, ...
    1, Nsmth, ...
    units_nflux_fert, stdLegend, 'N fertilization', title_suffix, lineWidth, fontSize, ...
    skip3rdColor)

clear tmp_*

%%%%%%%%%%%%%%%%%%%
% Plot timeseries: Irrigation
% Options %%%%%%%%%
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_this_ts_irrig_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, this_ts_irrig_bl, this_ts_irrig_yr, ignYrs, yearList_future) ;

subplot_tight(1,2,2,spacing)
tmp = this_ts_irrig_bl(yearList_baseline>=blYears(1) & yearList_baseline<=blYears(end)) ;
plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    tmp, [], tmp_this_ts_irrig_yr, ...
    1, Nsmth, ...
    units_irrig, stdLegend, 'Irrigation', title_suffix, lineWidth, fontSize, ...
    skip3rdColor)

clear tmp_*

if perArea
    file_suffix = [file_suffix '_perArea'] ;
end

if do_save
    export_fig([outDir_ts 'mgmt_inputs' file_suffix '.pdf'])
    close
end


%% Plot timeseries: N fertilizer on each crop

% Options %%%%%%%%%
thisVar = 'nflux_fert' ;
title_prefix = 'N applied' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = 1e-3*1e-6 ;   % kg to Mt
units = 'Mt N' ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
thisLegend = stdLegend ;

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames);for c = 1:length(cmds); eval(cmds{c}); end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, {}, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;
if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end

% Again, but with fertilizer according to crop0
if exist('crop0_area_YXBH','var')
    thisVar = 'nflux0_fert' ;
    cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    file_suffix = CFT_timeseries_plot(...
        cell_bl, cell_yr, {}, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
        theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
        'lineWidth',lineWidth, ...
        'fontSize',fontSize, ...
        'spacing',spacing, ...
        'ignYrs',ignYrs, ...
        'Nsmth',Nsmth, ...
        'conv_fact',conv_fact) ;
    if do_save
        export_fig([outDir_ts thisVar file_suffix '.pdf'])
        close
    end
end


%% Plot timeseries: N fertilizer PER HECTARE on each crop

% % Options %%%%%%%%%
% thisVar = 'nflux_fert' ;
% title_prefix = 'N applied per ha' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% ignYrs = 0 ;
% Nsmth = 1 ;
% conv_fact = 1e4 ;   % kg/m2 to kg/ha (kg/m2 * m2/ha)
% units = 'kg ha^{-1}' ;
% %%%%%%%%%%%%%%%%%%%
% 
% theseCFTnames = CFTnames ;
% thisLegend = stdLegend ;
% 
% clear cell_bl cell_yr
% cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_bl_nflux_fert = cell_bl ;
% cell_yr_nflux_fert = cell_yr ;
% clear cell_bl cell_yr
% cmds = get_cell_forPlot(whos, 'croparea', 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_bl_croparea = cell_bl ;
% clear cell_bl
% cmds = get_cell_forPlot(whos, 'croparea0', 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_bl_croparea0 = cell_bl ;
% cmds = get_cell_forPlot(whos, 'croparea', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_yr_croparea = cell_yr ;
% clear cell_bl cell_yr
% % Get kg/m2
% for i = 1:length(cell_bl_nflux_fert)
%     cell_bl{i} = cell_bl_nflux_fert{i} ./ cell_bl_croparea{i} ;
% %     cell_bl{i} = cell_bl_nflux_fert{i} ./ cell_bl_croparea0{i} ;
% %     cell_bl{i} = cell_bl_nflux_fert{i} ./ cell_bl_croparea{i} .* (cell_bl_croparea0{i} ./ cell_bl_croparea{i});
%     cell_yr{i} = cell_yr_nflux_fert{i} ./ cell_yr_croparea{i} ;
% end
% 
% file_suffix = CFT_timeseries_plot(...
%     cell_bl, cell_yr, {}, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'ignYrs',ignYrs, ...
%     'Nsmth',Nsmth, ...
%     'conv_fact',conv_fact) ;
% if do_save
%     export_fig([outDir_ts 'nfertPerHa' file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: Production of each crop

% Options %%%%%%%%%
thisVar = 'cropprod' ;
title_prefix = 'Production' ;
plum_area_adjustment = 1 - unhCropFrac ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = cf_kg2Mt ;
units = 'Mt DM' ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
thisLegend = stdLegend_plusFAO ;

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end

% Adjust for "technology" increase
cell_bl = adj_yield_tech(cell_bl, yearList_baseline) ;
cell_yr = adj_yield_tech(cell_yr, yearList_future) ;

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact, ...
    'plum_area_adjustment', plum_area_adjustment) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries: Production of each crop, EXPECTED

% % Options %%%%%%%%%
% % plum_area_adjustment = 1 ;
% % plum_area_adjustment = 1-0.28 ;
% plum_area_adjustment = 1 - unhCropFrac ;
% lpjg_area_adjustment = 1 ;
% % lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;
% thisVar = 'cropprodExp' ;
% title_prefix = 'Production (exp.)' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% ignYrs = 0 ;
% Nsmth = 1 ;
% conv_fact = cf_kg2Mt ;
% units = 'Mt DM' ;
% %%%%%%%%%%%%%%%%%%%
% 
% theseCFTnames = CFTnames ;
% thisLegend = stdLegend_plusFAO ;
% 
% clear cell_bl cell_yr
% cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% 
% file_suffix = CFT_timeseries_plot(...
%     cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'ignYrs',ignYrs, ...
%     'Nsmth',Nsmth, ...
%     'conv_fact',conv_fact, ...
%     'plum_area_adjustment', plum_area_adjustment, ...
%     'lpjg_area_adjustment', lpjg_area_adjustment) ;
% 
% if do_save
%     export_fig([outDir_ts thisVar file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: Production of each crop, ACTUAL/EXPECTED

% % Options %%%%%%%%%
% % title_prefix = 'Prod. diff: (Act.-Exp./Exp.)' ;
% title_prefix = 'Prod. diff' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% Nsmth = 1 ;
% units = '%' ;
% %%%%%%%%%%%%%%%%%%%
% 
% theseCFTnames = CFTnames ;
% thisLegend = stdLegend(2:end) ;
% 
% clear cell_yr cell_yr_act cell_yr_exp
% cmds = get_cell_forPlot(whos, 'cropprod', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_yr_act = cell_yr ; clear cell_yr
% cell_yr_act = adj_yield_tech(cell_yr_act, yearList_future) ;
% cmds = get_cell_forPlot(whos, 'cropprodExp', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_yr_exp = cell_yr ; clear cell_yr
% CFT_ActVsExp_plot(...
%     cell_yr_act, cell_yr_exp, yearList_future, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'Nsmth',Nsmth) ;
% 
% if do_save
%     export_fig([outDir_ts 'cropProdDiffFromExp.pdf'])
%     close
% end



%% Plot timeseries: Area of each crop

% Options %%%%%%%%%
plum_area_adjustment = 1 ;
% plum_area_adjustment = 1-0.28 ;
% plum_area_adjustment = 1 - unhCropFrac ;
lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;
thisVar = 'croparea' ;
title_prefix = 'Area' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = 1e-6*1e-6 ;
units = 'Million km^2' ;
%%%%%%%%%%%%%%%%%%%

if include_fao
    thisLegend = stdLegend_plusFAO ;
else
    thisLegend = stdLegend ;
end

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if include_fao
    cmds = get_cell_forPlot(whos, thisVar, 'fao', CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
else
    cell_fao = {} ;
    fao.tmp_fao_yearList = [] ;
end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    CFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact, ...
    'plum_area_adjustment', plum_area_adjustment, ...
    'lpjg_area_adjustment', lpjg_area_adjustment) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end



%% Plot timeseries: Yield of each crop

% Options %%%%%%%%%
thisVar = 'yield' ;
title_prefix = 'Yield' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = cf_kgPm2_to_tonsPha ;
units = 't ha^{-1}' ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
if include_fao
    thisLegend = stdLegend_plusFAO ;
else
    thisLegend = stdLegend ;
end

clear cell_bl cell_yr
cmds = get_cell_forPlot(whos, thisVar, 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
if include_fao
    cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
else
    cell_fao = {} ;
    fao.tmp_fao_yearList = [] ;
end

% Adjust for "technology" increase (do not include FAO!)
if do_adjYieldTech
    title_prefix = [title_prefix ' (techAdj)'] ;
end

file_suffix = CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact) ;

if do_adjYieldTech
    file_suffix = [file_suffix '_techAdj'] ;
end

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries: Yield of each crop (EXPECTED)

% % Options %%%%%%%%%
% thisVar = 'yieldExp' ;
% title_prefix = 'Yield (PLUM-exp)' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% ignYrs = 0 ;
% Nsmth = 1 ;
% conv_fact = cf_kgPm2_to_tonsPha ;
% units = 't ha^{-1}' ;
% %%%%%%%%%%%%%%%%%%%
% 
% theseCFTnames = CFTnames ;
% if include_fao
%     thisLegend = stdLegend_plusFAO ;
% else
%     thisLegend = stdLegend ;
% end
% 
% clear cell_bl cell_yr
% cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% if include_fao
%     cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% else
%     cell_fao = {} ;
%     fao.tmp_fao_yearList = [] ;
% end
% 
% file_suffix = CFT_timeseries_plot(...
%     cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'ignYrs',ignYrs, ...
%     'Nsmth',Nsmth, ...
%     'conv_fact',conv_fact) ;
% 
% if do_save
%     export_fig([outDir_ts thisVar file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: Yield of each crop, ACTUAL/EXPECTED

% % Options %%%%%%%%%
% title_prefix = 'Yield diff' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% Nsmth = 1 ;
% units = '%' ;
% %%%%%%%%%%%%%%%%%%%
% 
% theseCFTnames = CFTnames ;
% thisLegend = stdLegend(2:end) ;
% 
% clear cell_yr cell_yr_act cell_yr_exp
% cmds = get_cell_forPlot(whos, 'yield', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_yr_act = cell_yr ; clear cell_yr
% cmds = get_cell_forPlot(whos, 'yieldExp', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cell_yr_exp = cell_yr ; clear cell_yr
% 
% % Adjust for technology change, if doing so
% if do_adjYieldTech
%     file_suffix = [file_suffix '_techAdj'] ;
% end
% 
% CFT_ActVsExp_plot(...
%     cell_yr_act, cell_yr_exp, yearList_future, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'Nsmth',Nsmth) ;
% 
% if do_save
%     export_fig([outDir_ts 'yieldDiffFromExp' file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: total area (FOR TROUBLESHOOTING)

% error('Make this work with SI units!')
% figure ;
% plot(yearList_baseline,movmean(ts_LUarea_crop_bl+ts_LUarea_past_bl+ts_LUarea_ntrl_bl+ts_LUarea_bare_bl,Nsmth),'-k','LineWidth',lineWidth)
% hold on
% plot(yearList_future,ts_LUarea_crop_yr+ts_LUarea_past_yr+ts_LUarea_ntrl_yr+ts_LUarea_bare_yr,'LineWidth',lineWidth)
% hold off
% legend(stdLegend, ...
%        'Location','NorthEastOutside') ;
% set(gca,'FontSize',fontSize)
% xlabel('Year')
% ylabel('Million km^2')
% title('All land area')


%% Plot timeseries: N loss, total gaseous+dissolved

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
ignYrs = 0 ;
Nsmth = 1 ;
spacing = [0.1 0.1] ;   % [vert, horz]
%%%%%%%%%%%%%%%%%%%

tmp_ts_nflux_total_bl = ts_nflux_flux_bl + ts_nflux_leach_bl ;

[tmp_ts_nflux_total_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ...
    tmp_ts_nflux_total_bl, ...
    ts_nflux_flux_yr + ts_nflux_leach_yr, ...
    ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    tmp_ts_nflux_total_bl, [], tmp_ts_nflux_total_yr, ...
    cf_kg2Tg, Nsmth, ...
    'TgN', stdLegend, ['N losses (gaseous + dissolved)' title_suffix], '', lineWidth, fontSize, ...
    skip3rdColor) ;

clear tmp_*

if do_save
    export_fig([outDir_ts 'Nloss_gasdis' file_suffix '.pdf'])
    close
end


%% Plot timeseries: N loss, separate gaseous and dissolved

% % Options %%%%%%%%%
% lineWidth = 3 ;
% fontSize = 24 ;
% ignYrs = 0 ;
% Nsmth = 5 ;
% spacing = [0.1 0.1] ;   % [vert, horz]
% %%%%%%%%%%%%%%%%%%%
% 
% % % ts_nloss_bl = ts_nflux_flux_bl + ts_nflux_harvest_bl + ts_nflux_leach_bl + ts_nflux_LUch_bl ;
% % % ts_nloss_yr = ts_nflux_flux_yr + ts_nflux_harvest_yr + ts_nflux_leach_yr + ts_nflux_LUch_yr ;
% % ts_nloss_bl = ts_nflux_flux_bl + ts_nflux_leach_bl ;
% % ts_nloss_yr = ts_nflux_flux_yr + ts_nflux_leach_yr ;
% % 
% [tmp_ts_nflux_flux_yr, title_suffix, file_suffix] = ...
%     rebase_future2baseline(rebase, Nsmth, ts_nflux_flux_bl, ts_nflux_flux_yr, ignYrs, yearList_future) ;
% [tmp_ts_nflux_leach_yr, ~, ~] = ...
%     rebase_future2baseline(rebase, Nsmth, ts_nflux_leach_bl, ts_nflux_leach_yr, ignYrs, yearList_future) ;
% 
% figure('Position',figurePos,'Color','w') ;
% 
% subplot_tight(1,2,1,spacing)
% ht = plot_timeseries(...
%     yearList_baseline, [], yearList_future, ...
%     ts_nflux_flux_bl, [], tmp_ts_nflux_flux_yr, ...
%     cf_kg2Tg, Nsmth, ...
%     'TgN', stdLegend, 'N loss: Gaseous', '', lineWidth, fontSize, ...
%     skip3rdColor) ;
% letterlabel_align0('A',ht,do_caps) ;
% 
% subplot_tight(1,2,2,spacing)
% ht = plot_timeseries(...
%     yearList_baseline, [], yearList_future, ...
%     ts_nflux_leach_bl, [], tmp_ts_nflux_leach_yr, ...
%     cf_kg2Tg, Nsmth, ...
%     'TgN', stdLegend, 'N loss: Dissolved', '', lineWidth, fontSize, ...
%     skip3rdColor) ;
% letterlabel_align0('B',ht,do_caps) ;
% 
% clear tmp_*
% 
% if do_save
%     export_fig([outDir_ts 'Nloss' file_suffix '.pdf'])
%     close
% end


%% Bar graph: N loss

% years_endh = 2000:2009 ;
% years_begf = 2011:2020 ;
% years_endf = 2090:2099 ;
% 
% % Name, code, conversion factor, formatSpec mean, formatSpec SEM
% rowInfo = {'N loss (TgN)', 'nloss', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'N loss: Gaseous (TgN)', 'nflux_flux', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'N loss: Dissolved (TgN)', 'nflux_leach', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            } ;


%% Plot timeseries: Albedo

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_albedo1_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_albedo1_bl, ts_albedo1_yr, ignYrs, yearList_future) ;
[tmp_ts_albedo7_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_albedo7_bl, ts_albedo7_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

plot(yearList_baseline,movmean(ts_albedo1_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
if skip3rdColor
    plot(yearList_future,tmp_ts_albedo1_yr(:,1),'LineWidth',lineWidth)
    h = plot(yearList_future,tmp_ts_albedo1_yr(:,2),'LineWidth',lineWidth) ;
    set(h,'Color',[255 159 56]/255) ;
    set(gca,'ColorOrderIndex',4) ;
    plot(yearList_future,tmp_ts_albedo1_yr(:,3),'LineWidth',lineWidth)
else
    plot(yearList_future,tmp_ts_albedo1_yr,'LineWidth',lineWidth)
end
ax = gca;
ax.ColorOrderIndex = 1;
plot(yearList_baseline,movmean(ts_albedo7_bl,Nsmth),'--k','LineWidth',lineWidth)
% plot(yearList_future,tmp_ts_albedo7_yr,'--','LineWidth',lineWidth)
if skip3rdColor
    plot(yearList_future,tmp_ts_albedo7_yr(:,1),'--','LineWidth',lineWidth)
    h = plot(yearList_future,tmp_ts_albedo7_yr(:,2),'--','LineWidth',lineWidth) ;
    set(h,'Color',[255 159 56]/255) ;
    set(gca,'ColorOrderIndex',4) ;
    plot(yearList_future,tmp_ts_albedo7_yr(:,3),'--','LineWidth',lineWidth)
else
    plot(yearList_future,tmp_ts_albedo7_yr,'--','LineWidth',lineWidth)
end
hold off
legend([strcat(stdLegend,', Jan.') strcat(stdLegend,', Jul.')],...
       'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('Albedo')
title(['Albedo' title_suffix])

clear tmp_*

if do_save
    export_fig([outDir_ts 'albedo' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Albedo, boreal forest

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_albedo1_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_albedo1_borfor_bl, ts_albedo1_borfor_yr, ignYrs, yearList_future) ;
[tmp_ts_albedo7_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_albedo7_borfor_bl, ts_albedo7_borfor_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

plot(yearList_baseline,movmean(ts_albedo1_borfor_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
if skip3rdColor
    plot(yearList_future,tmp_ts_albedo1_yr(:,1),'LineWidth',lineWidth)
    h = plot(yearList_future,tmp_ts_albedo1_yr(:,2),'LineWidth',lineWidth) ;
    set(h,'Color',[255 159 56]/255) ;
    set(gca,'ColorOrderIndex',4) ;
    plot(yearList_future,tmp_ts_albedo1_yr(:,3),'LineWidth',lineWidth)
else
    plot(yearList_future,tmp_ts_albedo1_yr,'LineWidth',lineWidth)
end
ax = gca;
ax.ColorOrderIndex = 1;
plot(yearList_baseline,movmean(ts_albedo7_borfor_bl,Nsmth),'--k','LineWidth',lineWidth)
% plot(yearList_future,tmp_ts_albedo7_yr,'--','LineWidth',lineWidth)
if skip3rdColor
    plot(yearList_future,tmp_ts_albedo7_yr(:,1),'--','LineWidth',lineWidth)
    h = plot(yearList_future,tmp_ts_albedo7_yr(:,2),'--','LineWidth',lineWidth) ;
    set(h,'Color',[255 159 56]/255) ;
    set(gca,'ColorOrderIndex',4) ;
    plot(yearList_future,tmp_ts_albedo7_yr(:,3),'--','LineWidth',lineWidth)
else
    plot(yearList_future,tmp_ts_albedo7_yr,'--','LineWidth',lineWidth)
end
hold off
legend([strcat(stdLegend,', Jan.') strcat(stdLegend,', Jul.')],...
       'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('Albedo')
title(['Albedo (boreal forest)' title_suffix])

clear tmp_*

if do_save
    export_fig([outDir_ts 'albedo_borfor' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Albedo, tundra

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_albedo1_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_albedo1_tundra_bl, ts_albedo1_tundra_yr, ignYrs, yearList_future) ;
[tmp_ts_albedo7_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_albedo7_tundra_bl, ts_albedo7_tundra_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

plot(yearList_baseline,movmean(ts_albedo1_tundra_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on
if skip3rdColor
    plot(yearList_future,tmp_ts_albedo1_yr(:,1),'LineWidth',lineWidth)
    h = plot(yearList_future,tmp_ts_albedo1_yr(:,2),'LineWidth',lineWidth) ;
    set(h,'Color',[255 159 56]/255) ;
    set(gca,'ColorOrderIndex',4) ;
    plot(yearList_future,tmp_ts_albedo1_yr(:,3),'LineWidth',lineWidth)
else
    plot(yearList_future,tmp_ts_albedo1_yr,'LineWidth',lineWidth)
end
ax = gca;
ax.ColorOrderIndex = 1;
plot(yearList_baseline,movmean(ts_albedo7_tundra_bl,Nsmth),'--k','LineWidth',lineWidth)
% plot(yearList_future,tmp_ts_albedo7_yr,'--','LineWidth',lineWidth)
if skip3rdColor
    plot(yearList_future,tmp_ts_albedo7_yr(:,1),'--','LineWidth',lineWidth)
    h = plot(yearList_future,tmp_ts_albedo7_yr(:,2),'--','LineWidth',lineWidth) ;
    set(h,'Color',[255 159 56]/255) ;
    set(gca,'ColorOrderIndex',4) ;
    plot(yearList_future,tmp_ts_albedo7_yr(:,3),'--','LineWidth',lineWidth)
else
    plot(yearList_future,tmp_ts_albedo7_yr,'--','LineWidth',lineWidth)
end
hold off
legend([strcat(stdLegend,', Jan.') strcat(stdLegend,', Jul.')],...
       'Location','NorthEastOutside') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel('Albedo')
title(['Albedo (tundra)' title_suffix])

clear tmp_*

if do_save
    export_fig([outDir_ts 'albedo_tundra' file_suffix '.pdf'])
    close
end


%% Plot timeseries: C pools

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_cpool_VegC_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_cpool_VegC_bl, ts_cpool_VegC_yr, ignYrs, yearList_future) ;
[tmp_ts_cpool_Total_yr, ~, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_cpool_Total_bl, ts_cpool_Total_yr, ignYrs, yearList_future) ;


figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_cpool_VegC_bl, [], tmp_ts_cpool_VegC_yr, ...
    cf_kg2Pg, Nsmth, ...
    'PgC', stdLegend, 'Vegetation C', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('A',ht,do_caps) ;

subplot_tight(1,2,2,spacing)
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_cpool_Total_bl, [], tmp_ts_cpool_Total_yr, ...
    cf_kg2Pg, Nsmth, ...
    'PgC', stdLegend, 'Total C', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('B',ht,do_caps) ;

clear tmp_*

if do_save
    export_fig([outDir_ts 'cpools' file_suffix '.pdf'])
    close
end


%% Plot timeseries: BVOCs

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_amon_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_amon_bl, ts_amon_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing) ;
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_amon_bl, [], tmp_ts_amon_yr, ...
    cf_kg2Pg, Nsmth, ...
    'TgC', stdLegend, 'Monoterpene emissions', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('A',ht,do_caps) ;
clear tmp_*

[tmp_ts_aiso_yr, title_suffix, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_aiso_bl, ts_aiso_yr, ignYrs, yearList_future) ;

subplot_tight(1,2,2,spacing) ;
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_aiso_bl, [], tmp_ts_aiso_yr, ...
    cf_kg2Pg, Nsmth, ...
    'TgC', stdLegend, 'Isoprene emissions', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('B',ht,do_caps) ;
set(gca,'FontSize',fontSize)

clear tmp_*

if do_save
    export_fig([outDir_ts 'bvoc' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Evapotranspiration and runoff

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 5 ;
figurePosition = [1 376 1440 429] ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_aevapaaet_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_aevapaaet_bl, ts_aevapaaet_yr, ignYrs, yearList_future) ;

figure('Position',figurePosition,'Color','w') ;
subplot_tight(1,2,1,spacing) ;
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_aevapaaet_bl, [], tmp_ts_aevapaaet_yr, ...
    cf_m3_to_km3*1e-3, Nsmth, ...
    '1000 km^3', stdLegend, 'Evapotranspiration', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('A',ht,do_caps) ;
clear tmp_*

[tmp_ts_tot_runoff_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_tot_runoff_bl, ts_tot_runoff_yr, ignYrs, yearList_future) ;
subplot_tight(1,2,2,spacing) ;
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_tot_runoff_bl, [], tmp_ts_tot_runoff_yr, ...
    cf_m3_to_km3*1e-3, Nsmth, ...
    '1000 km^3', stdLegend, 'Runoff', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('B',ht,do_caps) ;
clear tmp_*

if do_save
    export_fig([outDir_ts 'evapotranspiration_runoff' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Calories

% Options %%%%%%%%%
plum_area_adjustment = 1 ;
% plum_area_adjustment = 1-0.28 ;

lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;

lineWidth = 2 ;
fontSize = 14 ;
ignYrs = 0 ;
Nsmth = 1 ;
figurePosition = [1 376 1440 429] ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_kcal_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_kcal_bl, ts_kcal_yr, ignYrs, yearList_future) ;

% Adjust for technology change, if doing so (do not include FAO!)
tmp_ts_kcal_bl = ts_kcal_bl ;
if do_adjYieldTech
    file_suffix = [file_suffix '_techAdj'] ;
    title_suffix = [title_suffix ' (techAdj)'] ;
end

% Adjust for unhandled crops, if doing so (DO include FAO)
tmp_ts_kcal_fao = ts_kcal_fao ;
if ~isequal(lpjg_area_adjustment,1)
    tmp_ts_kcal_bl = tmp_ts_kcal_bl .* lpjg_area_adjustment ;
    [~,IA] = intersect(yearList_baseline,fao.tmp_fao_yearList) ;
    tmp_ts_kcal_fao = tmp_ts_kcal_fao .* lpjg_area_adjustment(IA) ;
    title_suffix = [title_suffix ' (blAdj)'] ;
    file_suffix = [file_suffix '_blAdj'] ;
end
if plum_area_adjustment ~= 1
    tmp_ts_kcal_yr = tmp_ts_kcal_yr * plum_area_adjustment ;
    title_suffix = [title_suffix ' (PLUM\times' num2str(plum_area_adjustment) ')'] ;
    file_suffix = [file_suffix '_plumAdj' num2str(plum_area_adjustment)] ;
end

figure('Position',figurePosition,'Color','w') ;

plot_timeseries(...
    yearList_baseline, fao.tmp_fao_yearList, yearList_future, ...
    tmp_ts_kcal_bl, tmp_ts_kcal_fao, tmp_ts_kcal_yr, ...
    cf_kcalEcal, Nsmth, ...
    'Ecal', stdLegend_plusFAO, 'Caloric production', ...
    title_suffix, lineWidth, fontSize, skip3rdColor) ;
clear tmp_*

if do_save
    export_fig([outDir_ts 'calories' file_suffix '.pdf'])
    close
end


%% Plot timeseries: Calories, EXPECTED

% % Options %%%%%%%%%
% % 
% plum_area_adjustment = 1 ;
% % plum_area_adjustment = 1-0.28 ;
% 
% % lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;
% 
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% ignYrs = 0 ;
% Nsmth = 1 ;
% figurePosition = [1 376 1440 429] ;
% %%%%%%%%%%%%%%%%%%%
% 
% [tmp_ts_kcalExp_yr, title_suffix, file_suffix] = ...
%     rebase_future2baseline(rebase, Nsmth, ts_kcal_bl, ts_kcalExp_yr, ignYrs, yearList_future) ;
% 
% % Adjust for technology change, if doing so (do not include FAO!)
% tmp_ts_kcal_bl = ts_kcal_bl ;
% if do_adjYieldTech
%     file_suffix = [file_suffix '_techAdj'] ;
%     title_suffix = [title_suffix ' (techAdj)'] ;
% end
% 
% % Adjust for unhandled crops, if doing so (DO include FAO)
% tmp_ts_kcal_fao = ts_kcal_fao ;
% if lpjg_area_adjustment ~= 1
%     tmp_ts_kcal_bl = tmp_ts_kcal_bl .* lpjg_area_adjustment ;
%     [~,IA] = intersect(yearList_baseline,fao.tmp_fao_yearList) ;
%     tmp_ts_kcal_fao = tmp_ts_kcal_fao .* lpjg_area_adjustment(IA) ;
%     title_suffix = [title_suffix ' (blAdj)'] ;
%     file_suffix = [file_suffix '_blAdj'] ;
% end
% if plum_area_adjustment ~= 1
%     tmp_ts_kcalExp_yr = tmp_ts_kcalExp_yr * plum_area_adjustment ;
%     title_suffix = [title_suffix ' (PLUM\times' num2str(plum_area_adjustment) ')'] ;
%     file_suffix = [file_suffix '_plumAdj' num2str(plum_area_adjustment)] ;
% end
% 
% figure('Position',figurePosition,'Color','w') ;
% plot(yearList_baseline,cf_kcalEcal*movmean(tmp_ts_kcal_bl,Nsmth),'-k','LineWidth',lineWidth)
% hold on
% plot(fao.tmp_fao_yearList,cf_kcalEcal*tmp_ts_kcal_fao,'--k','LineWidth',lineWidth)
% plot(yearList_future,cf_kcalEcal*tmp_ts_kcalExp_yr,'LineWidth',lineWidth)
% hold off
% legend(stdLegend_plusFAO, ...
%        'Location','NorthWest') ;
% set(gca,'FontSize',fontSize)
% xlabel('Year')
% ylabel('Ecal')
% ht = title(['Caloric production (PLUM-exp)' title_suffix]) ;
% clear tmp_*
% 
% if do_save
%     export_fig([outDir_ts 'caloriesExp' file_suffix '.pdf'])
%     close
% end


%% Plot timeseries: Calories, ACTUAL/EXPECTED

% % Options %%%%%%%%%
% lineWidth = 2 ;
% fontSize = 14 ;
% Nsmth = 1 ;
% figurePosition = [1 376 1440 429] ;
% %%%%%%%%%%%%%%%%%%%
% 
% file_suffix = '' ;
% title_suffix = '' ;
% if do_adjYieldTech
%     file_suffix = [file_suffix '_techAdj'] ;
%     title_suffix = [title_suffix ' (techAdj)'] ;
% end
% 
% figure('Position',figurePosition,'Color','w') ;
% tmp_ts_kcal_yr = ts_kcal_yr ;
% plot(yearList_future,movmean((tmp_ts_kcal_yr - ts_kcalExp_yr)./ts_kcalExp_yr*100,Nsmth,1),'-','LineWidth',lineWidth)
% hold on
% plot(get(gca,'XLim'),[0 0],'--k')
% hold off
% legend(stdLegend(2:end),'Location','SouthWest') ;
% set(gca,'FontSize',fontSize)
% xlabel('Year')
% ylabel('Ecal')
% title(['Calories diff.' title_suffix]) ;
% clear tmp_*
% 
% if do_save
%     export_fig([outDir_ts 'caloriesDiffFromExp' file_suffix '.pdf'])
%     close
% end



%% Map end-of-run mean cropfracs

% % Options %%%%%%%%%
% title_prefix = 'Frac. crop' ;
% whichCFTs = 'lpjg' ;
% lineWidth = 2 ;
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% y2include = 65:350 ;
% %%%%%%%%%%%%%%%%%%%
% 
% error('Make this work as PLUMnative!')
% 
% Ny = 2 ;
% Nmaps = length(runList) + 1 ;
% if strcmp(whichCFTs,'lpjg')
%     Ncrops = Ncrops_lpjg ;
%     theseCrops = CFTnames ;
%     Nx = 3 ;
%     theseMaps = maps_cropfracs ;
%     figure_position = figurePos ;
% elseif strcmp(whichCFTs,'plum')
%     Nx = 4 ;
%     theseMaps = maps_cropfracs_plum7 ;
%     figure_position = figurePos ;
% else
%     error(['whichCFTs (' whichCFTs ') not recognized!']) ;
% end
% 
% 
% for c = 1:Ncrops
%     figure('Position',figure_position,'Color','w') ;
%     thisCrop = theseCrops{c} ;
%     for p = 1:Nmaps
%         subplot_tight(Ny,Nx,p,spacing)
%         if p==1
%             tmp = theseMaps.maps_YXvyB(:,:,strcmp(theseMaps.varNames,thisCrop),:) ;
%             tmp(maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'CROPLAND'),:)==0) = 0 ;
%             tmp = mean(tmp,4) ;
%             pcolor(tmp(y2include,:)) ;
%             thisTitle = 'Baseline' ;
%         else
%             pcolor(mean(theseMaps.maps_YXvyr(y2include,:,strcmp(theseMaps.varNames,thisCrop),:,p-1),4)) ;
%             thisTitle = runList{p-1} ;
%         end
%         shading flat ; axis equal tight off
%         caxis([0 1]) ; colorbar('SouthOutside')
%         title([thisCrop ' (' thisTitle ')'])
%         set(gca,'FontSize',fontSize) ;
%     end
%     
%     if do_save
%         export_fig([outDir_maps 'maps_cropfracs_' thisCrop '.png'],['-r' num2str(pngres)])
%         close
%     end
% end
% 
% clear theseMaps Ncrops Nmaps Nx Ny
% 
% 
% %% Map changes in BD hotspot area: CI
% 
% % Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edgecolor = 0.6*ones(3,1) ;
% latlim = [-60,80];
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% textX = 25 ;
% textY_1 = 50 ;
% textY_2 = 20 ;
% cbarOrient = 'SouthOutside' ;
% lineWidth = 0.25 ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% error('Make this work with SI units!')
% 
% if strcmp(version('-release'),'2014b')
%     
%     % Biodiversity hotspots
%     hotspot_shp = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/hotspots_clipByGridlist.shp' ;
%     hotspot_YX = dlmread('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/hotspots_raster.txt',...
%         ' ', 6, 0) ;
%     hotspot_YX(hotspot_YX==-9999) = NaN ;
%     hotspot_YX(:,721) = [] ;
%     hotspot_YX = flipud(hotspot_YX) ;
%     hotspot_YX(nanmask) = NaN ;
%     hotspot_YX(~nanmask & isnan(hotspot_YX)) = 0 ;
%     hotspot_area_YX = hotspot_YX.*land_area_YX ;
%     hotspot_area_YXB = hotspot_area_YX .* maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'NATURAL'),end) ;
%     % hotspot_area_YXB = hotspot_area_YX ;
%     hotspot_area_YXr = repmat(hotspot_area_YX,[1 1 Nruns]) .* squeeze(maps_LU.maps_YXvyr(:,:,strcmp(maps_LU.varNames,'NATURAL'),end,:)) ;
%     
%     hotspot_diff_YXr = hotspot_area_YXr - repmat(hotspot_area_YXB,[1 1 Nruns]) ;
%     
%     map_hotspot_diffs(...
%         hotspot_area_YXB, hotspot_diff_YXr, hotspot_YX, hotspot_shp, ...
%         spacing, latlim, edgecolor, cbarOrient, fontSize, ...
%         textX, textY_1, textY_2, ssp_plot_index, lineWidth, ...
%         yearList_baseline, yearList_future, runList)
%     
%     if do_save
%         export_fig([outDir_maps 'areaDiff_BDhotspots_CI.png'],['-r' num2str(pngres)])
%         close
%     end
%     
% else
%     warning('Skipping hotspots')
% end


%% Map changes in BD hotspot area: glob200

% % Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edgecolor = 0.6*ones(3,1) ;
% latlim = [-60,80];
% fontSize = 14 ;
% spacing = [0.1 0.05] ;   % [vert, horz]
% textX = 25 ;
% textY_1 = 50 ;
% textY_2 = 20 ;
% cbarOrient = 'SouthOutside' ;
% lineWidth = 0.25 ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% error('Make this work with SI units!')
% 
% if strcmp(version('-release'),'2014b')
% 
%     % Biodiversity hotspots
%     hotspot_YX = imread('/Users/sam/Geodata/global200ecoregions/g200_terr_raster0.5deg.tif') ;
%     hotspot_shp = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work/g200_terr_from0.5raster.shp' ;
%     hotspot_YX = flipud(hotspot_YX) ;
%     hotspot_area_YX = hotspot_YX.*land_area_YX ;
%     hotspot_area_YXB = hotspot_area_YX .* maps_LU.maps_YXvyB(:,:,strcmp(maps_LU.varNames,'NATURAL'),end) ;
%     % hotspot_area_YXB = hotspot_area_YX ;
%     hotspot_area_YXr = repmat(hotspot_area_YX,[1 1 Nruns]) .* squeeze(maps_LU.maps_YXvyr(:,:,strcmp(maps_LU.varNames,'NATURAL'),end,:)) ;
%     
%     hotspot_diff_YXr = hotspot_area_YXr - repmat(hotspot_area_YXB,[1 1 Nruns]) ;
%     
%     map_hotspot_diffs(...
%         hotspot_area_YXB, hotspot_diff_YXr, hotspot_YX, hotspot_shp, ...
%         spacing, latlim, edgecolor, cbarOrient, fontSize, ...
%         textX, textY_1, textY_2, ssp_plot_index, lineWidth, ...
%         yearList_baseline, yearList_future, runList)
%     
%     if do_save
%         export_fig([outDir_maps 'areaDiff_BDhotspots_glob200.png'],['-r' num2str(pngres)])
%         close
%     end
% 
% else
%     warning('Skipping hotspots')
% end


%% Table after Krause et al. (2017) Table 2

% disp('Making table...')
% 
% years_endh = 2000:2009 ;
% years_begf = 2011:2020 ;
% years_endf = 2090:2099 ;
% 
% % Name, code, conversion factor, formatSpec mean, formatSpec SEM
% rowInfo = {'Vegetation C (GtC)', 'cpool_VegC', cf_kg2Pg, '%d', '%d' ;
%            'Soil and litter C (GtC)', 'cpool_LitterSoilC', cf_kg2Pg, '%d', '%d' ;
% %            'Product C (GtC)', 'cpool_HarvSlowC', cf_kg2Pg, '%0.1f', '%0.1f' ;
%            'Total C (GtC)', 'cpool_Total', cf_kg2Pg, '%d', '%d' ;
%            'January albedo', 'albedo1', 1, '%0.3f', '%0.3f' ;
%            'July albedo', 'albedo7', 1, '%0.3f', '%0.3f' ;
%            'Evapotranspiration (1000 km^3)', 'aevapaaet', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
%            'Runoff (1000 km^3)', 'tot_runoff', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
% %            'Peak monthly runoff (1000 km^3)', 'pkrunoff', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
%            'Crop production (Ecal)', 'kcal', cf_kcalEcal, '%0.1f', '%0.1f' ;
%            'N loss (TgN)', 'nloss', cf_kg2Tg, '%0.1f', '%0.1f' ;
% %            'N loss: Gaseous (TgN)', 'nflux_flux', cf_kg2Tg, '%0.1f', '%0.1f' ;
% %            'N loss: Dissolved (TgN)', 'nflux_leach', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'Isoprene emissions (TgC)', 'aiso', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'Monoterpene emissions (TgC)', 'amon', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'January albedo, boreal forest', 'albedo1_borfor', 1, '%0.3f', '%0.3f' ;
%            'January albedo,  tundra', 'albedo1_tundra', 1, '%0.3f', '%0.3f' ;
%            } ;
% 
% Nvars = size(rowInfo,1) ;
% mean_endh_v = nan(Nvars,1) ;
% mean_begf = nan(Nvars,Nruns) ;
% mean_endf_vr = nan(Nvars,Nruns) ;
% sem_endh_v = nan(Nvars,1) ;
% sem_begf = nan(Nvars,Nruns) ;
% sem_endf_vr = nan(Nvars,Nruns) ;
% string_endh = cell(Nvars,1) ;
% string_begf = cell(Nvars,Nruns) ;
% string_endf = cell(Nvars,Nruns) ;
% for c = 1:Nvars
%     
%     % Get values
%     thisVar = rowInfo{c,2} ;
%     thisConv = rowInfo{c,3} ;
%     mean_endh_v(c) = thisConv*eval(['mean(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
%     sem_endh_v(c) = thisConv*eval(['std(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
%     mean_begf(c,:) = thisConv*eval(['mean(ts_' thisVar '_yr(yearList_future>=min(years_begf) & yearList_future<=max(years_begf),:))']) ;
%     sem_begf(c,:) = thisConv*eval(['std(ts_' thisVar '_yr(yearList_future>=min(years_begf) & yearList_future<=max(years_begf),:))']) ;
%     mean_endf_vr(c,:) = thisConv*eval(['mean(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
%     sem_endf_vr(c,:) = thisConv*eval(['std(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
%     
%     % Turn into strings
%     if strcmp(rowInfo{c,4},'%d')
%         thisMean = round(mean_endh_v(c)) ;
%     else
%         thisMean = mean_endh_v(c) ;
%     end
%     if strcmp(rowInfo{c,4},'%d')
%         thisSD = round(sem_endh_v(c)) ;
%     else
%         thisSD = sem_endh_v(c) ;
%     end
%     string_endh{c} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean thisSD]) ;
%     for r = 1:Nruns
%         if strcmp(rowInfo{c,4},'%d')
%             thisMean_begf = round(mean_begf(c,r)) ;
%             thisMean_endf = round(mean_endf_vr(c,r)) ;
%         else
%             thisMean_begf = mean_begf(c,r) ;
%             thisMean_endf = mean_endf_vr(c,r) ;
%         end
%         if strcmp(rowInfo{c,4},'%d')
%             thisSD_begf = round(sem_begf(c,r)) ;
%             thisSD_endf = round(sem_endf_vr(c,r)) ;
%         else
%             thisSD_begf = sem_begf(c,r) ;
%             thisSD_endf = sem_endf_vr(c,r) ;
%         end
%         string_begf{c,r} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean_begf thisSD_begf]) ;
%         string_endf{c,r} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean_endf thisSD_endf]) ;
%     end
% 
% end
% 
% table_out = table(collate_empties(rowInfo(:,1)),...
%                   collate_empties(string_endh)) ;
% for r = 1:Nruns
%     table_out = [table_out collate_twocells(string_begf(:,r),string_endf(:,r))] ;
% end
% table_out.Properties.VariableNames = [{'Ecosystem_function','Baseline'} runColNames] ;
% 
% if do_save
%     writetable(table_out,[outDir_base 'summary_table.xlsx'],'Sheet',1) ;
% end
% disp('Done making table.')


%% Table after Krause et al. (2017) Table 2 but without first fut. decade

% disp('Making table...')
% 
% years_endh = 2000:2009 ;
% years_endf = 2090:2099 ;
% 
% % Name, code, conversion factor, formatSpec mean, formatSpec SEM
% rowInfo = {'Vegetation C (GtC)', 'cpool_VegC', cf_kg2Pg, '%d', '%d' ;
%            'Soil and litter C (GtC)', 'cpool_LitterSoilC', cf_kg2Pg, '%d', '%d' ;
% %            'Product C (GtC)', 'cpool_HarvSlowC', cf_kg2Pg, '%0.1f', '%0.1f' ;
%            'Total C (GtC)', 'cpool_Total', cf_kg2Pg, '%d', '%d' ;
%            'January albedo', 'albedo1', 1, '%0.3f', '%0.3f' ;
%            'July albedo', 'albedo7', 1, '%0.3f', '%0.3f' ;
%            'Evapotranspiration (1000 km^3)', 'aevapaaet', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
%            'Runoff (1000 km^3)', 'tot_runoff', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
% %            'Peak monthly runoff (1000 km^3)', 'pkrunoff', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
%            'Crop production (Ecal)', 'kcal', cf_kcalEcal, '%0.1f', '%0.1f' ;
%            'N loss (TgN)', 'nloss', cf_kg2Tg, '%0.1f', '%0.1f' ;
% %            'N loss: Gaseous (TgN)', 'nflux_flux', cf_kg2Tg, '%0.1f', '%0.1f' ;
% %            'N loss: Dissolved (TgN)', 'nflux_leach', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'Isoprene emissions (TgC)', 'aiso', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'Monoterpene emissions (TgC)', 'amon', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'January albedo, boreal forest', 'albedo1_borfor', 1, '%0.3f', '%0.3f' ;
%            'January albedo,  tundra', 'albedo1_tundra', 1, '%0.3f', '%0.3f' ;
%            } ;
% 
% Nvars = size(rowInfo,1) ;
% mean_endh_v = nan(Nvars,1) ;
% mean_endf_vr = nan(Nvars,Nruns) ;
% sem_endh_v = nan(Nvars,1) ;
% sem_endf_vr = nan(Nvars,Nruns) ;
% string_endh = cell(Nvars,1) ;
% string_endf = cell(Nvars,Nruns) ;
% for c = 1:Nvars
%     
%     % Get values
%     thisVar = rowInfo{c,2} ;
%     thisConv = rowInfo{c,3} ;
%     mean_endh_v(c) = thisConv*eval(['mean(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
%     sem_endh_v(c) = thisConv*eval(['std(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
%     mean_endf_vr(c,:) = thisConv*eval(['mean(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
%     sem_endf_vr(c,:) = thisConv*eval(['std(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
%     
%     % Turn into strings
%     if strcmp(rowInfo{c,4},'%d')
%         thisMean = round(mean_endh_v(c)) ;
%     else
%         thisMean = mean_endh_v(c) ;
%     end
%     if strcmp(rowInfo{c,4},'%d')
%         thisSD = round(sem_endh_v(c)) ;
%     else
%         thisSD = sem_endh_v(c) ;
%     end
%     string_endh{c} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean thisSD]) ;
%     for r = 1:Nruns
%         if strcmp(rowInfo{c,4},'%d')
%             thisMean_endf = round(mean_endf_vr(c,r)) ;
%         else
%             thisMean_endf = mean_endf_vr(c,r) ;
%         end
%         if strcmp(rowInfo{c,4},'%d')
%             thisSD_endf = round(sem_endf_vr(c,r)) ;
%         else
%             thisSD_endf = sem_endf_vr(c,r) ;
%         end
%         string_endf{c,r} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean_endf thisSD_endf]) ;
%     end
% 
% end
% 
% % table_out = table(collate_empties(rowInfo(:,1)),...
% %                   collate_empties(string_endh)) ;
% table_out = table(rowInfo(:,1),...
%                   string_endh) ;
%               
% for r = 1:Nruns
% %     table_out = [table_out collate_twocells(string_begf(:,r),string_endf(:,r))] ;
%     table_out = [table_out string_endf(:,r)] ;
% end
% table_out.Properties.VariableNames = [{'Ecosystem_function','Baseline'} runColNames] ;
% 
% if do_save
%     writetable(table_out,[outDir_base 'summary_table_nobegf.xlsx'],'Sheet',1) ;
% end
% disp('Done making table.')


%%

disp('All done!')



