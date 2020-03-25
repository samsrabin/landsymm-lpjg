%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLUM2LPJG figures: Multiple runs: %%%
%%% PLUM-style native %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_figs = false ;

if false
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
end

% thisVer = 'harm3' ;
thisVer = 'harm3_constLU' ;
% thisVer = 'harm3_constClim' ;
% thisVer = 'harm3_constCO2' ;
% thisVer = 'harm3_constClimCO2' ;
% thisVer = 'harm3_onlyCO2' ;
% thisVer = 'harm3_onlyClim' ;
% thisVer = 'harm3_S1R4.5_attr' ;
% thisVer = 'harm3_S3R6.0_attr' ;
% thisVer = 'harm3_S4R6.0_attr' ;
% thisVer = 'harm3_S5R8.5_attr' ;

do_adjYieldTech = true ; % Apply annual tech. change increase to yields?

unhCropFrac = 0 ; % Set to zero for previous behavior. v10 = 0.177

% ignored_crops = {'CC3G','CC4G'} ;
% ignored_crops = {'CC3G','CC4G','Miscanthus'} ;
ignored_crops = {'CC3G','CC4G','ExtraCrop'} ;

do_save = false ;
rebase = false ;
pngres = 300 ;

do_caps = -1 ;

years_endh = 2001:2010 ;
years_endf = 2091:2100 ;

       
%% Setup and import

tic
run('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/MATLAB_work/plum2lpjg_figs_setup_import_garr.m') ;
% addpath('/Users/sam/Documents/Dropbox/FireMIP/Stijn_files_2019/')
% addpath('/Users/Shared/PLUM/crop_calib_code/')
toc


%% Big bar graph: Drivers

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
orientation = 'h' ; % v(ertical) or h(orizontal)
% sd_or_sem = 'st. dev.' ;
% sd_or_sem = 'SEM' ;
sd_or_sem = '' ;
% errbar_color = 'k' ;
errbar_color = 0.5*ones(3,1) ;
fontSize = 14 ;
% figure_position = [1    33   846   772] ;
% figure_position = [1    33   720   772] ;
figure_position = [2041    -539   720   972] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = groot ;
mp = r.MonitorPositions ;
if max(mp(:,1) + mp(:,3)) < figure_position(1)
    error('This figure is going to be off the screen!')
end

% Name, code, conversion factor, formatSpec mean, formatSpec SEM, units
where2sep = [0.5 4.5 7.5 10.5] ;
sep_labels = {...
    'Exogenous forcing' ;
    'Commodity demand' ;
    'Land use areas' ;
    'Management inputs' ;
    } ;
rowInfo = { ...
           % Exogenous inputs
           'Population', 'pop', 1e-9, '%0.1f', '%0.1f', 'billion' ;
           '[CO_2]', 'co2', 1, '%0.0f', '%0.0f', 'ppm' ;
%            'Temperature', 'temp', 1, '%0.1f', '%0.1f', 'K' ;
           'Temperature', 'temp', 1, '%0.1f', '%0.1f', [char(176) 'C'] ;
%            'Precipitation', 'prec', 1, '%0.1f', '%0.1f', 'mm' ;
           'Precipitation', 'prec', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f', 'Kkm^3' ;
% Demand and production
           'Crops', 'Demand.crops', 1e-3*1e-6, '%.0f', '%.0f', 'Mt' ;
           'Ruminants', 'Demand.ruminants', 1e-6*1e-6, '%.1f', '%.1f', 'Gt' ;
           'Monogastrics', 'Demand.monogastrics', 1e-3*1e-6, '%.0f', '%.0f', 'Mt' ;
%            'Crop prod. (wt.)', 'cropprod', 1e-3*1e-6, '%.0f', '%.0f', 'Mt' ;
%            'Crop demand (kcal)', 'Demand_kcal.crops', cf_kcalEcal, '%.0f', '%.0f', 'Ecal' ;
%            'Crop prod. (cal)', 'kcal', cf_kcalEcal, '%.0f', '%.0f', 'Ecal' ;
%            'Ruminant demand', 'DemandPC.ruminants', 1, ' %.0f', '%.0f', 'kg person^{-1} yr^{-1}' ;
%            'Monogastric demand', 'DemandPC.monogastrics', 1, '%.0f', '%.0f', 'kg person^{-1} yr^{-1}' ;
%            'Crop demand', 'DemandPC.crops', 1, '%.0f', '%.0f', 'kg person^{-1} yr^{-1}' ;
%            'Crop prod.', 'kcalPC', 1/365, '%.0f', '%.0f', 'kcal person^{-1} day^{-1}' ;

           % Land use areas
           'Cropland', 'LUarea_crop', 1e-6*1e-6, '%0.1f', '%0.1f', 'Mkm^2' ;
           'Pasture', 'LUarea_past', 1e-6*1e-6, '%0.1f', '%0.1f', 'Mkm^2' ;
%            'Agriculture', 'LUarea_crop+LUarea_past', 1e-6*1e-6, '%0.1f', '%0.1f', 'Mkm^2' ;
           'Non-agri.', 'LUarea_ntrl', 1e-6*1e-6, '%.0f', '%.0f', 'Mkm^2' ;
           
           % Management inputs
           'Fertilizer', 'nflux_fert', -1e-9, '%.0f', '%.0f', 'TgN' ;
           'Irrigation', 'irrig', cf_m3_to_km3, '%.0f', '%.0f', 'km^3' ;
           } ;

       
%%%%%%%%%%%%%%%%%%%
%%% Make figure %%%
%%%%%%%%%%%%%%%%%%%

run('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/MATLAB_work/make_big_bar_graph_garr.m') ;

ht = title(sprintf('Change in land use and drivers, %d-%d to %d-%d', ...
    min(years_endh), max(years_endh), min(years_endf), max(years_endf))) ;
ht.Position = [100.0003 0.3 0] ;

%%%%%%%%%%%%
%%% Save %%%
%%%%%%%%%%%%

if do_save
%     set(gcf,'Renderer', 'Painters') % Ensures vector graphics; https://www.mathworks.com/matlabcentral/answers/2755-printing-figure-to-pdf-produces-bitmap-instead-of-vector
    export_fig([outDir_base 'bargraph.drivers.pdf'],['-r' num2str(300)])
    close
end


%% Big bar graph: Indicators

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
orientation = 'h' ; % v(ertical) or h(orizontal)
% sd_or_sem = 'st. dev.' ;
% sd_or_sem = 'SEM' ;
sd_or_sem = '' ;
% errbar_color = 'k' ;
errbar_color = 0.5*ones(3,1) ;
fontSize = 14 ;
figure_position = [1    33   720   772] ;
% figure_position = [1    33   846   772] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Name, code, conversion factor, formatSpec mean, formatSpec SEM, units
where2sep = [] ;
% where2sep = [0.5 5.5 6.5] ;
% sep_labels = {...
%     'Beneficial' ; ...
%     'Detrimental' ; ...
%     'Neutral/varying' ;
%     } ;
rowInfo = { ... 
           'Veg. C', 'cpool_VegC', cf_kg2Pg, '%.0f', '%.0f', 'GtC' ;
%            'Soil/litter C', 'cpool_LitterSoilC', cf_kg2Pg, '%d', '%d', 'GtC' ;
           'Total C', 'cpool_Total', cf_kg2Pg, '%.0f', '%.0f', 'GtC' ;
%            'Jan. albedo', 'albedo1', 1, '%0.3f', '%0.3f', '' ;
%            'Jan. albedo, borfor+tundra', 'albedo1_borfor+albedo1_tundra', 1, '%0.3f', '%0.3f', '' ;
           'Hotspot area: CI', 'hotspot_area', 1e-6*1e-6, '%0.1f', '%0.1f', 'Mkm^2' ;
           'Hotspot area: CI+CSLF', 'hotspotCSLF_area', 1e-6*1e-6, '%0.1f', '%0.1f', 'Mkm^2' ;
%            'Hotspot area: CI+CSLF (no Misc.)', 'hotspotCSLF_area_nobioenergy', 1e-6*1e-6, '%0.1f', '%0.1f', 'Mkm^2' ;
%            'ET', 'aevapaaet', cf_m3_to_km3, '%.0f', '%.0f', 'km^3' ;
           'Runoff', 'tot_runoff', cf_m3_to_km3*1e-3, '%.0f', '%.1f', 'Kkm^3' ;
           'N loss', 'nloss', cf_kg2Tg, '%.0f', '%.0f', 'TgN' ;
           'BVOC emis.', 'aiso+amon', cf_kg2Tg, '%.0f', '%.0f', 'TgC' ;
%            'Iso. emis.', 'aiso', cf_kg2Tg, '%.0f', '%.0f', 'TgC' ;
%            'Mon. emis.', 'amon', cf_kg2Tg, '%.0f', '%.0f', 'TgC' ;
           } ;

%%%%%%%%%%%%%%%%%%%
%%% Make figure %%%
%%%%%%%%%%%%%%%%%%%

run('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/MATLAB_work/make_big_bar_graph_garr.m') ;


title(sprintf('Change in ecosystem service indicators, %d-%d to %d-%d', ...
    min(years_endh), max(years_endh), min(years_endf), max(years_endf)))


%%%%%%%%%%%%
%%% Save %%%
%%%%%%%%%%%%

if do_save
    export_fig([outDir_base 'bargraph.indicators.pdf'],['-r' num2str(300)])
    close
end


%% Big table

disp('Making table...')

years_endh = 2001:2010 ;
years_endf = 2091:2100 ;

% Name, code, conversion factor, formatSpec mean, formatSpec SD
rowInfo = {'Area: Crop.', 'LUarea_crop', 1e-6*1e-6, '%0.1f', '%0.2f' ;
           'Area: Past.', 'LUarea_past', 1e-6*1e-6, '%0.1f', '%0.2f' ;
           'Area: Non-agri.', 'LUarea_ntrl', 1e-6*1e-6, '%0.1f', '%0.2f' ;
           'Vegetation C (GtC)', 'cpool_VegC', cf_kg2Pg, '%d', '%d' ;
           'Soil and litter C (GtC)', 'cpool_LitterSoilC', cf_kg2Pg, '%d', '%d' ;
%            'Product C (GtC)', 'cpool_HarvSlowC', cf_kg2Pg, '%0.1f', '%0.1f' ;
           'Total C (GtC)', 'cpool_Total', cf_kg2Pg, '%d', '%d' ;
           'January albedo', 'albedo1', 1, '%0.3f', '%0.3f' ;
           'July albedo', 'albedo7', 1, '%0.3f', '%0.3f' ;
           'Evapotranspiration (1000 km^3)', 'aevapaaet', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
           'Runoff (1000 km^3)', 'tot_runoff', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
%            'Peak monthly runoff (1000 km^3)', 'pkrunoff', cf_m3_to_km3*1e-3, '%0.1f', '%0.1f' ;
           'Crop production (Ecal)', 'kcal', cf_kcalEcal, '%0.1f', '%0.1f' ;
           'N loss (TgN)', 'nloss', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'N loss: Gaseous (TgN)', 'nflux_flux', cf_kg2Tg, '%0.1f', '%0.1f' ;
%            'N loss: Dissolved (TgN)', 'nflux_leach', cf_kg2Tg, '%0.1f', '%0.1f' ;
           'Isoprene emissions (TgC)', 'aiso', cf_kg2Tg, '%0.1f', '%0.1f' ;
           'Monoterpene emissions (TgC)', 'amon', cf_kg2Tg, '%0.1f', '%0.1f' ;
           'January albedo, boreal forest', 'albedo1_borfor', 1, '%0.3f', '%0.3f' ;
           'January albedo,  tundra', 'albedo1_tundra', 1, '%0.3f', '%0.3f' ;
           } ;

Nvars = size(rowInfo,1) ;
mean_endh_v = nan(Nvars,1) ;
mean_endf_vr = nan(Nvars,Nruns) ;
errb_endh_v = nan(Nvars,1) ;
errb_endf_vr = nan(Nvars,Nruns) ;
string_endh = cell(Nvars,1) ;
string_endf = cell(Nvars,Nruns) ;
for c = 1:Nvars
    
    % Get values
    thisVar = rowInfo{c,2} ;
    thisConv = rowInfo{c,3} ;
    mean_endh_v(c) = thisConv*eval(['mean(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
    errb_endh_v(c) = thisConv*eval(['std(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
    mean_endf_vr(c,:) = thisConv*eval(['mean(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
    errb_endf_vr(c,:) = thisConv*eval(['std(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
    
    % Turn into strings
    if strcmp(rowInfo{c,4},'%d')
        thisMean = round(mean_endh_v(c)) ;
    else
        thisMean = mean_endh_v(c) ;
    end
    if strcmp(rowInfo{c,4},'%d')
        thisSD = round(errb_endh_v(c)) ;
    else
        thisSD = errb_endh_v(c) ;
    end
    string_endh{c} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean thisSD]) ;
    for r = 1:Nruns
        if strcmp(rowInfo{c,4},'%d')
            thisMean_endf = round(mean_endf_vr(c,r)) ;
        else
            thisMean_endf = mean_endf_vr(c,r) ;
        end
        if strcmp(rowInfo{c,4},'%d')
            thisSD_endf = round(errb_endf_vr(c,r)) ;
        else
            thisSD_endf = errb_endf_vr(c,r) ;
        end
        string_endf{c,r} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean_endf thisSD_endf]) ;
    end

end

% table_out = table(collate_empties(rowInfo(:,1)),...
%                   collate_empties(string_endh)) ;
table_out = table(rowInfo(:,1),...
                  string_endh) ;
              
for r = 1:Nruns
%     table_out = [table_out collate_twocells(string_begf(:,r),string_endf(:,r))] ;
    table_out = [table_out string_endf(:,r)] ;
end
table_out.Properties.VariableNames = [{'Ecosystem_function','Baseline'} runColNames] ;

if do_save
    writetable(table_out,[outDir_base 'summary_table_nobegf.xlsx'],'Sheet',1) ;
end
disp('Done making table.')


%% Map changes in each crop yield: End-Historical to End-Future

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
conv_fact_map = cf_kgPm2_to_tonsPha ;   % m2 to km2
conv_fact_total = [] ;   % Do not show total
units_map = 'tons ha^{-1}' ;
units_total = '' ;
only1bl = true ;
thisStat = 'yield' ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;
tf = [true false] ;
for c = 1:Ncrops
    for as_pct_change = tf
        make_LUdiff_fig_v2p1_garr(...
            mean(garr_yield_d9.garr_xvyB(:,c,:),3), squeeze(mean(garr_yield_d9.garr_xvyr(:,c,:,:),3)), ...
            thisY1, thisYN, garr_yield_d9.varNames{c}, runList, ...
            spacing, fontSize, textX, textY_1, textY_2, ...
            nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
            Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps, ...
            thisStat, as_pct_change, map_size, list2map) ;
                stop
        if do_save
            if as_pct_change
                thisForm = '%s/yieldPctDiff_%ds-%ds_%s.png' ;
            else
                thisForm = '%s/yieldDiff_%ds-%ds_%s.png' ;
            end
            this_filename = sprintf(thisForm, ...
                removeslashifneeded(outDir_maps), thisY1, thisYN, garr_yield_d9.varNames{c}) ;
            export_fig(this_filename, ['-r' num2str(pngres)])
            close
        end
    end
end


%% Map changes in each crop fertilizer application: End-Historical to End-Future

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
conv_fact_map = 1e4 ;
conv_fact_total = [] ;   % Do not show total
units_map = 'kg N ha^{-1}' ;
units_total = '' ;
only1bl = true ;
thisStat = 'fert. appl.' ;
as_pct_change = true ;
%%%%%%%%%%%%%%%%%%%

garr_tmp_d9_xvyB = garr_Nfert_d9.garr_xvyB .* garr_cropareas_d9.garr_xvyB ;
garr_tmp_d9_xvyr = garr_Nfert_d9.garr_xvyr .* garr_cropareas_d9.garr_xvyr ;

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;
for c = 1:Ncrops
    make_LUdiff_fig_v2p1_garr(...
        mean(garr_tmp_d9_xvyB(:,c,:),3), squeeze(mean(garr_tmp_d9_xvyr(:,c,:,:),3)), ...
        thisY1, thisYN, garr_Nfert_d9.varNames{c}, runList, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
        Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps, ...
        thisStat, as_pct_change, map_size, list2map) ;
    stop
    if do_save
        if as_pct_change
            thisForm = '%s/NfertTotalPctDiff_%ds-%ds_%s.png' ;
        else
            thisForm = '%s/NfertTotalDiff_%ds-%ds_%s.png' ;
        end
        this_filename = sprintf(thisForm, ...
                removeslashifneeded(outDir_maps), thisY1, thisYN, garr_Nfert_d9.varNames{c}) ;
        export_fig(this_filename, ['-r' num2str(pngres)])
        close
    end
end


%% Map changes in each crop irrigation use: End-Historical to End-Future

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
% conv_fact_map = cf_m3_to_km3 ;
% conv_fact_total = [] ;   % Do not show total
% units_map = 'km^{-3}' ;
% units_total = '' ;
% only1bl = true ;
% thisStat = 'irrig.' ;
% as_pct_change = false ;
% %%%%%%%%%%%%%%%%%%%
% 
% error('There''s something wrong with gsirr maps! No baseline rice irrigation, for example')
% warning('Not tested after maps->garr conversion.')
% 
% [~,IA] = intersect(garr_cropareas_d9.varNames, garr_gsirrig_d9.varNames, 'stable') ;
% garr_tmp_d9_xvyB = garr_gsirrig_d9.garr_xvyB(:,IA,:) .* garr_cropareas_d9.garr_xvyB ;
% garr_tmp_d9_xvyr = garr_gsirrig_d9.garr_xvyr(:,IA,:,:) .* garr_cropareas_d9.garr_xvyr ;
% 
% thisY1 = yearList_baseline(end) ;
% thisYN = yearList_future(end) ;
% for c = 1:Ncrops
%     make_LUdiff_fig_v2p1_garr(...
%         mean(garr_tmp_d9_xvyB(:,c,:),3), squeeze(mean(garr_tmp_d9_xvyr(:,c,:,:),3)), ...
%         thisY1, thisYN, garr_gsirrig_d9.varNames{c}, runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps, ...
%         thisStat, as_pct_change, map_size, list2map) ;
%     
%     if do_save
%         if as_pct_change
%             thisForm = '%s/IrrigTotalPctDiff_%ds-%ds_%s.png' ;
%         else
%             thisForm = '%s/IrrigTotalDiff_%ds-%ds_%s.png' ;
%         end
%         this_filename = sprintf(thisForm, ...
%                 removeslashifneeded(outDir_maps), thisY1, thisYN, garr_Nfert_d9.varNames{c}) ;
%         export_fig(this_filename, ['-r' num2str(pngres)])
%         close
%     end
% end


%% Map areas with changing drought/flood risk
% Based on Asadieh & Krakauer (2017). Their DROUGHT RISK indicator was the
% 5th percentile of ANNUAL streamflow, whereas their FLOODING indicator was
% the 95th percentile of DAILY streamflow. I don't have daily runoff data,
% so I'll have to use MONTHLY for FLOOD RISK.

% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_aggregate_basins = false ;
spacing = [0.03, 0.03] ;   % v, h
norm_ticks = [1.5 2 3 5 10 Inf] ; % Tick marks: Multiply q20 by this to get q21
fontSize = 14 ;
fontSize_text = 14 ;
Ystart = 69 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup
file_suffix = '' ;
if do_aggregate_basins
    file_suffix = '.basins' ;
end

do_norm = false ;

% Change in 10-year lowest annual runoff
pctDiff_xr = make_runoffFigs_Asadieh_garr( ...
    garr_mon_runoff_last30, garr_awater_last30, runList, 'drought', land_area_unmasked_weights_x, ...
    do_norm, spacing, norm_ticks, fontSize, fontSize_text, Ystart, ...
    basins_x, do_aggregate_basins, map_size, list2map) ;
if do_save
    thisFile = sprintf('%s/pkRunoff_drought_last3decs_20th-21st%s.png', ...
        outDir_maps, file_suffix) ;
    export_fig(thisFile, ['-r' num2str(pngres)])
    for r = 1:Nruns
        thisYX = lpjgu_vector2map(pctDiff_xr(:,r), map_size, list2map) ;
        thisFile = sprintf('%s/pkRunoff_drought_last3decs_20th-21st.%s%s.tif', ...
            outDir_gtif, runList{r}, file_suffix) ;
        geotiffwrite_ssr(thisFile, thisYX, R, gtif_missing) ;
    end
    close
end

% Map land use contribution to 10-year lowest annual runoff (percentage points)
%%% Calculated as fully-varying minus constLU. This means that LU legacy
%%% effects are NOT included in these numbers---only NEW land use changes.
if contains(thisVer, 'attr')
     LUCcont_drought_YX = ...
        map_contribution_garr(pctDiff_xr, continents_shp, ...
        {runList{1}, 'constLU'}, runList, ...
        'thisTitle', 'LUC contribution to declining 10-year lowest annual runoff', ...
        'units_map', 'Percentage points', ...
        'nanmask_x', pctDiff_xr(:,1)>0, ...
        'latlim', [-60 80], ...
        'thisColormap', 'rdbu_ssr', ...
        'fontSize', 14, ...
        'edgeColor', 0.6*ones(3,1), ...
        'lineWidth', 1, ...
        'cbarOrient', 'SouthOutside', ...
        'caxis_lims', 100*[-1 1]) ;
    if do_save
        thisFile = sprintf('%s/pkRunoff_drought_last3decs_20th-21st.LUCcont%s.png', ...
            outDir_maps, file_suffix) ;
        export_fig(thisFile, ['-r' num2str(pngres)])
        thisFile = sprintf('%s/pkRunoff_drought_last3decs_20th-21st.LUCcont%s.tif', ...
            outDir_gtif, file_suffix) ;
        geotiffwrite_ssr(thisFile, LUCcont_drought_YX, R, gtif_missing) ;
        close
    end
end

% Change in 10-year highest monthly runoff
pctDiff_xr = make_runoffFigs_Asadieh_garr( ...
    garr_mon_runoff_last30, garr_awater_last30, runList, 'flood', land_area_unmasked_weights_x, ...
    do_norm, spacing, norm_ticks, fontSize, fontSize_text, Ystart, ...
    basins_x, do_aggregate_basins, map_size, list2map) ;
if do_save
    thisFile = sprintf('%s/pkRunoff_flood_last3decs_20th-21st%s.png', ...
        outDir_maps, file_suffix) ;
    export_fig(thisFile,['-r' num2str(pngres)])
    for r = 1:Nruns
        thisYX = lpjgu_vector2map(pctDiff_xr(:,r), map_size, list2map) ;
        thisFile = sprintf('%s/pkRunoff_flood_last3decs_20th-21st.%s%s.tif', ...
            outDir_gtif, runList{r}, file_suffix) ;
        geotiffwrite_ssr(thisFile, thisYX, R, gtif_missing) ;
    end
    close
end

% Map land use contribution to 10-year highest monthly runoff (percentage points)
%%% Calculated as fully-varying minus constLU. This means that LU legacy
%%% effects are NOT included in these numbers---only NEW land use changes.
if contains(thisVer, 'attr')
    warning('Not tested after maps->garr conversion')
    LUCcont_flood_YX = ...
        map_contribution_garr(pctDiff_xr, continents_shp, ...
        {runList{1}, 'constLU'}, runList, ...
        'thisTitle', 'LUC contribution to increasing 10-year highest monthly runoff', ...
        'units_map', 'Percentage points', ...
        'nanmask_x', pctDiff_xr(:,1)<0, ...
        'latlim', [-60 80], ...
        'thisColormap', 'rdbu_ssr', ...
        'fontSize', 14, ...
        'edgeColor', 0.6*ones(3,1), ...
        'lineWidth', 1, ...
        'cbarOrient', 'SouthOutside', ...
        'caxis_lims', 100*[-1 1]) ;
    if do_save
        thisFile = sprintf('%s/pkRunoff_flood_last3decs_20th-21st.LUCcont%s.png', ...
            outDir_maps, file_suffix) ;
        export_fig(thisFile,['-r' num2str(pngres)])
        thisFile = sprintf('%s/pkRunoff_flood_last3decs_20th-21st.LUCcont%s.tif', ...
            outDir_gtif, file_suffix) ;
        geotiffwrite_ssr(thisFile, LUCcont_flood_YX, R, gtif_missing) ;
        close
    end
end

do_norm = true ;

% Normalized change in 10-year lowest annual runoff
make_runoffFigs_Asadieh_garr( ...
    garr_mon_runoff_last30, garr_awater_last30, runList, 'drought', land_area_unmasked_weights_x, ...
    do_norm, spacing, norm_ticks, fontSize, fontSize_text, Ystart, ...
    basins_x, do_aggregate_basins, map_size, list2map) ;
if do_save
    thisFile = sprintf('%s/pkRunoff_drought_last3decs_20th-21st.norm%s.png', ...
        outDir_maps, file_suffix) ;
    export_fig(thisFile,['-r' num2str(pngres)])
    close
end

% Normalized change in 10-year highest monthly runoff
make_runoffFigs_Asadieh_garr( ...
    garr_mon_runoff_last30, garr_awater_last30, runList, 'flood', land_area_unmasked_weights_x, ...
    do_norm, spacing, norm_ticks, fontSize, fontSize_text, Ystart, ...
    basins_x, do_aggregate_basins, map_size, list2map) ;
if do_save
    thisFile = sprintf('%s/pkRunoff_flood_last3decs_20th-21st.norm%s.png', ...
        outDir_maps, file_suffix) ;
    export_fig(thisFile,['-r' num2str(pngres)])
    close
end


%% Map differences from end of historical to end of future: Total N loss

% Options %%%%%%%%%
filename_base = [outDir_maps 'nloss'] ;
title_text = 'total N loss' ;
sumvars = {'flux','leach'} ;
pct_clim = 100 ;
equalize_cbars = true ;
% conv_fact_map = 1 ;   % kgN/m2
% units_map = 'kgN m^{-2} yr ^{-1}' ;
conv_fact_map = gcel_area_YX*cf_kg2Tg*1e3 ;   % kgN/m2 to GgN
units_map = 'GgN yr ^{-1}' ;
conv_fact_total = cf_kg2Tg ;   % kgN to TgN
units_total = 'TgN yr ^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
prctile_clim = 99 ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, garr_nflux_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;


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
prctile_clim = 99 ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, garr_nflux_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;


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
prctile_clim = 99 ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, garr_nflux_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;


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
prctile_clim = [] ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, garr_aiso_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;


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
prctile_clim = [] ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, garr_amon_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;


%% Map differences from end of historical to end of future: January albedo

% Options %%%%%%%%%
filename_base = [outDir_maps 'albedo_jan'] ;
title_text = 'January albedo' ;
sumvars = 'January' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % unitless
units_map = 'unitless' ;
conv_fact_total = 1 ;   % unitless
units_total = '' ;
this_land_area_map = lpjgu_vector2map( ...
    vegd_area_YXmean/nansum(nansum(vegd_area_YXmean)), ...
    map_size, list2map) ; % Set to [] if not needed to calculate total
prctile_clim = [] ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, garr_albedo_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;


%% Map differences from end of historical to end of future: July albedo

% Options %%%%%%%%%
filename_base = [outDir_maps 'albedo_jul'] ;
title_text = 'July albedo' ;
sumvars = 'July' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1 ;   % unitless
units_map = 'unitless' ;
conv_fact_total = 1 ;   % unitless
units_total = '' ;
this_land_area_map = lpjgu_vector2map( ...
    vegd_area_YXmean/nansum(nansum(vegd_area_YXmean)), ...
    map_size, list2map) ; % Set to [] if not needed to calculate total
prctile_clim = [] ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, garr_albedo_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;


%% Map differences from end of historical to end of future: Vegetation C

% Options %%%%%%%%%
filename_base = [outDir_maps 'cpool_veg'] ;
title_text = 'vegetation C' ;
sumvars = 'VegC' ;
conv_fact_map = 1e-3*1e4 ;   % kgC/m2 to tonsC/ha
units_map = 'tons C ha^{-1}' ;
conv_fact_total = cf_kg2Pg ;   % kgC to PgC
units_total = 'PgC' ;
precision_total = 1 ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
%%%%%%%%%%%%%%%%%%%

% pct_clim = 100 ;
% equalize_cbars = true ;
% prctile_clim = [] ;
% textX = 25/720 ; textY_1 = 50/360 ; textY_2 = 20/360 ; 
% shiftup = 10/360 ; textY_1 = textY_1 + shiftup ; textY_2 = textY_2 + shiftup ; 
% % fontSize = 14 ;
% fontSize = 16 ;
% colorBarLoc = 'SouthOutside' ; 
% nx = 2 ; ny = 2 ;
% % nx = 1 ; ny = 4 ;
% % thisPos = figurePos ;
% thisPos = [1         130        1320         675] ;
% % spacing = [0.1 0.05] ;% [v h]
% spacing = [0.05 0.05] ;% [v h]
% do_map_run_diffs_fromEndHist(do_save, garr_cpool_d9, sumvars, title_text, filename_base, ...
%     equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
%     thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
%     conv_fact_map, units_map, conv_fact_total, units_total, ...
%     pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;

thisPos = [1    33   407   772] ;
spacing = [0.05 0.05] ; % [v h]
fontSize = 11 ;
% textX = 0 ; textY_1 = 50/360 ; textY_2 = 20/360 ;
textX = 80/720 ; textY_1 = 90/360 ; textY_2 = 20/360 ;
bins_lowBnds = [-250:50:-50 -5 5 50:50:200] ;
map_run_diffs_fromEndHist_oneCol_v2(garr_cpool_d9, title_text, sumvars, ...
    fontSize, spacing, textX, textY_1, textY_2, ...
    thisPos, runList_titles, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    bins_lowBnds, 'BrBG', ...
    map_size, list2map, precision_total) ;

if do_save
    filename = [filename_base '_diff'] ;
    if exist('prctile_clim', 'var') && ~isempty(prctile_clim)
        if isint(prctile_clim)
            filename = sprintf('%s_limPrctile%df', filename, prctile_clim) ;
        else
            filename = sprintf('%s_limPrctile%0.1f', filename, prctile_clim) ;
        end
    end
    filename = [filename '_2000s-2090s_1col.png'] ;
    export_fig(filename,['-r' num2str(pngres)])
    close
end


%% Map differences from end of historical to end of future: Total C

% Options %%%%%%%%%
filename_base = [outDir_maps 'cpool_tot'] ;
title_text = 'total C' ;
sumvars = 'Total' ;
pct_clim = 100 ;
equalize_cbars = true ;
conv_fact_map = 1e-3*1e4 ;   % kgC/m2 to tonsC/ha
units_map = 'tons C ha^{-1}' ;
conv_fact_total = cf_kg2Pg ;   % kgC to PgC
units_total = 'PgC' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
prctile_clim = [] ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 16 ;
thisPos = [233 33 1208 772] ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.025] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, garr_cpool_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;


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
prctile_clim = [] ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

if all_figs
    do_map_run_diffs_fromEndHist(do_save, garr_awater_d9, sumvars, title_text, filename_base, ...
        equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
        thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
        conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;
end


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
units_total = 'Kkm^3 yr^{-1}' ;
this_land_area_map = gcel_area_YX ; % Set to [] if not needed to calculate total
prctile_clim = [] ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

do_map_run_diffs_fromEndHist(do_save, garr_awater_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;


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
prctile_clim = [] ;
%%%%%%%%%%%%%%%%%%%
textX = 25 ; textY_1 = 50 ; textY_2 = 20 ; fontSize = 14 ;
thisPos = figurePos ; colorBarLoc = 'SouthOutside' ; nx = 2 ; ny = 2 ;
spacing = [0.1 0.05] ;% [v h]
%%%%%%%%%%%%%%%%%%%

if all_figs 
    if exist('maps_pk_runoff_d9', 'var')
        do_map_run_diffs_fromEndHist(do_save, garr_pk_runoff_d9, sumvars, title_text, filename_base, ...
            equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
            thisPos, nx, ny, colorBarLoc, runList_titles, do_caps, this_land_area_map, ...
            conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, map_size, list2map) ;
    else
        warning('maps_pk_runoff_d9 does not exist (maybe just last30?, so skipping peak runoff maps')
    end
end

%% Map AREA changes in cropland and pasture area: End-historical to End-Future

% Options %%%%%%%%%
fontSize = 14 ;
spacing = [0.02 0.02] ;   % [vert, horz]
% textX = 25 ;
textX = 45 ;
textY_1 = 50 ;
textY_2 = 20 ;
% shiftup = 0 ; textY_1 = textY_1 + shiftup ; textY_2 = textY_2 + shiftup - shiftup/3 ; 
shiftup = 15 ; textY_1 = textY_1 + shiftup ; textY_2 = textY_2 + shiftup - shiftup/3 ; 
thisPos = [1    33   770   772] ;
nx = 2 ;
ny = 4 ;
colorBarLoc = 'EastOutside' ;
as_frac_land = false ;
only1bl = true ;
same_caxis = true ;
Nbins = 11 ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;

crop_area_YXB = lpjgu_vector2map(crop_area_xBH, map_size, list2map) ;
past_area_YXB = lpjgu_vector2map(past_area_xBH, map_size, list2map) ;
crop_diff_YXrH = nan([map_size Nruns]) ;
past_diff_YXrH = nan([map_size Nruns]) ;
for r = 1:Nruns
    crop_diff_YXrH(:,:,r) = lpjgu_vector2map(crop_diff_xrH(:,r), map_size, list2map) ;
    past_diff_YXrH(:,:,r) = lpjgu_vector2map(past_diff_xrH(:,r), map_size, list2map) ;
end

if as_frac_land
    crop_area_YXB = crop_area_YXB ./ gcel_area_YX ;
    past_area_YXB = past_area_YXB ./ gcel_area_YX ;
    crop_diff_YXrH = crop_diff_YXrH ./ repmat(gcel_area_YX, [1 1 Nruns]) ;
    past_diff_YXrH = past_diff_YXrH ./ repmat(gcel_area_YX, [1 1 Nruns]) ;
    conv_fact_map = 100 ;
    conv_fact_total = 0 ;
    units_map = '%' ;
    units_total = '???' ;
else
    conv_fact_map = 1e-6 ;   % m2 to km2
    conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
    units_map = 'km^2' ;
    units_total = 'Mkm^2' ;
end


[diff_crop_YXr, diff_past_YXr] = make_LUdiff_fig_v4(...
    crop_area_YXB, ...
    past_area_YXB, ...
    crop_diff_YXrH, past_diff_YXrH, ...
    thisY1, thisYN, 'Cropland', runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
    Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps, ...
    same_caxis, Nbins) ;
if do_save
    filename = [outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_croppast.png'] ;
    export_fig(filename,['-r' num2str(pngres)])
    diff_agri_YXr = diff_crop_YXr + diff_past_YXr ;
    filename_crop = strrep(filename,'croppast','crop') ;
    filename_past = strrep(filename,'croppast','past') ;
    filename_agri = strrep(filename,'croppast','agri') ;
    save_geotiffs(diff_crop_YXr, filename_crop, runList, R, gtif_missing) ;
    save_geotiffs(diff_past_YXr, filename_past, runList, R, gtif_missing) ;
    save_geotiffs(diff_agri_YXr, filename_agri, runList, R, gtif_missing) ;
    clear diff_agri_YXr filename*
    close
end
clear diff_crop_YXr diff_past_YXr


%% Map %CELL changes in cropland and pasture area: End-historical to End-Future

% Options %%%%%%%%%
fontSize = 14 ;
% spacing = [0.02 0.02] - 0.0025*8 ;   % [vert, horz]
spacing = 0 ;
textX = 0.115 ;
textY_1 = 50/360 ;
textY_2 = 20/360 ;
% shiftup = 0 ; textY_1 = textY_1 + shiftup ; textY_2 = textY_2 + shiftup - shiftup/3 ; 
shiftup = 15/360 ; textY_1 = textY_1 + shiftup ; textY_2 = textY_2 + shiftup - shiftup/3 ; 
thisPos = [1    33   770   772] ;
nx = 2 ;
ny = 4 ;
as_frac_land = true ;
bins_lowBnds = [-100:20:-20 -3 3 20:20:80] ;
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
this_colormap_name = '-PiYG_ssr' ;
lines_overlay = 'landareas.shp' ;
% lines_overlay = '/Users/sam/Geodata/General/continents_from_countries/continents_from_countries.shp' ;
%%%%%%%%%%%%%%%%%%%

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;
col_titles = { ...
    sprintf('%s cropland area, %d%s%d','\Delta',thisY1,char(8211),thisYN), ...
    sprintf('%s pasture area, %d%s%d','\Delta',thisY1,char(8211),thisYN) ...
    } ;

area_crop_bl = sum(crop_area_xBH)*conv_fact_total ;
area_past_bl = sum(past_area_xBH)*conv_fact_total ;
total_cropDiff_r = sum(crop_diff_xrH, 1)*conv_fact_total ;
total_pastDiff_r = sum(past_diff_xrH, 1)*conv_fact_total ;

% Get difference (%)
crop_diff_YXrH = nan([map_size Nruns]) ;
past_diff_YXrH = nan([map_size Nruns]) ;
for r = 1:Nruns
    crop_diff_YXrH(:,:,r) = lpjgu_vector2map(100*crop_diff_xrH(:,r)./gcel_area_x, map_size, list2map) ;
    past_diff_YXrH(:,:,r) = lpjgu_vector2map(100*past_diff_xrH(:,r)./gcel_area_x, map_size, list2map) ;
end

units_map = '%' ;
units_total = 'Mkm^2' ;

if ~as_frac_land
    error('This only works with as_frac_land TRUE')
end
[diff_crop_YXr, diff_past_YXr] = make_LUdiff_fig_v5(...
    area_crop_bl, area_past_bl, total_cropDiff_r, total_pastDiff_r, ...
    crop_diff_YXrH, past_diff_YXrH, ...
    thisY1, thisYN, runList, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, ...
    Nruns, thisPos, units_map, units_total, do_caps, ...
    bins_lowBnds, col_titles) ;
if do_save
    filename = [outDir_maps 'areaPctDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_croppast.png'] ;
    export_fig(filename,['-r' num2str(pngres)])
    diff_agri_YXr = diff_crop_YXr + diff_past_YXr ;
    filename_crop = strrep(filename,'croppast','crop') ;
    filename_past = strrep(filename,'croppast','past') ;
    filename_agri = strrep(filename,'croppast','agri') ;
    save_geotiffs(diff_crop_YXr, filename_crop, runList, R, gtif_missing) ;
    save_geotiffs(diff_past_YXr, filename_past, runList, R, gtif_missing) ;
    save_geotiffs(diff_agri_YXr, filename_agri, runList, R, gtif_missing) ;
    clear diff_agri_YXr filename*
    close
end
clear diff_crop_YXr diff_past_YXr


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

% this_ssp_plot_index = ssp_plot_index ;
this_ssp_plot_index = sort(ssp_plot_index) ;

thisY1 = yearList_baseline(end) ;
thisYN = yearList_future(end) ;
for c = 1:Ncrops
    make_LUdiff_fig_v2(...
        lpjgu_vector2map(garr_cropareas_xvBH(:,c), map_size, list2map), ...
        lpjgu_xz_to_YXz(squeeze(garr_cropareasDiffs_xvrH(:,c,:)), map_size, list2map), ...
        thisY1, thisYN, CFTnames_garr{c}, runList, ...
        spacing, fontSize, textX, textY_1, textY_2, ...
        nx, ny, 1, colorBarLoc, this_ssp_plot_index, only1bl, ...
        Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
    stop
    if do_save
        export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_garr{c} '.png'],['-r' num2str(pngres)])
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
%         thisY1, thisYN, CFTnames_garr{c}, runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_' CFTnames_garr{c} '.png'],['-r' num2str(pngres)])
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
% Nsmth = 1 ;
Nsmth = 5 ; warning('SMOOTHING WINDOW >1 (%d)', Nsmth)
spacing = [0.15 0.1] ;   % vert, horiz
blYears = yearList_baseline ;
% blYears = 1960:yearList_baseline(end) ;
perArea = true ;
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
    skip3rdColor) ;

clear tmp_*

set(gca,'XLim',[1950 2100])

%%%%%%%%%%%%%%%%%%%
% Plot timeseries: Irrigation
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
    skip3rdColor) ;

clear tmp_*

set(gca,'XLim',[1950 2100])

if perArea
    file_suffix = [file_suffix '_perArea'] ;
end

if do_save
    export_fig([outDir_ts 'mgmt_inputs' file_suffix '.1950-2100.pdf'])
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
if exist('crop0_area_xBH','var')
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


%% Plot timeseries: Irrigation on each crop

% Options %%%%%%%%%
thisVar = 'gsirrig' ;
title_prefix = 'Irrigation' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = 1e-3*1e-6 ;   % kg to Mt
units = 'units?' ;
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

% Options %%%%%%%%%
% plum_area_adjustment = 1 ;
% plum_area_adjustment = 1-0.28 ;
plum_area_adjustment = 1 - unhCropFrac ;
lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;
thisVar = 'cropprodExp' ;
title_prefix = 'Production (exp.)' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = cf_kg2Mt ;
units = 'Mt DM' ;
%%%%%%%%%%%%%%%%%%%

if all_figs && false
    
    theseCFTnames = CFTnames ;
    thisLegend = stdLegend_plusFAO ;
    
    clear cell_bl cell_yr
    cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    
    file_suffix = CFT_timeseries_plot(...
        cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, fao.tmp_fao_yearList, rebase, ...
        theseCFTnames, units, title_prefix, thisLegend, do_caps, skip3rdColor, ...
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
    
end


%% Plot timeseries: Production of each crop, ACTUAL/EXPECTED

% Options %%%%%%%%%
% title_prefix = 'Prod. diff: (Act.-Exp./Exp.)' ;
title_prefix = 'Prod. diff' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
Nsmth = 1 ;
units = '%' ;
%%%%%%%%%%%%%%%%%%%

if all_figs && false
    
    theseCFTnames = CFTnames ;
    thisLegend = stdLegend(2:end) ;
    
    clear cell_yr cell_yr_act cell_yr_exp
    cmds = get_cell_forPlot(whos, 'cropprod', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    cell_yr_act = cell_yr ; clear cell_yr
    cell_yr_act = adj_yield_tech(cell_yr_act, yearList_future) ;
    cmds = get_cell_forPlot(whos, 'cropprodExp', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    cell_yr_exp = cell_yr ; clear cell_yr
    CFT_ActVsExp_plot(...
        cell_yr_act, cell_yr_exp, yearList_future, ...
        theseCFTnames, units, title_prefix, thisLegend, do_caps, ...
        'lineWidth',lineWidth, ...
        'fontSize',fontSize, ...
        'spacing',spacing, ...
        'Nsmth',Nsmth) ;
    
    if do_save
        export_fig([outDir_ts 'cropProdDiffFromExp.pdf'])
        close
    end
    
end


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

if all_figs
    
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
    
end


%% Plot timeseries: Yield of each crop (EXPECTED)

% Options %%%%%%%%%
thisVar = 'yieldExp' ;
title_prefix = 'Yield (PLUM-exp)' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = cf_kgPm2_to_tonsPha ;
units = 't ha^{-1}' ;
%%%%%%%%%%%%%%%%%%%

if all_figs && false
    
    theseCFTnames = CFTnames ;
    if include_fao
        thisLegend = stdLegend_plusFAO ;
    else
        thisLegend = stdLegend ;
    end
    
    clear cell_bl cell_yr
    cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'bl', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    if include_fao
        cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    else
        cell_fao = {} ;
        fao.tmp_fao_yearList = [] ;
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
    
    if do_save
        export_fig([outDir_ts thisVar file_suffix '.pdf'])
        close
    end

end


%% Plot timeseries: Yield of each crop, ACTUAL/EXPECTED

% Options %%%%%%%%%
title_prefix = 'Yield diff' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
Nsmth = 1 ;
units = '%' ;
%%%%%%%%%%%%%%%%%%%

if all_figs && false
    
    theseCFTnames = CFTnames ;
    thisLegend = stdLegend(2:end) ;
    
    clear cell_yr cell_yr_act cell_yr_exp
    cmds = get_cell_forPlot(whos, 'yield', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    cell_yr_act = cell_yr ; clear cell_yr
    cmds = get_cell_forPlot(whos, 'yieldExp', 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
    cell_yr_exp = cell_yr ; clear cell_yr
    
    % Adjust for technology change, if doing so
    if do_adjYieldTech
        file_suffix = [file_suffix '_techAdj'] ;
    end
    
    CFT_ActVsExp_plot(...
        cell_yr_act, cell_yr_exp, yearList_future, ...
        theseCFTnames, units, title_prefix, thisLegend, do_caps, ...
        'lineWidth',lineWidth, ...
        'fontSize',fontSize, ...
        'spacing',spacing, ...
        'Nsmth',Nsmth) ;
    
    if do_save
        export_fig([outDir_ts 'yieldDiffFromExp' file_suffix '.pdf'])
        close
    end

end


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

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
ignYrs = 0 ;
Nsmth = 5 ;
spacing = [0.1 0.1] ;   % [vert, horz]
%%%%%%%%%%%%%%%%%%%

if all_figs
    
    % % ts_nloss_bl = ts_nflux_flux_bl + ts_nflux_harvest_bl + ts_nflux_leach_bl + ts_nflux_LUch_bl ;
    % % ts_nloss_yr = ts_nflux_flux_yr + ts_nflux_harvest_yr + ts_nflux_leach_yr + ts_nflux_LUch_yr ;
    % ts_nloss_bl = ts_nflux_flux_bl + ts_nflux_leach_bl ;
    % ts_nloss_yr = ts_nflux_flux_yr + ts_nflux_leach_yr ;
    %
    [tmp_ts_nflux_flux_yr, title_suffix, file_suffix] = ...
        rebase_future2baseline(rebase, Nsmth, ts_nflux_flux_bl, ts_nflux_flux_yr, ignYrs, yearList_future) ;
    [tmp_ts_nflux_leach_yr, ~, ~] = ...
        rebase_future2baseline(rebase, Nsmth, ts_nflux_leach_bl, ts_nflux_leach_yr, ignYrs, yearList_future) ;
    
    figure('Position',figurePos,'Color','w') ;
    
    subplot_tight(1,2,1,spacing)
    ht = plot_timeseries(...
        yearList_baseline, [], yearList_future, ...
        ts_nflux_flux_bl, [], tmp_ts_nflux_flux_yr, ...
        cf_kg2Tg, Nsmth, ...
        'TgN', stdLegend, 'N loss: Gaseous', '', lineWidth, fontSize, ...
        skip3rdColor) ;
    letterlabel_align0('A',ht,do_caps) ;
    
    subplot_tight(1,2,2,spacing)
    ht = plot_timeseries(...
        yearList_baseline, [], yearList_future, ...
        ts_nflux_leach_bl, [], tmp_ts_nflux_leach_yr, ...
        cf_kg2Tg, Nsmth, ...
        'TgN', stdLegend, 'N loss: Dissolved', '', lineWidth, fontSize, ...
        skip3rdColor) ;
    letterlabel_align0('B',ht,do_caps) ;
    
    clear tmp_*
    
    if do_save
        export_fig([outDir_ts 'Nloss' file_suffix '.pdf'])
        close
    end

end


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
thisCF = cf_kg2Tg ;
units = 'TgC' ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_amon_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_amon_bl, ts_amon_yr, ignYrs, yearList_future) ;

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing) ;
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_amon_bl, [], tmp_ts_amon_yr, ...
    thisCF, Nsmth, ...
    units, stdLegend, 'Monoterpene emissions', title_suffix, lineWidth, fontSize, ...
    skip3rdColor) ;
letterlabel_align0('A',ht,do_caps) ;
clear tmp_*

[tmp_ts_aiso_yr, title_suffix, ~] = ...
    rebase_future2baseline(rebase, Nsmth, ts_aiso_bl, ts_aiso_yr, ignYrs, yearList_future) ;

subplot_tight(1,2,2,spacing) ;
ht = plot_timeseries(...
    yearList_baseline, [], yearList_future, ...
    ts_aiso_bl, [], tmp_ts_aiso_yr, ...
    thisCF, Nsmth, ...
    units, stdLegend, 'Isoprene emissions', title_suffix, lineWidth, fontSize, ...
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


%% Plot timeseries: Production (weight)

% Options %%%%%%%%%
plum_area_adjustment = 1 ;
% plum_area_adjustment = 1-0.28 ;

lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;

lineWidth = 2 ;
fontSize = 22 ;
ignYrs = 0 ;
Nsmth = 1 ;
figurePosition = [1 376 1440 429] ;
%%%%%%%%%%%%%%%%%%%

[tmp_ts_cropprod_yr, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_cropprod_bl, ts_cropprod_yr, ignYrs, yearList_future) ;

% Adjust for technology change, if doing so (do not include FAO!)
tmp_ts_cropprod_bl = ts_cropprod_bl ;
if do_adjYieldTech
    file_suffix = [file_suffix '_techAdj'] ;
%     title_suffix = [title_suffix ' (techAdj)'] ;
end

% Adjust for unhandled crops, if doing so (DO include FAO)
tmp_ts_cropprod_fao = ts_cropprod_fao ;
if ~isequal(lpjg_area_adjustment,1)
    tmp_ts_cropprod_bl = tmp_ts_cropprod_bl .* lpjg_area_adjustment ;
    [~,IA] = intersect(yearList_baseline,fao.tmp_fao_yearList) ;
    tmp_ts_cropprod_fao = tmp_ts_cropprod_fao .* lpjg_area_adjustment(IA) ;
    title_suffix = [title_suffix ' (blAdj)'] ;
    file_suffix = [file_suffix '_blAdj'] ;
end
if plum_area_adjustment ~= 1
    tmp_ts_cropprod_yr = tmp_ts_cropprod_yr * plum_area_adjustment ;
    title_suffix = [title_suffix ' (PLUM\times' num2str(plum_area_adjustment) ')'] ;
    file_suffix = [file_suffix '_plumAdj' num2str(plum_area_adjustment)] ;
end

figure('Position',figurePosition,'Color','w') ;

plot_timeseries(...
    yearList_baseline, fao.tmp_fao_yearList, yearList_future, ...
    tmp_ts_cropprod_bl, tmp_ts_cropprod_fao, tmp_ts_cropprod_yr, ...
    cf_kg2Pg, Nsmth, ...
    'Dry matter (Gt)', stdLegend_plusFAO, 'Crop production', ...
    title_suffix, lineWidth, fontSize, skip3rdColor) ;
clear tmp_*

if do_save
    export_fig([outDir_ts 'cropprod_total' file_suffix '.pdf'])
    close
end



%% Plot timeseries: Calories

% Options %%%%%%%%%
plum_area_adjustment = 1 ;
% plum_area_adjustment = 1-0.28 ;

lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;

lineWidth = 2 ;
fontSize = 22 ;
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
%     title_suffix = [title_suffix ' (techAdj)'] ;
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


%% Stop

stop


%%

mean_2010s_vr = squeeze(mean(ts_commodDemand_yvr(yearList_future>=2011 & yearList_future<=2020,:,:),1)) ;
mean_2050s_vr = squeeze(mean(ts_commodDemand_yvr(yearList_future>=2051 & yearList_future<=2060,:,:),1)) ;
mean_2090s_vr = squeeze(mean(ts_commodDemand_yvr(yearList_future>=2091 & yearList_future<=2100,:,:),1)) ;

disp('Demand')
(mean_2050s_vr - mean_2010s_vr) ./ mean_2010s_vr
(mean_2090s_vr - mean_2050s_vr) ./ mean_2050s_vr
(mean_2090s_vr - mean_2010s_vr) ./ mean_2010s_vr


mean_2010s_r = squeeze(mean(ts_nflux_fert_yr(yearList_future>=2011 & yearList_future<=2020,:),1)) ;
mean_2050s_r = squeeze(mean(ts_nflux_fert_yr(yearList_future>=2051 & yearList_future<=2060,:),1)) ;
mean_2090s_r = squeeze(mean(ts_nflux_fert_yr(yearList_future>=2091 & yearList_future<=2100,:),1)) ;

disp('Nfert')
(mean_2050s_r - mean_2010s_r) ./ mean_2010s_r
(mean_2090s_r - mean_2050s_r) ./ mean_2050s_r
(mean_2090s_r - mean_2010s_r) ./ mean_2010s_r


mean_2010s_r = squeeze(mean(ts_LUarea_crop_yr(yearList_future>=2011 & yearList_future<=2020,:),1)) ;
mean_2050s_r = squeeze(mean(ts_LUarea_crop_yr(yearList_future>=2051 & yearList_future<=2060,:),1)) ;
mean_2090s_r = squeeze(mean(ts_LUarea_crop_yr(yearList_future>=2091 & yearList_future<=2100,:),1)) ;

disp('Crop area')
(mean_2050s_r - mean_2010s_r) ./ mean_2010s_r
(mean_2090s_r - mean_2050s_r) ./ mean_2050s_r
(mean_2090s_r - mean_2010s_r) ./ mean_2010s_r


%% Plot timeseries: Global demand

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
conv_fact = 1e-3*1e-6 ;
units = 'Mt' ;
thisPos = figurePos ;
%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',thisPos) ;
for m = 1:Ncommods
    h = subplot_tight(2,4,m,spacing) ;
    ts_tmp_dmnd_yr = conv_fact*squeeze(ts_commodDemand_yvr(:,m,:)) ;
    plot(yearList_PLUMout,ts_tmp_dmnd_yr,'-','LineWidth',lineWidth) ;
    title(get_commods_title(commods{m}))
    ylabel(sprintf('Demand (%s)',units))
%     legend(runList,'Location','Best')
    legend(runList,'Location','Northwest')
    h.FontSize = fontSize ;
end

if do_save
    export_fig([outDir_ts 'commodity_demand.pdf'])
    close
end



%% Plot timeseries: Global per-capita demand

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
conv_fact = 1 ;
units = 'kg/person/yr' ;
thisPos = figurePos ;
%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',thisPos) ;
for m = 1:Ncommods
    h = subplot_tight(2,4,m,spacing) ;
    ts_tmp_dmnd_yr = conv_fact*squeeze(ts_commodDemandPC_yvr(:,m,:)) ;
    plot(yearList_PLUMout,ts_tmp_dmnd_yr,'-','LineWidth',lineWidth) ;
    title(get_commods_title(commods{m}))
    ylabel(sprintf('Demand (%s)',units))
    legend(runList,'Location','Best')
    h.FontSize = fontSize ;
end

if do_save
    export_fig([outDir_ts 'commodity_demand_perCapita.pdf'])
    close
end


%% Plot timeseries: Per-country (or group) production and exports

thisCountry_short = '' ;

%%%%%%%%%%%%%%%
% Options
thisYear = 2010 ;
% thisCountry = {'United States of America','Canada'} ; thisCountry_short = 'USA + Canada' ;
% thisCountry = {'United States of America'} ;
% thisCountry = {'Germany Austria & Switzerland'} ;
% thisCountry = {'Central Africa', 'Democratic Republic of the Congo', ...
%     'East Africa', 'Ethiopia', 'Kenya', 'Nigeria', 'South Africa', ...
%     'Southern Africa other', 'Sudan', 'Uganda', ...
%     'United Republic of Tanzania', 'West Africa'} ; thisCountry_short = 'Sub-Saharan Africa' ;
% thisCountry = 'Russian Federation' ;
% thisCountry = {'India  & Sri Lanka','Pakistan & Afghanistan', 'Bangladesh', 'Nepal & Butan'} ; thisCountry_short = 'South Asia' ;
thisCountry = 'China' ; thisCountry_short = 'China' ;
thisPos = figurePos ;
lineWidth = 2 ;
spacing = [0.1 0.05] ; % [v h]
fontSize = 14 ;
%%%%%%%%%%%%%%%

if ~exist('countryList_dem','var')
    import_country_demandEtc
end

if ischar(thisCountry)
    thisCountry = {thisCountry} ;
end

[~,IA] = intersect(countryList_dem,thisCountry) ;
if isempty(IA)
    error('No country match(es) found!')
elseif length(IA) < length(thisCountry)
    error('Only %d of %d countries found!', length(IA), length(thisCountry))
end
ts_tmp_prod_ymr = squeeze(sum(ts_countryProdnet_ymur(:,:,IA,:),3)) ;
ts_tmp_expt_ymr = squeeze(sum(-ts_countryImps_ymur(:,:,IA,:),3)) ;
ts_tmp_expt_ymr(ts_tmp_expt_ymr<0) = 0 ;

hf1 = figure('Color','w','Position',thisPos) ;
hf2 = figure('Color','w','Position',thisPos) ;
hf3 = figure('Color','w','Position',thisPos) ;
for m = 1:Ncommods
    
    ts_tmp_prod_yr = squeeze(ts_tmp_prod_ymr(:,m,:)) ;
    ts_tmp_expt_yr = squeeze(ts_tmp_expt_ymr(:,m,:)) ;

    % Production (net of seed/other waste) and export numbers
    set(0,'CurrentFigure',hf1) ;
    h = subplot_tight(2,4,m,spacing) ;
    plot(yearList_PLUMout,ts_tmp_prod_yr,'-','LineWidth',lineWidth) ;
    hold on
    h.ColorOrderIndex = 1 ;
    plot(yearList_PLUMout,ts_tmp_expt_yr,':','LineWidth',lineWidth) ;
    hold off
    ylabel('Prod. (dots=exports)')
    title(get_commods_title(commods{m}))
    h.FontSize = fontSize ;
    
    % Same, but as % change from thisYear
    ts_tmp_prod_yr_pct = 100*(-1 + ts_tmp_prod_yr / ts_tmp_prod_yr(yearList_PLUMout==thisYear)) ;
    ts_tmp_expt_yr_pct = 100*(-1 + ts_tmp_expt_yr / ts_tmp_expt_yr(yearList_PLUMout==thisYear)) ;
    set(0,'CurrentFigure',hf2) ;
    h = subplot_tight(2,4,m,spacing) ;
    plot(yearList_PLUMout,ts_tmp_prod_yr_pct,'-','LineWidth',lineWidth) ;
    hold on
    h.ColorOrderIndex = 1 ;
    plot(yearList_PLUMout,ts_tmp_expt_yr_pct,':','LineWidth',lineWidth) ;
    plot(h.XLim, [0 0], ':k')
    hold off
    ylabel(sprintf('Prod. %s from %d (%%, dots=exports)', '\Delta', thisYear))
    title(get_commods_title(commods{m}))
    h.FontSize = fontSize ;
    
    % Fraction of net production going to (net) exports
    set(0,'CurrentFigure',hf3) ;
    h = subplot_tight(2,4,m,spacing) ;
    plot(yearList_PLUMout,ts_tmp_expt_yr./ts_tmp_prod_yr,'-','LineWidth',lineWidth) ;
    ylabel('Fraction of production to exports')
    title(get_commods_title(commods{m}))
    h.FontSize = fontSize ;
    
    clear h
end; clear m

% Add figure titles
if ~isempty(thisCountry_short)
    thisSGtitle = thisCountry_short ;
else
    thisSGtitle = strjoin(thisCountry,' + ') ;
end
for hf = {hf1 hf2 hf3}
    sgtitle(hf{:}, thisSGtitle, ...
        'FontSize', fontSize+6, ...
        'FontWeight', 'Bold') ;
end

% Save figures
if do_save
    if isempty(thisCountry_short)
        error('Specify thisCountry_short to save figure')
    end
    thisCountry_short_filename = thisCountry_short_toFilename(thisCountry_short) ;
    set(0,'CurrentFigure',hf1) ;
    export_fig(sprintf('%s/production_%s.pdf', removeslashifneeded(outDir_ts), thisCountry_short_filename), ...
        '-r150') ;
    close(hf1)
    set(0,'CurrentFigure',hf2) ;
    export_fig(sprintf('%s/productionDelta_%s.pdf', removeslashifneeded(outDir_ts), thisCountry_short_filename), ...
        '-r150') ;
    close(hf2)
    set(0,'CurrentFigure',hf3) ;
    export_fig(sprintf('%s/exportfrac_%s.pdf', removeslashifneeded(outDir_ts), thisCountry_short_filename), ...
        '-r150') ;
    close(hf3)
end


clear hf* ts_tmp*


%% Plot timeseries: Per-country (or group) demand and domestic supply frac.

thisCountry_short = '' ;

%%%%%%%%%%%%%%%
% Options
thisYear = 2010 ;
% thisCountry = {'United States of America','Canada'} ; thisCountry_short = 'USA + Canada' ;
% thisCountry = {'United States of America'} ;
% thisCountry = {'Germany Austria & Switzerland'} ;
% thisCountry = {'India  & Sri Lanka','Pakistan & Afghanistan', 'Bangladesh', 'Nepal & Butan'} ; thisCountry_short = 'South Asia' ;
% thisCountry = {'Eastern Europe','France Netherlands & Benlex','Italy', ...
%     'Germany Austria & Switzerland', 'Other former USSR', 'Poland', ...
%     'Spain & Portugal', 'Ukraine', 'United Kingdom', 'ex-Yugoslavia'} ;
% thisCountry = 'Russian Federation' ;
thisCountry = 'China' ; thisCountry_short = 'China' ;
% thisCountry = {'Central Africa', 'Democratic Republic of the Congo', ...
%     'East Africa', 'Ethiopia', 'Kenya', 'Nigeria', 'South Africa', ...
%     'Southern Africa other', 'Sudan', 'Uganda', ...
%     'United Republic of Tanzania', 'West Africa'} ; thisCountry_short = 'Sub-Saharan Africa' ;
% thisCountry = {'India  & Sri Lanka','Pakistan & Afghanistan', 'Bangladesh', ...
%     'Nepal & Butan', 'China', 'Indonesia'} ; thisCountry_short = 'China, South Asia, Indonesia' ;
thisPos = figurePos ;
lineWidth = 2 ;
spacing = [0.1 0.05] ; % [v h]
fontSize = 14 ;
% units = 'kg' ; conv_fact = 1 ;
units = 'Mt' ; conv_fact = 1e-3*1e-6 ;
%%%%%%%%%%%%%%%

if ~exist('countryList_dem','var')
    import_country_demandEtc
end

if ischar(thisCountry)
    thisCountry = {thisCountry} ;
end

[~,IA] = intersect(countryList_dem,thisCountry) ;
if isempty(IA)
    error('No country match(es) found!')
elseif length(IA) < length(thisCountry)
    error('Only %d of %d countries found!', length(IA), length(thisCountry))
end
ts_tmp_dmnd_ymr = squeeze(sum(ts_countryDemand_ymur(:,:,IA,:),3)) ;
ts_tmp_dmndF_ymr = squeeze(sum(ts_countryDemandWithFeed_ymur(:,:,IA,:),3)) ;
ts_tmp_self_ymr = ...
    squeeze(sum(ts_countryProdnet_ymur(:,:,IA,:),3) ...
    ./ sum(ts_countryDemandWithFeed_ymur(:,:,IA,:),3)) ;
ts_tmp_self_ymr(ts_tmp_dmnd_ymr==0) = NaN ;
ts_tmp_self_ymr(ts_tmp_dmndF_ymr==0) = NaN ;
% ts_tmp_self_ymr(ts_tmp_self_ymr>1) = 1 ;

ts_tmp_dmnd_ymr = ts_tmp_dmnd_ymr * conv_fact ;
ts_tmp_dmndF_ymr = ts_tmp_dmndF_ymr * conv_fact ;

% figure('Color','w','Position',[350   268   720   350]) ;
% ydata = ts_tmp_dmndF_ymr(:,strcmp(commods, 'monogastrics'),:) ...
%     ./ sum(ts_tmp_dmndF_ymr(:,contains(commods, {'monogastrics','ruminants'}),:),2) ;
% ydata = squeeze(ydata) ;
% plot(yearList_PLUMout, ydata)
% title(thisCountry_short)
% ylabel('Monogastric fraction')
% set(gca, 'FontSize', fontSize) ;
% stop

hf1 = figure('Color','w','Position',thisPos) ;
hf2 = figure('Color','w','Position',thisPos) ;
hf3 = figure('Color','w','Position',thisPos) ;
for m = 1:Ncommods
    
    ts_tmp_dmndF_yr = squeeze(ts_tmp_dmndF_ymr(:,m,:)) ;
    ts_tmp_dmnd_yr = squeeze(ts_tmp_dmnd_ymr(:,m,:)) ;
    
    % Demand
    set(0,'CurrentFigure',hf1) ;
    h = subplot_tight(2,4,m,spacing) ;
    plot(yearList_PLUMout,ts_tmp_dmndF_yr,'-','LineWidth',lineWidth) ;
    ylabel(units)
    title(get_commods_title(commods{m}))
%     hold on
%     h.ColorOrderIndex = 1 ;
%     plot(yearList_PLUMout,ts_tmp_dmnd_yr,':','LineWidth',lineWidth) ;
%     hold off
%     ylabel(sprintf('%s (dots=w/o feed)', units))
    title(get_commods_title(commods{m}))
    legend(runList,'Location','Best')
    h.FontSize = fontSize ;
    
    % Same, but as % change from thisYear
    ts_tmp_dmndF_yr_pct = 100*(-1 + ts_tmp_dmndF_yr / ts_tmp_dmndF_yr(yearList_PLUMout==thisYear)) ;
    ts_tmp_dmnd_yr_pct = 100*(-1 + ts_tmp_dmnd_yr / ts_tmp_dmnd_yr(yearList_PLUMout==thisYear)) ;
    set(0,'CurrentFigure',hf2) ;
    h = subplot_tight(2,4,m,spacing) ;
    plot(yearList_PLUMout,ts_tmp_dmndF_yr_pct,'-','LineWidth',lineWidth) ;
    hold on
    plot(h.XLim, [0 0], ':k')
    h.ColorOrderIndex = 1 ;
    plot(yearList_PLUMout,ts_tmp_dmnd_yr_pct,':','LineWidth',lineWidth) ;
    ylabel(sprintf('%% %s from %d (dots=w/o feed)', '\Delta', thisYear))
    title(get_commods_title(commods{m}))
    hold off
    legend(runList,'Location','Best')
    h.FontSize = fontSize ;
    
    % Fraction of demand satisfied by domestic production
    set(0,'CurrentFigure',hf3) ;
    h = subplot_tight(2,4,m,spacing) ;
    ts_tmp_self_yr = squeeze(ts_tmp_self_ymr(:,m,:)) ;
    plot(yearList_PLUMout,ts_tmp_self_yr,'-','LineWidth',lineWidth) ;
    title(get_commods_title(commods{m}))
    legend(runList,'Location','Best')
    ylabel('"Self-sufficiency"')
    legend(runList,'Location','Best')
    h.FontSize = fontSize ;
    if max(h.YLim) - min(h.YLim) < 1e-9
        for y = 1:length(h.YTick)
            h.YTickLabel{y} = sprintf('%0.1f...',h.YTick(y)) ;
        end
    end
    clear h
end; clear m

% Add figure titles
if ~isempty(thisCountry_short)
    thisSGtitle = sprintf('%s: Demand', thisCountry_short) ;
else
    thisSGtitle = strjoin(thisCountry,' + ') ;
end
for hf = {hf1 hf2}
    sgtitle(hf{:}, thisSGtitle, ...
        'FontSize', fontSize+6, ...
        'FontWeight', 'Bold') ;
end
if ~isempty(thisCountry_short)
    thisSGtitle = thisCountry_short ;
else
    thisSGtitle = strjoin(thisCountry,' + ') ;
end
for hf = {hf3}
    sgtitle(hf{:}, thisSGtitle, ...
        'FontSize', fontSize+6, ...
        'FontWeight', 'Bold') ;
end

% Save figures
thisCountry_short_filename = thisCountry_short_toFilename(thisCountry_short) ;
if do_save
    if isempty(thisCountry_short)
        error('Specify thisCountry_short to save figure')
    end
    set(0,'CurrentFigure',hf1) ;
    export_fig(sprintf('%s/demand_%s.pdf', removeslashifneeded(outDir_ts), thisCountry_short_filename), ...
        '-r150') ;
    close(hf1)
    set(0,'CurrentFigure',hf2) ;
    export_fig(sprintf('%s/demandDelta_%s.pdf', removeslashifneeded(outDir_ts), thisCountry_short_filename), ...
        '-r150') ;
    close(hf2)
    set(0,'CurrentFigure',hf3) ;
    export_fig(sprintf('%s/selfsuff_%s.pdf', removeslashifneeded(outDir_ts), thisCountry_short_filename), ...
        '-r150') ;
    close(hf3)
    disp('Done')
    clear hf* ts_tmp*
end




%% Plot timeseries: Per-country (or group) demand and production

thisCountry_short = '' ;

%%%%%%%%%%%%%%%
% Options
thisYear = 2010 ;
% thisCountry = {'United States of America','Canada'} ; thisCountry_short = 'USA + Canada' ;
% thisCountry = {'United States of America'} ;
% thisCountry = {'Germany Austria & Switzerland'} ;
% thisCountry = {'India  & Sri Lanka','Pakistan & Afghanistan', 'Bangladesh', 'Nepal & Butan'} ; thisCountry_short = 'South Asia' ;
% thisCountry = {'Eastern Europe','France Netherlands & Benlex','Italy', ...
%     'Germany Austria & Switzerland', 'Other former USSR', 'Poland', ...
%     'Spain & Portugal', 'Ukraine', 'United Kingdom', 'ex-Yugoslavia'} ;
% thisCountry = 'Russian Federation' ;
thisCountry = 'China' ; thisCountry_short = 'China' ;
% thisCountry = {'Central Africa', 'Democratic Republic of the Congo', ...
%     'East Africa', 'Ethiopia', 'Kenya', 'Nigeria', 'South Africa', ...
%     'Southern Africa other', 'Sudan', 'Uganda', ...
%     'United Republic of Tanzania', 'West Africa'} ; thisCountry_short = 'Sub-Saharan Africa' ;
% thisCountry = {'India  & Sri Lanka','Pakistan & Afghanistan', 'Bangladesh', ...
%     'Nepal & Butan', 'China', 'Indonesia'} ; thisCountry_short = 'China, South Asia, Indonesia' ;
thisPos = figurePos ;
lineWidth = 2 ;
spacing = [0.1 0.05] ; % [v h]
fontSize = 14 ;
% units = 'kg' ; conv_fact = 1 ;
units = 'Mt' ; conv_fact = 1e-3*1e-6 ;
%%%%%%%%%%%%%%%

if ~exist('countryList_dem','var')
    import_country_demandEtc
end

if ischar(thisCountry)
    thisCountry = {thisCountry} ;
end

[~,IA] = intersect(countryList_dem,thisCountry) ;
if isempty(IA)
    error('No country match(es) found!')
elseif length(IA) < length(thisCountry)
    error('Only %d of %d countries found!', length(IA), length(thisCountry))
end
ts_tmp_dmndF_ymr = squeeze(sum(ts_countryDemandWithFeed_ymur(:,:,IA,:),3)) ;
ts_tmp_prod_ymr = squeeze(sum(ts_countryProdnet_ymur(:,:,IA,:),3)) ;

ts_tmp_dmndF_ymr = ts_tmp_dmndF_ymr * conv_fact ;
ts_tmp_prod_ymr = ts_tmp_prod_ymr * conv_fact ;

% figure('Color','w','Position',[350   268   720   350]) ;
% ydata = ts_tmp_dmndF_ymr(:,strcmp(commods, 'monogastrics'),:) ...
%     ./ sum(ts_tmp_dmndF_ymr(:,contains(commods, {'monogastrics','ruminants'}),:),2) ;
% ydata = squeeze(ydata) ;
% plot(yearList_PLUMout, ydata)
% title(thisCountry_short)
% ylabel('Monogastric fraction')
% set(gca, 'FontSize', fontSize) ;
% stop

hf1 = figure('Color','w','Position',thisPos) ;
for m = 1:Ncommods
    
    ts_tmp_dmndF_yr = squeeze(ts_tmp_dmndF_ymr(:,m,:)) ;
    ts_tmp_prod_yr = squeeze(ts_tmp_prod_ymr(:,m,:)) ;
    
    % Demand
    set(0,'CurrentFigure',hf1) ;
    h = subplot_tight(2,4,m,spacing) ;
    plot(yearList_PLUMout,ts_tmp_dmndF_yr,'-','LineWidth',lineWidth) ;
    
    % Production
    hold on
    h.ColorOrderIndex = 1 ;
    plot(yearList_PLUMout,ts_tmp_prod_yr,':','LineWidth',lineWidth) ;
    hold off
    
    % Finish
    legend(runList,'Location','Best')
    title(get_commods_title(commods{m}))
    ylabel(units)
    h.FontSize = fontSize ;
    
    clear h
end; clear m

% Add figure titles
if ~isempty(thisCountry_short)
    thisSGtitle = sprintf('%s: Demand and production', thisCountry_short) ;
else
    thisSGtitle = strjoin(thisCountry,' + ') ;
end
for hf = {hf1}
    sgtitle(hf{:}, thisSGtitle, ...
        'FontSize', fontSize+6, ...
        'FontWeight', 'Bold') ;
end

% Save figures
if do_save
    if isempty(thisCountry_short)
        error('Specify thisCountry_short to save figure')
    end
    set(0,'CurrentFigure',hf1) ;
    thisCountry_short_filename = thisCountry_short_toFilename(thisCountry_short) ;    set(0,'CurrentFigure',hf1) ;
    export_fig(sprintf('%s/demandVsProd_%s.pdf', removeslashifneeded(outDir_ts), thisCountry_short_filename), ...
        '-r150') ;
    close(hf1)
    disp('Done')
end

clear hf* ts_tmp*


%% Plot timeseries: Calories, EXPECTED

if all_figs && false
    
    % Options %%%%%%%%%
    %
    plum_area_adjustment = 1 ;
    % plum_area_adjustment = 1-0.28 ;
    
    % lpjg_area_adjustment = 1 ;
    lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;
    
    lineWidth = 2 ;
    fontSize = 14 ;
    spacing = [0.1 0.05] ;   % [vert, horz]
    ignYrs = 0 ;
    Nsmth = 1 ;
    figurePosition = [1 376 1440 429] ;
    %%%%%%%%%%%%%%%%%%%
    
    [tmp_ts_kcalExp_yr, title_suffix, file_suffix] = ...
        rebase_future2baseline(rebase, Nsmth, ts_kcal_bl, ts_kcalExp_yr, ignYrs, yearList_future) ;
    
    % Adjust for technology change, if doing so (do not include FAO!)
    tmp_ts_kcal_bl = ts_kcal_bl ;
    if do_adjYieldTech
        file_suffix = [file_suffix '_techAdj'] ;
        title_suffix = [title_suffix ' (techAdj)'] ;
    end
    
    % Adjust for unhandled crops, if doing so (DO include FAO)
    tmp_ts_kcal_fao = ts_kcal_fao ;
    if lpjg_area_adjustment ~= 1
        tmp_ts_kcal_bl = tmp_ts_kcal_bl .* lpjg_area_adjustment ;
        [~,IA] = intersect(yearList_baseline,fao.tmp_fao_yearList) ;
        tmp_ts_kcal_fao = tmp_ts_kcal_fao .* lpjg_area_adjustment(IA) ;
        title_suffix = [title_suffix ' (blAdj)'] ;
        file_suffix = [file_suffix '_blAdj'] ;
    end
    if plum_area_adjustment ~= 1
        tmp_ts_kcalExp_yr = tmp_ts_kcalExp_yr * plum_area_adjustment ;
        title_suffix = [title_suffix ' (PLUM\times' num2str(plum_area_adjustment) ')'] ;
        file_suffix = [file_suffix '_plumAdj' num2str(plum_area_adjustment)] ;
    end
    
    figure('Position',figurePosition,'Color','w') ;
    plot(yearList_baseline,cf_kcalEcal*movmean(tmp_ts_kcal_bl,Nsmth),'-k','LineWidth',lineWidth)
    hold on
    plot(fao.tmp_fao_yearList,cf_kcalEcal*tmp_ts_kcal_fao,'--k','LineWidth',lineWidth)
    plot(yearList_future,cf_kcalEcal*tmp_ts_kcalExp_yr,'LineWidth',lineWidth)
    hold off
    legend(stdLegend_plusFAO, ...
        'Location','NorthWest') ;
    set(gca,'FontSize',fontSize)
    xlabel('Year')
    ylabel('Ecal')
    ht = title(['Caloric production (PLUM-exp)' title_suffix]) ;
    clear tmp_*
    
    if do_save
        export_fig([outDir_ts 'caloriesExp' file_suffix '.pdf'])
        close
    end
    
end

stop


%% Plot timeseries: Calories, ACTUAL/EXPECTED

% Options %%%%%%%%%
lineWidth = 2 ;
fontSize = 14 ;
Nsmth = 1 ;
figurePosition = [1 376 1440 429] ;
%%%%%%%%%%%%%%%%%%%

if all_figs
    
    file_suffix = '' ;
    title_suffix = '' ;
    if do_adjYieldTech
        file_suffix = [file_suffix '_techAdj'] ;
        title_suffix = [title_suffix ' (techAdj)'] ;
    end
    
    figure('Position',figurePosition,'Color','w') ;
    tmp_ts_kcal_yr = ts_kcal_yr ;
    plot(yearList_future,movmean((tmp_ts_kcal_yr - ts_kcalExp_yr)./ts_kcalExp_yr*100,Nsmth,1),'-','LineWidth',lineWidth)
    hold on
    plot(get(gca,'XLim'),[0 0],'--k')
    hold off
    legend(stdLegend(2:end),'Location','SouthWest') ;
    set(gca,'FontSize',fontSize)
    xlabel('Year')
    ylabel('Ecal')
    title(['Calories diff.' title_suffix]) ;
    clear tmp_*
    
    if do_save
        export_fig([outDir_ts 'caloriesDiffFromExp' file_suffix '.pdf'])
        close
    end

end



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
%     theseGarrs = garr_cropfracs ;
%     figure_position = figurePos ;
% elseif strcmp(whichCFTs,'plum')
%     Nx = 4 ;
%     theseGarrs = garr_cropfracs_plum7 ;
%     figure_position = figurePos ;
% else
%     error(['whichCFTs (' whichCFTs ') not recognized!']) ;
% end
% 
% warning('Not tested after maps->garr conversion')
% for c = 1:Ncrops
%     figure('Position',figure_position,'Color','w') ;
%     thisCrop = theseCrops{c} ;
%     for p = 1:Nmaps
%         subplot_tight(Ny,Nx,p,spacing)
%         if p==1
%             tmp = theseGarrs.garr_xvyB(:,strcmp(theseGarrs.varNames,thisCrop),:) ;
%             tmp(garr_LU.garr_xvyB(:,strcmp(garr_LU.varNames,'CROPLAND'),:)==0) = 0 ;
%             tmp = mean(tmp,3) ;
%             tmp = lpjgu_vector2map(tmp, map_size, list2map) ;
%             pcolor(tmp(y2include,:)) ;
%             thisTitle = 'Baseline' ;
%         else
%             tmp = mean(theseGarrs.garr_xvyr(:,strcmp(theseGarrs.varNames,thisCrop),:,p-1),3) ;
%             tmp = lpjgu_vector2map(tmp, map_size, list2map) ;
%             pcolor(tmp(y2include,:)) ;
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
%% Map changes in BD hotspot area: CI

% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edgecolor = 0.6*ones(3,1) ;
% latlim = [-60,80];
latlim = [-60,50];
fontSize = 14 ;
spacing = [0.07 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
cbarOrient = 'SouthOutside' ;
lineWidth = 0.25 ;
conv_fact_map = 1e-6 ;   % m2 to km2
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
units_map = 'km^2' ;
units_total = 'Mkm^2' ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Biodiversity hotspots
hotspot_area_YXB = lpjgu_vector2map( ...
    mean(repmat(hotspot_area_x, [1 length(years_endh)]) .* permute(garr_LU_d9.garr_xvyB(:,strcmp(garr_LU_d9.varNames,'NATURAL'),:),[1 3 2]), 2), ...
    map_size, list2map) ;
hotspot_area_YXr = lpjgu_xz_to_YXz( ...
    squeeze(mean(repmat(hotspot_area_x, [1 length(years_endh) Nruns]) .* permute(garr_LU_d9.garr_xvyr(:,strcmp(garr_LU_d9.varNames,'NATURAL'),:,:),[1 3 4 2]), 2)), ...
    map_size, list2map) ;
hotspot_diff_YXr = hotspot_area_YXr - repmat(hotspot_area_YXB, [1 1 Nruns]) ;

map_hotspot_diffs(...
    hotspot_area_YXB, hotspot_diff_YXr, hotspot_YX, hotspot_shp, [], ...
    spacing, latlim, edgecolor, cbarOrient, fontSize, ...
    textX, textY_1, textY_2, ssp_plot_index, lineWidth, ...
    yearList_baseline, yearList_future, runList, ...
    conv_fact_map, conv_fact_total, units_map, units_total, do_caps)

if do_save
    export_fig( ...
        sprintf('%s/areaDiff_BDhotspots_CI.%d-%d.png',removeslashifneeded(outDir_maps),yearList_baseline(end),yearList_future(end)), ...
        ['-r' num2str(pngres)])
    close
end


%% Map changes in BD hotspot area: CI + CSLF

lu_source = 'plum' ;
% lu_source = 'luh1' ;

% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edgecolor = 0.6*ones(3,1) ;
% latlim = [-60,80];
latlim = [-60,50];
fontSize = 14 ;
% spacing = [0.07 0.05] ;   % [vert, horz]
spacing = [0.03 0.05] ;   % [vert, horz]
textX = 25 ;
textY_1 = 50 ;
textY_2 = 20 ;
cbarOrient = 'SouthOutside' ;
lineWidth = 0.5 ;
conv_fact_map = 1e-6 ;   % m2 to km2
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
units_map = 'km^2' ;
units_total = 'Mkm^2' ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Biodiversity hotspots
if strcmp(lu_source, 'plum')
    tmp_lu_d9_xyB = permute(garr_LU_d9.garr_xvyB(:,strcmp(garr_LU_d9.varNames,'NATURAL'),:),[1 3 2]) ;
    tmp_lu_d9_xyr = permute(garr_LU_d9.garr_xvyr(:,strcmp(garr_LU_d9.varNames,'NATURAL'),:,:),[1 3 4 2]) ;
    tmp_hotspotCSLF_area_YXB = lpjgu_vector2map( ...
        mean(repmat(hotspotCSLF_area_x, [1 length(years_endh)]) .* tmp_lu_d9_xyB, 2), ...
        map_size, list2map) ;
    tmp_hotspotCSLF_area_YXr = lpjgu_xz_to_YXz( ...
        squeeze(mean(repmat(hotspotCSLF_area_x, [1 length(years_endh) Nruns]) .* tmp_lu_d9_xyr, 2)), ...
        map_size, list2map) ;
    tmp_hotspotCSLF_diff_YXr = tmp_hotspotCSLF_area_YXr - repmat(tmp_hotspotCSLF_area_YXB, [1 1 Nruns]) ;
elseif strcmp(lu_source, 'luh1')
    import_luh1_NTRLfrac
    tmp_hotspotCSLF_area_YXrB = lpjgu_xz_to_YXz( ...
        repmat(hotspotCSLF_area_x, [1 Nruns]) .* garr_NTRLfrac_luh1_d9.garr_xrB, ...
        map_size, list2map) ;
    tmp_hotspotCSLF_area_YXB = tmp_hotspotCSLF_area_YXrB ;
    tmp_hotspotCSLF_area_YXrF = lpjgu_xz_to_YXz( ...
        repmat(hotspotCSLF_area_x, [1 Nruns]) .* garr_NTRLfrac_luh1_d9.garr_xrF, ...
        map_size, list2map) ;
    tmp_hotspotCSLF_diff_YXr = tmp_hotspotCSLF_area_YXrF - tmp_hotspotCSLF_area_YXrB ;
else
    error('lu_source (%s) not recognized', lu_source)
end

map_hotspot_diffs(...
    tmp_hotspotCSLF_area_YXB, tmp_hotspotCSLF_diff_YXr, (hotspot_YX | cslf_YX), hotspot_shp, cslf_shp, ...
    spacing, latlim, edgecolor, cbarOrient, fontSize, ...
    textX, textY_1, textY_2, ssp_plot_index, lineWidth, ...
    yearList_baseline, yearList_future, runList, ...
    conv_fact_map, conv_fact_total, units_map, units_total, do_caps)

if do_save
    export_fig( ...
        sprintf('%s/areaDiff_BDhotspots_CI_CSLF.%d-%d.png',removeslashifneeded(outDir_maps),yearList_baseline(end),yearList_future(end)), ...
        ['-r' num2str(pngres*2)])
    close
end

clear tmp_hotspot* tmp_lu*
    


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
% error('Make this work with SI units and garr (instead of maps)!')
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
%         hotspot_area_YXB, hotspot_diff_YXr, hotspot_YX, hotspot_shp, [], ...
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
% errb_endh_v = nan(Nvars,1) ;
% sem_begf = nan(Nvars,Nruns) ;
% errb_endf_vr = nan(Nvars,Nruns) ;
% string_endh = cell(Nvars,1) ;
% string_begf = cell(Nvars,Nruns) ;
% string_endf = cell(Nvars,Nruns) ;
% for c = 1:Nvars
%     
%     % Get values
%     thisVar = rowInfo{c,2} ;
%     thisConv = rowInfo{c,3} ;
%     mean_endh_v(c) = thisConv*eval(['mean(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
%     errb_endh_v(c) = thisConv*eval(['std(ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
%     mean_begf(c,:) = thisConv*eval(['mean(ts_' thisVar '_yr(yearList_future>=min(years_begf) & yearList_future<=max(years_begf),:))']) ;
%     sem_begf(c,:) = thisConv*eval(['std(ts_' thisVar '_yr(yearList_future>=min(years_begf) & yearList_future<=max(years_begf),:))']) ;
%     mean_endf_vr(c,:) = thisConv*eval(['mean(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
%     errb_endf_vr(c,:) = thisConv*eval(['std(ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
%     
%     % Turn into strings
%     if strcmp(rowInfo{c,4},'%d')
%         thisMean = round(mean_endh_v(c)) ;
%     else
%         thisMean = mean_endh_v(c) ;
%     end
%     if strcmp(rowInfo{c,4},'%d')
%         thisSD = round(errb_endh_v(c)) ;
%     else
%         thisSD = errb_endh_v(c) ;
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
%             thisSD_endf = round(errb_endf_vr(c,r)) ;
%         else
%             thisSD_begf = sem_begf(c,r) ;
%             thisSD_endf = errb_endf_vr(c,r) ;
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

if all_figs
    
    thisY1 = yearList_baseline(end) ;
    thisYN = yearList_future(end) ;
        
    % Natural area
    make_LUdiff_fig_v2(...
        lpjgu_vector2map(ntrl_area_xBH, map_size, list2map), ...
        lpjgu_xz_to_YXz(ntrl_diff_xrH, map_size, list2map), ...
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
        lpjgu_vector2map(crop_area_xBH, map_size, list2map), ...
        lpjgu_xz_to_YXz(crop_diff_xrH, map_size, list2map), ...
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
        lpjgu_vector2map(past_area_xBH, map_size, list2map), ...
        lpjgu_xz_to_YXz(past_diff_xrH, map_size, list2map), ...
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
        lpjgu_vector2map(agri_area_xBH, map_size, list2map), ...
        lpjgu_xz_to_YXz(agri_diff_xrH, map_size, list2map), ...
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
    if exist('crop0_area_xBH','var')
        make_LUdiff_fig_v2(...
            lpjgu_vector2map(crop0_area_xBH, map_size, list2map), ...
            lpjgu_xz_to_YXz(crop0_diff_xrH, map_size, list2map), ...
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
    if exist('past0_area_xBH','var')
        make_LUdiff_fig_v2(...
            lpjgu_vector2map(past0_area_xBH, map_size, list2map), ...
            lpjgu_xz_to_YXz(past0_diff_xrH, map_size, list2map), ...
            thisY1, thisYN, '"Pasture0"', runList, ...
            spacing, fontSize, textX, textY_1, textY_2, ...
            nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
            Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
        if do_save
            export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past0.png'],['-r' num2str(pngres)])
            close
        end
    end
    
end


%% Vegetation C if we had used LUH1
% Assume that all vegC is on natural land. This is thus an upper limit for
% the difference.

import_luh1_NTRLfrac

% Get natural veg C (assuming cropland and pasture veg C = 0)
% is_vegC = strcmp(garr_cpool_d9.varNames, 'VegC') ;
is_vegC = strcmp(garr_cpool_d9.varNames, 'Total') ;
garr_NTRLc_d9.garr_xrB = repmat(mean(sum(garr_cpool_d9.garr_xvyB(:,is_vegC,6:end),2),3), [1 Nruns]) ...
    ./ garr_NTRLfrac_plum_d9.garr_xrB ;
garr_NTRLc_d9.garr_xrF = squeeze(mean(sum(garr_cpool_d9.garr_xvyr(:,is_othr,6:end,:),2),3)) ...
    ./ garr_NTRLfrac_plum_d9.garr_xrF ;

% Get total veg C (PLUM)
garr_NTRLc_plum_d9.garr_xrB = garr_NTRLfrac_plum_d9.garr_xrB .* garr_NTRLc_d9.garr_xrB .* repmat(land_area_x, [1 Nruns]) ;
garr_NTRLc_plum_d9.garr_xrF = garr_NTRLfrac_plum_d9.garr_xrF .* garr_NTRLc_d9.garr_xrF .* repmat(land_area_x, [1 Nruns]) ;

% Convert PLUM total veg C to LUH1 version
garr_NTRLc_luh1_d9.garr_xrB = garr_NTRLc_plum_d9.garr_xrB .* garr_NTRLfrac_luh1_d9.garr_xrB ./ garr_NTRLfrac_plum_d9.garr_xrB ;
garr_NTRLc_luh1_d9.garr_xrF = garr_NTRLc_plum_d9.garr_xrF .* garr_NTRLfrac_luh1_d9.garr_xrF ./ garr_NTRLfrac_plum_d9.garr_xrF ;

% Ignore where PLUM had little natural area
garr_NTRLc_plum_d9.garr_xrB(garr_NTRLfrac_plum_d9.garr_xrB<0.001 & garr_NTRLfrac_luh1_d9.garr_xrB==0) = 0 ;
garr_NTRLc_plum_d9.garr_xrB(garr_NTRLfrac_plum_d9.garr_xrB<0.001 & garr_NTRLfrac_luh1_d9.garr_xrB>0) = NaN ;
garr_NTRLc_plum_d9.garr_xrF(garr_NTRLfrac_plum_d9.garr_xrF<0.001 & garr_NTRLfrac_luh1_d9.garr_xrF==0) = 0 ;
garr_NTRLc_plum_d9.garr_xrF(garr_NTRLfrac_plum_d9.garr_xrF<0.001 & garr_NTRLfrac_luh1_d9.garr_xrF>0) = NaN ;
garr_NTRLc_luh1_d9.garr_xrB(garr_NTRLfrac_plum_d9.garr_xrB<0.001 & garr_NTRLfrac_luh1_d9.garr_xrB==0) = 0 ;
garr_NTRLc_luh1_d9.garr_xrB(garr_NTRLfrac_plum_d9.garr_xrB<0.001 & garr_NTRLfrac_luh1_d9.garr_xrB>0) = NaN ;
garr_NTRLc_luh1_d9.garr_xrF(garr_NTRLfrac_plum_d9.garr_xrF<0.001 & garr_NTRLfrac_luh1_d9.garr_xrF==0) = 0 ;
garr_NTRLc_luh1_d9.garr_xrF(garr_NTRLfrac_plum_d9.garr_xrF<0.001 & garr_NTRLfrac_luh1_d9.garr_xrF>0) = NaN ;

for r = 1:Nruns
    disp(runList{r})
    fprintf('     \t Baseline \t Future \t Delta\n')
    plum_b = nansum(garr_NTRLc_plum_d9.garr_xrB(:,r)) ;
    plum_f = nansum(garr_NTRLc_plum_d9.garr_xrF(:,r)) ;
    plum_b2f = plum_f - plum_b ;
    luh1_b = nansum(garr_NTRLc_luh1_d9.garr_xrB(:,r)) ;
    luh1_f = nansum(garr_NTRLc_luh1_d9.garr_xrF(:,r)) ;
    luh1_b2f = luh1_f - luh1_b ;
    diff_b = (luh1_b - plum_b) / plum_b * 100 ;
    diff_f = (luh1_f - plum_f) / plum_f * 100 ;
    diff_b2f = (luh1_b2f - plum_b2f) / abs(plum_b2f) * 100 ;
    fprintf('PLUM:\t %0.3g \t %0.3g \t %0.3g\n', plum_b, plum_f, plum_b2f) ;
    fprintf('LUH1:\t %0.3g \t %0.3g \t %0.3g\n', luh1_b, luh1_f, luh1_b2f) ;
    fprintf('diff:\t %0.1f%% \t\t %0.1f%% \t\t %0.1f%%\n', diff_b, diff_f, diff_b2f) ;
    disp(' ')
end


%% Map lost primary land

if ~exist('dPrim_YX', 'var')
    % fileName = '/Users/sam/Geodata/LUH2/v2h/states.1850-2015.nc' ;
    % Ntimes = 166 ;
    % date1 = [1850 1 1] ;
    % dateN = [2010 1 1] ;
    fileName = '/Volumes/WDMPP_Storage/ExternalGeodata/LUH2/v2h/states.nc' ;
    Ntimes = 1166 ;
    date1 = [850 1 1] ;
    dateN = [2010 1 1] ;
    
    Nlons = 1440 ;
    Nlats = 720 ;
    verbose = false ;
    do_byPFT = false ;
    
    varName = 'primf' ;
    primf_t1_YX = flipud(read_and_arrange_selectTime(fileName, ...
        Nlons, Nlats, Ntimes, varName, verbose, do_byPFT, ...
        date1)) ;
    primf_tN_YX = flipud(read_and_arrange_selectTime(fileName, ...
        Nlons, Nlats, Ntimes, varName, verbose, do_byPFT, ...
        dateN)) ;
    
    varName = 'primn' ;
    primn_t1_YX = flipud(read_and_arrange_selectTime(fileName, ...
        Nlons, Nlats, Ntimes, varName, verbose, do_byPFT, ...
        date1)) ;
    primn_tN_YX = flipud(read_and_arrange_selectTime(fileName, ...
        Nlons, Nlats, Ntimes, varName, verbose, do_byPFT, ...
        dateN)) ;
    
    prim_t1_YX = primf_t1_YX + primn_t1_YX ;
    prim_tN_YX = primf_tN_YX + primn_tN_YX ;
    carea_YX = 1e6*flipud(transpose(ncread( ...
        '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc', ...
        'carea'))) ;
    prim_t1_YX = prim_t1_YX .* carea_YX ;
    prim_tN_YX = prim_tN_YX .* carea_YX ;
    
    prim_t1_YX(isnan(prim_t1_YX)) = 0 ;
    tmp = prim_t1_YX(:,1:2:1440) + prim_t1_YX(:,2:2:1440) ;
    prim_t1_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
    clear tmp
    
    prim_tN_YX(isnan(prim_tN_YX)) = 0 ;
    tmp = prim_tN_YX(:,1:2:1440) + prim_tN_YX(:,2:2:1440) ;
    prim_tN_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
    clear tmp
    
    dPrim_YX = prim_tN_YX - prim_t1_YX ;
end

if ~exist('biomeID_x', 'var')
    biomeID_YX = flipud(imread( ...
        '/Users/sam/Geodata/General/WWF terrestrial ecosystems/wwf_terr_ecos_dissolveBiome_halfDeg_id.tif')) ;
    biomeID_YX(biomeID_YX<0) = NaN ;
    biomeID_x = biomeID_YX(list2map) ;
end

disp('Done')

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
% if all_figs
%     
%     thisY1 = yearList_baseline(end) ;
%     thisYN = yearList_future(1) ;
%     
%     % Natural area
%     make_LUdiff_fig_v2(...
%         lpjgu_vector2map(ntrl_area_xBH, map_size, list2map), ntrl_area_YXBFr - repmat(lpjgu_vector2map(ntrl_area_xBH, map_size, list2map),[1 1 Nruns]), ...
%         thisY1, thisYN, '"Natural"', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
%         close
%     end
%     
%     % Cropland area
%     make_LUdiff_fig_v2(...
%         lpjgu_vector2map(crop_area_xBH, map_size, list2map), crop_area_YXBFr - repmat(lpjgu_vector2map(crop_area_xBH, map_size, list2map),[1 1 Nruns]), ...
%         thisY1, thisYN, 'Cropland', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
%         close
%     end
%     
%     % "Pasture" area
%     make_LUdiff_fig_v2(...
%         lpjgu_vector2map(past_area_xBH, map_size, list2map), past_area_YXBFr - repmat(lpjgu_vector2map(past_area_xBH, map_size, list2map),[1 1 Nruns]), ...
%         thisY1, thisYN, '"Pasture"', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
%         close
%     end
%     
%     % Agricultural area
%     make_LUdiff_fig_v2(...
%         lpjgu_vector2map(agri_area_xBH, map_size, list2map), agri_area_YXBFr - repmat(lpjgu_vector2map(agri_area_xBH, map_size, list2map),[1 1 Nruns]), ...
%         thisY1, thisYN, 'Agricultural', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_agri.png'],['-r' num2str(pngres)])
%         close
%     end
%     
%     % Cropland0 area
%     if exist('crop0_area_xBH','var')
%         make_LUdiff_fig_v2(...
%             lpjgu_vector2map(crop0_area_xBH, map_size, list2map), crop_area_YXBFr - repmat(lpjgu_vector2map(crop0_area_xBH, map_size, list2map),[1 1 Nruns]), ...
%             thisY1, thisYN, 'Cropland0', runList, ...
%             spacing, fontSize, textX, textY_1, textY_2, ...
%             nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%             Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%         if do_save
%             export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop0.png'],['-r' num2str(pngres)])
%             close
%         end
%     end
%     
%     % "Pasture0" area
%     if exist('past0_area_xBH','var')
%         make_LUdiff_fig_v2(...
%             lpjgu_vector2map(past0_area_xBH, map_size, list2map), past_area_YXBFr - repmat(lpjgu_vector2map(past0_area_xBH, map_size, list2map),[1 1 Nruns]), ...
%             thisY1, thisYN, '"Pasture0"', runList, ...
%             spacing, fontSize, textX, textY_1, textY_2, ...
%             nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%             Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%         if do_save
%             export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past0.png'],['-r' num2str(pngres)])
%             close
%         end
%     end
%     
% end


%% Map changes in LU area: Begin-Future to End-Future

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
% if all_figs
%     
%     thisY1 = yearList_future(1) ;
%     thisYN = yearList_future(end) ;
%     
%     % Natural area
%     make_LUdiff_fig_v2(...
%         ntrl_area_YXBFr, ntrl_diff_YXrF, ...
%         thisY1, thisYN, '"Natural"', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_ntrl.png'],['-r' num2str(pngres)])
%         close
%     end
%     
%     % Cropland area
%     make_LUdiff_fig_v2(...
%         crop_area_YXBFr, crop_diff_YXrF, ...
%         thisY1, thisYN, 'Cropland', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_crop.png'],['-r' num2str(pngres)])
%         close
%     end
%     
%     % "Pasture" area
%     make_LUdiff_fig_v2(...
%         past_area_YXBFr, past_diff_YXrF, ...
%         thisY1, thisYN, '"Pasture"', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_past.png'],['-r' num2str(pngres)])
%         close
%     end
%     
%     % Agricultural area
%     make_LUdiff_fig_v2(...
%         agri_area_YXBFr, agri_diff_YXrF, ...
%         thisY1, thisYN, 'Agricultural', runList, ...
%         spacing, fontSize, textX, textY_1, textY_2, ...
%         nx, ny, 1, colorBarLoc, ssp_plot_index, only1bl, ...
%         Nruns, thisPos, conv_fact_map, conv_fact_total, units_map, units_total, do_caps) ;
%     if do_save
%         export_fig([outDir_maps 'areaDiff_' num2str(thisY1) '-' num2str(thisYN) '_LU_agri.png'],['-r' num2str(pngres)])
%         close
%     end
% 
% end


%% Save outputs

if false
    disp('Saving outputs...')
    save(['/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/fig_script_outputs/' thisVer '.mat'], ...
        'agri_*', 'bare_*', 'bl_*', 'crop*', 'firstdec_tmp', 'lastdec_tmp', ...
        'land_area_YX', 'gcel_area_YX', 'maps_*', 'nanmask', 'ntrl_*', 'past*', ...
        'ts_*_bl','ts_*_yr')
    disp('Done.')
end


%%

disp('All done!')



