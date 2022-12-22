%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLUM2LPJG figures: Compare experiments %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expt_names = {'v4s1_v20180426','v4s1_v20180426_asPLUMout2011'} ;
expt_titles = {'Default','PLUM 2011'} ;

% expt_names = {'v4s1_v20180426_asPLUMout2011','v4s1_v20180426_asLUH2_2010'} ;
% expt_titles = {'PLUM 2011','LUH2 2010'} ;

do_save = true ;
rebase = false ;
do_pct = false ;
dashed_zero_line = true ;


%% Setup

% Where are outputs saved?
% fig_script_outputs_dir = addslashifneeded('/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/fig_script_outputs') ;
fig_script_outputs_dir = addslashifneeded('/Volumes/Crucial 480GB SSD/PLUM/trunk_runs/fig_script_outputs') ;

% Where to save figures?
outDir_ts = addslashifneeded(...
    ['/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/' ...
    'figs_comp_' expt_names{1} '_VS_' expt_names{2}]) ;
if ~exist(outDir_ts,'dir')
    unix(['mkdir -p ' outDir_ts]) ;
end

% Set up path
addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work/')

% Build runList
runList = {'SSP1 (R45)','SSP3 (R60)','SSP4 (R60)','SSP5 (R85)'} ;
% runList = {} ;
% for e = 1:2
%     thisName = expt_names{e} ;
%     if strcmp(thisName,'v4s1_v20180426')
%         runList = [runList {'SSP1 (R45)','SSP3 (R60)','SSP4 (R60)','SSP5 (R85)'}] ;
%     elseif strcmp(thisName,'v4s1_v20180426_asPLUMout2011')
%         runList = [runList {'SSP1constPLUM2011 (R45)','SSP3constPLUM2011 (R60)','SSP4constPLUM2011 (R60)','SSP5constPLUM2011 (R85)'}] ;
%     else
%         error([thisName ' not recognized in building of runList!'])
%     end
% end

% Get year info
yearList_baseline = 1850:2010 ;
yearList_future = 2011:2100 ;

% Conversion factors
cf_kg2Mt = 1e-3*1e-6 ;
cf_t2kg = 1e3 ;   % For FAO only
cf_ha2m2 = 1e4 ; % For FAO only
cf_kgPm2_to_tonsPha = 1e-3*1e4 ;
cf_kg2Tg = 1e-9 ;
cf_kg2Pg = 1e-12 ;
cf_m3_to_km3 = (1e-3)^3 ;

% Capitalize subplot labels?
do_caps = -1 ;


%% Import

disp('Importing time series: experiment 1...')
ts1 = load([fig_script_outputs_dir expt_names{1} '.mat'],'ts_*_yr') ;

disp('Importing time series: experiment 2...')
ts2 = load([fig_script_outputs_dir expt_names{2} '.mat'],'ts_*_yr') ;

% Get CFT names
tmp = fieldnames(ts1) ;
CFTnames = strrep(strrep(tmp(contains(tmp,'ts_yield_')),'ts_yield_',''),'_yr','') ;
clear tmp

disp('Done importing.')


%% Plot timeseries of diff: Area of each crop

% Options %%%%%%%%%
thisVar = 'croparea' ;
title_prefix = 'Area' ;
lineWidth = 2 ;
fontSize = 12 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = 1e-6*1e-6 ;
units = 'Million km^2' ;
legend_loc = 'best' ;
%%%%%%%%%%%%%%%%%%%

clear cell_1 cell_2
cmds = get_cell_forPlot_fromStruct(ts1, thisVar, 1, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot_fromStruct(ts2, thisVar, 2, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cell_diff = cell(length(cell_1),1) ;
for c = 1:length(cell_1)
    cell_diff{c} = cell_1{c} - cell_2{c} ;
end

file_suffix = CFT_timeseries_plot(...
    [], cell_diff, [], [], yearList_future, [], false, ...
    CFTnames, units, title_prefix, runList, do_caps, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact, ...
    'legend_loc',legend_loc, ...
    'dashed_zero_line',dashed_zero_line) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries of diff: Production of each crop

% Options %%%%%%%%%
thisVar = 'cropprod' ;
title_prefix = 'Production' ;
lineWidth = 2 ;
fontSize = 12 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = cf_kg2Mt ;
units = 'Mt DM' ;
legend_loc = 'best' ;
%%%%%%%%%%%%%%%%%%%

clear cell_1 cell_2
cmds = get_cell_forPlot_fromStruct(ts1, thisVar, 1, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot_fromStruct(ts2, thisVar, 2, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cell_diff = cell(length(cell_1),1) ;
for c = 1:length(cell_1)
    cell_diff{c} = cell_1{c} - cell_2{c} ;
end

file_suffix = CFT_timeseries_plot(...
    [], cell_diff, [], [], yearList_future, [], false, ...
    CFTnames, units, title_prefix, runList, do_caps, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact, ...
    'legend_loc',legend_loc, ...
    'dashed_zero_line',dashed_zero_line) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries of diff: Yield of each crop

% Options %%%%%%%%%
thisVar = 'yield' ;
title_prefix = 'Yield' ;
lineWidth = 2 ;
fontSize = 12 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = cf_kgPm2_to_tonsPha ;
units = 't ha^{-1}' ;
legend_loc = 'best' ;
%%%%%%%%%%%%%%%%%%%

clear cell_1 cell_2
cmds = get_cell_forPlot_fromStruct(ts1, thisVar, 1, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot_fromStruct(ts2, thisVar, 2, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cell_diff = cell(length(cell_1),1) ;
for c = 1:length(cell_1)
    cell_diff{c} = cell_1{c} - cell_2{c} ;
end

file_suffix = CFT_timeseries_plot(...
    [], cell_diff, [], [], yearList_future, [], false, ...
    CFTnames, units, title_prefix, runList, do_caps, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact, ...
    'legend_loc',legend_loc, ...
    'dashed_zero_line',dashed_zero_line) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries of diff: Yield of each crop (%)

% Options %%%%%%%%%
thisVar = 'yield' ;
title_prefix = 'Yield' ;
lineWidth = 2 ;
fontSize = 12 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = 1 ;
units = '%' ;
legend_loc = 'best' ;
%%%%%%%%%%%%%%%%%%%

clear cell_1 cell_2
cmds = get_cell_forPlot_fromStruct(ts1, thisVar, 1, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot_fromStruct(ts2, thisVar, 2, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cell_diff = cell(length(cell_1),1) ;
for c = 1:length(cell_1)
    cell_diff{c} = (cell_1{c} - cell_2{c}) ./ cell_2{c} * 100 ;
end

file_suffix = CFT_timeseries_plot(...
    [], cell_diff, [], [], yearList_future, [], false, ...
    CFTnames, units, title_prefix, runList, do_caps, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact, ...
    'legend_loc',legend_loc, ...
    'dashed_zero_line',dashed_zero_line) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries of diff: Nfert per ha on each crop

% Options %%%%%%%%%
thisVar = 'nflux_fert' ;
title_prefix = 'N applied per ha' ;
lineWidth = 2 ;
fontSize = 12 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = 1e4 ;   % kg/m2 to kg/ha (kg/m2 * m2/ha)
units = 'kg ha^{-1}' ;
legend_loc = 'best' ;
%%%%%%%%%%%%%%%%%%%

% clear cell_1 cell_2
% cmds = get_cell_forPlot_fromStruct(ts1, thisVar, 1, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cmds = get_cell_forPlot_fromStruct(ts2, thisVar, 2, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end


cmds = get_cell_forPlot_fromStruct(ts1, thisVar, 1, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot_fromStruct(ts2, thisVar, 2, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cell_1_nflux_fert = cell_1 ;
cell_2_nflux_fert = cell_2 ;
clear cell_1 cell_2
cmds = get_cell_forPlot_fromStruct(ts1, 'croparea', 1, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cell_1_croparea = cell_1 ;
clear cell_1
cmds = get_cell_forPlot_fromStruct(ts2, 'croparea', 2, CFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cell_2_croparea = cell_2 ;
clear cell_2
% Get kg/m2
cell_1 = cell(length(cell_1_nflux_fert),1) ;
cell_2 = cell(length(cell_1_nflux_fert),1) ;
for i = 1:length(cell_1_nflux_fert)
    cell_1{i} = cell_1_nflux_fert{i} ./ cell_1_croparea{i} ;
    cell_2{i} = cell_2_nflux_fert{i} ./ cell_2_croparea{i} ;
end

cell_diff = cell(length(cell_1),1) ;
for c = 1:length(cell_1)
    cell_diff{c} = cell_1{c} - cell_2{c} ;
end

file_suffix = CFT_timeseries_plot(...
    [], cell_diff, [], [], yearList_future, [], false, ...
    CFTnames, units, title_prefix, runList, do_caps, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact, ...
    'legend_loc',legend_loc, ...
    'dashed_zero_line',dashed_zero_line) ;

if do_save
    export_fig([outDir_ts thisVar file_suffix '.pdf'])
    close
end


%% Plot timeseries of diff: N loss

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
Nsmth = 1 ;
units = 'TgN' ;
conv_fact = cf_kg2Tg ;
legend_loc = 'best' ;
filename = 'Nloss' ;
%%%%%%%%%%%%%%%%%%%

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
plot_scen_diff_timeseries(...
    'ts_nflux_flux_yr', ts1, ts2, yearList_future, conv_fact, ...
    Nsmth, lineWidth, runList, legend_loc, fontSize, units, ...
    'N loss: Gaseous', do_caps, '', dashed_zero_line, do_pct, ...
    expt_titles) ;

subplot_tight(1,2,2,spacing)
plot_scen_diff_timeseries(...
    'ts_nflux_leach_yr', ts1, ts2, yearList_future, conv_fact, ...
    Nsmth, lineWidth, runList, legend_loc, fontSize, units, ...
    'N loss: Dissolved', do_caps, '', dashed_zero_line, do_pct, ...
    expt_titles) ;

if do_save
    file_suffix = '' ;
    if do_pct
        file_suffix = [file_suffix '_pct'] ;
    end
    export_fig([outDir_ts filename file_suffix '.pdf'])
    close
end



%% Plot timeseries of diff: Albedo

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
Nsmth = 1 ;
units = '' ;
legend_loc = 'best' ;
% filename = 'cpools' ;
%%%%%%%%%%%%%%%%%%%

figure('Position',figurePos,'Color','w') ;

if do_pct
    ydata1 = movmean((ts1.ts_albedo1_yr-ts2.ts_albedo1_yr)./ts2.ts_albedo1_yr*100,Nsmth) ;
    ydata7 = movmean((ts1.ts_albedo7_yr-ts2.ts_albedo7_yr)./ts2.ts_albedo7_yr*100,Nsmth) ;
    units = ['% (rel. to "' expt_titles{2} '")'] ;
else
    ydata1 = movmean(ts1.ts_albedo1_yr-ts2.ts_albedo1_yr,Nsmth) ;
    ydata7 = movmean(ts1.ts_albedo7_yr-ts2.ts_albedo7_yr,Nsmth) ;
end

plot(yearList_future,ydata1,'LineWidth',lineWidth)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(yearList_future,ydata7,'--','LineWidth',lineWidth)
hold off
clear ydata1 ydata7
legend([strcat(runList,' (Jan.)') strcat(runList,' (Jul.)')],...
    'Location',legend_loc,'AutoUpdate','off') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel(units)

title_suffix = '' ;
if Nsmth>1
    title_suffix = [title_suffix ' (Nsmth=' num2str(Nsmth)] ;
end

ht = title([title_prefix ': "' expt_titles{1} '" - "' expt_titles{2} '"' title_suffix]) ;

if dashed_zero_line
    ylims = get(gca,'YLim') ;
    if min(ylims)<0 && max(ylims)>0
        hold on
        xlims = get(gca,'XLim') ;
        plot(xlims,[0 0],'--k') ;
        hold off
    end
end

if do_save
    export_fig([outDir_ts 'albedo' file_suffix '.pdf'])
    close
end


%% Plot timeseries of diff: C pools

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
units = 'Pg' ;
conv_fact = cf_kg2Pg ;
legend_loc = 'best' ;
filename = 'cpools' ;
%%%%%%%%%%%%%%%%%%%

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
plot_scen_diff_timeseries(...
    'ts_cpool_VegC_yr', ts1, ts2, yearList_future, conv_fact, ...
    Nsmth, lineWidth, runList, legend_loc, fontSize, units, ...
    'Vegetation C', do_caps, '', dashed_zero_line, do_pct, ...
    expt_titles) ;

subplot_tight(1,2,2,spacing)
plot_scen_diff_timeseries(...
    'ts_cpool_Total_yr', ts1, ts2, yearList_future, conv_fact, ...
    Nsmth, lineWidth, runList, legend_loc, fontSize, units, ...
    'Total C', do_caps, '', dashed_zero_line, do_pct, ...
    expt_titles) ;

if do_save
    file_suffix = '' ;
    if do_pct
        file_suffix = [file_suffix '_pct'] ;
    end
    export_fig([outDir_ts filename file_suffix '.pdf'])
    close
end


%% Plot timeseries of diff: BVOCs

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
units = 'TgC' ;
conv_fact = cf_kg2Tg ;
legend_loc = 'best' ;
filename = 'bvocs' ;
%%%%%%%%%%%%%%%%%%%

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
plot_scen_diff_timeseries(...
    'ts_amon_yr', ts1, ts2, yearList_future, conv_fact, ...
    Nsmth, lineWidth, runList, legend_loc, fontSize, units, ...
    'BVOCs: Monoterpenes', do_caps, '', dashed_zero_line, do_pct, ...
    expt_titles) ;

subplot_tight(1,2,2,spacing)
plot_scen_diff_timeseries(...
    'ts_aiso_yr', ts1, ts2, yearList_future, conv_fact, ...
    Nsmth, lineWidth, runList, legend_loc, fontSize, units, ...
    'BVOCs: Isoprene', do_caps, '', dashed_zero_line, do_pct, ...
    expt_titles) ;

if do_save
    file_suffix = '' ;
    if do_pct
        file_suffix = [file_suffix '_pct'] ;
    end
    export_fig([outDir_ts filename file_suffix '.pdf'])
    close
end


%% Plot timeseries of diff: Evapotranspiration & runoff

% Options %%%%%%%%%
lineWidth = 3 ;
fontSize = 24 ;
spacing = [0.1 0.1] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
units = '1000 km^3' ;
conv_fact = cf_m3_to_km3*1e-3 ;
legend_loc = 'best' ;
filename = 'evapotranspiration_runoff' ;
%%%%%%%%%%%%%%%%%%%

figure('Position',figurePos,'Color','w') ;

subplot_tight(1,2,1,spacing)
plot_scen_diff_timeseries(...
    {'ts_aevap_yr','ts_aaet_yr'}, ts1, ts2, yearList_future, conv_fact, ...
    Nsmth, lineWidth, runList, legend_loc, fontSize, units, ...
    'Evapotranspiration', do_caps, '', dashed_zero_line, do_pct, ...
    expt_titles) ;

subplot_tight(1,2,2,spacing)
plot_scen_diff_timeseries(...
    'ts_tot_runoff_yr', ts1, ts2, yearList_future, conv_fact, ...
    Nsmth, lineWidth, runList, legend_loc, fontSize, units, ...
    'Runoff', do_caps, '', dashed_zero_line, do_pct, ...
    expt_titles) ;

if do_save
    file_suffix = '' ;
    if do_pct
        file_suffix = [file_suffix '_pct'] ;
    end
    export_fig([outDir_ts filename file_suffix '.pdf'])
    close
end


