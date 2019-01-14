%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare FAO data with PLUM-expected values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plum_ver = 7 ;
% plum_ver = 8 ;
% plum_ver = 9 ;
plum_ver = 10 ;

faoCommBalElements = {'Production'} ;
% faoCommBalElements = {'Feed'} ;
% faoCommBalElements = {'Waste','Processing','Food','Other uses'} ;
% faoCommBalElements = {'Food'} ;

do_save = true ;
rebase = false ;
pngres = 150 ;

do_caps = -1 ;


%% Setup

CFTnames = {'Wheat','Maize','Rice','Oilcrops','Pulses','StarchyRoots'} ;
Ncrops = length(CFTnames) ;

yearList_plum = shiftdim(2010:2100) ;
Nyears_plum = length(yearList_plum) ;

outDir_ts = addslashifneeded(['/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/'...
                           'figures_v' num2str(plum_ver) '_FAOvsPLUMfbs']) ;
if ~exist(outDir_ts,'dir')
    mkdir(outDir_ts)
end

% Conversion factors
cf_kg2Mt = 1e-3*1e-6 ;
cf_t2kg = 1e3 ;
cf_ha2m2 = 1e4 ;
cf_kcal2Ecal = 1e-15 ;

addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work/')


%% Import FAO data

for f = 1:length(faoCommBalElements)
    
    faoCommBalElement = faoCommBalElements{f} ;
    disp(faoCommBalElement)
    fao = load(['/Users/Shared/PLUM/crop_calib_data/fao/FAOdata_1961-2010_calibVer16_' faoCommBalElement '.mat']) ;
    
    yearList_fao = shiftdim(fao.tmp_fao_yearList) ;
    
    % Adjust FAO data
    if Ncrops>0
        if any(strcmp(fao.listCrops_fa2o,'Starchy roots'))
            fao.listCrops_fa2o{strcmp(fao.listCrops_fa2o,'Starchy roots')} = 'StarchyRoots' ;
        end
        if any(contains(CFTnames,'_plum'))
            fao.listCrops_fa2o = strcat(fao.listCrops_fa2o,'_plum') ;
        end
        for c = 1:Ncrops
            thisCrop = CFTnames{c} ;
            if isfield(fao.listCrops_fa2o,'a') || isfield(fao.listCrops_fa2o,'p')
                if ~(isfield(fao.listCrops_fa2o,'a') && isfield(fao.listCrops_fa2o,'p'))
                    error('One but not both present of fao.listCrops_fa2o.a and .b')
                end
            elseif length(find(strcmp(fao.listCrops_fa2o,thisCrop)))==1
                if f==1
                    eval(['ts_cropprod_' thisCrop '_fao = single(cf_t2kg*squeeze(nansum(fao.tmp_total_fa2_Ccy(:,strcmp(fao.listCrops_fa2o,thisCrop),:),1))) ;']) ;
                else
                    eval(['ts_cropprod_' thisCrop '_fao = ts_cropprod_' thisCrop '_fao + single(cf_t2kg*squeeze(nansum(fao.tmp_total_fa2_Ccy(:,strcmp(fao.listCrops_fa2o,thisCrop),:),1))) ;']) ;
                end
            else
                warning(['Assuming zeros for FAO data for ' thisCrop])
                if f==1
                    eval(['ts_cropprod_' thisCrop '_fao = zeros(length(fao.tmp_fao_yearList),1,''single'') ;']) ;
                end
            end
        end
    end
    
end
disp('Done')


%% Import PLUM FBS data

%%%%%%%%%%%%
% Options
top_dir = '/Volumes/WDMPP_Storage/Shared/PLUM/PLUM_outputs_for_LPJG' ;
%%%%%%%%%%%%

% Decide which fbs.txt column to use
if isequal(faoCommBalElements,{'Waste','Processing','Food','Other uses'})
    thisCol = {'SeedAndOtherLosses','FoodAnd1stGen'} ;
elseif strcmp(faoCommBalElement,'Production')
    thisCol = 'Production' ;
    % thisCol = 'FoodAnd1stGen' ; calories = true ;
elseif strcmp(faoCommBalElement,'Feed')
    thisCol = {'MonogastricsFeed','RuminantsFeed'} ;
elseif strcmp(faoCommBalElement,'Food')
    thisCol = {'FoodAnd1stGen'} ;
else
    error(['This faoCommBalElement (' faoCommBalElement ') not recognized.'])
end

% Get available directories
top_dir = addslashifneeded(top_dir) ;
dirList = dir([top_dir 'SSP*.v' num2str(plum_ver) '.s1']) ;
Nruns = length(dirList) ;

% Set up empties
for c = 1:Ncrops
    thisCrop = CFTnames{c} ;
    eval(['ts_cropprodExp_' thisCrop '_yr = nan(Nyears_plum,Nruns) ;']) ;
end

for r = 1:Nruns
    file_in = [dirList(r).folder '/' dirList(r).name '/fbs.txt'] ;
    table_in = readtable(file_in) ;
    
    for c = 1:Ncrops
        thisCrop = CFTnames{c} ;
%         eval(['ts_cropprodExp_' thisCrop '_yr(:,r) = cf_t2kg*1e6*table2array(table_in(strcmpi(table_in.Crop,thisCrop),strcmp(table_in.Properties.VariableNames,thisCol))) ;']) ;
        eval(['ts_cropprodExp_' thisCrop '_yr(:,r) = cf_t2kg*1e6*sum(table2array(table_in(strcmpi(table_in.Crop,thisCrop),contains(table_in.Properties.VariableNames,thisCol))),2) ;']) ;
    end
end

% Calculate calorie production: FAO
Nyears_fao = length(ts_cropprod_Oilcrops_fao) ;
tmp = whos('ts_cropprod_*') ;
tmp_name = {tmp.name}' ;
ts_kcal_fao = zeros(size(Nyears_fao,1)) ;
for i = 1:length(tmp_name)
    thisCrop = strrep(strrep(strrep(strrep(tmp_name{i},'ts_cropprod_',''),'_bl',''),'_yr',''),'_fao','') ;
    if strcmp(thisCrop,'Miscanthus')
        continue
    end
    thisSuffix = strrep(tmp_name{i},['ts_cropprod_' thisCrop '_'],'') ;
    kcal_per_g = get_kcalDensity(strrep(strrep(thisCrop,'Maize','CerealsC4'),'Wheat','CerealsC3')) ;
    kcal_per_kg = 1e3 * kcal_per_g ;
    eval(['ts_kcal_' thisSuffix ' = ts_kcal_' thisSuffix ' + kcal_per_kg * eval(tmp_name{i}) ;']) ;
end ; clear i

% Calculate calorie production: PLUM
Nyears_plum = length(ts_cropprodExp_Maize_yr) ;
tmp = whos('ts_cropprodExp_*') ;
tmp_name = {tmp.name}' ;
ts_kcal_yr = zeros(size(Nyears_plum,Nruns)) ;
for i = 1:length(tmp_name)
    thisCrop = strrep(strrep(strrep(strrep(tmp_name{i},'ts_cropprodExp_',''),'_bl',''),'_yr',''),'_fao','') ;
    if strcmp(thisCrop,'Miscanthus')
        continue
    end
    thisSuffix = strrep(tmp_name{i},['ts_cropprodExp_' thisCrop '_'],'') ;
    kcal_per_g = get_kcalDensity(strrep(strrep(thisCrop,'Maize','CerealsC4'),'Wheat','CerealsC3')) ;
    kcal_per_kg = 1e3 * kcal_per_g ;
    eval(['ts_kcal_' thisSuffix ' = ts_kcal_' thisSuffix ' + kcal_per_kg * eval(tmp_name{i}) ;']) ;
end


%% Plot production (kg)

%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot Options %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
plum_area_adjustment = 1 ;
% plum_area_adjustment = 1-0.28 ;
lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;
thisVar = 'cropprodExp' ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = cf_kg2Mt ;
units = 'Mt DM' ;
%%%%%%%%%%%%%%%%%%%


if iscell(thisCol)
    if isequal(thisCol,{'SeedAndOtherLosses','FoodAnd1stGen'})
        title_prefix = 'All but feed' ;
    else
        title_prefix = '' ;
        Ncols = length(thisCol) ;
        for i = 1:Ncols
            title_prefix = [title_prefix thisCol{i}] ;
            if i < Ncols
                title_prefix = [title_prefix '+'] ;
            end
        end
    end
else
    title_prefix = thisCol ;
end

theseCFTnames = CFTnames ;
thisLegend = {'FAO','S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;

cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
cmds = get_cell_forPlot(whos, strrep(thisVar,'Exp',''), 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end

file_suffix = CFT_timeseries_plot(...
    {}, cell_yr, cell_fao, yearList_fao, yearList_plum, yearList_fao, rebase, ...
    theseCFTnames, units, title_prefix, thisLegend, do_caps, ...
    'lineWidth',lineWidth, ...
    'fontSize',fontSize, ...
    'spacing',spacing, ...
    'ignYrs',ignYrs, ...
    'Nsmth',Nsmth, ...
    'conv_fact',conv_fact, ...
    'plum_area_adjustment', plum_area_adjustment, ...
    'lpjg_area_adjustment', lpjg_area_adjustment, ...
    'fao_linestyle','-') ;

file_suffix = ['_fbs_' strrep(title_prefix,' ','') file_suffix] ;

if do_save
    export_fig([outDir_ts 'ssp' num2str(plum_ver) '_' thisVar file_suffix '.pdf'])
    close
end


%% Plot production (kcal)

%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot Options %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
plum_area_adjustment = 1 ;
% plum_area_adjustment = 1-0.28 ;
lpjg_area_adjustment = 1 ;
% lpjg_area_adjustment = ts_LUarea_crop0_bl ./ ts_LUarea_crop_bl ;
thisVar = 'kcal' ;
title_prefix = thisCol ;
lineWidth = 2 ;
fontSize = 14 ;
spacing = [0.1 0.05] ;   % [vert, horz]
ignYrs = 0 ;
Nsmth = 1 ;
conv_fact = cf_kcal2Ecal ;
units = 'Ecal' ;
%%%%%%%%%%%%%%%%%%%

theseCFTnames = CFTnames ;
thisLegend = {'FAO','S1/R4.5','S3/R6.0','S4/R6.0','S5/R8.5'} ;

% cmds = get_cell_forPlot(whos, thisVar, 'yr', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% cmds = get_cell_forPlot(whos, thisVar, 'fao', theseCFTnames); for c = 1:length(cmds); eval(cmds{c}); end
% 
% file_suffix = CFT_timeseries_plot(...
%     {}, cell_yr, cell_fao, yearList_fao, yearList_plum, yearList_fao, rebase, ...
%     theseCFTnames, units, title_prefix, thisLegend, do_caps, ...
%     'lineWidth',lineWidth, ...
%     'fontSize',fontSize, ...
%     'spacing',spacing, ...
%     'ignYrs',ignYrs, ...
%     'Nsmth',Nsmth, ...
%     'conv_fact',conv_fact, ...
%     'plum_area_adjustment', plum_area_adjustment, ...
%     'lpjg_area_adjustment', lpjg_area_adjustment) ;

figure('Color','w') ;
plot(yearList_fao,ts_kcal_fao,'-k','LineWidth',lineWidth)
hold on
plot(yearList_plum,ts_kcal_yr,'LineWidth',lineWidth)
hold off
legend(thisLegend,'Location','Best')
set(gca,'FontSize',fontSize)
ylabel(units)
title('Crop production (calories)')

% file_suffix = ['_fbs_' thisCol file_suffix] ;
% 
% if do_save
%     export_fig([outDir_ts thisVar file_suffix '.pdf'])
%     close
% end





