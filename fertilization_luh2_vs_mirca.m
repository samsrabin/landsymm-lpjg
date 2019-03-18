%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exploring LUH2 data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

yearList_toRead = 1960:2010 ;


%% Import LUH2 data

cd '/Users/sam/Geodata/LUH2/'

crop_names = {'c3ann','c4ann','c3per','c4per','c3nfx'} ;
Ncrops = length(crop_names) ;

% lu_names = [crop_names ...
%             {'primf','primn','secdf','secdn','urban','pastr','range'}] ;
lu_names = crop_names ;
Nlu = length(lu_names) ;

Nyears = length(yearList_toRead) ;

yearList_inFile = 1850:2010 ;
[~,yearIndices,~] = intersect(yearList_inFile,yearList_toRead) ;
starts = [1 1 min(yearIndices)] ;
counts = [Inf Inf Nyears] ;

% Import cell area (km2)
carea_YX = flipud(transpose(ncread('supporting/staticData_quarterdeg.nc','carea'))) ;
icwtr_YX = flipud(transpose(ncread('supporting/staticData_quarterdeg.nc','icwtr'))) ;
larea_YX = (1-icwtr_YX) .* carea_YX ;

% Land use fractions
luh2_LUfrac_YXyv = nan([size(carea_YX) Nyears Nlu]) ;
for i = 1:Nlu
    thisVar = lu_names{i} ;
    fprintf('Reading %s...\n', thisVar) ;
    tmp = ncread('v2h/states.1850-2015.nc',thisVar,starts,counts) ;
    luh2_LUfrac_YXyv(:,:,:,i) = flipud(permute(tmp,[2 1 3])) ;
end
disp('Done reading land uses.')

% Crop fractions and areas (km2)
[~,IA,~] = intersect(lu_names,crop_names,'stable') ;
luh2_cropfrac_YXyv = luh2_LUfrac_YXyv(:,:,:,IA) ;
luh2_croparea_YXyv = luh2_cropfrac_YXyv .* repmat(carea_YX,[1 1 Nyears Ncrops]) ;

% Fertilization (kg/ha)
luh2_nfert_YXyv = nan([size(carea_YX) Nyears Ncrops]) ;
for i = 1:Ncrops
    thisVar = ['fertl_' crop_names{i}] ;
    fprintf('Reading %s...\n', thisVar) ;
    tmp = ncread('v2h/management.1850-2015.nc',thisVar,starts,counts) ;
    luh2_nfert_YXyv(:,:,:,i) = flipud(permute(tmp,[2 1 3])) ;
end
disp('Done reading fertilization.')

% Fertilization (kg)
luh2_nfertTot_YXyv = luh2_nfert_YXyv .* (1e2*luh2_croparea_YXyv) ;


%% Import MIRCA data

addpath(genpath('/Users/sam/Geodata/MIRCA/harvested_area_grids_26crops_30mn/MATLAB work'))

% Set up for area import
cd '/Users/sam/Geodata/MIRCA/harvested_area_grids_26crops_30mn'
list_crops_MIRCA = {'Wheat' ; 'Maize' ; 'Rice' ; 'Barley' ; 'Rye' ;
                    'Millet' ; 'Sorghum' ; 'Soybeans' ; 'Sunflower' ;
                    'Potatoes' ; 'Cassava' ;'Sugarcane' ; 'Sugarbeet' ;
                    'Oilpalm' ; 'RapeseedCanola' ; 'GroundnutsPeanuts' ;
                    'Pulses' ; 'Citrus' ; 'Datepalm' ; 'GrapesVine' ;
                    'Cotton' ; 'Cocoa' ; 'Coffee' ; 'OtherPerennials' ;
                    'FodderGrasses' ; 'OtherAnnuals'} ;
Ncrops_MIRCA = length(list_crops_MIRCA) ;

% Import MIRCA area (originally ha; convert to km2)
maps_ir_YXc = nan(360,720,Ncrops_MIRCA) ;
maps_rf_YXc = nan(360,720,Ncrops_MIRCA) ;
disp('Importing MIRCA files:')
for c = 1:Ncrops_MIRCA
    thisCrop = list_crops_MIRCA{c} ;
    disp(['   ' thisCrop ' (' num2str(c) ' of ' num2str(Ncrops_MIRCA) ')'])
    % Note that import_mirca does flipud() to make MATLAB plots look right.
    maps_ir_YXc(:,:,c) = import_mirca(['originals/annual_area_harvested_irc_crop' sprintf('%02.0f',c) '_ha_30mn.asc']) ;
    maps_rf_YXc(:,:,c) = import_mirca(['originals/annual_area_harvested_rfc_crop' sprintf('%02.0f',c) '_ha_30mn.asc']) ;
end
mirca_croparea_YXv = 1e-2*(maps_ir_YXc + maps_rf_YXc) ;
% clear maps_*_YXc

% Set up for Nfert import
cd '/Users/Shared/unifying_gridlist/AgGRID_nutrient_input_v1.1'
list_crops_agmip = {'cassava' ; 'groundnut' ; 'maize' ; 'millet' ;
                    'potato' ; 'rapeseed' ; 'rice' ; 'sorghum' ;
                    'soybean' ; 'sugarbeet' ; 'sunflower' ; 'wheat'} ;
list_crops_agmip_asMIRCA = {'Cassava' ; 'GroundnutsPeanuts' ; 'Maize' ;
                            'Millet' ; 'Potatoes' ; 'RapeseedCanola' ;
                            'Rice' ; 'Sorghum' ; 'Soybeans' ; 'Sugarbeet' ;
                            'Sunflower' ; 'Wheat'} ;
list_crops_agmip_ignore = {'Coffee';'Cotton';'Sugarcane'} ;
Ncrops_agmip = length(list_crops_agmip) ;

% Sanity check
if Ncrops_agmip ~= length(list_crops_agmip_asMIRCA)
    error('Lengths of list_crops_agmip and list_crops_agmip_asMIRCA do not match.')
end
    
% Import AgMIP files
agmip_nfert_YXv = nan(360,720,Ncrops_agmip) ;
disp('Importing AgMIP files:')
for c = 1:Ncrops_agmip
    thisCrop = list_crops_agmip{c} ;
    disp(['   ' thisCrop ' (' num2str(c) ' of ' num2str(Ncrops_agmip) ')'])
    filename = ['agmip_' thisCrop '_apprate_fill_NPK_0.5.nc4'] ;
    % Note that flipud(transpose()) to make MATLAB plots look right.
    agmip_nfert_YXv(:,:,c) = flipud(transpose(ncread(filename,'Napprate'))) ;
end
disp('Done.')

% Get crops included in AgMIP
[~,IA,IB] = intersect(list_crops_MIRCA, list_crops_agmip_asMIRCA,'stable') ;
agmip_croparea_YXv = mirca_croparea_YXv(:,:,IA) ;
agmip_nfert_YXv = agmip_nfert_YXv(:,:,IB) ;

% Get total N applied (kg)
agmip_nfertTot_YXv = (1e2*agmip_croparea_YXv) .* agmip_nfert_YXv ;


%% Figure: Time series: Global cropland area
%%%%%%%%%%%%%%%%
% Options
lineWidth = 3 ;
markerSize = 10 ;
fontSize = 18 ;
%%%%%%%%%%%%%%%%

luh2_crop_km2_y = squeeze(nansum(nansum(nansum(luh2_croparea_YXyv,2),1),4)) ;
mirca_crop_km2 = nansum(mirca_croparea_YXv(:)) ;
agmip_crop_km2 = nansum(agmip_croparea_YXv(:)) ;

figure('Color','w','Position',figurePos) ;
plot(yearList_toRead,luh2_crop_km2_y,'LineWidth',lineWidth) ;
hold on
plot(2000,mirca_crop_km2,'.k','MarkerSize',4*markerSize)
plot(2000,agmip_crop_km2,'ok','MarkerSize',markerSize)
hold off
xlabel('Year')
ylabel('Area (km^2)')
legend({'LUH2','MIRCA','AgMIP'},'Location','West') ;
title('Global cropland area')
set(gca,'FontSize',fontSize,'YLim',[0 max(get(gca,'YLim'))]) ;


%% Figure: Time series: Global N applied
%%%%%%%%%%%%%%%%
% Options
lineWidth = 3 ;
markerSize = 10 ;
fontSize = 18 ;
%%%%%%%%%%%%%%%%

luh2_nfertTot_Mt_yv = (1e-3*1e-6)*squeeze(nansum(nansum(luh2_nfertTot_YXyv,1),2)) ;
agmip_nfertTot_Mt = (1e-3*1e-6)*nansum(agmip_nfertTot_YXv(:)) ;

figure('Color','w','Position',figurePos) ;
plot(yearList_toRead, luh2_nfertTot_Mt_yv,'LineWidth',lineWidth) ;
hold on
plot(yearList_toRead, sum(luh2_nfertTot_Mt_yv,2), '-k','LineWidth',lineWidth) ;
plot(2000,agmip_nfertTot_Mt,'ok','MarkerSize',markerSize) ;
hold off
xlabel('Year')
ylabel('N applied (Mt)')
legend([crop_names {'LUH2 total','AgMIP total'}],'Location','NorthWest')
set(gca,'FontSize',18)






