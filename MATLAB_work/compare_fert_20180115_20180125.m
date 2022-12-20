%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare baseline fertilization 20180115 vs 20180125 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup

file1_fert = '/Users/Shared/PLUM/input/Nfert/LUH2/Nfert_fromLUH2.min5kgha.1901-2015.NfertIrrFactEmpties.64493.20161221.txt' ;
file2_fert = '/Users/Shared/PLUM/input/Nfert/LUH2/Nfert_fromLUH2.min5kgha.1850-2010.64493.20170116.assignWWorSW.txt' ;

file1_land = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
file2_land = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;

file1_crop = '/Users/Shared/remapping_MIRCA/crop_fractions_fromMIRCA_PLUM_1yr_20180105b.MP.out' ;
file2_crop = '/Users/Shared/PLUM/input/LU/crop_fractions_fromMIRCA_PLUM7_1850-2010_20180105a.assignWWorSW.out' ;

% Import land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calibration_for_PLUM/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd


%% Import

disp('Importing...')

fert1 = lpjgu_matlab_readTable_then2map(file1_fert,'force_mat_save',true) ;
fert2 = lpjgu_matlab_readTable_then2map(file2_fert,'force_mat_save',true) ;

land1 = lpjgu_matlab_readTable_then2map(file1_land,'force_mat_save',true) ;
if ~strcmp(file1_land,file2_land)
    land2 = lpjgu_matlab_readTable_then2map(file2_land,'force_mat_save',true) ;
else
    land2 = land1 ;
end
land1.maps_YXvy = land1.maps_YXvy .* repmat(land_area_YX,[1 1 size(land1.maps_YXvy,3) size(land1.maps_YXvy,4)]) ;
land2.maps_YXvy = land2.maps_YXvy .* repmat(land_area_YX,[1 1 size(land2.maps_YXvy,3) size(land2.maps_YXvy,4)]) ;

crop1 = lpjgu_matlab_readTable_then2map(file1_crop,'force_mat_save',true) ;
crop2 = lpjgu_matlab_readTable_then2map(file2_crop,'force_mat_save',true) ;

disp('Done.')


%% Match years

disp('Matching years...')

if ~isequal(fert1.yearList,fert2.yearList)
    [yearList,IA,IB] = intersect(fert1.yearList,fert2.yearList) ;
    
    fert1.maps_YXvy = fert1.maps_YXvy(:,:,:,IA) ;
    fert2.maps_YXvy = fert2.maps_YXvy(:,:,:,IB) ;
    
    fert1.yearList = yearList ;
    fert2.yearList = yearList ;
end
yearList = fert1.yearList ;
Nyears = length(yearList) ;

if ~isfield(crop1,'yearList')
    tmp = crop1.maps_YXv ;
    crop1 = rmfield(crop1,'maps_YXv') ;
    crop1.yearList = yearList ;
    crop1.maps_YXvy = repmat(tmp,[1 1 1 Nyears]) ;
end

if ~isfield(crop2,'yearList')
    tmp = crop2.maps_YXv ;
    crop2 = rmfield(crop2,'maps_YXv') ;
    crop2.yearList = yearList ;
    crop2.maps_YXvy = repmat(tmp,[1 1 1 Nyears]) ;
end

if ~isequal(yearList,crop1.yearList)
    [~,IA] = intersect(yearList,crop1.yearList) ;
    crop1.maps_YXvy = crop1.maps_YXvy(:,:,:,IA) ;
    crop1.yearList = yearList ;
end

if ~isequal(yearList,crop2.yearList)
    [~,IA] = intersect(yearList,crop2.yearList) ;
    crop2.maps_YXvy = crop2.maps_YXvy(:,:,:,IA) ;
    crop2.yearList = yearList ;
end

if ~isequal(yearList,land1.yearList)
    [~,IA] = intersect(yearList,land1.yearList) ;
    land1.maps_YXvy = land1.maps_YXvy(:,:,:,IA) ;
    land1.yearList = yearList ;
end

if ~isequal(yearList,land2.yearList)
    [~,IA] = intersect(yearList,land2.yearList) ;
    land2.maps_YXvy = land2.maps_YXvy(:,:,:,IA) ;
    land2.yearList = yearList ;
end

disp('Done.')


%% Get ancillary maps

disp('Getting ancillary maps...')

map1_croparea_YXy = squeeze(land1.maps_YXvy(:,:,strcmp(land1.varNames,'CROPLAND'),:)) ;
map2_croparea_YXy = squeeze(land2.maps_YXvy(:,:,strcmp(land2.varNames,'CROPLAND'),:)) ;

rice_rfFrac_1 = strcmp(crop1.varNames,'TrRi') | strcmp(crop1.varNames,'Rice') ;
rice_rfFrac_2 = strcmp(crop2.varNames,'TrRi') | strcmp(crop2.varNames,'Rice') ;
if length(find(rice_rfFrac_1))>1
    error('length(find(rice_rfFrac_1))>1')
elseif length(find(rice_rfFrac_2))>1
    error('length(find(rice_rfFrac_2))>1')
elseif isempty(find(rice_rfFrac_1,1))
    error('length(find(rice_rfFrac_1))==0')
elseif isempty(find(rice_rfFrac_2,1))
    error('length(find(rice_rfFrac_2))==0')
end
map1_ricefrac_rf_YXy = squeeze(crop1.maps_YXvy(:,:,rice_rfFrac_1,:)) ;
map2_ricefrac_rf_YXy = squeeze(crop2.maps_YXvy(:,:,rice_rfFrac_2,:)) ;
map1_ricearea_rf_YXy = map1_croparea_YXy .* map1_ricefrac_rf_YXy ;
map2_ricearea_rf_YXy = map2_croparea_YXy .* map2_ricefrac_rf_YXy ;

rice_irFrac_1 = strcmp(crop1.varNames,'TrRii') | strcmp(crop1.varNames,'Ricei') ;
rice_irFrac_2 = strcmp(crop2.varNames,'TrRii') | strcmp(crop2.varNames,'Ricei') ;
if length(find(rice_irFrac_1))>1
    error('length(find(rice_irFrac_1))>1')
elseif length(find(rice_irFrac_2))>1
    error('length(find(rice_irFrac_2))>1')
elseif isempty(find(rice_irFrac_1,1))
    error('length(find(rice_irFrac_1))==0')
elseif isempty(find(rice_irFrac_2,1))
    error('length(find(rice_irFrac_2))==0')
end
map1_ricefrac_ir_YXy = squeeze(crop1.maps_YXvy(:,:,rice_irFrac_1,:)) ;
map2_ricefrac_ir_YXy = squeeze(crop2.maps_YXvy(:,:,rice_irFrac_2,:)) ;
map1_ricearea_ir_YXy = map1_croparea_YXy .* map1_ricefrac_ir_YXy ;
map2_ricearea_ir_YXy = map2_croparea_YXy .* map2_ricefrac_ir_YXy ;

map1_ricearea_YXy = map1_ricearea_rf_YXy + map1_ricearea_ir_YXy ;
map2_ricearea_YXy = map2_ricearea_rf_YXy + map2_ricearea_ir_YXy ;

rice_rfFert_1 = strcmp(fert1.varNames,'TrRi') | strcmp(fert1.varNames,'Rice') ;
rice_rfFert_2 = strcmp(fert2.varNames,'TrRi') | strcmp(fert2.varNames,'Rice') ;
if length(find(rice_rfFert_1))>1
    error('length(find(rice_rfFert_1))>1')
elseif length(find(rice_rfFert_2))>1
    error('length(find(rice_rfFert_2))>1')
elseif isempty(find(rice_rfFert_1,1))
    error('length(find(rice_rfFert_1))==0')
elseif isempty(find(rice_rfFert_2,1))
    error('length(find(rice_rfFert_2))==0')
end
map1_ricefert_rf_YXy = map1_croparea_YXy .* squeeze(fert1.maps_YXvy(:,:,rice_rfFert_1,:)) ;
map2_ricefert_rf_YXy = map2_croparea_YXy .* squeeze(fert2.maps_YXvy(:,:,rice_rfFert_2,:)) ;

rice_irFert_1 = strcmp(fert1.varNames,'TrRii') | strcmp(fert1.varNames,'Ricei') ;
rice_irFert_2 = strcmp(fert2.varNames,'TrRii') | strcmp(fert2.varNames,'Ricei') ;
if length(find(rice_irFert_1))>1
    error('length(find(rice_irFert_1))>1')
elseif length(find(rice_irFert_2))>1
    error('length(find(rice_irFert_2))>1')
elseif isempty(find(rice_irFert_1,1))
    error('length(find(rice_irFert_1))==0')
elseif isempty(find(rice_irFert_2,1))
    error('length(find(rice_irFert_2))==0')
end
map1_ricefert_ir_YXy = map1_croparea_YXy .* squeeze(fert1.maps_YXvy(:,:,rice_irFert_1,:)) ;
map2_ricefert_ir_YXy = map2_croparea_YXy .* squeeze(fert2.maps_YXvy(:,:,rice_irFert_2,:)) ;

map1_ricefert_YXy = map1_ricefert_rf_YXy + map1_ricefert_ir_YXy ;
map2_ricefert_YXy = map2_ricefert_rf_YXy + map2_ricefert_ir_YXy ;

disp('Done.')


%% Compare cropland area

tmp1 = squeeze(nansum(nansum(map1_croparea_YXy,2),1)) ;
tmp2 = squeeze(nansum(nansum(map2_croparea_YXy,2),1)) ;

figure('Color','w')
plot(yearList,tmp1,'-b','LineWidth',3)
hold on
plot(yearList,tmp2,'--r','LineWidth',3)
hold off
legend('20180115','20180125')
xlabel('Year')
ylabel('Area (km^2)')
title('Total cropland area')
set(gca,'FontSize',14)


%% Compare rice area

figure('Color','w','Position',figurePos)

subplot(2,2,1)
tmp1 = squeeze(nansum(nansum(map1_ricearea_rf_YXy,2),1)) ;
tmp2 = squeeze(nansum(nansum(map2_ricearea_rf_YXy,2),1)) ;
plot(yearList,tmp1,'-b','LineWidth',3)
hold on
plot(yearList,tmp2,'--r','LineWidth',3)
hold off
legend('20180115','20180125','Location','Best')
xlabel('Year')
ylabel('Area (km^2)')
title('Rice area (rainfed)')
set(gca,'FontSize',14)

subplot(2,2,2)
tmp1 = squeeze(nansum(nansum(map1_ricearea_ir_YXy,2),1)) ;
tmp2 = squeeze(nansum(nansum(map2_ricearea_ir_YXy,2),1)) ;
plot(yearList,tmp1,'-b','LineWidth',3)
hold on
plot(yearList,tmp2,'--r','LineWidth',3)
hold off
legend('20180115','20180125','Location','Best')
xlabel('Year')
ylabel('Area (km^2)')
title('Rice area (irrigated)')
set(gca,'FontSize',14)

subplot(2,2,3)
tmp1 = squeeze(nansum(nansum(map1_ricearea_rf_YXy+map1_ricearea_ir_YXy,2),1)) ;
tmp2 = squeeze(nansum(nansum(map2_ricearea_rf_YXy+map2_ricearea_ir_YXy,2),1)) ;
plot(yearList,tmp1,'-b','LineWidth',3)
hold on
plot(yearList,tmp2,'--r','LineWidth',3)
hold off
legend('20180115','20180125','Location','Best')
xlabel('Year')
ylabel('Area (km^2)')
title('Rice area (rf+irr)')
set(gca,'FontSize',14)


%% Compare rice fertilization

figure('Color','w','Position',figurePos)

subplot(2,2,1)
tmp1 = squeeze(nansum(nansum(map1_ricefert_rf_YXy,2),1)) ;
tmp2 = squeeze(nansum(nansum(map2_ricefert_rf_YXy,2),1)) ;
plot(yearList,tmp1,'-b','LineWidth',3)
hold on
plot(yearList,tmp2,'--r','LineWidth',3)
hold off
legend('20180115','20180125','Location','Best')
xlabel('Year')
ylabel('???')
title('Rice fert (rainfed)')
set(gca,'FontSize',14)

subplot(2,2,2)
tmp1 = squeeze(nansum(nansum(map1_ricefert_ir_YXy,2),1)) ;
tmp2 = squeeze(nansum(nansum(map2_ricefert_ir_YXy,2),1)) ;
plot(yearList,tmp1,'-b','LineWidth',3)
hold on
plot(yearList,tmp2,'--r','LineWidth',3)
hold off
legend('20180115','20180125','Location','Best')
xlabel('Year')
ylabel('???')
title('Rice fert (irrigated)')
set(gca,'FontSize',14)

subplot(2,2,3)
tmp1 = squeeze(nansum(nansum(map1_ricefert_rf_YXy+map1_ricefert_ir_YXy,2),1)) ;
tmp2 = squeeze(nansum(nansum(map2_ricefert_rf_YXy+map2_ricefert_ir_YXy,2),1)) ;
plot(yearList,tmp1,'-b','LineWidth',3)
hold on
plot(yearList,tmp2,'--r','LineWidth',3)
hold off
legend('20180115','20180125','Location','Best')
xlabel('Year')
ylabel('???')
title('Rice fert (rf+irr)')
set(gca,'FontSize',14)









