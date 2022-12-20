%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare PLUM trajectory with LUH1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_luh1 = '/Users/sam/Geodata/LU_HurttEA2011/future/LUHa.v1_message.v1.updated_states.mat' ;
file_plum = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1.harm.forLPJG/landcover.txt.maps.mat' ;


%% Setup

% Map LUH1 types to PLUM types
plum_from_luh1 = {{'past'};
                  {'crop'};
                  {'secd','othr'};
                  {'icwr'}} ;


%% Import

% PLUM land use
load(file_plum) ; LU_plum = out_struct ; clear out_struct

% LUH1 land use
load(file_luh1) ;

% Get equal yearLists
[yearList,IA,IB] = intersect(LU_luh1.yearList,LU_plum.yearList) ;
fprintf('New yearList is %d-%d\n', min(yearList), max(yearList)) ;
stop
if ~isequal(yearList,LU_luh1.yearList)
    LU_luh1.maps_YXvy = LU_luh1.maps_YXvy(:,:,:,IA) ;
    LU_luh1.yearList = yearList ;
end
if ~isequal(yearList,LU_plum.yearList)
    LU_plum.maps_YXvy = LU_plum.maps_YXvy(:,:,:,IB) ;
    LU_plum.yearList = yearList ;
end

% Get equal varNames
tmp_YXvy = nan(size(LU_plum.maps_YXvy)) ;
for k = 1:length(plum_from_luh1)
    [~,IA] = intersect(LU_luh1.varNames,plum_from_luh1{k}) ;
    tmp_YXvy(:,:,k,:) = sum(LU_luh1.maps_YXvy(:,:,IA,:),3) ;
end
LU_luh1.varNames = LU_plum.varNames ;
LU_luh1.maps_YXvy = tmp_YXvy ;
clear tmp_YXvy

% Import land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = gcel_area_YXqd(:,1:2:1440) + gcel_area_YXqd(:,2:2:1440) ;
gcel_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
land_area_unmasked_YX = land_area_YX ;
%%% Convert to m2
land_area_YX = land_area_YX*1e6 ;
gcel_area_YX = gcel_area_YX*1e6 ;
clear tmp land_frac_YXqd land_area_YXqd

% Convert fractions to area
gcel_area_YXvy = repmat(gcel_area_YX,[1 1 length(LU_plum.varNames) length(yearList)]) ;
LU_luh1.maps_YXvy = LU_luh1.maps_YXvy .* gcel_area_YXvy ;
LU_plum.maps_YXvy = LU_plum.maps_YXvy .* gcel_area_YXvy ;


%% Make figure: Time series

%%%%%%%%%%%%%%%
% Options %%%%%
pct_change_baseyear = 2010 ;   % Leave empty for actual areas
abs_change_baseyear = [] ;   % Leave empty for actual areas
units = 'Mkm^2' ;
conv_fact = 1e-12 ;
thisPos = figurePos ;
spacing = [0.1 0.1] ;   % v, h
lineWidth = 3 ;
fontSize = 14 ;
%%%%%%%%%%%%%%%


ts_luh1_vy = squeeze(nansum(nansum(LU_luh1.maps_YXvy,1),2)) ;
ts_plum_vy = squeeze(nansum(nansum(LU_plum.maps_YXvy,1),2)) ;

do_zero_line = false ;
if ~isempty(pct_change_baseyear)
    base_vy = repmat(ts_luh1_vy(:,yearList==pct_change_baseyear), [1 size(ts_luh1_vy,2)]) ;
    ts_luh1_vy = 100 * (ts_luh1_vy - base_vy) ./ base_vy ;
    base_vy = repmat(ts_plum_vy(:,yearList==pct_change_baseyear), [1 size(ts_luh1_vy,2)]) ;
    ts_plum_vy = 100 * (ts_plum_vy - base_vy) ./ base_vy ;
    units = sprintf('Pct. change from %d', pct_change_baseyear) ;
    conv_fact = 1 ;
    do_zero_line = true ;
elseif ~isempty(abs_change_baseyear)
    base_vy = repmat(ts_luh1_vy(:,yearList==abs_change_baseyear), [1 size(ts_luh1_vy,2)]) ;
    ts_luh1_vy = ts_luh1_vy - base_vy ;
    base_vy = repmat(ts_plum_vy(:,yearList==abs_change_baseyear), [1 size(ts_luh1_vy,2)]) ;
    ts_plum_vy = ts_plum_vy - base_vy ;
    units = sprintf('Change from %d (%s)', abs_change_baseyear, units) ;
    do_zero_line = true ;
end

ts_luh1_vy = ts_luh1_vy * conv_fact ;
ts_plum_vy = ts_plum_vy * conv_fact ;

figure('Color','w','Position',thisPos) ;
for v = 1:length(LU_plum.varNames)
    h = subplot_tight(2,2,v,spacing) ;
    plot(yearList,[ts_luh1_vy(v,:) ; ts_plum_vy(v,:)], ...
        'LineWidth', lineWidth)
    
    ylims = h.YLim ;
    if do_zero_line && ylims(1)<0 && ylims(2)>0
        hold on
        plot(h.XLim, [0 0], ':k')
        hold off
    end
    
    legend({'LUH1','PLUM'},'Location','Best')
    title(LU_plum.varNames{v}) ;
    ylabel(units)
    set(h, 'FontSize',fontSize)
end


%% Make figure: Maps

thisLU = 'CROPLAND' ;

%%%%%%%%%%%%%%%
% Options %%%%%
baseyear = 2010 ;
units = 'km^2' ;
conv_fact = 1e-6 ;
thisPos = figurePos ;
spacing = [0.1 0.05] ;   % v, h
lineWidth = 3 ;
fontSize = 14 ;
yrange = 65:360 ;
%%%%%%%%%%%%%%%

k = find(strcmp(LU_plum.varNames,thisLU)) ;
diff_luh1_YX = conv_fact*(LU_luh1.maps_YXvy(yrange,:,k,yearList==2100) - LU_luh1.maps_YXvy(yrange,:,k,yearList==baseyear)) ;
diff_plum_YX = conv_fact*(LU_plum.maps_YXvy(yrange,:,k,yearList==2100) - LU_plum.maps_YXvy(yrange,:,k,yearList==baseyear)) ;

clims = [-1 1]*max(max(max(abs(diff_luh1_YX))),max(max(abs(diff_plum_YX)))) ;

figure('Color','w','Position',thisPos) ;

h = subplot_tight(1,2,1,spacing) ;
pcolor(diff_luh1_YX); shading flat; axis equal tight off
caxis(clims); colorbar('Location','SouthOutside')
title(sprintf('LUH1 change %d to 2100: %s', baseyear, thisLU))

h = subplot_tight(1,2,2,spacing) ;
pcolor(diff_plum_YX); shading flat; axis equal tight off
caxis(clims); colorbar('Location','SouthOutside')
title(sprintf('PLUM change %d to 2100: %s', baseyear, thisLU))





