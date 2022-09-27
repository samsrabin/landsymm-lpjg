%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figures for Almut 2022-07-27 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

thisDir = '/Users/Shared/landsymm_forestry/uc_work/landsymm/runs-forestonly/runs-2022-05/ukesm_actual_hist2015soc_default/potential/hist/1850pot_1848-2014/output-2022-06-30-201751' ;


% Make figure with time series and map

% thisFile = 'cmass_wood_sts.out.gz' ;
% thisTitle = 'Woody biomass' ;
% thisUnit = 'PgC' ;
% thisConv = 1e-12 ;
% thisUnit_map = 'kgC m^{-2}' ;
% thisConv_map = 1 ;

thisFile = 'cmass_sts.out.gz' ;
thisTitle = 'Live biomass' ;
thisUnit = 'PgC' ;
thisConv = 1e-12 ;
thisUnit_map = 'kgC m^{-2}' ;
thisConv_map = 1 ;

% thisFile = 'lai_natural.out.gz' ;
% thisTitle = 'Woody LAI' ;
% thisUnit = 'm^{2} m^{-2}' ;
% thisUnit_map = thisUnit ;
% thisConv = 1 ;
% thisConv_map = 1 ;

% thisFile = 'anpp_natural.out.gz' ;
% thisTitle = 'Woody NPP' ;
% thisUnit = 'PgC yr^{-1}' ;
% thisUnit_map = 'kgC m^{-2} yr^{-1}' ;
% thisConv = 1e-12 ;
% thisConv_map = 1 ;

data = lpjgu_matlab_read2geoArray(sprintf('%s/%s', thisDir, thisFile)) ;

% Import land area (km2)
landarea_file = '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc' ;
% Use gridcell area instead of land area
land_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
% Coarsen to half-degree
addpath '/Users/Shared/PLUM/crop_calib_code'
land_area_YX = aggregate_land_area(land_area_YXqd,0.5,0.5) ;
clear land_area_YXqd
% Convert to m2
land_area_YX = land_area_YX * 1e6 ;

if any(strcmp(data.varNames, 'ntrl'))
    data.garr_xvy = data.garr_xvy(:,strcmp(data.varNames, 'ntrl'),:) ;
    data.varNames = data.varNames(strcmp(data.varNames, 'ntrl')) ;
else
    [~, IA] = intersect(data.varNames, ...
        {'BNE', 'BINE', 'BNS', 'TeNE', 'TeBS', 'IBS', 'TeBE', 'TrBE', 'TrIBE', 'TrBR'}) ;
    data.garr_xvy = data.garr_xvy(:,IA,:) ;
    data.varNames = data.varNames(IA) ;
end

Nyears = length(data.yearList) ;
Ncells = length(data.list2map) ;


% Make figure

spacing = [0.08 0.06] ; % v h

f = figure('Color', 'w', 'Position', [1    41   947   960]) ;

subplot_tight(2,1,1,spacing) ;
land_area_x = land_area_YX(data.list2map) ;
data_y = squeeze(sum(sum(data.garr_xvy .* repmat(land_area_x, [1 1 Nyears]), 1), 2)) * thisConv ;
plot(data.yearList, data_y, ...
    'LineWidth', 3)
title(sprintf('%s (global total)', thisTitle))
ylabel(thisUnit)
set(gca, ...
    'FontSize', 14)

subplot_tight(2,1,2,spacing) ;
endYears = 10 ;
data_x = squeeze(mean(sum(data.garr_xvy(:,:,end-(endYears-1):end), 2), 3)) * thisConv_map;
data_YX = lpjgu_vector2map(data_x, size(land_area_YX), data.list2map) ;
pcolor(data_YX); shading flat; axis equal tight off
y1 = data.yearList(end-(endYears-1)) ;
yN = data.yearList(end) ;
title(sprintf('%s, %d-%d (%s)', thisTitle, y1, yN, thisUnit_map))
set(gca, ...
    'FontSize', 14)
hcb = colorbar(gca, 'Location', 'South') ;
xlabel(hcb, thisUnit_map)

export_fig(sprintf('%s/%s.png', thisDir, strrep(strrep(thisFile, '.gz', ''), '.out', '')), '-r150')
close





