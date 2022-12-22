%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Is biomass decreasing?? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bart noticed that aboveground tree biomass was decreasing, e.g., -104.75,
% 59.25.

% thisDir = '/Users/Shared/landsymm_forestry/uc_work/landsymm/2022-03-08/1970past.forPLUM.20220305002322/1850-2014' ;
% thisDir  = '/Users/Shared/landsymm_forestry/uc_work/landsymm/2022-03-08/1970past_trunkmerge20220309/1850-2014' ;
thisDir  = '/Users/Shared/landsymm_forestry/uc_work/landsymm/2022-03-08/1850past.forPLUM.20220302184520/1850-2014' ;

%% Import

% ctree_ag_sts = lpjgu_matlab_read2geoArray(sprintf('%s/landsymm_ctree_ag_sts.out.gz', thisDir)) ;
% ctree_sts = lpjgu_matlab_read2geoArray(sprintf('%s/landsymm_ctree_sts.out.gz', thisDir)) ;
% ctotal = lpjgu_matlab_read2geoArray(sprintf('%s/ctotal_sts.out.gz', thisDir)) ;
% cmass_sts = lpjgu_matlab_read2geoArray(sprintf('%s/cmass_sts.out.gz', thisDir)) ;
% cmass = lpjgu_matlab_read2geoArray(sprintf('%s/cmass.out.gz', thisDir)) ;
% cmass_wood_sts = lpjgu_matlab_read2geoArray(sprintf('%s/cmass_wood_sts.out.gz', thisDir)) ;
plutW_from_forC = lpjgu_matlab_read2geoArray(sprintf('%s/landsymm_plutW_from_forC.out.gz', thisDir)) ;

% Bug fix (not needed as of commit 55b91c)
if strcmp(thisDir, '/Users/Shared/landsymm_forestry/uc_work/landsymm/2022-03-08/1970past.forPLUM.20220305002322/1850-2014')
    ctree_ag_sts.garr_xvy = ctree_ag_sts.garr_xvy / 20 ;
end

yearList = plutW_from_forC.yearList ;
Nyears = length(yearList) ;

% Import cell area (km2)
cell_area_YXqd = transpose(ncread( ...
    '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc', ...
    'carea')) ;
addpath(genpath(landsymm_lpjg_path()))
cell_area_YX = aggregate_land_area(cell_area_YXqd,0.5,0.5) ;
cell_area_x = cell_area_YX(plutW_from_forC.list2map) ;
clear cell_area_YX*


%% One cell

thisLU = 'forC' ;
lineWidth = 2 ;
thisLon = 13.75 ;
thisLat = 55.75 ;

is_thisCell = plutW_from_forC.lonlats(:,1)==thisLon & plutW_from_forC.lonlats(:,2)==thisLat ;
is_thisCol = strcmp(plutW_from_forC.varNames, 'to_forC') ;

% cmass_sts_thisCell_y = squeeze(cmass_sts.garr_xvy(is_thisCell, is_thisCol, :)) ;
plutW_from_forC_thisCell_y = squeeze(plutW_from_forC.garr_xvy(is_thisCell, is_thisCol, :)) ;

figure('Color', 'w', 'Position', [865 706 883 412]) ;
% plot(yearList, cmass_sts_thisCell_y, 'LineWidth', lineWidth)
% hold on
plot(yearList, plutW_from_forC_thisCell_y, 'LineWidth', lineWidth)
% hold off
% legend('cmass\_sts', 'cmass\_wood\_potharv\_sts', ...
%     'Location', 'Best')
ylabel('kgC m^{-2}')
title(sprintf('lon %0.2f, lat %0.2f', thisLon, thisLat))
set(gca, 'FontSize', 14)


%% Global

thisLU = 'ntrl' ;
lineWidth = 2 ;

ctree_sts_global_y = area_weighted_mean(ctree_sts.garr_xvy, ctree_sts.varNames, thisLU, cell_area_x) ;
ctree_ag_sts_global_y = area_weighted_mean(ctree_ag_sts.garr_xvy, ctree_sts.varNames, thisLU, cell_area_x) ;
% ctotal_sts_global_y = area_weighted_mean(ctotal.garr_xvy, ctree_sts.varNames, thisLU, cell_area_x) ;
cmass_sts_global_y = area_weighted_mean(cmass_sts.garr_xvy, ctree_sts.varNames, thisLU, cell_area_x) ;
cmass_global_total_y = area_weighted_mean(cmass.garr_xvy, cmass.varNames, 'Total', cell_area_x) ;
cmass_global_ntrl_y = area_weighted_mean(cmass.garr_xvy, cmass.varNames, 'Natural_sum', cell_area_x) ;
cmass_wood_sts_global_y = area_weighted_mean(cmass_wood_sts.garr_xvy, ctree_sts.varNames, thisLU, cell_area_x) ;

figure('Color', 'w', 'Position', [865 706 883 412]) ;
plot(ctree_sts.yearList, ctree_sts_global_y, 'LineWidth', lineWidth)
hold on
plot(ctree_ag_sts.yearList, ctree_ag_sts_global_y, 'LineWidth', lineWidth)
% plot(ctree_ag.yearList, ctotal_sts_global_y, 'LineWidth', lineWidth)
plot(ctree_ag_sts.yearList, cmass_sts_global_y, 'LineWidth', lineWidth)
plot(ctree_ag_sts.yearList, cmass_wood_sts_global_y, 'LineWidth', lineWidth)
plot(cmass.yearList, cmass_global_total_y, 'LineWidth', lineWidth)
plot(cmass.yearList, cmass_global_ntrl_y, '--', 'LineWidth', lineWidth)
hold off
% legend('ctree\_sts', 'ctree\_ag\_sts', 'ctotal\_sts', 'cmass\_sts', 'cmass', 'cmass\_wood\_sts', ...
%     'Location', 'West')
legend('ctree\_sts', 'ctree\_ag\_sts', 'cmass\_sts', 'cmass\_wood\_sts', 'cmass\_total', 'cmass\_ntrl', ...
    'Location', 'Best')
ylabel('kgC m^{-2}')
set(gca, 'FontSize', 14)


%% FUNCTIONS

function out = area_weighted_mean(in_xvy, varNames, thisLU, cell_area_x)

Nyears = size(in_xvy,3) ;

out = squeeze(sum(in_xvy(:,strcmp(varNames, thisLU), :) ...
    .* repmat(cell_area_x/sum(cell_area_x), [1 1 Nyears]),1)) ;

end