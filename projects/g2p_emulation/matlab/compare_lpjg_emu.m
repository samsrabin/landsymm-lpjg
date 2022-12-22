%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare LPJ-GUESS and emulator outputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% thisEmu = 'LPJ-GUESS' ;
thisEmu = 'LPJmL' ;
% thisEmu = 'pDSSAT' ;


%% Read files

yield_lpj1 = lpjgu_matlab_read2geoArray('/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/LPJGPLUM_2001-2100_remap6p7_forPotYields_rcp45/LPJGPLUM_2001-2100_remap6p7_forPotYields_rcp45_forED_20190225101539/2011-2015/yield.out.gz') ;
yield_lpj2 = lpjgu_matlab_read2geoArray('/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/LPJGPLUM_2001-2100_remap6p7_forPotYields_rcp45/LPJGPLUM_2001-2100_remap6p7_forPotYields_rcp45_forED_20190225101539/2016-2020/yield.out.gz') ;
yield_lpj = yield_lpj1 ;
yield_lpj.garr_xv = nanmean(cat(3,yield_lpj1.garr_xv,yield_lpj2.garr_xv),3) ;
clear yield_lpj1 yield_lpj2

yield_emu = lpjgu_matlab_read2geoArray(...
    sprintf(...
    '/Users/Shared/GGCMI2PLUM/emulator/Sam/outputs_20190830/IPSL-CM5A-MR/%s/rcp45/2011-2020/yield.out', ...
    thisEmu)) ;


%% Compare for a given crop

thisVar_lpj = 'CerealsC40' ;
thisVar_emu = 'CerealsC4010' ;
% thisVar_lpj = 'CerealsC40200' ;
% thisVar_emu = 'CerealsC4200' ;
% thisVar_lpj = 'Rice0' ;
% thisVar_emu = 'Rice010' ;
% thisVar_lpj = 'Rice0200' ;
% thisVar_emu = 'Rice200' ;
% thisVar_lpj = 'Ricei0200' ;
% thisVar_emu = 'Ricei200' ;


%%%% Options
spacing = [0.025 0.025] ; % v h
yrange = 70:360 ;
thisPos = [1    34   720   771] ;
%%%%

[cf_lpj, cf_emu] = get_cf(thisVar_lpj, thisEmu) ;

tmp_lpj = cf_lpj * yield_lpj.garr_xv(:,strcmp(yield_lpj.varNames,thisVar_lpj)) ;
tmp_emu = cf_emu * yield_emu.garr_xv(:,strcmp(yield_emu.varNames,thisVar_emu)) ;
new_caxis = [0 max(max(tmp_lpj),max(tmp_emu))] ;

tmp_lpj_YX = nan(360,720) ;
tmp_lpj_YX(yield_lpj.list2map) = tmp_lpj ;
tmp_emu_YX = nan(360,720) ;
tmp_emu_YX(yield_emu.list2map) = tmp_emu ;

figure('Color','w','Position',thisPos)

subplot_tight(2,1,1,spacing)
pcolor(tmp_lpj_YX(yrange,:)); shading flat; axis equal tight off
caxis(new_caxis) ; colormap('jet'); colorbar('Location','SouthOutside')
title(sprintf('LPJ-GUESS: %s', thisVar_lpj))

subplot_tight(2,1,2,spacing)
pcolor(tmp_emu_YX(yrange,:)); shading flat; axis equal tight off
caxis(new_caxis) ; colormap('jet'); colorbar('Location','SouthOutside')
title(sprintf('%s emulator: %s', thisEmu, thisVar_emu))


%% Map N sensitivity

thisCrop = 'Pulsesi' ;
%%%% Options
spacing = [0.025 0.025] ; % v h
yrange = 70:360 ;
thisPos = [1    34   720   771] ;
% new_caxis = [] ;
% new_caxis = 100*[-1 1] ;
new_caxis = [-100 250] ;
%%%%

tmp_lpj_lo = yield_lpj.garr_xv(:,strcmp(yield_lpj.varNames,[thisCrop '0'])) ;
tmp_lpj_hi = yield_lpj.garr_xv(:,strcmp(yield_lpj.varNames,[thisCrop '0200'])) ;
tmp_emu_lo = yield_emu.garr_xv(:,strcmp(yield_emu.varNames,[thisCrop '010'])) ;
tmp_emu_hi = yield_emu.garr_xv(:,strcmp(yield_emu.varNames,[thisCrop '200'])) ;

tmp_lpj = 100*(tmp_lpj_hi-tmp_lpj_lo)./tmp_lpj_lo ; tmp_lpj(tmp_lpj_lo==0) = NaN ;
tmp_emu = 100*(tmp_emu_hi-tmp_emu_lo)./tmp_emu_lo ; tmp_emu(tmp_emu_lo==0) = NaN ;

tmp_lpj_YX = get_map(tmp_lpj, yield_lpj.list2map) ;
tmp_emu_YX = get_map(tmp_emu, yield_emu.list2map) ;

figure('Color','w','Position',thisPos)

subplot_tight(2,1,1,spacing)
pcolor(tmp_lpj_YX(yrange,:)); shading flat; axis equal tight off
if ~isempty(new_caxis); caxis(new_caxis) ; end
colormap('jet'); colorbar('Location','SouthOutside')
title(sprintf('LPJ-GUESS: %s, 0 to 200 kgN/ha', thisCrop))

subplot_tight(2,1,2,spacing)
pcolor(tmp_emu_YX(yrange,:)); shading flat; axis equal tight off
if ~isempty(new_caxis); caxis(new_caxis) ; end 
colormap('jet'); colorbar('Location','SouthOutside')
title(sprintf('%s emulator: %s 10 to 200 kgN/ha', thisEmu, thisCrop))





%% FUNCTIONS

function [cf_lpj, cf_emu] = get_cf(thisVar_lpj, thisEmu)

thisCrop = thisVar_lpj ;

xtra = regexp(thisCrop,'(i|0|1|2)+$') ;
if ~isempty(xtra)
    thisCrop(xtra:end) = [] ;
end

crop_list = {'CerealsC3','CerealsC4','Rice','Oilcrops','Pulses','StarchyRoots'} ;
thisCrop_i = find(strcmp(crop_list, thisCrop)) ;
if length(thisCrop_i) ~= 1
    error('Error finding index of thisCrop')
end

cf_lpj_list = [1.046 0.654 0.972 0.578 0.686 4.560] ;
switch lower(thisEmu)
    case 'lpj-guess'; cf_emu_list = [1.078 0.874 1.815 0.711 1.084 5.748] ;
    case 'lpjml';     cf_emu_list = [0.682 0.716 1.420 0.592 0.642 4.071] ;
    case 'pdssat';    cf_emu_list = [0.710 0.562 0.483 0.402 0.425 4.579] ;
    otherwise ; error('%s not recognized in code to get calib factors', thisEmu)
end

if length(cf_lpj_list) ~= length(crop_list)
    error('Mismatch in length of crop_list and cf_lpj_list')
elseif length(cf_emu_list) ~= length(crop_list)
    error('Mismatch in length of crop_list and cf_emu_list')
end

cf_lpj = cf_lpj_list(thisCrop_i) ;
cf_emu = cf_emu_list(thisCrop_i) ;

end


function out_YX = get_map(in_x, list2map)

out_YX = nan(360,720) ;
out_YX(list2map) = in_x ;

end
