%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Are we getting the deltas we expect? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% thisVer = 'orig.' ;
thisVer = '' ;
% thisVer = '2deg.' ;

year0 = 2048 ;


%% Setup

topDir_orig = addslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1') ;
topDir_harm = addslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm') ;

addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;

year1 = year0 + 1 ;
year2 = year1 + 1 ;
year0s = num2str(year0) ;
year1s = num2str(year1) ;
year2s = num2str(year2) ;


%% Import LUH2

disp('Importing LUH2...')

% Get LUH2 land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
landArea_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
% Harmonize LUH2 mask and PLUM mask
file_in = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/2011/LandCoverFract.txt' ;
S = lpjgu_matlab_readTable_then2map(file_in,'verboseIfNoMat',false,'force_mat_nosave',true) ;
mask_YX = isnan(S.maps_YXv(:,:,1)) | landArea_YX==0 ;
landArea_YX(mask_YX) = 0 ;
clear S
% Convert to 2-degree
clear tmp
tmp = landArea_YX(:,1:4:end) ...
    + landArea_YX(:,2:4:end) ...
    + landArea_YX(:,3:4:end) ...
    + landArea_YX(:,4:4:end) ;
landArea_2deg_YX = tmp(1:4:end,:) ...
                 + tmp(2:4:end,:) ...
                 + tmp(3:4:end,:) ...
                 + tmp(4:4:end,:) ;
clear tmp

% Import LUH2
luh2_file = '/Users/Shared/PLUM/input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
luh2 = lpjgu_matlab_readTable_then2map(luh2_file) ;
luh2.maps_YXvy = luh2.maps_YXvy(:,:,~contains(luh2.varNames,{'URBAN','PEATLAND'}),:) ;
luh2.varNames = luh2.varNames(~contains(luh2.varNames,{'URBAN','PEATLAND'})) ;
LUnames = luh2.varNames ;
Nlu = length(LUnames) ;
landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;
landArea_2deg_YXv = repmat(landArea_2deg_YX,[1 1 Nlu]) ;

% Extract what's needed; aggregate to 2deg
luh2_2010_YXv = luh2.maps_YXvy(:,:,:,luh2.yearList==2010) ;
luh2_2010_YXv(isnan(luh2_2010_YXv)) = 0 ;
luh2_2010_YXv = luh2_2010_YXv .* landArea_YXv ;
tmp = luh2_2010_YXv(:,1:4:end,:) ...
    + luh2_2010_YXv(:,2:4:end,:) ...
    + luh2_2010_YXv(:,3:4:end,:) ...
    + luh2_2010_YXv(:,4:4:end,:) ;
luh2_2010_2deg_YXv = tmp(1:4:end,:,:) ...
               + tmp(2:4:end,:,:) ...
               + tmp(3:4:end,:,:) ...
               + tmp(4:4:end,:,:) ;
clear tmp

disp('Done importing LUH2.')


%% Import PLUM

disp('Importing PLUM...')

is2deg = strcmp(thisVer,'2deg.') ;

% Originals
if is2deg
    [~,PLUM_orig_year0] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year0s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, []) ;
    [~,PLUM_orig_year1] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year1s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, []) ;
    [~,PLUM_orig_year2] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year2s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, []) ;
else
    [PLUM_orig_year0,~] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year0s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, []) ;
    [PLUM_orig_year1,~] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year1s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, []) ;
    [PLUM_orig_year2,~] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/' year2s '/LandCoverFract.txt'],...
        landArea_YX, LUnames, []) ;
end

% Harmonized
if is2deg
    [~,PLUM_harm_year1] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm/' year1s '/LandCoverFract.base2010.' thisVer 'txt'],...
        landArea_2deg_YX, LUnames, []) ;
    [~,PLUM_harm_year2] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm/' year2s '/LandCoverFract.base2010.' thisVer 'txt'],...
        landArea_2deg_YX, LUnames, []) ;
else
    [PLUM_harm_year1,~] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm/' year1s '/LandCoverFract.base2010.' thisVer 'txt'],...
        landArea_YX, LUnames, []) ;
    [PLUM_harm_year2,~] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm/' year2s '/LandCoverFract.base2010.' thisVer 'txt'],...
        landArea_YX, LUnames, []) ;
end

if is2deg
    thisLandArea_YX = landArea_2deg_YX ;
    Yincl = 15:90 ;
    thisluh2_2010_YXv = luh2_2010_2deg_YXv ;
else
    thisLandArea_YX = landArea_YX ;
    thisluh2_2010_YXv = luh2_2010_YXv ;
    Yincl = 60:360 ;
end
thisLandArea_YXv = repmat(thisLandArea_YX,[1 1 Nlu]) ;

disp('Done importing PLUM.')


%% Map transitions for 2010 to 2011, including original L2010-->P2011
spacing = 0.025 ;
thisPos = figurePos ;
do_debug = true ;
thisCell = [67 42] ; thisV = 2 ; % 2-degree, cropland gain where previously none

if year0~=2010
    error('This doesn''t make sense to do without year0==2010.')
end

figure('Color','w','Position',thisPos) ;

h = cell(12,1) ;
for v = 1:4
    h{v} = subplot_tight(3,4,v,spacing) ;
    thisDiff_YX = PLUM_orig_year1.maps_YXv(:,:,v) ...
                - PLUM_orig_year0.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
%     pcolor(thisDiff_YX(Yincl,:))
    pcolor(thisDiff_YX)
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Original year1-year0: ' LUnames{v}])
    set(gca,'FontSize',14)
    if do_debug && v==thisV
        disp(['PLUM_orig_year0 = ' num2str(PLUM_orig_year0.maps_YXv(thisCell(1),thisCell(2),v))])
        disp(['PLUM_orig_year1 = ' num2str(PLUM_orig_year1.maps_YXv(thisCell(1),thisCell(2),v))])
        disp(['thisDiff = ' num2str(thisDiff_YX(thisCell(1),thisCell(2)))])
    end
end
for v = 1:4
    h{4+v} = subplot_tight(3,4,4+v,spacing) ;
    thisDiff_YX = PLUM_harm_year1.maps_YXv(:,:,v) ...
                - thisluh2_2010_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
%     pcolor(thisDiff_YX(Yincl,:))
    pcolor(thisDiff_YX)
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized year1-year0: ' LUnames{v}])
    set(gca,'FontSize',14)
    if do_debug && v==thisV
        disp(['luh2_year0      = ' num2str(thisluh2_2010_YXv(thisCell(1),thisCell(2),v))])
        disp(['PLUM_harm_year1 = ' num2str(PLUM_harm_year1.maps_YXv(thisCell(1),thisCell(2),v))])
        disp(['thisDiff = ' num2str(thisDiff_YX(thisCell(1),thisCell(2)))])
    end
end
for v = 1:4
    h{8+v} = subplot_tight(3,4,8+v,spacing) ;
    thisDiff_YX = PLUM_orig_year1.maps_YXv(:,:,v) ...
                - thisluh2_2010_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Oyear1-Lyear0: ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    new_caxis = get_combined_caxis(h{v},h{4+v}) ;
    caxis(h{v},new_caxis) ;
    caxis(h{4+v},new_caxis) ;
end



%% Map transitions for year0 to year1
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

h = cell(8,1) ;
for v = 1:4
    h{v} = subplot_tight(2,4,v,spacing) ;
    thisDiff_YX = PLUM_orig_year1.maps_YXv(:,:,v) ...
                - PLUM_orig_year0.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Original ' year1s '-' year0s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    h{4+v} = subplot_tight(2,4,4+v,spacing) ;
    thisDiff_YX = PLUM_harm_year1.maps_YXv(:,:,v) ...
                - thisluh2_2010_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized ' year1s '-' year0s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
% for v = 1:4
%     new_caxis = get_combined_caxis(h{v},h{4+v}) ;
%     caxis(h{v},new_caxis) ;
%     caxis(h{4+v},new_caxis) ;
% end

%% Map difference in transitions for year0-year1
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

for v = 1:4
    subplot_tight(1,4,v,spacing)
    thisDiffO_YX = PLUM_orig_year1.maps_YXv(:,:,v) ...
                 - PLUM_orig_year0.maps_YXv(:,:,v) ;
    thisDiffH_YX = PLUM_harm_year1.maps_YXv(:,:,v) ...
                 - thisluh2_2010_YXv(:,:,v) ;
    thisDiff_YX = thisDiffH_YX - thisDiffO_YX ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['(H' year1s '-Lyear0)-(O' year1s '-Oyear0): ' LUnames{v}])
    set(gca,'FontSize',14)
end


%% Map transitions for year1 to year2
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

h = cell(8,1) ;
for v = 1:4
    h{v} = subplot_tight(2,4,v,spacing) ;
    thisDiff_YX = PLUM_orig_year2.maps_YXv(:,:,v) ...
                - PLUM_orig_year1.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Original ' year2s '-' year1s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    h{4+v} = subplot_tight(2,4,4+v,spacing) ;
    thisDiff_YX = PLUM_harm_year2.maps_YXv(:,:,v) ...
                - PLUM_harm_year1.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized ' year2s '-' year1s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    new_caxis = get_combined_caxis(h{v},h{4+v}) ;
    caxis(h{v},new_caxis) ;
    caxis(h{4+v},new_caxis) ;
end


%% Map difference in transitions for year1-year2
spacing = 0.025 ;

thisPos = figurePos ;

figure('Color','w','Position',thisPos) ;

for v = 1:4
    subplot_tight(2,2,v,spacing)
    thisDiffO_YX = PLUM_orig_year2.maps_YXv(:,:,v) ...
                 - PLUM_orig_year1.maps_YXv(:,:,v) ;
    thisDiffH_YX = PLUM_harm_year2.maps_YXv(:,:,v) ...
                 - PLUM_harm_year1.maps_YXv(:,:,v) ;
    thisDiff_YX = thisDiffH_YX - thisDiffO_YX ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['(H' year2s '-H' year1s ')-(O' year2s '-O' year1s '): ' LUnames{v}])
    set(gca,'FontSize',14)
end


%% Maps in year1
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

h = cell(8,1) ;
for v = 1:4
    h{v} = subplot_tight(2,4,v,spacing) ;
    thisDiff_YX = PLUM_orig_year1.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Original ' year1s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    h{4+v} = subplot_tight(2,4,4+v,spacing) ;
    thisDiff_YX = PLUM_harm_year1.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized ' year1s ': ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    new_caxis = get_combined_caxis(h{v},h{4+v}) ;
    caxis(h{v},new_caxis) ;
    caxis(h{4+v},new_caxis) ;
end


%% Maps in year2
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

h = cell(8,1) ;
for v = 1:4
    h{v} = subplot_tight(2,4,v,spacing) ;
    thisDiff_YX = PLUM_orig_year2.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Original year2: ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    h{4+v} = subplot_tight(2,4,4+v,spacing) ;
    thisDiff_YX = PLUM_harm_year2.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized year2: ' LUnames{v}])
    set(gca,'FontSize',14)
end
for v = 1:4
    new_caxis = get_combined_caxis(h{v},h{4+v}) ;
    caxis(h{v},new_caxis) ;
    caxis(h{4+v},new_caxis) ;
end

%% Map differences in year1
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

for v = 1:4
    subplot_tight(1,4,v,spacing)
    thisDiff_YX = PLUM_harm_year1.maps_YXv(:,:,v) ...
                - PLUM_orig_year1.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized-Original year1: ' LUnames{v}])
    set(gca,'FontSize',14)
end


%% Map differences in year2
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

for v = 1:4
    subplot_tight(1,4,v,spacing)
    thisDiff_YX = PLUM_harm_year2.maps_YXv(:,:,v) ...
                - PLUM_orig_year2.maps_YXv(:,:,v) ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['Harmonized-Original year2: ' LUnames{v}])
    set(gca,'FontSize',14)
end


%% Map "difference in differences" of year1 vs year2
% How persistent are the differences induced by the harmonization to year0?
spacing = 0.025 ;

thisPos = [1 313 1440 492] ;

figure('Color','w','Position',thisPos) ;

for v = 1:4
    subplot_tight(1,4,v,spacing)
    thisDiff1_YX = PLUM_harm_year1.maps_YXv(:,:,v) ...
                 - PLUM_orig_year1.maps_YXv(:,:,v) ;
    thisDiff2_YX = PLUM_harm_year2.maps_YXv(:,:,v) ...
                 - PLUM_orig_year2.maps_YXv(:,:,v) ;
    thisDiff_YX = thisDiff2_YX - thisDiff1_YX ;
    thisDiff_YX(thisLandArea_YX==0) = NaN ;
    pcolor(thisDiff_YX(Yincl,:))
    shading flat
    axis equal tight off
    colorbar('Location','SouthOutside')
    title(['(H-O year2) - (H-O year1): ' LUnames{v}])
    set(gca,'FontSize',14)
end

