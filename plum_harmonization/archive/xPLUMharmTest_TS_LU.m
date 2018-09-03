%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing harmonized PLUM land use trajectory %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisVer = '' ;
% thisVer = 'v2.' ;


%% Setup

topDir_orig = addslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1') ;
topDir_harm = addslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm') ;


%% Import LUH2

yearList = 2011:2050 ;
Nyears = length(yearList) ;

% Get LUH2 land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
landArea_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
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


%% Import PLUM (original + harmonized)

for y = 1:Nyears
    thisYear = yearList(y) ;
    disp(['Reading ' num2str(thisYear) '...'])
    inDir_orig = addslashifneeded([topDir_orig num2str(thisYear)]) ;
    inDir_harm = addslashifneeded([topDir_harm num2str(thisYear)]) ;
    
    % Original
    LU_orig = lpjgu_matlab_readTable_then2map([inDir_orig 'LandCoverFract.txt'],'verboseIfNoMat',false,'force_mat_nosave',true) ;
    LU_orig.maps_YXv(:,:,strcmp(LU_orig.varNames,'BARREN')) = ...
        sum(LU_orig.maps_YXv(:,:,contains(LU_orig.varNames,{'BARREN','URBAN'})),3) ;
    LU_orig.maps_YXv(:,:,strcmp(LU_orig.varNames,'URBAN')) = [] ;
    LU_orig.varNames(strcmp(LU_orig.varNames,'URBAN')) = [] ;
    
    % Harmonized
    LU_harm = lpjgu_matlab_readTable_then2map([inDir_harm 'LandCoverFract.base2010.' thisVer 'txt'],'verboseIfNoMat',false,'force_mat_nosave',true) ;
    
    % Setup
    if y==1
        LU_orig_YXvy = nan([size(LU_orig.maps_YXv) Nyears]) ;
        LU_harm_YXvy = nan([size(LU_harm.maps_YXv) Nyears]) ;
        landArea_YXv = repmat(landArea_YX,[1 1 size(LU_harm.maps_YXv,3)]) ;
        mask_orig_YX = isnan(LU_orig.maps_YXv(:,:,1)) ...
            | sum(LU_orig.maps_YXv(:,:,contains(LU_orig.varNames,{'CROPLAND','PASTURE','NATURAL'})),3)==0 ;
        mask_harm_YX = isnan(LU_harm.maps_YXv(:,:,1)) ...
            | sum(LU_harm.maps_YXv(:,:,contains(LU_harm.varNames,{'CROPLAND','PASTURE','NATURAL'})),3)==0 ;
    end
    
    % Get area (reordered)
    [~,~,IB_orig] = intersect(LUnames,LU_orig.varNames,'stable') ;
    [~,~,IB_harm] = intersect(LUnames,LU_harm.varNames,'stable') ;
    LU_orig_YXvy(:,:,:,y) = LU_orig.maps_YXv(:,:,IB_orig) .* landArea_YX ;
    LU_harm_YXvy(:,:,:,y) = LU_harm.maps_YXv(:,:,IB_harm) .* landArea_YX ;
    clear LU_orig LU_harm
end

ts_orig = squeeze(nansum(nansum(LU_orig_YXvy,1),2)) ;
ts_harm = squeeze(nansum(nansum(LU_harm_YXvy,1),2)) ;

disp('Done reading PLUM.')


%% Get LUH2 area

landArea_YX(mask_orig_YX) = 0 ;
luh2_area_YXvy = luh2.maps_YXvy .* repmat(landArea_YX,[1 1 size(luh2.maps_YXvy,3) size(luh2.maps_YXvy,4)]) ;
ts_luh2 = squeeze(nansum(nansum(luh2_area_YXvy,1),2)) ;


%% Plot
spacing = 0.1 ;

figure('Color','w','Position',figurePos)

for v = 1:4
    subplot_tight(2,2,v,spacing) ;
    plot(1850:2010,ts_luh2(v,1:161)*1e-6,'-k','LineWidth',2) ;
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList,ts_orig(v,:)*1e-6,'-','LineWidth',2)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList,ts_harm(v,:)*1e-6,'--','LineWidth',2)
    hold off
    title(LUnames{v})
    set(gca,'FontSize',14)
    ylabel('Million km2')
    legend({'LUH2','PLUM orig','PLUM harm'},'Location','NorthWest')
end


%% Map difference
spacing = 0.025 ;

ntrl_luh2_2010_YX = luh2_area_YXvy(:,:,strcmp(LUnames,'NATURAL'),luh2.yearList==2010) ;
ntrl_plum_2011_orig_YX = LU_orig_YXvy(:,:,strcmp(LUnames,'NATURAL'),1) ;
ntrl_plum_2011_harm_YX = LU_harm_YXvy(:,:,strcmp(LUnames,'NATURAL'),1) ;

% Aggregate LUH2 and PLUM_orig to 2deg
ntrl_luh2_2010_YX(isnan(ntrl_luh2_2010_YX)) = 0 ;
tmp = ntrl_luh2_2010_YX(:,1:4:end) ...
    + ntrl_luh2_2010_YX(:,2:4:end) ...
    + ntrl_luh2_2010_YX(:,3:4:end) ...
    + ntrl_luh2_2010_YX(:,4:4:end) ;
ntrl_luh2_2010_YX = tmp(1:4:end,:) ...
                  + tmp(2:4:end,:) ...
                  + tmp(3:4:end,:) ...
                  + tmp(4:4:end,:) ;
ntrl_plum_2011_orig_YX(isnan(ntrl_plum_2011_orig_YX)) = 0 ;
tmp = ntrl_plum_2011_orig_YX(:,1:4:end) ...
    + ntrl_plum_2011_orig_YX(:,2:4:end) ...
    + ntrl_plum_2011_orig_YX(:,3:4:end) ...
    + ntrl_plum_2011_orig_YX(:,4:4:end) ;
ntrl_plum_2011_orig_YX = tmp(1:4:end,:) ...
                  + tmp(2:4:end,:) ...
                  + tmp(3:4:end,:) ...
                  + tmp(4:4:end,:) ;
              
% Get % diffs
diff_orig_YX = 100 * (ntrl_plum_2011_orig_YX - ntrl_luh2_2010_YX)./ntrl_luh2_2010_YX ;
diff_orig_YX(ntrl_luh2_2010_YX==0 & ntrl_plum_2011_orig_YX==0) = 0 ;
diff_orig_YX(ntrl_luh2_2010_YX==0 & ntrl_plum_2011_orig_YX>0) = Inf ;
diff_harm_YX = 100 * (ntrl_plum_2011_harm_YX - ntrl_luh2_2010_YX)./ntrl_luh2_2010_YX ;
diff_harm_YX(ntrl_luh2_2010_YX==0 & ntrl_plum_2011_harm_YX==0) = 0 ;
diff_harm_YX(ntrl_luh2_2010_YX==0 & ntrl_plum_2011_harm_YX>0) = Inf ;
diff_orig_YX(isnan(ntrl_plum_2011_harm_YX)) = NaN ;

% Make figure
figure ;
subplot_tight(2,1,1,spacing)
pcolor(diff_orig_YX); shading flat; axis equal tight off
caxis([-100 100])
colorbar
subplot_tight(2,1,2,spacing)
pcolor(diff_harm_YX); shading flat; axis equal tight off
caxis([-100 100])
colorbar


