

%% Setup

year0 = 2010 ;

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
% % Convert to 2-degree
% clear tmp
% tmp = landArea_YX(:,1:4:end) ...
%     + landArea_YX(:,2:4:end) ...
%     + landArea_YX(:,3:4:end) ...
%     + landArea_YX(:,4:4:end) ;
% landArea_2deg_YX = tmp(1:4:end,:) ...
%                  + tmp(2:4:end,:) ...
%                  + tmp(3:4:end,:) ...
%                  + tmp(4:4:end,:) ;
% clear tmp

% Import LUH2
luh2_file = '/Users/Shared/PLUM/input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
luh2 = lpjgu_matlab_readTable_then2map(luh2_file) ;
luh2.maps_YXvy = luh2.maps_YXvy(:,:,~contains(luh2.varNames,{'URBAN','PEATLAND'}),:) ;
luh2.varNames = luh2.varNames(~contains(luh2.varNames,{'URBAN','PEATLAND'})) ;
LUnames = luh2.varNames ;
Nlu = length(LUnames) ;
% landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;
% landArea_2deg_YXv = repmat(landArea_2deg_YX,[1 1 Nlu]) ;

% % Extract what's needed; aggregate to 2deg
% luh2_2010_YXv = luh2.maps_YXvy(:,:,:,luh2.yearList==2010) ;
% luh2_2010_YXv(isnan(luh2_2010_YXv)) = 0 ;
% luh2_2010_YXv = luh2_2010_YXv .* landArea_YXv ;
% tmp = luh2_2010_YXv(:,1:4:end,:) ...
%     + luh2_2010_YXv(:,2:4:end,:) ...
%     + luh2_2010_YXv(:,3:4:end,:) ...
%     + luh2_2010_YXv(:,4:4:end,:) ;
% luh2_2010_2deg_YXv = tmp(1:4:end,:,:) ...
%                + tmp(2:4:end,:,:) ...
%                + tmp(3:4:end,:,:) ...
%                + tmp(4:4:end,:,:) ;
clear tmp

disp('Done importing LUH2.')

%%
[test1,~] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm/2011/LandCoverFract.base2010.orig.txt'],...
        landArea_YX, LUnames, []) ;
    
[test2,~] = PLUMharm_processPLUMin(...
        ['/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1.harm/2011/LandCoverFract.base2010.txt'],...
        landArea_YX, LUnames, []) ;
    
minmax_ssr(test2.maps_YXv - test1.maps_YXv)


%%

spacing = 0.03 ;

figure('Color','w','Position',figurePos) ;

for v = 1:Nlu
    
    subplot_tight(2,2,v,spacing) ;
    tmp = test2.maps_YXv(:,:,v) - test1.maps_YXv(:,:,v) ;
    tmp(landArea_YX==0) = NaN ;
    pcolor(tmp)
    shading flat; axis equal tight off
    colorbar ;
    title(LUnames{v})
    
end


