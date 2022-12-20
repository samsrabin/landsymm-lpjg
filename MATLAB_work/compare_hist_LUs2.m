%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare two landuse input data sets %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lu_file_1 = '/Users/Shared/PLUM/input/remaps_v6p3/LU.remapv6p3.txt' ;
% cf_file_1 = '/Users/Shared/PLUM/input/remaps_v6p3/cropfracs.remapv6p3.txt' ;
% nf_file_1 = '/Users/Shared/PLUM/input/remaps_v6p3/nfert.remapv6p3.txt' ;
% lu_file_2 = '/Users/Shared/PLUM/input/remaps_v6p6/LU.remapv6p6.txt' ;
% cf_file_2 = '/Users/Shared/PLUM/input/remaps_v6p6/cropfracs.remapv6p6.txt' ;
% nf_file_2 = '/Users/Shared/PLUM/input/remaps_v6p6/nfert.remapv6p6.txt' ;
% outDir = '/Users/Shared/PLUM/input/remaps_v6p3_versus_v6p6/' ;

lu_file_1 = '/Users/Shared/PLUM/input/remaps_v7a/LU.remapv7a.txt' ;
cf_file_1 = '/Users/Shared/PLUM/input/remaps_v7a/cropfracs.remapv7a.txt' ;
nf_file_1 = '/Users/Shared/PLUM/input/remaps_v7a/nfert.remapv7a.txt' ;
lu_file_2 = '/Users/Shared/PLUM/input/remaps_v7b/LU.remapv7b.txt' ;
cf_file_2 = '/Users/Shared/PLUM/input/remaps_v7b/cropfracs.remapv7b.txt' ;
nf_file_2 = '/Users/Shared/PLUM/input/remaps_v7b/nfert.remapv7b.txt' ;
outDir = '/Users/Shared/PLUM/input/remaps_v7a_versus_v7b/' ;


%% Setup

if ~exist(outDir,'dir')
    mkdir(outDir) ;
end

addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;

% Get LUH2 land area (m2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = 1e6*transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
landArea_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear *_YXqd tmp


%% Import

% Land use
disp('Importing land use...')
lu1 = lpjgu_matlab_readTable(lu_file_1) ;
lu2 = lpjgu_matlab_readTable(lu_file_2) ;

% Check for identical gridlists
latlons1 = table2array(unique(lu1(:,[1 2]),'rows')) ;
latlons2 = table2array(unique(lu2(:,[1 2]),'rows')) ;
if ~isequal(latlons1,latlons2)
    error('Inputs 1 and 2 have different gridlists!')
end

% Get gridlist
lon_i = (latlons1(:,1)-0.25)*2+361 ;
lat_i = (latlons1(:,2)-0.25)*2+181 ;
I = sub2ind(size(landArea_YX),lat_i,lon_i) ;
Ncells = length(I) ;
landArea_x = landArea_YX(I) ;
landArea_11x = permute(landArea_x,[2 3 1]) ;

% Check for identical yearLists
yearList = unique(lu1.Year) ;
if ~isequal(yearList,unique(lu2.Year))
    error('Deal with non-identical yearLists!')
end
Nyears = length(yearList) ;

% Check for identical luLists
luList = lu1.Properties.VariableNames ;
luList = luList(4:end) ;
luList2 = lu2.Properties.VariableNames ;
luList2 = luList2(4:end) ;
if ~isequal(luList, luList2)
    error('Inputs have different luLists!')
end
clear luList2
Nlus = length(luList) ;

% Reorganize land use
lu1_yvx = lpjgu_table2yvx(lu1,Nyears,Ncells) ; clear lu1
lu2_yvx = lpjgu_table2yvx(lu2,Nyears,Ncells) ; clear lu2

% Crop fractions
disp('Importing crop fractions...')
cf1 = lpjgu_matlab_readTable(cf_file_1) ;
cf2 = lpjgu_matlab_readTable(cf_file_2) ;

% Check for identical cropLists
cropList = cf1.Properties.VariableNames ;
cropList = cropList(4:end) ;
cropList2 = cf2.Properties.VariableNames ;
cropList2 = cropList2(4:end) ;
if ~isequal(cropList, cropList2)
    warning('Inputs have different cropLists!')
end
clear cropList2
Ncrops = length(cropList) ;

% Reorganize crops
cf1_yvx = lpjgu_table2yvx(cf1,Nyears,Ncells) ; clear cf1
cf2_yvx = lpjgu_table2yvx(cf2,Nyears,Ncells) ; clear cf2

% Nfert
disp('Importing N fert...')
nf1 = lpjgu_matlab_readTable(nf_file_1) ;
nf2 = lpjgu_matlab_readTable(nf_file_2) ;
nf1_yvx = lpjgu_table2yvx(nf1,Nyears,Ncells) ; clear nf1
nf2_yvx = lpjgu_table2yvx(nf2,Nyears,Ncells) ; clear nf2

% Getting ancillaries
disp('Calculating ancillary arrays...')
luArea1_yvx = lu1_yvx .* repmat(landArea_11x,[Nyears Nlus 1]) ;
luArea2_yvx = lu2_yvx .* repmat(landArea_11x,[Nyears Nlus 1]) ;
cropArea1_yvx = ...
    repmat(lu1_yvx(:,strcmp(luList,'CROPLAND'),:),[1 Ncrops 1]) ...
    .* cf1_yvx ...
    .* repmat(landArea_11x,[Nyears Ncrops 1]) ;
cropArea2_yvx = ...
    repmat(lu2_yvx(:,strcmp(luList,'CROPLAND'),:),[1 Ncrops 1]) ...
    .* cf2_yvx ...
    .* repmat(landArea_11x,[Nyears Ncrops 1]) ;
nfertTot1_yvx = nf1_yvx .* cropArea1_yvx ;
nfertTot2_yvx = nf2_yvx .* cropArea2_yvx ;

disp('Done importing.')



%% Plot N fert

figure ;
plot(yearList,1e-9*[sum(sum(nfertTot1_yvx,3),2) sum(sum(nfertTot2_yvx,3),2)])
ylabel('TgN')
legend('Old','New')






