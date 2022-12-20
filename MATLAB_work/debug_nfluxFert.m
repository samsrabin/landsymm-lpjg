%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Diagnosing nflux-fert problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisLon = 118.25 ;
thisLat = 35.75 ;

% file_lu = '/Users/Shared/PLUM/input/remaps_v6p2/LU.remapv6p2.txt' ;
% file_cf = '/Users/Shared/PLUM/input/remaps_v6p2/cropfracs.remapv6p2.txt' ;
% file_nf = '/Users/Shared/PLUM/input/remaps_v6p2/nfert.remapv6p2.txt' ;
file_lu = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_remap6/test2/test_lu.txt' ;
file_cf = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_remap6/test2/test_cf.txt' ;
file_nf = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_remap6/test2/test_nf.txt' ;


%% Setup

lons = -179.75:0.5:179.75 ;
lats = (-89.75:0.5:89.75)' ;
lons_map = repmat(lons,[360 1]) ;
lats_map = repmat(lats,[1 720]) ;
% i = (thisLat - -89.75)/0.5 + 1 ;
% j = (thisLon - -179.75)/0.5 + 1 ;

% Import land area (km2 to m2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
tmp = gcel_area_YXqd(:,1:2:1440) + gcel_area_YXqd(:,2:2:1440) ;
gcel_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
% Convert to m2
land_area_YX = land_area_YX*1e6 ;
gcel_area_YX = gcel_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd


%% Import

disp('Importing...')

lu = lpjgu_matlab_readTable(file_lu) ;
lu = lu(lu.Lon==thisLon & lu.Lat==thisLat,:) ;
if any(strcmp(lu.Properties.VariableNames,'Year'))
    lu_yearList = lu.Year ;
end

cf = lpjgu_matlab_readTable(file_cf) ;
cf = cf(cf.Lon==thisLon & cf.Lat==thisLat,:) ;
cropList = cf.Properties.VariableNames ;
cropList = cropList(~contains(cropList,{'Lon','Lat','Year'})) ;
cf_yc = table2array(cf(:,~contains(cf.Properties.VariableNames,{'Lon','Lat','Year'}))) ;
if any(strcmp(lu.Properties.VariableNames,'Year'))
    cf_yearList = cf.Year ;
end
clear cf

nf = lpjgu_matlab_readTable(file_nf) ;
nf = nf(nf.Lon==thisLon & nf.Lat==thisLat,:) ;
[~,~,IB] = intersect(cropList,nf.Properties.VariableNames,'stable') ;
if ~isequal(cropList,nf.Properties.VariableNames(IB))
    error('Incompatible crop lists?')
end
nf_yc = table2array(nf(:,IB)) ;
if any(strcmp(lu.Properties.VariableNames,'Year'))
    nf_yearList = nf.Year ;
end
clear nf

% Align yearLists
if any(strcmp(lu.Properties.VariableNames,'Year'))
    [yearList,~,~] = intersect(lu_yearList, cf_yearList) ;
    [yearList,~,~] = intersect(yearList, nf_yearList) ;
    y1 = min(yearList) ;
    yN = max(yearList) ;
    lu = lu(lu.Year>=y1 & lu.Year<=yN, :) ;
    cf = cf_yc(cf_yearList>=y1 & cf_yearList<=yN, :) ;
    nf = nf_yc(nf_yearList>=y1 & nf_yearList<=yN, :) ;
end

disp('Done.')


%% Calculate

thisYear = 2010 ;

thisLandArea = land_area_YX(lons_map==thisLon & lats_map==thisLat) ;
expected = lu.CROPLAND .* sum(cf_yc .* nf_yc, 2) ;
% Convert from kgN/m2 to kgN/ha
expected = expected .* 1e4 ;
if exist('yearList','var')
    expected = expected(yearList==thisYear) ;
end

fprintf('Cell %0.2f %0.2f in %d should have had nflux(fert) = %0.1e\n', ...
    thisLon, thisLat, thisYear, -expected) ;







