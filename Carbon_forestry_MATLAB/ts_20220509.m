%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check discontinuity at 50°N %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hook://email/AM8PR05MB8034FE54542443718939642A8CC39%40AM8PR05MB8034.eurprd05.prod.outlook.com

topDir = '/Users/Shared/landsymm_forestry/1970past/1850-2014/output-2022-05-02-083616' ;
% topDir = '/Users/Shared/landsymm_forestry/1970past/1850-2014/output-2022-05-11-000026' ;

newLims_lon = [-150 -60] ;
newLims_lat = [48 52] ;


%% Import

thisFile = 'cmass_natural.out' ;
data = lpjgu_matlab_read2geoArray(sprintf('%s/%s', topDir, thisFile)) ;
disp('Done importing')


%% Make map

thisVar = 'Total' ;
thisYear = 2014 ;

thisMap_YX = lpjgu_vector2map(data.garr_xvy(:,strcmp(data.varNames, thisVar),data.yearList==thisYear), ...
    [360 720], data.list2map) ;

lons_map = repmat(-179.75:0.5:179.75, [360 1]) ;
lats_map = repmat(transpose(-89.75:0.5:89.75), [1 720]) ;
ok_lon = lons_map > newLims_lon(1) & lons_map < newLims_lon(2) ;
ok_lat = lats_map > newLims_lat(1) & lats_map < newLims_lat(2) ;
ok_map = ok_lon & ok_lat ;
thisMap_YX(~ok_map) = NaN ;

shademap(thisMap_YX) ;
title(sprintf('%s %d', thisVar, thisYear))

set(gca, 'XLim', [98.8919  248.1563], ...
         'YLim', [270.8888  289.7167])
set(gcf, 'Position', [1000        1136         560         201], ...
         'Color', 'w')
axis off


%% Make new gridlist just for N America

gl_in = sprintf('%s/gridlist_ggcmi_v1.1.gapfilled.lpjg.remapv12_g2p.soilmap_center_interpolated.txt', topDir) ;
gridlist_orig = lpjgu_matlab_readTable(gl_in) ;
gridlist = gridlist_orig ;
ok_lon = gridlist.Lon > newLims_lon(1) & gridlist.Lon < newLims_lon(2) ;
ok_lat = gridlist.Lat > newLims_lat(1) & gridlist.Lat < newLims_lat(2) ;
ok = ok_lon & ok_lat ;
gridlist = gridlist(ok,:) ;

gl_out = sprintf('%s/gridlist_NAm.txt', topDir) ;
save_file(gridlist.Lon, gridlist.Lat, gl_out)



%% FUNCTIONS

function save_file(lons, lats, file_out)

fid = fopen(file_out,'w') ;
fprintf(fid,'%0.2f %0.2f\n',[lons lats]') ;
fclose(fid) ;

end