load ('~/Desktop/mycells.mat')
figure ;
worldmap(mycells.map,mycells.R);
geoshow(mycells.map,mycells.R,'DisplayType','surface','ZData',zeros(size(mycells.map)),'CData',mycells.map)

%%

figure ;
tmp = mycells.map ;
% tmp(isnan(tmp)) = 1e9 ;
tmp(40:41,70:75) = 1e9 ;
% figure ;
worldmap(tmp,mycells.R);
h = geoshow(tmp,mycells.R,'DisplayType','surface','ZData',zeros(size(tmp)),'CData',tmp) ;
set(h,'AlphaData',isnan(tmp))


%%

[i,j] = find(~isnan(mycells.map)) ;
lats = -89 + (i-1)*2 ;
lons = -178.75 + (j-1)*2.5 ;
worldmap(tmp,mycells.R);
geoshow(lats,lons,'DisplayType','point')
