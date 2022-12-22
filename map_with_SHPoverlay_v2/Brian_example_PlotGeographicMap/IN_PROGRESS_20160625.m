R = georasterref('RasterSize',[90 144],...
     'Latlim',[-90 90], 'Lonlim', [-180 180]);

% tmp = mean(BF_obs_totl_YXy(15:90,:,:),3) ;
tmp = mean(BF_obs_totl_YXy,3) ;

figure ;
worldmap(tmp,R);
geoshow(tmp,R,'DisplayType','surface')


%%

[X,Y] = meshgrid(-178.75:2.5:178.75,-89:2:89) ;

%%

figure ;
worldmap(tmp,R);
geoshow(Y,X,tmp,'DisplayType','surface')

%%

fontSize = 14 ;

map_data = mean(BF_obs_totl_YXy,3) ;

[~,indices] = histc(map_data,BF_binEdges);
indices(map_data>max(BF_binEdges)) = max(indices(:)) ;
indices(isnan(map_data)) = NaN ;
indices(isinf(map_data)) = NaN ;

if ~exist('land','var')
    land = shaperead('landareas','UseGeoCoords', true);
end
if ~exist('countries','var')
%     countries = shaperead('world_borders','UseGeoCoords', true);
    countries = shaperead('/Users/sam/Geodata/General/countriesEtc_NaturalEarth/ne_110m_admin_0_countries_lakes/ne_110m_admin_0_countries_lakes.shp','UseGeoCoords', true);
end

% pad = false ;
pad = true ;

lonlim = [-180,180];
latlim = [-60,80];
% latlim = [-90,90];
if ~pad
    lat = -89:2:89 ;
    lon = -178.75:2.5:178.75 ;
else
%     lat = -89:2:89 ;
    lat = -89:2:89+2 ;
    lon = -178.75:2.5:178.75+2.5 ;
    indices = [indices nan(90,1)] ;
%     indices = [nan(1,145) ; indices] ;
    indices = [indices ; nan(1,145)] ;
end

figure ;

worldmap(latlim,lonlim);
mlabel('off'); % ,'fontsize',5)
plabel('off'); % ,'fontsize',5)
framem off
% geoshow(land,'FaceColor',[0 0 0],'EdgeColor',[0.1,0.1,0.1],'linewidth',0.1);
% pcolorm(lat,lon,map_data) ;
% pcolorm(lat,lon,indices-1) ;
pcolorm(lat-1,lon-1.25,indices-1) ;
% pcolorm(lat,lon-1.25,indices-1) ;
% pcolorm(lat+1,lon+1.25,indices-1) ;
% pcolorm(latlim,lonlim,indices-1) ;

% [cmap, lims, ticks, bfncol, ctable] = ...
% cptcmap('/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr.BF_colormap.cpt',...
%         'ncol',length(BF_binEdges)-1);%,...
        %'mapping','direct') ;
cptcmap('/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr.BF_colormap.indices.cpt',...
        'ncol',length(BF_binEdges)-1);%,...
        %'mapping','direct') ;
        
% geoshow(land,'FaceColor','none','EdgeColor',[0.1,0.1,0.1],'linewidth',0.1);
% geoshow(countries,'FaceColor','none','EdgeColor',[1 1 1],'linewidth',3);
geoshow(countries,'FaceColor','none','EdgeColor',[0.1,0.1,0.1],'linewidth',0.1);

caxis([0 length(BF_binEdges)-1])
h = colorbar;
% set(h,'FontSize',fontSize,'Ticks',0:(length(BF_binEdges)-1),'TickLabels',BF_labels)
set(h,'FontSize',fontSize,'Ticks',0:(length(BF_binEdges)-1),'TickLabels',BF_labels)

%%

lonlim = [-180,180] ;
latlim = [-60,80] ;
colormap = 'ssr_BFbins' ;
fontSize = 14 ;
% edgecolor = [0.1 0.1 0.1] ;
edgecolor = 0.4*ones(3,1) ;
% linewidth = 0.1 ;
linewidth = 0.25 ;

figure ;
map_with_SHPoverlay(mean(BF_obs_totl_YXy,3),countries,lonlim,...
                    latlim,colormap,fontSize, ...
                    edgecolor,linewidth)
                
set(gcf,'Position',[440 292 948 506])
%%
export_fig '~/Desktop/test.png' -r300
disp('Done.')
                


%%

Example_PlottingGeographicMap_ssr()


%%
