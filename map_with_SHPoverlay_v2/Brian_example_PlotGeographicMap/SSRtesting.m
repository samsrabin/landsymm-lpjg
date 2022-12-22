%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing Brian's mapping code %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup

addpath(genpath('/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap'))

%% Do it
Example_PlottingGeographicMap()

%% Test

if ~exist('countries','var')
    states = shaperead('usastatehi','UseGeoCoords', true);
    land = shaperead('landareas','UseGeoCoords', true);
    countries = shaperead('world_borders','UseGeoCoords', true);
end

latlonlim = [] ;
savefig = 1;
plotFont = 7;
figure;
clf reset;
set(gcf,'Color','white');
if isempty(latlonlim);
    % full bounds
    lonlim = [-180,180];
    latlim = [-60,80];
    % % %         lonlim = [-130,-65];
    % % %         latlim = [24,50];
else
    latlim = latlonlim(1,:);
    lonlim = latlonlim(2,:);
end;  % if loop
% usamap(latlim,lonlim);
worldmap(latlim,lonlim);
mlabel('off'); % ,'fontsize',5)
plabel('off'); % ,'fontsize',5)
framem off
geoshow(states,'FaceColor',[0.9,0.9,0.9],'EdgeColor',[0.1,0.1,0.1],'linewidth',0.1);
geoshow(land,'FaceColor',[0.9,0.9,0.9],'linewidth',0.1);
% geoshow(lakes,'FaceColor',[0.8,0.8,1]);
geoshow(countries,'FaceColor',[0.9,0.9,0.9],'EdgeColor',[0.1,0.1,0.1],'linewidth',0.1);
% pcolorm(lat,lon,data);