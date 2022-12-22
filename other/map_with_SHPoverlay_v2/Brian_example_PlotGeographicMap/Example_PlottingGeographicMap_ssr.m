function [] = Example_PlottingGeographicMap();
%
%  Uses a matlab central function (again, by Kelly Kearney, who i think i
%  played softball with once or twice) cptcmap.m at
%  http://www.mathworks.com/matlabcentral/fileexchange/28943-color-palette-tables---cpt--for-matlab 
%  or https://github.com/kakearney/cptcmap-pkg
%
if nargin==0;
    output.script = 'Example_PlottingGeographicMap.m';
    output.processingDate = date;
    [output] = GetDataSets(output);
else
    output = input;
end;  % if nargin loop
disp('Loading shapefiles, which causes my win7 machine to take about 40-60 seconds to finish.  not sure why.');
GetShapefilesForPlots;
doPlotMonthlyMapsOfFire = 1;
PlotMonthlyMapsOfFire(output,doPlotMonthlyMapsOfFire);
end % main

function [output] = GetDataSets(input);
%
%
%
output = input;
output.dataSetsFire = {'MYDIS'};
% cycle over fireDataSets
for ii=1:length(output.dataSetsFire);
    dataSetFire = output.dataSetsFire{ii};
    fileNameFire = cat(2,'firecounts.',dataSetFire,'.0.5x0.5.2011-2011.nc');
    fieldNameFire = 'CloudCorrFirePix';
    % specify output into a single structure field
    dataFire.fileName{ii} = fileNameFire;
    dataFire.fieldName{ii} = fieldNameFire;
    % get the lat lon
    dataFire.lat = double(ncread(cat(2,fileNameFire),'lat'));
    dataFire.lon = double(ncread(cat(2,fileNameFire),'lon'));
    % get the fire count data
    dataFire.(dataSetFire) = double(permute(ncread(cat(2,fileNameFire),fieldNameFire),[3,2,1]));
end;  % for ii loop
output.dataFire = dataFire;
end % function

function [] = PlotMonthlyMapsOfFire(input,doThis);
%
%
%
output = input;
if doThis;
    dataFire = output.dataFire;
    % this is the colormap as a *.cpt file in the cptfiles directory per
    % the cptcmap.m script.  i have a bunch in the cptfiles directory that
    % i fiddle around with.  tmeperature is a useful general colorscale.
    % let me know if you find something better!
    colormapCPT = 'temperature';
    latlonlimInput = [];
    totalRangeRelative = [0:0.1:1];
    doPlotGeoMap = 1;
    plotMonthsStart = 1;
    plotMonthsStop = 1;
    % cycle over months
    for kk=plotMonthsStart:plotMonthsStop;
        % delete old figures every month
        close all
        % cycle over fireDataSets
        for ii=1:length(output.dataSetsFire);
            dataSetFire = output.dataSetsFire{ii};
            plotData = permute(dataFire.(dataSetFire)(kk,:,:),[2,3,1]);
            plotData(plotData<0) = NaN;
            % normalize the monthly output to display the maximum fire or lightning
            % month as 1 and other months scaled
            normValueFire = permute(max(output.dataFire.(dataSetFire),[],1),[2,3,1]);
            plotData = plotData./normValueFire;
            lat = dataFire.lat;
            lon = dataFire.lon;
            if kk<10;
                savePlotString = cat(2,'month','0',num2str(kk),'-fire-',dataSetFire);
            else
                savePlotString = cat(2,'month',num2str(kk),'-fire-',dataSetFire);
            end;  % if loop
            PlotGeoMap(lat,lon,plotData,totalRangeRelative,latlonlimInput,savePlotString,colormapCPT,doPlotGeoMap);
        end;  % for ii loop
    end;  % for kk loop
end;  % if doThis loop
end % function

function [] = PlotGeoMap(lat,lon,data,caxisValues,latlonlim,figureSuffix,colormapCPT,doThis);
%
%
%
global states land countries
if doThis;
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
% % %     geoshow(states,'FaceColor',[0.9,0.9,0.9],'EdgeColor',[0.1,0.1,0.1],'linewidth',0.1);
% % %     geoshow(land,'FaceColor',[0.9,0.9,0.9],'linewidth',0.1);
    % geoshow(lakes,'FaceColor',[0.8,0.8,1]);
    geoshow(countries,'FaceColor',[0.9,0.9,0.9],'EdgeColor',[0.1,0.1,0.1],'linewidth',0.1);
% % %     geoshow(countries,'FaceColor',[0 0 0],'EdgeColor',[0.1,0.1,0.1],'linewidth',0.1);
    pcolorm(lat,lon,data);
    if ~isempty(caxisValues);
        caxis([min(caxisValues),max(caxisValues)]);
        % custom color bar script
        cptcmap(colormapCPT,'ncol',length(caxisValues)-1);
    end;  % if loop
    Hc = colorbar('YTick',caxisValues,'YTickLabel',caxisValues);
    set(Hc,'YTickMode','manual','fontsize',plotFont-1);
    Ht = title(figureSuffix);
    set(Ht,'fontsize',plotFont);
    tightmap
    if savefig;
        ChangePlotDimensions(5,3);
% % %         if isunix;
% % %             error('uh oh!');
% % %         else
% % %             figpath = [];
% % %         end;  % if isunix loop
        figpath = [];
        figpre = cat(2,'maps-',figureSuffix,'.');
        figname = cat(2,figpath,figpre,'png');
        disp(cat(2,'** Overwriting fig to ',figname));
        print('-dpng','-r600',figname);
    end;  % if savefig loop
end;  % if doThis loop
end % function

function [] = GetShapefilesForPlots;
%
%  many shapefiles are built in to standard matlab paths, but
%  world borders path needs to be customized or built into custom default
%  paths (if applicable) or in the current directory
%
global states land countries
states = shaperead('usastatehi','UseGeoCoords', true);
land = shaperead('landareas','UseGeoCoords', true);
% lakes = shaperead('worldlakes', 'UseGeoCoords', true);
countries = shaperead('world_borders','UseGeoCoords', true);
end % function

function [] = ChangePlotDimensions(width,height,orientation);
%
%  Change dimensions of a plot and (optionally) the paper/layout
%  orientation.  Run as:
%
%  Universal_ChangePlotDimensions(width,height,orientation);
%  Universal_ChangePlotDimensions(6,9);
%  Universal_ChangePlotDimensions(6,9,'portrait');
%  Universal_ChangePlotDimensions(9,6,'landscape');
%
if nargin==2;
    orientation = 'portrait';
end;  % if nargin loop
set(gcf,'PaperPositionMode','manual');
set(gcf, 'PaperUnits', 'inches');
switch orientation
    case 'portrait'
        figPos = [0.5,0.5];
    case 'landscape'
        figPos = [0.5,0.5];
    otherwise
        error('!! Unknown paper orientation, dude');
end;  % switch loop
% PaperPosition is left, bottom, width, height but i think this flips in
% landscape mode
set(gcf,'PaperPosition',[figPos(1),figPos(2),width,height]);
set(gcf,'PaperOrientation',orientation);
end  % function
