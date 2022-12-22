function [h,hcb] = map_with_SHPoverlay_v2(map_data,shapedata,varargin)
% Optional input arguments:
%   lonlim:     2-element vector, default [-180 180]
%   latlim:     2-element vector, default [-90 90]
%   thisColormap:   string
%   fontSize:   number
%   edgeColor:  3-element vector, default [0.1,0.1,0.1]
%   lineWidth:  number, default 0.1
%   cbarOrient: string, default 'eastoutside'
%   bground:    2-d array
%   do_log10x_ticks: 1x1 logical

% Set up & parse input arguments
p = inputParser ;
isnumarray = @(x) isnumeric(x) && ismatrix(x) && min(size(x))>1 ;
is2elemvec = @(x) (length(x)==2 && ismatrix(x) && min(size(x))==1) || isempty(x) ;
is3elemvec01 = @(x) length(x)==3 && ismatrix(x) && min(size(x))==1 && min(x)>=0 && max(x)<=1 ;
is1x1num = @(x) isequal(size(x),[1 1]) ;
addRequired(p,'map_data',isnumarray) ;
addRequired(p,'shapedata') ;
addOptional(p,'lonlim',[-180 180],is2elemvec) ;
addOptional(p,'latlim',[-90 90],is2elemvec) ;
addOptional(p,'thisColormap',[],@isstr) ;
addOptional(p,'fontSize',[],is1x1num) ;
addOptional(p,'edgeColor',[0.1,0.1,0.1],is3elemvec01) ;
addOptional(p,'lineWidth',0.1,is1x1num) ;
addOptional(p,'cbarOrient','eastoutside',@isstr) ;
addOptional(p,'bground',[]) ;
addOptional(p,'caxis_lims',[],is2elemvec) ;
addOptional(p,'do_log10x_ticks',false,is1x1num) ;
addOptional(p,'do_log10_ticks',false,is1x1num) ;
addOptional(p,'flip',false,is1x1num) ;
addOptional(p,'units_map','',@isstr) ;
addOptional(p,'shapedata2',[]) ;
addOptional(p,'shiftup',0,is1x1num) ;
% Options for SSE_diffs
addOptional(p,'clim_int',1,@isint) ;
addOptional(p,'only_optd',true,is1x1num) ;
addOptional(p,'cbar_minmax',[-Inf Inf],is2elemvec) ;

parse(p,map_data,shapedata,varargin{:});
lonlim = p.Results.lonlim ;
latlim = p.Results.latlim ;
thisColormap = p.Results.thisColormap ;
fontSize = p.Results.fontSize ;
edgeColor = p.Results.edgeColor ;
lineWidth = p.Results.lineWidth ;
cbarOrient = p.Results.cbarOrient ;
bground = p.Results.bground ;
caxis_lims = p.Results.caxis_lims ;
do_log10x_ticks = p.Results.do_log10x_ticks ;
do_log10_ticks = p.Results.do_log10_ticks ;
clim_int = p.Results.clim_int ;
only_optd = p.Results.only_optd ;
cbar_minmax = p.Results.cbar_minmax ;
flip = p.Results.flip ;
units_map = p.Results.units_map ;
shapedata2 = p.Results.shapedata2 ;
shiftup = p.Results.shiftup ;


%%
% Calculate cell size: Assumes that this is a global geographic map
cellsize_lat = 180 / size(map_data,1) ;
cellsize_lon = 360 / size(map_data,2) ;
lat = (-90 + 0.5*cellsize_lat):cellsize_lat:(90 + 0.5*cellsize_lat) ;
lon = (-180 + 0.5*cellsize_lon):cellsize_lon:(180 + 0.5*cellsize_lon) ;

worldmap(latlim,lonlim) ;
mlabel('off'); % ,'fontsize',5)
plabel('off'); % ,'fontsize',5)
framem off
gridm off

if strcmp(thisColormap,'ssr_BFbins') || strcmp(thisColormap,'ssr_CemitBins') ...
   || strcmp(thisColormap,'ssr_BFbins_v2') || strcmp(thisColormap,'ssr_CemitBins_v2') ...
   || strcmp(thisColormap,'ssr_BFbins_v3') || strcmp(thisColormap,'ssr_CemitBins_v3')
    if strcmp(thisColormap,'ssr_BFbins') || strcmp(thisColormap,'ssr_BFbins_v2') || strcmp(thisColormap,'ssr_BFbins_v3')
        binEdges = [0 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1] ;
        cAxisLabels = textscan(num2str(binEdges*100),'%s') ;
        cAxisLabels = cAxisLabels{1} ;
        cAxisLabels{length(cAxisLabels)} = [cAxisLabels{length(cAxisLabels)} '%'] ;
    elseif strcmp(thisColormap,'ssr_CemitBins') || strcmp(thisColormap,'ssr_CemitBins_v2') || strcmp(thisColormap,'ssr_CemitBins_v3')
%         binEdges = [0 10^10 5*10^10 10^11 5*10^11 10^12 5*10^12 10^13 5*10^13 10^14] ;
%         binEdges = [0 5*10^10 10^11 5*10^11 10^12 5*10^12 10^13 5*10^13 10^14] ;
        binEdges = [0 10^10 5*10^10 10^11 5*10^11 10^12 5*10^12 10^13 5*10^13] ;
        cAxisLabels = {'0','0.01','0.05','0.1','0.5','1','5','10','50'} ;
    end
    
    [~,indices] = histc(map_data,binEdges);
    indices(map_data>max(binEdges)) = max(indices(:)) ;
    indices(isnan(map_data)) = NaN ;
    indices(isinf(map_data)) = NaN ;
%     indices = [indices nan(90,1)] ;
%     indices = [indices ; nan(1,145)] ;
    indices = [indices nan(size(indices,1),1)] ;
    indices = [indices ; nan(1,size(indices,2))] ;
    
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,indices-1) ;
    
    if strcmp(thisColormap,'ssr_BFbins') || strcmp(thisColormap,'ssr_CemitBins')
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr.BF_colormap.indices.cpt' ;
    elseif strcmp(thisColormap,'ssr_BFbins_v2') || strcmp(thisColormap,'ssr_CemitBins_v2')
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr.BF_colormap.indices_v2.cpt' ;
    elseif strcmp(thisColormap,'ssr_BFbins_v3') || strcmp(thisColormap,'ssr_CemitBins_v3')
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr.BF_colormap.indices_v3.cpt' ;
    else
        error(['No cpt_file specified for thisColormap = ' thisColormap '.'])
    end
    cptcmap(cpt_file,'ncol',length(binEdges)-1);
    caxis([0 length(binEdges)-1])
%     hcb = colorbar(cbarOrient) ;
%     if ~isempty(fontSize)
%         set(hcb,'FontSize',fontSize)
%     end
%     set(hcb,'Ticks',0:(length(binEdges)-1),'TickLabels',cAxisLabels)
    hcb = do_colorbar(cbarOrient,fontSize,0:(length(binEdges)-1),cAxisLabels) ;


elseif strcmp(thisColormap,'rdbu_ssr')
    map_data = [map_data nan(size(map_data,1),1)] ;
    map_data = [map_data ; nan(1,size(map_data,2))] ;
    
    if ~isempty(bground)
        bground = [bground nan(size(bground,1),1)] ;
        bground = [bground ; nan(1,size(bground,2))] ;
        
        % Get initial colorbar
%         pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
        htmp = gca ;
        pcolor(map_data); shading flat; axis equal tight
%         cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/rdbu_ssr.cpt' ;
%         cptcmap(cpt_file,'ncol',65);
%         old_cmap = colormap(gca,brewermap(64,'rdbu_ssr')) ;
        old_cmap = brewermap(64,'rdbu_ssr') ;
        colormap(htmp,old_cmap) ;
        if isempty(caxis_lims)
            caxis_lims = [-max(abs(caxis)) max(abs(caxis))] ;
        end
        caxis(caxis_lims)
        cb = colorbar(cbarOrient) ;
        cb_lims = [cb.Limits(1) cb.Limits(2)] ;
        negVal = (cb.Limits(2) - cb.Limits(1)) / size(old_cmap,1) ;
        map_data_2 = map_data ;
        map_data_2(bground==1 & isnan(map_data)) = cb_lims(1)-negVal ;
        
        colorbar('off')
        
        pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data_2) ;
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/rdbu_ssr_plusBGgray_dark.cpt' ;
        cptcmap(cpt_file,htmp,'ncol',65,'flip',flip);
        hcb = colorbar(cbarOrient) ;
        caxis([cb_lims(1)-negVal cb_lims(2)])
        set(hcb,'YLim',cb_lims) ;
        if ~isempty(fontSize)
            set(hcb,'FontSize',fontSize)
        end
        if do_log10x_ticks
            gca_pos = get(gca,'Position') ;
            hcb_pos = get(hcb,'Position') ;
            ticks = get(hcb,'Ticks') ;
            ticklabels = get(hcb,'TickLabels') ;
            for tt = 1:length(ticks)
                ticklabels{tt} = [num2str(10^ticks(tt)) 'x'] ;
%                 ticklabels{tt} = ['10^{' num2str(ticks(tt)) '}x'] ;
            end
            set(hcb,'TickLabels',ticklabels)
            set(gca,'Position',gca_pos)
            set(hcb,'Position',hcb_pos)
        end
        if ~isempty(units_map)
            hcb.Label.String = units_map ;
        end
        
        colorbar('off')
        if shiftup > 0
            set(gca,'Position',get(gca,'Position') + [0 shiftup 0 0])
        end
        clear hcb
        hcb = cb_lims ;
        h = gca ;
                
    else
        hpcolorm = pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/rdbu_ssr.cpt' ;
        cptcmap(cpt_file,'ncol',64);
        if isempty(caxis_lims)
            caxis_lims = [-max(abs(caxis)) max(abs(caxis))] ;
        end
        if do_log10x_ticks
            caxis_lims_new = caxis_lims ;
            caxis_lims_new(1) = floor(caxis_lims_new(1)) ;
            caxis_lims_new(2) = ceil(caxis_lims_new(2)) ;
            caxis_lims = caxis_lims_new ;
        end
        caxis(caxis_lims)
        %     hcb = colorbar(cbarOrient) ;
        %     if ~isempty(fontSize)
        %         set(hcb,'FontSize',fontSize)
        %     end
        hcb = do_colorbar(cbarOrient,fontSize) ;
    end
    
    
    
    
    
elseif strcmp(thisColormap,'parula')
    map_data = [map_data nan(size(map_data,1),1)] ;
    map_data = [map_data ; nan(1,size(map_data,2))] ;
    
    if ~isempty(bground)
        bground = [bground nan(size(bground,1),1)] ;
        bground = [bground ; nan(1,size(bground,2))] ;
        
        % Get initial colorbar
%         pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
        pcolor(map_data); shading flat; axis equal tight
%         cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/rdbu_ssr.cpt' ;
%         cptcmap(cpt_file,'ncol',65);
%         old_cmap = colormap(gca,brewermap(64,'rdbu_ssr')) ;
        old_cmap = brewermap(64,'Reds') ;
        colormap(old_cmap) ;
%         if isempty(caxis_lims)
%             caxis_lims = [-max(abs(caxis)) max(abs(caxis))] ;
%         end
%         caxis(caxis_lims)
        cb = colorbar(cbarOrient) ;
        cb_lims = [cb.Limits(1) cb.Limits(2)] ;
        negVal = (cb.Limits(2) - cb.Limits(1)) / size(old_cmap,1) ;
        map_data_2 = map_data ;
        map_data_2(bground==1 & isnan(map_data)) = cb_lims(1)-negVal ;
        
        pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data_2) ;
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr_parula_plusBGgray_light.cpt' ;
        cptcmap(cpt_file,'ncol',65); % good thru here
        hcb = colorbar(cbarOrient) ; % good thru here
%         caxis([cb_lims(1)-negVal cb_lims(2)])
        set(hcb,'YLim',cb_lims) ; % good thru here
        if ~isempty(fontSize)
            set(hcb,'FontSize',fontSize)
        end
        if do_log10x_ticks
            gca_pos = get(gca,'Position') ;
            hcb_pos = get(hcb,'Position') ;
            ticks = get(hcb,'Ticks') ;
            ticklabels = get(hcb,'TickLabels') ;
            for tt = 1:length(ticks)
                ticklabels{tt} = [num2str(10^ticks(tt)) 'x'] ;
%                 ticklabels{tt} = ['10^{' num2str(ticks(tt)) '}x'] ;
            end
            set(hcb,'TickLabels',ticklabels)
            set(gca,'Position',gca_pos)
            set(hcb,'Position',hcb_pos)
        end
    else
        pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/rdbu_ssr.cpt' ;
        cptcmap(cpt_file,'ncol',64);
%         if isempty(caxis_lims)
%             caxis_lims = [-max(abs(caxis)) max(abs(caxis))] ;
%         end
        if do_log10x_ticks
            caxis_lims_new = caxis_lims ;
            caxis_lims_new(1) = floor(caxis_lims_new(1)) ;
            caxis_lims_new(2) = ceil(caxis_lims_new(2)) ;
            caxis_lims = caxis_lims_new ;
        end
%         caxis(caxis_lims)
        %     hcb = colorbar(cbarOrient) ;
        %     if ~isempty(fontSize)
        %         set(hcb,'FontSize',fontSize)
        %     end
        hcb = do_colorbar(cbarOrient,fontSize) ;
    end
    
    if ~isempty(units_map)
        hcb.Label.String = units_map ;
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
elseif strcmp(thisColormap,'rdbu_ssr_neg1to1')
    map_data = [map_data nan(90,1)] ;
    map_data = [map_data ; nan(1,145)] ;
    
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
    cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/rdbu_ssr.cpt' ;
    cptcmap(cpt_file,'ncol',64);
    caxis([-1 1])
%     hcb = colorbar(cbarOrient) ;
%     if ~isempty(fontSize)
%         set(hcb,'FontSize',fontSize)
%     end
    hcb = do_colorbar(cbarOrient,fontSize) ;
elseif strcmp(thisColormap,'reds_ssr_withgray_10_0to1')
    map_data = [map_data nan(90,1)] ;
    map_data = [map_data ; nan(1,145)] ;
    
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
    cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/reds_ssr_withgray_10.cpt' ;
    cptcmap(cpt_file,'ncol',10);
    caxis([0 1])
%     caxis([0 max(caxis)])
%     hcb = colorbar(cbarOrient) ;
%     if ~isempty(fontSize)
%         set(hcb,'FontSize',fontSize)
%     end
    hcb = do_colorbar(cbarOrient,fontSize) ;
elseif strcmp(thisColormap,'blues_ssr_withgray_10_0to1')
    map_data = [map_data nan(90,1)] ;
    map_data = [map_data ; nan(1,145)] ;
    
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
    cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/blues_ssr_withgray_10.cpt' ;
    cptcmap(cpt_file,'ncol',10);
    caxis([-1 0])
%     caxis([min(caxis) 0])
%     hcb = colorbar(cbarOrient) ;
%     if ~isempty(fontSize)
%         set(hcb,'FontSize',fontSize)
%     end
    hcb = do_colorbar(cbarOrient,fontSize) ;
elseif strcmp(thisColormap,'reds_ssr_withgray_10_0toMax')
    map_data = [map_data nan(90,1)] ;
    map_data = [map_data ; nan(1,145)] ;
    
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
    cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/reds_ssr_withgray_10.cpt' ;
    cptcmap(cpt_file,'ncol',10);
    caxis([0 max(caxis)])
%     hcb = colorbar(cbarOrient) ;
%     if ~isempty(fontSize)
%         set(hcb,'FontSize',fontSize)
%     end
    hcb = do_colorbar(cbarOrient,fontSize) ;
elseif strcmp(thisColormap,'blues_ssr_withgray_10_minTo0')
    map_data = [map_data nan(90,1)] ;
    map_data = [map_data ; nan(1,145)] ;
    
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
    cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/blues_ssr_withgray_10.cpt' ;
    cptcmap(cpt_file,'ncol',10);
    caxis([min(caxis) 0])
%     hcb = colorbar(cbarOrient) ;
%     if ~isempty(fontSize)
%         set(hcb,'FontSize',fontSize)
%     end
    hcb = do_colorbar(cbarOrient,fontSize) ;
elseif strcmp(thisColormap,'ssrLC_orange_5') || strcmp(thisColormap,'ssrLC_green_5') || strcmp(thisColormap,'ssrLC_purple_5')
    map_data = [map_data nan(90,1)] ;
    map_data = [map_data ; nan(1,145)] ;
    bground = [bground nan(90,1)] ;
    bground = [bground ; nan(1,145)] ;
    
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,bground) ;
    cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr_allgray.cpt' ;
    cptcmap(cpt_file,'ncol',1);
    hold on
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
    hold off
    cpt_file = ['/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/' thisColormap] ;
    cptcmap(cpt_file,'ncol',5);
    caxis([0 1])
%     hcb = colorbar(cbarOrient) ;
%     if ~isempty(fontSize)
%         set(hcb,'FontSize',fontSize)
%     end
    hcb = do_colorbar(cbarOrient,fontSize) ;
elseif strcmp(thisColormap,'ssrLC_orange_5_pluslo') || strcmp(thisColormap,'ssrLC_green_5_pluslo') || strcmp(thisColormap,'ssrLC_purple_5_pluslo')
%     map_data = [map_data nan(90,1)] ;
%     map_data = [map_data ; nan(1,145)] ;
    
    binEdges = [0 eps 0.2:0.2:1] ;
    cAxisLabels = textscan(num2str(binEdges*100),'%s') ;
    cAxisLabels = cAxisLabels{1} ;
    cAxisLabels{length(cAxisLabels)} = [cAxisLabels{length(cAxisLabels)} '%'] ;
    cAxisLabels{2} = '>0' ;
    [~,indices] = histc(map_data,binEdges);
    indices(map_data>max(binEdges)) = max(indices(:)) ;
    indices(isnan(map_data)) = NaN ;
    indices(isinf(map_data)) = NaN ;
    indices = [indices nan(90,1)] ;
    indices = [indices ; nan(1,145)] ;
    
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,indices-1) ;
    
%     bground = [bground nan(90,1)] ;
%     bground = [bground ; nan(1,145)] ;
%     
%     pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,bground) ;
%     cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr_allgray.cpt' ;
%     cptcmap(cpt_file,'ncol',1);
%     hold on
%     pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
%     hold off
    cpt_file = ['/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/' thisColormap] ;
    cptcmap(cpt_file,'ncol',length(binEdges)-1);
    caxis([0 length(binEdges)-1])
%     hcb = colorbar(cbarOrient) ;
%     if ~isempty(fontSize)
%         set(hcb,'FontSize',fontSize)
%     end
%     set(hcb,'Ticks',0:(length(binEdges)-1),'TickLabels',cAxisLabels)
    hcb = do_colorbar(cbarOrient,fontSize,0:(length(binEdges)-1),cAxisLabels) ;
elseif strcmp(thisColormap,'heterogeneity_maps')
    map_data = [map_data nan(size(map_data,1),1)] ;
    map_data = [map_data ; nan(1,size(map_data,2))] ;
    if ~isempty(bground)
        bground = [bground nan(size(bground,1),1)] ;
        bground = [bground ; nan(1,size(bground,2))] ;
        pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,bground) ;
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr_allgray.cpt' ;
        cptcmap(cpt_file,'ncol',1);
        hold on
    end
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
    hold off
    cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr_parula.cpt' ;
    cptcmap(cpt_file,'ncol',64);
    if ~isempty(caxis_lims)
        caxis(caxis_lims) ;
    end
    cbarticks = log10([0.1 0.2 0.5 1 2 5 10 20 50 100]) ;
    cbarticklabels = {} ;
    for t = 1:length(cbarticks)
        cbarticklabels{t} = num2str(10^cbarticks(t)) ;
    end
%     hcb = colorbar(cbarOrient) ;
%     if ~isempty(fontSize)
%         set(hcb,'FontSize',fontSize)
%     end
%     set(hcb,'Ticks',cbarticks,'TickLabels',cbarticklabels)
    hcb = do_colorbar(cbarOrient,fontSize,cbarticks,cbarticklabels) ;
elseif strcmp(thisColormap,'months')
    map_data = [map_data nan(90,1)] ;
    map_data = [map_data ; nan(1,145)] ;
    
    cmap_tmp = [0 0 255 ;    % Blue
                0 128 64 ;   % Dark green
                255 255 0 ;  % Yellow
                255 128 0 ;  % Orange
                255 0 255 ;  % Magenta
                108 0 110 ;  % Purple
                0 0 255 ...  % Blue
                ]/255;
    cmap_tmp = brighten(cmap_tmp,0.5) ;
    months_colorscheme = interp1(1:size(cmap_tmp,1),cmap_tmp,1:((size(cmap_tmp,1)-1)/364):size(cmap_tmp,1),'pchip') ;
%     months_colorscheme = interp1(1:size(cmap_tmp,1),cmap_tmp,1:((size(cmap_tmp,1)-1)/11):size(cmap_tmp,1),'pchip') ;
    months_ticks = ([1 32 60 91 121 152 182 213 244 274 305 335]+14) ;
    months_labels = {'\begin{tabular}{c} Jan. \\ 1 \end{tabular}';
                    '\begin{tabular} Feb.\\1\end{tabular}';
                    '\begin{tabular} Mar.\\1\end{tabular}';
                    '\begin{tabular} Apr.\\1\end{tabular}';
                    '\begin{tabular} May\\1\end{tabular}';
                    '\begin{tabular} Jun.\\1\end{tabular}';
                    '\begin{tabular} Jul.\\1\end{tabular}';
                    '\begin{tabular} Aug.\\1\end{tabular}';
                    '\begin{tabular} Sep.\\1\end{tabular}';
                    '\begin{tabular} Oct.\\1\end{tabular}';
                    '\begin{tabular} Nov.\\1\end{tabular}';
                    '\begin{tabular} Dec.\\1\end{tabular}'} ;
     
    pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
    cpt_file = '/Users/sam/Documents/Dropbox/FireMIP/FireMIP_paper2_results/maps/cptfiles/ssr.months_colormap.cpt' ;
    cptcmap(cpt_file,'ncol',255);
    if ~isempty(caxis_lims)
        caxis(caxis_lims)
    else
        caxis([-max(abs(caxis)) max(abs(caxis))])
    end
%     hcb = colorbar(cbarOrient) ;
%     if ~isempty(fontSize)
%         set(hcb,'FontSize',fontSize)
%     end
% % %     hcb = do_colorbar(cbarOrient,fontSize) ;
    hcb = do_colorbar(cbarOrient,fontSize,months_ticks,months_labels) ;
elseif strcmp(thisColormap,'SSEs_optd')
    
    
    map_data = [map_data nan(90,1)] ;
    map_data = [map_data ; nan(1,145)] ;
    if ~isempty(bground)
        bground = [bground nan(90,1)] ;
        bground = [bground ; nan(1,145)] ;
        
        % Get initial colorbar
        pcolor(map_data); shading flat; axis equal tight
        old_cmap = colormap(gca,'parula') ;
        colormap(old_cmap) ;
        cb = colorbar(gca) ;
        cb_lims = [cb.Limits(1) cb.Limits(2)] ;
        negVal = (cb.Limits(2) - cb.Limits(1)) / size(old_cmap,1) ;
        map_data_2 = map_data ;
        map_data_2(bground==1 & isnan(map_data)) = cb_lims(1)-negVal ;
        delete(cb)
        
        pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data_2) ;
%         cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr_parula_plusBGgray_dark.cpt' ;
%         cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr_parula_plusBGgray_light.cpt' ;
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr_parulaBright0.5_plusBGgray_dark.cpt' ;
        cptcmap(cpt_file,'ncol',65);
        hcb = colorbar(cbarOrient) ;
        set(hcb,'YLim',cb_lims) ;
        if ~isempty(fontSize)
            set(hcb,'FontSize',fontSize)
        end
        if do_log10_ticks
            gca_pos = get(gca,'Position') ;
            hcb_pos = get(hcb,'Position') ;
            ticks = get(hcb,'Ticks') ;
            ticklabels = get(hcb,'TickLabels') ;
            for tt = 1:length(ticks)
%                 ticklabels{tt} = num2str(10^ticks(tt)) ;
                ticklabels{tt} = ['10^{' num2str(ticks(tt)) '}'] ;
            end
            set(hcb,'TickLabels',ticklabels)
            set(gca,'Position',gca_pos)
            set(hcb,'Position',hcb_pos)
        end
    else
        pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/ssr_parula.cpt' ;
        cptcmap(cpt_file,'ncol',64);
        if ~isempty(caxis_lims)
            caxis(caxis_lims) ;
        end
        hcb = colorbar(cbarOrient) ;
        if ~isempty(fontSize)
            set(hcb,'FontSize',fontSize)
        end
    end
elseif strcmp(thisColormap,'SSE_diffs')
    map_data = [map_data nan(90,1)] ;
    map_data = [map_data ; nan(1,145)] ;
    if ~isempty(bground)
        bground = [bground nan(90,1)] ;
        bground = [bground ; nan(1,145)] ;
        
        pcolorm(lat-0.5*cellsize_lat,lon-0.5*cellsize_lon,map_data) ;
        cpt_file = '/Users/sam/Documents/Dropbox/Dissertation/MATLAB work/Brian_example_PlotGeographicMap/cptfiles/tmp.cpt' ;
        ncols = linecount(cpt_file) - 1 ;
        if only_optd
            ncols_eachSide = (ncols-2)/2 ;
        else
            ncols_eachSide = (ncols-1)/2 ;
        end
        cptcmap(cpt_file)
        if ~isempty(caxis_lims)
            caxis(caxis_lims) ;
        end
        hcb1 = colorbar(cbarOrient) ;
        if ~isempty(fontSize)
            set(hcb1,'FontSize',fontSize)
        end
        newTickInt = (ncols-1)/ncols ;
        if only_optd
            newTicks = 0:newTickInt:ncols-1 ;
        else
            newTicks = 1:newTickInt:ncols ;
        end
        set(hcb1,'Ticks',newTicks)
        
        if isinf(cbar_minmax(1))
            cbar_bottom = '-Inf' ;
        elseif isnan(cbar_minmax(1))
            cbar_bottom = ' ' ;
        else
            cbar_bottom = ['-10^{' num2str(floor(log10(cbar_minmax(1)))) '} \times ' num2str(round(cbar_minmax(1)/10^floor(log10(cbar_minmax(1))),1))] ;
        end
        
        if isinf(cbar_minmax(2))
            cbar_top = 'Inf' ;
        elseif isnan(cbar_minmax(2))
            cbar_top = ' ' ;
        else
            cbar_top = ['10^{' num2str(floor(log10(cbar_minmax(2)))) '} \times ' num2str(round(cbar_minmax(2)/10^floor(log10(cbar_minmax(2))),1))] ;
        end
        
        if only_optd
            newTickLabels = {'Missing',cbar_bottom} ;
        else
            newTickLabels = {cbar_bottom} ;
        end
        for t = -(clim_int*ncols_eachSide):clim_int:-clim_int
            if t==-1
                newTickLabels = [newTickLabels '-10'] ;
            else
                newTickLabels = [newTickLabels ['-10^{' num2str(-t) '}']] ;
            end
        end
        for t = clim_int:clim_int:(clim_int*ncols_eachSide)
            if t==1
                newTickLabels = [newTickLabels '10'] ;
            else
                newTickLabels = [newTickLabels ['10^{' num2str(t) '}']] ;
            end
        end
        newTickLabels = [newTickLabels cbar_top] ;
        set(hcb1,'TickLabels',newTickLabels)
        if only_optd
            hcb2 = colorbar ;
            if ~isempty(fontSize)
                set(hcb2,'FontSize',fontSize)
            end
%             gca_pos = get(gca,'Position') ;
%             hcb_pos = get(hcb2,'Position') ;
            set(hcb2, ...
                     'Location', cbarOrient,...
                     'Ticks',newTicks(2:end), ...
                     'TickLabels', newTickLabels(2:end), ...
                     'YLim', newTicks([2 end])) ;
            hcb = hcb2 ;
            delete(hcb1) ;
%             set(gca,'Position',gca_pos)
%             set(hcb,'Position',hcb_pos)
        else
            hcb = hcb1 ;
        end
    end
else
    error('Add code to map this "colorbar"!')
end
% good thru here
if ~isempty(shapedata)
    geoshow(shapedata,'FaceColor','none','EdgeColor',edgeColor,'lineWidth',lineWidth);
    if ~isempty(shapedata2)
        geoshow(shapedata2,'FaceColor','none','EdgeColor','m',...edgeColor, ...
            'lineWidth',lineWidth, 'lineStyle','-');
    end
end

h = gca ;

set(h,'FontSize',fontSize)


    function hcb = do_colorbar(cbarOrient,fontSize,varargin)
        if ~isempty(varargin)
            cbarticks2 = varargin{1} ;
            cbarticklabels2 = varargin{2} ;
        end
        if ~strcmp(cbarOrient,'none')
            hcb = colorbar(cbarOrient) ;
            if ~isempty(fontSize)
                set(hcb,'FontSize',fontSize)
            end
            if ~isempty(varargin)
                set(hcb,'Ticks',cbarticks2,'TickLabels',cbarticklabels2)
            end
        else
            hcb = [] ;
        end
        if do_log10x_ticks
            gca_pos = get(gca,'Position') ;
            hcb_pos = get(hcb,'Position') ;
            ticks = get(hcb,'Ticks') ;
            ticklabels = get(hcb,'TickLabels') ;
            for tt = 1:length(ticks)
                if ticks(tt) ~= round(ticks(tt))
                    ticks(tt) = NaN ;
                elseif ticks(tt) == 0
                    ticklabels{tt} = '1' ;
                else
                    ticklabels{tt} = num2str(10^ticks(tt)) ;
                end
                if ~isnan(ticks(tt))
                    ticklabels{tt} = [ticklabels{tt} 'x'] ;
                end
            end
            ticklabels = {ticklabels{~isnan(ticks)}} ;
            ticks = ticks(~isnan(ticks)) ;
            ticklabels{1} = ['<=' ticklabels{1}] ;
            ticklabels{length(ticks)} = ['>=' ticklabels{length(ticks)}] ;
            set(hcb,'Ticks',ticks,'TickLabels',ticklabels)
            set(gca,'Position',gca_pos)
            set(hcb,'Position',hcb_pos)
        end
    end

    function n = linecount(filename)
        fid = fopen(filename,'r') ;
        n = 0;
        tline = fgetl(fid);
        while ischar(tline)
          tline = fgetl(fid);
          n = n+1;
        end
        fclose(fid) ;
    end


end % function

