function ESintxns_subplots_biomes( ...
    xdata_YX, ydata_YX, biome_names, biome_maps_YXv, ...
    varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parsing arguments %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

is1x1num = @(x) numel(x)==1 & isnumeric(x) ;
is1x1int = @(x) numel(x)==1 & isint(x) ;
% is1x1logical = @(x) numel(x)==1 & (x==0 | x==1) ;
% iscellorempty = @(x) iscell(x) | isempty(x) ;
iscolor = @(x) (isnumeric(x) & isequal(size(shiftdim(x)),[1 3])) ...
    | (ischar(x) & length(x)==1) ;
isAlpha = @(x) is1x1num(x) & x>=0 & x<=1 ;
isSpacing = @(x) (numel(x)==1 | numel(x)==2) & isnumeric(x) ;

p = inputParser ;
addRequired(p,'xdata_YX',@ismatrix) ;
addRequired(p,'ydata_YX',@ismatrix) ;
addRequired(p,'biome_names',@iscellstr) ;
addRequired(p,'biome_maps_YXv',@islogical) ;

% Subplots
addParameter(p,'spacing',[0.1 0.05], isSpacing) ; % [vertical horizontal]

% Text
addParameter(p,'fontSize',14,is1x1int) ;
addParameter(p,'xlab','',@ischar) ;
addParameter(p,'ylab','',@ischar) ;
addParameter(p,'titl_prefix','',@ischar) ;

% Points
addParameter(p,'marker','.',@ischar) ;
addParameter(p,'markerSize',30,is1x1int) ; % Size in one-d points
addParameter(p,'markerColor','b',iscolor) ;
addParameter(p,'markerAlpha',1,isAlpha) ;

% Best-fit line
addParameter(p,'bestFitType', 'none', @ischar)
addParameter(p,'lineStyle','-',@ischar) ;
addParameter(p,'lineWidth',2,is1x1int) ; % Size in one-d points
addParameter(p,'lineColor','r',iscolor) ;

% Parse inputs
parse(p,...
    xdata_YX, ydata_YX, biome_names, biome_maps_YXv, ...
    varargin{:});
pFields = fieldnames(p.Results) ;
Nfields = length(pFields) ;
for f = 1:Nfields
    thisField = pFields{f} ;
    if ~exist(thisField,'var')
        eval([thisField ' = p.Results.' thisField ' ;']) ;
    end
    clear thisField
end ; clear f
clear p


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The actual function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up subplot arrangement
Nbiomes = length(biome_names) ;
if Nbiomes > 9
    error('Set up nx and ny for Nbiomes > 9.')
elseif Nbiomes > 6
    ny = 3 ;
    nx = 3 ;
elseif Nbiomes > 4
    ny = 2 ;
    nx = 3 ;
elseif Nbiomes == 4
    ny = 2 ;
    nx = 2 ;
elseif Nbiomes >= 1
    ny = 1 ;
    nx = Nbiomes ;
else
    error('No biomes??')
end

% Harmonize NaNs
isBad_YX = ~isequal(isnan(xdata_YX), isnan(ydata_YX)) ;
if any(isBad_YX)
    warning('Removing %d cells to harmonize xdata_YX and ydata_YX.', length(find(isBad_YX))) ;
    xdata_YX(isBad_YX) = NaN ;
    ydata_YX(isBad_YX) = NaN ;
end
biome_maps_YXv(repmat(isnan(xdata_YX), [1 1 Nbiomes])) = false ;


for b = 1:Nbiomes
    
    subplot_tight(ny, nx, b, spacing) ;
    
    % Setup
    thisBiome = biome_names{b} ;
    thisBiome_YX = biome_maps_YXv(:,:,b) ;

    % Points
    xdata = xdata_YX(thisBiome_YX) ;
    ydata = ydata_YX(thisBiome_YX) ;
    ESintxns_scatter( ...
        xdata, ydata, ...
        'marker', marker, ...
        'markerSize', markerSize, ...
        'markerColor', markerColor, ...
        'markerAlpha', markerAlpha) ;
    
    % Text
    set(gca,'FontSize',fontSize) ;
    if ~isempty(xlab)
        xlabel(xlab) ;
    end
    if ~isempty(ylab)
        ylabel(ylab) ;
    end
    title(thisBiome) ;
    
    % Best-fit line
    if ~strcmp(bestFitType,'none')
        ESintxns_bestFit(xdata, ydata, bestFitType) ;
    end

end

subplot_tight(ny, nx, Nbiomes+1, spacing) ;
imshow('/Users/sam/Geodata/KoeppenGeiger/Koeppen-Geiger_v1.1/koeppen-geiger_0.5_QGISprint.png')

if ~isempty(titl_prefix)
    ht = sgtitle(titl_prefix,'FontSize',round(fontSize*1.5),'FontWeight','bold') ;
end


end