function thisMap_YX = map_contribution(maps_YXr, shapefile_path, theseRuns, runList, varargin)

% Set up & parse input arguments
p = inputParser ;
is2elemvec = @(x) (length(x)==2 && ismatrix(x) && min(size(x))==1) || isempty(x) ;
is3elemvec01 = @(x) length(x)==3 && ismatrix(x) && min(size(x))==1 && min(x)>=0 && max(x)<=1 ;
is1x1num = @(x) isequal(size(x), [1 1]) ;


% Required input arguments
addRequired(p, 'maps_YXr', @(x) isnumeric(x) && ndims(x)==3) ;
addRequired(p, 'shapefile_path', @ischar) ;
addRequired(p, 'theseRuns', @(x) (iscellstr(x) && length(x)<=2) || ischar(x)) ;
addRequired(p, 'runList', @iscellstr) ;

% Optional input arguments
addOptional(p, 'lonlim', [-180 180], is2elemvec) ;
addOptional(p, 'latlim', [-90 90], is2elemvec) ;
addOptional(p, 'thisColormap', 'rdbu_ssr', @isstr) ;
addOptional(p, 'fontSize', [], is1x1num) ;
addOptional(p, 'edgeColor', [0.1, 0.1, 0.1], is3elemvec01) ;
addOptional(p, 'lineWidth', 0.1, is1x1num) ;
addOptional(p, 'cbarOrient', 'eastoutside', @isstr) ;
addOptional(p, 'caxis_lims', [], is2elemvec) ;
addOptional(p, 'units_map', '', @isstr) ;
addOptional(p, 'thisTitle', '', @isstr) ;
addOptional(p, 'nanmask_YX', [], @(x) ismatrix(x) || isempty(x)) ;


% Parse inputs
parse(p, maps_YXr, shapefile_path, theseRuns, runList, varargin{:});
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

% Get map to plot
if length(theseRuns)==2
    run1 = find(contains(runList, theseRuns{1})) ;
    run2 = find(contains(runList, theseRuns{2})) ;
    if length(run1)~=1
        error('length(run1)~=1')
    elseif length(run2)~=1
        error('length(run2)~=1')
    end
    thisMap_YX = maps_YXr(:, :, run1) - maps_YXr(:, :, run2) ;
else
    run1 = find(contains(runList, theseRuns)) ;
    if length(run1)~=1
        error('length(run1)~=1')
    end
    thisMap_YX = maps_YXr(:,:,run1) ;
end
if ~isempty(nanmask_YX)
    thisMap_YX(nanmask_YX) = NaN ;
end

% Make figure
figure('Color', 'w', 'Position', figurePos);
[~, hcb] = map_with_SHPoverlay_v2(thisMap_YX, shapefile_path, ...
    'lonlim', lonlim, ...
    'latlim', latlim, ...
    'thisColormap', thisColormap, ...
    'fontSize', fontSize, ...
    'edgeColor', edgeColor, ...
    'lineWidth', lineWidth, ...
    'cbarOrient', cbarOrient, ...
    'caxis_lims', caxis_lims) ;

% Edit figure
if strcmpi(thisColormap, 'RdBu_ssr')
    colormap(brighten(brewermap(64, 'RdBu_ssr'), -0.3))
end
if ~isempty(thisTitle)
    title(thisTitle)
end
hcb.Label.String = units_map ;


end