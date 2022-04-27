%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check patchiness Bart found %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hook://email/AM8PR05MB80346516E607A885B18661678CF89%40AM8PR05MB8034.eurprd05.prod.outlook.com

% topDir = '/Users/Shared/landsymm_forestry/1970past/1850-2014/output-2022-04-15-182232' ;
% topDir = '/Users/Shared/landsymm_forestry/1970past/1850-2014_germany/output-2022-04-26-233030' ;
% topDir = '/Users/Shared/landsymm_forestry/1970past/1850-2014_germany_nohistluc/output-2022-04-27-002417' ;
topDir = '/Users/Shared/landsymm_forestry/1970past/1850-2014_germany_nohistluc_oldprodpart/output-2022-04-27-004851' ;
estYr = 1970 ;

LUlist = {'forC', 'ntrl'} ;
LUlist_long = {'forest', 'natural'} ;


%% Setup

Nlu = length(LUlist) ;


%% Import

% plutW
thisPattern = sprintf('%s/landsymm_plutW_from_%%s.out', topDir) ;
plutW = read_LUlist(thisPattern, LUlist) ;

% ctree
ctree = lpjgu_matlab_read2geoArray(sprintf('%s/landsymm_ctree_sts.out', topDir)) ;

% cmass_wood
cmass_wood = lpjgu_matlab_read2geoArray(sprintf('%s/cmass_wood_sts.out', topDir)) ;

% cmass_wood_potharv
cmass_wood_potharv = lpjgu_matlab_read2geoArray(sprintf('%s/cmass_wood_potharv_sts.out', topDir)) ;

disp('Done importing.')


%% Make maps of 45-year yield

thisYr = estYr + 44 ;
nx = 3 ;

% Country bounding boxes: https://gist.github.com/graydon/11198540
% bbox = [] ; % Use all cells
bbox = [5.98865807458, 47.3024876979, 15.0169958839, 54.983104153] ; bbox_name = 'Germany' ;
thisPos = [1   484   957   517] ; spacing = [0.1 0.05] ;

ny = 2 ;
fontSize = 14 ;

figure('Color', 'w', 'Position', thisPos) ;

for L = 1:Nlu
    thisLU = LUlist{L} ;

    y = find(plutW.yearList == thisYr) ;
    if length(y) ~= 1
        error('Expected to find 1 match of %d in plutW.yearList; found %d', thisYr, length(y))
    end

    make_map(ny, nx, (L-1)*nx+1, plutW, thisLU, L, y, 'plutW', fontSize, bbox, spacing)

    make_map(ny, nx, (L-1)*nx+2, ctree, thisLU, L, y, 'ctree', fontSize, bbox, spacing)

%     make_map(ny, nx, (L-1)*nx+3, cmass_wood, thisLU, L, y, 'cmass wood', fontSize, bbox, spacing)
    make_map(ny, nx, (L-1)*nx+3, cmass_wood_potharv, thisLU, L, y, 'cmass wood potharv', fontSize, bbox, spacing)

    sgtitle(bbox_name, 'FontSize', fontSize+4, 'FontWeight', 'bold')

end


%% Scatter plot for a bounding box, same-year same-LU

thisYr = estYr + 44 ;
thisLU = 'ntrl' ;

% Country bounding boxes: https://gist.github.com/graydon/11198540
% bbox = [] ; % Use all cells
bbox = [5.98865807458, 47.3024876979, 15.0169958839, 54.983104153] ; bbox_name = 'Germany' ;

xdata = get_scatter_data(cmass_wood, bbox, thisYr, thisLU) ; xlab = 'cmass wood' ;

ydata = get_scatter_data(plutW, bbox, thisYr, thisLU, LUlist, 'to_forC') ; ylab = 'PLUT wood to-forC'  ;
% ydata = get_scatter_data(cmass_wood_potharv, bbox, thisYr, thisLU) ; ylab = 'cmass wood potharv' ;


make_scatter_plot(xdata, ydata, bbox_name, thisYr, thisLU, xlab, ylab) ;


%% Save bbox_lonlats as gridlist

file_out = sprintf('%s/../../../gridlist_bbox_%s.txt', topDir, bbox_name) ;
save_file(bbox_lonlats(:,1), bbox_lonlats(:,2), file_out) ;


%% FUNCTIONS

function sdata = get_scatter_data(S, bbox, thisYr, thisLU, varargin)

if isfield(S, 'garr_xvyL')
    LUlist = varargin{1} ;
    toLU = varargin{2} ;
elseif ~isempty(varargin)
    warning('Ignoring optional arguments because S does not have garr_xvyL')
end

y = find(S.yearList == thisYr) ;
if length(y) ~= 1
    error('Expected to find 1 match of %d in S.yearList; found %d', thisYr, length(y))
end

sdata = extract_bbox(S, bbox) ;

if isfield(S, 'garr_xvyL')
    v = find(strcmp(S.varNames, toLU)) ;
    if length(v) ~= 1
        error('Expected to find 1 match of toLU (%s) in S.varNames; found %d', toLU, length(v))
    end
    L = find(strcmp(LUlist, thisLU)) ;
    if length(L) ~= 1
        error('Expected to find 1 match of thisLU (%s) in LUlist; found %d', thisLU, length(L))
    end
    sdata = sdata.garr_xvyL(:,v,y,L) ;
else
    v = find(strcmp(S.varNames, thisLU)) ;
    if length(v) ~= 1
        error('Expected to find 1 match of thisLU (%s) in S.varNames; found %d', thisLU, length(v))
    end
    sdata = sdata.garr_xvy(:,v,y) ;
end

end


function make_map(ny, nx, n, S, thisLU, L, y, title_var, fontSize, bbox, spacing)

subplot_tight(ny, nx, n, spacing) ;

if isfield(S, 'garr_xvy')
    v = find(strcmp(S.varNames, thisLU)) ;
    thisMap = lpjgu_vector2map(S.garr_xvy(:,v,y), [360 720], S.list2map) ;
else
    v = find(strcmp(S.varNames, 'to_forC')) ;
    if length(v) ~= 1
        error('Expected to find 1 match of to_forC in S.varNames; found %d', length(v))
    end
    thisMap = lpjgu_vector2map(S.garr_xvyL(:,v,y,L), [360 720], S.list2map) ;
end
thisMap = crop_to_bbox(thisMap, bbox) ;

pcolor(thisMap); shading flat; axis equal tight
colorbar
title(sprintf('%s %s', thisLU, title_var))
set(gca, 'FontSize', fontSize)

end


function thisMap = crop_to_bbox(thisMap, bbox)

lon_to_x = @(x) (x+180)*2 ;
lat_to_y = @(x) (x+90)*2 ;

if ~isempty(bbox)
    xlims = [max(1,floor(lon_to_x(bbox(1)))-1) min(720,ceil(lon_to_x(bbox(3)))+1)] ;
    ylims = [max(1,floor(lat_to_y(bbox(2)))-1) min(360,ceil(lat_to_y(bbox(4)))+1)] ;
    thisMap = thisMap(ylims(1):ylims(2),xlims(1):xlims(2)) ;
end

end

function save_file(lons, lats, file_out)

fid = fopen(file_out,'w') ;
fprintf(fid,'%0.2f %0.2f\n',[lons lats]') ;
fclose(fid) ;

end


function S = extract_bbox(S, bbox)

if isempty(bbox)
    return
end

lons = S.lonlats(:,1) ;
lats = S.lonlats(:,2) ;

lon_ok = lons >= bbox(1) & lons <= bbox(3) ;
lat_ok = lats >= bbox(2) & lats <= bbox(4) ;
incl = lon_ok & lat_ok ;

if isfield(S, 'garr_xvy')
    S.garr_xvy = S.garr_xvy(incl,:,:) ;
elseif isfield(S, 'garr_xvyL')
    S.garr_xvyL = S.garr_xvyL(lon_ok & lat_ok,:,:,:) ;
else
    error('garr field not found')
end
S.lonlats = S.lonlats(incl,:) ;
S.list2map = S.list2map(incl) ;

end

function S = read_LUlist(thisPattern, LUlist)

Nlu = length(LUlist) ;

for L = 1:Nlu
    thisLU = LUlist{L} ;
    thisFile= sprintf(thisPattern, thisLU) ;
    tmp = lpjgu_matlab_read2geoArray(thisFile) ;
    
    if L == 1
        S = rmfield(tmp, "garr_xvy") ;
        S.garr_xvyL = nan([size(tmp.garr_xvy) Nlu]) ;
    else
        if ~isequal(tmp.list2map, S.list2map)
            error('Gridlist mismatch')
        elseif ~isequal(tmp.yearList, S.yearList)
            error('Yearlist mismatch')
        end
    end

    S.garr_xvyL(:,:,:,L) = tmp.garr_xvy ;
    clear tmp
end

end

function make_scatter_plot(xdata, ydata, bbox_name, thisYr, thisLU, xlab, ylab)
figure('Color', 'w') ;
plot(xdata, ydata, '.')
axis equal
hold on
if min(get(gca, 'XLim')) < 0
    set(gca, 'XLim', [0, max(get(gca, 'XLim'))])
    set(gca, 'YLim', [0, max(get(gca, 'YLim'))])
end
if min(get(gca, 'YLim')) < 0
    set(gca, 'YLim', [0, max(get(gca, 'YLim'))])
    set(gca, 'XLim', [0, max(get(gca, 'XLim'))])
end
min11 = min(min(get(gca, 'XLim')), min(get(gca, 'YLim'))) ;
max11 = max(max(get(gca, 'XLim')), max(get(gca, 'YLim'))) ;
plot([min11, max11], [min11, max11], '--k')
hold off
title(sprintf('%s: %s %d', bbox_name, thisLU, thisYr))
xlabel(xlab)
ylabel(ylab)
set(gca, 'FontSize', 14)
end