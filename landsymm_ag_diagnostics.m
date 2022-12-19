%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make some diagnostic LandSyMM-ag figures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Settings for individual crop runs %%%%%

% Directory from which we'll be importing
inDir = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851' ;

% File we'll be importing
%   - If *_st.out present, use that instead of *.out
% thisFile = 'yield.out' ;
thisFile = 'gsirrigation.out' ;

% Land use info
filename_landuse = '/Users/Shared/PLUM/input/remaps_v6p7/LU.remapv6p7.txt' ;
filename_cropfrac = '/Users/Shared/PLUM/input/remaps_v6p7/cropfracs.remapv6p7.txt' ;


%%%%% General settings %%%%%

filename_staticData_quarterdeg = '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc' ;


%% Setup

% Output directory
outDir = sprintf('%s/figs', inDir) ;
if ~exist(outDir, 'dir')
    mkdir(outDir)
end


%% Import

% Import LPJ-GUESS output
fprintf('Importing %s...\n', thisFile)
S = lpjgu_matlab_read2geoArray(sprintf('%s/%s', inDir, thisFile)) ;

% Import land use and crop fractions, ensuring match to LPJ-GUESS gridlist
disp('Importing land use fractions...')
lu_frac = lpjgu_matlab_read2geoArray(filename_landuse) ;
if ~isequal(S.list2map, lu_frac.list2map)
    lu_frac = harmonize_gridlists(S, lu_frac) ;
end
lu_frac = harmonize_yearLists(S, lu_frac) ;
disp('Importing crop fractions...')
crop_frac = lpjgu_matlab_read2geoArray(filename_cropfrac) ;
if ~isequal(S.list2map, crop_frac.list2map)
    crop_frac = harmonize_gridlists(S, crop_frac) ;
end
crop_frac = harmonize_yearLists(S, crop_frac) ;

% Import cell area
disp('Importing cell area...')
cell_area_YXqd = transpose(ncread(filename_staticData_quarterdeg, 'carea')) ;
cell_area_YX = aggregate_cell_area(cell_area_YXqd, 0.5, 0.5) ;
% Convert km2 to m2
cell_area_YX = cell_area_YX * 1e6 ;
% Vectorize
cell_area_x = cell_area_YX(S.list2map) ;

% Get land use and crop area
disp('Getting land use areas...')
lu_area = get_lu_area(lu_frac, cell_area_x) ;
crop_area = get_lu_area(crop_frac, squeeze(lu_area.garr_xvy(:,strcmp(lu_frac.varNames, 'CROPLAND'),:))) ;

disp('Done importing.')


%% Plot area-weighted global time series for each crop stand type
%%%%%%%%%%%%%%%%%%
% Plot per-area values (e.g., yield in tons/ha)? If not, will plot totals
% (e.g., production in Mt).
perarea = true ;

% Other options
figure_window_position = [56 78 1139 563] ; % From get(gcf, 'Position')
fontSize = 14 ;
lineWidth = 2 ;
%%%%%%%%%%%%%%%%%%

% Process units
if perarea
    if contains(thisFile, 'yield')
        units = 'tons dry matter ha^{-1}' ;
        conversion_factor = 10 ; % From kg/m2
        titleName = 'yield' ;
    elseif contains(thisFile, 'gsirrig')
        units = 'mm' ;
        conversion_factor = 1 ; % Native LPJ-GUESS output unit
        titleName = 'average irrigation' ;
    else
        error('What units etc. should be used for %s per-area?', thisFile)
    end
else
    if contains(thisFile, 'yield')
        units = 'Mt DM' ;
        conversion_factor = 1e-9 ; % From kg
        titleName = 'production' ;
    elseif contains(thisFile, 'gsirrig')
        units = 'km^3' ;
        conversion_factor = 1e-12 ; % From mm*m2
        titleName = 'irrigation volume' ;
    else
        error('What units etc. should be used for %s totals?', thisFile)
    end
end

% Get crop list
is_irrigated = @(x) strcmp(x(end), 'i') ;
cropList = crop_frac.varNames(~cellfun(is_irrigated, crop_frac.varNames)) ;

% Sanity check
if ~isfield(S, 'yearList')
    error('You''re trying to make a time series, but S appears to be for just one timestep')
end

disp('Plotting area-weighted global time series:')
for c = 1:length(cropList)

    % Find stands of this crop
    thisCrop_rf = cropList{c} ;
    fprintf('    %s...\n', thisCrop_rf)
    thisCrop_ir = [thisCrop_rf 'i'] ;
    [theseVars, ~, IB] = intersect({thisCrop_rf, thisCrop_ir}, S.varNames, 'stable') ;
    if isempty(theseVars)
        error('No %s stands found in LPJ-GUESS output')
    end

    % Get LPJ-GUESS output data
    perarea_xvy = S.garr_xvy(:,IB,:) ;

    % Get crop areas
    [C, ~, IB] = intersect(theseVars, crop_area.varNames, 'stable') ;
    if ~isequal(C, theseVars)
        error('Some crop(s) in LPJ-GUESS output data not found in crop area input')
    end
    area_xvy = crop_area.garr_xvy(:,IB,:) ;

    % Get totals
    totals_xvy = perarea_xvy .* area_xvy ;
    totals_xvy(:,end+1,:) = sum(totals_xvy, 2) ; %#ok<SAGROW> 
    theseVars{end+1} = 'Combined' ; %#ok<SAGROW>

    % Get global sums
    data_vy = squeeze(sum(totals_xvy, 1)) ;
    
    % Convert back to per-area value, if needed
    if perarea
        area_xvy(:,end+1,:) = sum(area_xvy, 2) ; %#ok<UNRCH,SAGROW> 
        area_vy = squeeze(sum(area_xvy, 1)) ;
        data_vy = data_vy ./ area_vy ;
        data_vy(area_vy == 0) = NaN ;
    end

    % Get legend
    thisLegend = theseVars ;
    if any(strcmp(thisLegend, thisCrop_rf))
        thisLegend{strcmp(thisLegend, thisCrop_rf)} = 'Rainfed' ;
    end
    if any(strcmp(thisLegend, thisCrop_ir))
        thisLegend{strcmp(thisLegend, thisCrop_ir)} = 'Irrigated' ;
    end

    % Plot
    figure('Color', 'w', 'Position', figure_window_position)
    data_vy = data_vy * conversion_factor ;
    for v = 1:length(theseVars)
        thisVar = thisLegend{v} ;
        switch thisVar
            case 'Rainfed' ; lineStyle = '--r' ;
            case 'Irrigated' ; lineStyle = ':b' ;
            case 'Combined' ; lineStyle = '-k' ;
            otherwise
                error('%s not recognized for lineStyle', thisVar)
        end
        hold on
        plot(S.yearList, data_vy(v,:), lineStyle, ...
            'LineWidth', lineWidth)
        hold off
    end
    legend(thisLegend, 'Location', 'Best')
    xlabel('Year')
    ylabel(units)
    set(gca, 'FontSize', fontSize)
    thisTitle = sprintf('%s global %s', thisCrop_rf, titleName) ;
    title(thisTitle)
    
    % Save figure
    outFile = sprintf('%s/%s.pdf', outDir, thisTitle) ;
    export_fig(outFile, '-r300')
    close
end
disp('Done.')


%% FUNCTIONS

function S = get_lu_area(S, area_x)

ndims_area = length(find(size(area_x) > 1)) ;

if isfield(S, 'garr_xvy')
    if ndims_area == 2
        if size(area_x, 2) ~= size(S.garr_xvy, 3)
            error('Mismatch between area_x and S.garr_xvy on Year dimension')
        end
        area_x1y = permute(area_x, [1 3 2]) ;
        S.garr_xvy = S.garr_xvy .* repmat(area_x1y, [1 size(S.garr_xvy, 2) 1]) ;
    elseif ndims_area == 1
        S.garr_xvy = S.garr_xvy .* repmat(area_x, [1 size(S.garr_xvy, 2:3)]) ;
    else
        error('Too many dimensions on area_x')
    end
else
    if ndims_area > 1
        error('Too many dimensions on area_x for one-timstep S')
    end
    S.garr_xv = S.garr_xv .* repmat(area_x, [1 size(S.garr_xvy, 2)]) ;
end

end

% Ensure that two geo-array structures have the same yearlist
function S_newyears = harmonize_yearLists(S_target, S_newyears)

if ~isfield(S_target, 'yearList') && ~isfield(S_newyears, 'yearList')
    return
elseif ~isfield(S_target, 'yearList')
    error('S_target has no yearList but S_newyears does')
end

disp('    Harmonizing year lists...')

if ~isfield(S_newyears, 'yearList')
    garr_xv = S_newyears.garr_xv ;
    S_newyears = rmfield(S_newyears, 'garr_xv') ;
    S_newyears.yearList = S_target.yearList ;
    S_newyears.garr_xvy = repmat(garr_xv, [1 1 length(S_target.yearList)]) ;
else
    [C, ~, IB] = intersect(S_target.yearList, S_newyears.yearList, 'stable') ;
    if ~isequal(C, S_target.yearList)
        error('yearList mismatch')
    elseif length(C) < length(S_newyears.yearList)
        fprintf('    (Trimming %d years not present in S_target.yearList)\n', length(S_newyears.yearList) - length(C))
    end
    S_newyears.garr_xvy = S_newyears.garr_xvy(:,:,IB) ;
    S_newyears.yearList = C ;
end

end

% Ensure that two geo-array structures have the same gridlist
function S_rearr = harmonize_gridlists(S_target, S_rearr)

if ~isequal(S_target.list2map, S_rearr.list2map)
    % Test whether rearranging will fix it
    [C, ~, IB] = intersect(S_target.list2map, S_rearr.list2map, 'stable') ;
    if ~isequal(S_target.list2map, C)
        error('Gridlist mismatch')
    end
    % Rearrange to fix
    disp('    Rearranging gridlist to match target...')
    if isfield(S_rearr, 'garr_xv')
        S_rearr.garr_xv = S_rearr.garr_xv(IB,:) ;
    elseif isfield(S_rearr, 'garr_xvy')
        S_rearr.garr_xvy = S_rearr.garr_xvy(IB,:,:) ;
    else
        error('No apparent geo-array field in S')
    end
    S_rearr.lonlats = S_rearr.lonlats(IB,:) ;
    S_rearr.list2map = C ;
end

end