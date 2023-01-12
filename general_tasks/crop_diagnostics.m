%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make some diagnostic LPJ-GUESS crop figures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSR 2023-01-12
% - Makes area-weighted global timeseries as well as time-averaged maps.
% - 

addpath(genpath(landsymm_lpjg_path()))

% crop_diagnostic_options.m must be somewhere on your path.
% There, specify the following variables:
%
% SETTINGS RELATED TO LPJ-GUESS RUN
%    inDir: Directory containing the LPJ-GUESS outputs.
%    thisFile: File we'll be importing. If *_st.out is present, use that instead of *.out.
%              E.g.: thisFile = 'gsirrigation.out' ;
%    filename_guess_landuse: The land use file to be used for area-weighted averaging in
%                            calibration. This is usually the same as what was used in
%                            the actual runs, but doesn't need to be. It just needs to
%                            have columns for all the crops in the LPJ-GUESS output.
%                            E.g.: filename_guess_landuse = '/Users/Shared/LandSyMM/inputs/LU/remaps_v10_old_62892_gL/LU.remapv10_old_62892_gL.txt' ;
%    filename_guess_cropfrac: As filename_guess_landuse, but for crop fractions instead
%                             of general land use fractions.
%                             E.g.: filename_guess_cropfrac = '/Users/Shared/LandSyMM/inputs/LU/remaps_v10_old_62892_gL/cropfracs.remapv10_old_62892_gL.txt' ;
%    xres and yres: The longitude and latitude resolution.
%                   E.g.: xres = 0.5; yres = 0.5;
%
% OTHER INPUT FILES
%    filename_staticData_quarterdeg: Path to staticData_quarterdeg.nc.
%
% SETTINGS FOR TIMESERIES FIGURES
% These are all fields in the timeseries_opts structure.
%    perarea: Plot per-area values (e.g., yield in tons/ha)? If not, will plot totals
%             (e.g., production in Mt).
%             E.g.: timeseries_opts.perarea = true ;
%    figure_window_position: For figure() 'Position' property (location and size of 
%                            drawable area). To get these values, make an empty figure
%                            window with figure(), then move and resize the window to
%                            your liking. Then do get('gcf', 'Position'), and that will
%                            give you the numbers you need.
%                            E.g.: timeseries_opts.figure_window_position = [56 78 1139 563] ;
%    fontSize: Font size for plots.
%              E.g.: timeseries_opts.fontSize = 14 ;
%    lineWidth: Line width for plots.
%               E.g.: timeseries_opts.lineWidth = 2 ;
%
% SETTINGS FOR MAPS FIGURES
% These are all fields in the map_opts structure.
%    perarea: As for "SETTINGS FOR TIMESERIES FIGURES" above.
%    figure_window_position: As for "SETTINGS FOR TIMESERIES FIGURES" above.
%    fontSize: As for "SETTINGS FOR TIMESERIES FIGURES" above.
%    incl_years: Years to take the mean over.
%                E.g.: map_opts.incl_years = 1995:2005 ;
%    lineWidth: Width of line used to draw coastlines on maps.
%               E.g.: map_opts.lineWidth = 0.5 ;
%    subplot_spacing: Spacing between subplots, [vertical horizontal].
%                     E.g.: map_opts.subplot_spacing = [0.075 0.05] ;

crop_diagnostics_options


%% Setup

% Output directory
outDir = fullfile(inDir, 'figs') ;
fprintf('Saving figures to %s\n', outDir)
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

% Map info
map_opts.lons = (-180+xres/2):xres:(180-xres/2) ;
map_opts.lats = (-90+yres/2):yres:(90-yres/2) ;
map_opts.mapSize = [180/yres 360/xres] ;


%% Import

% Import LPJ-GUESS output
fprintf('Importing %s...\n', thisFile)
S = lpjgu_matlab_read2geoArray(sprintf('%s/%s', inDir, thisFile)) ;

% Import land use and crop fractions, ensuring match to LPJ-GUESS gridlist
disp('Importing land use fractions...')
lu_frac = lpjgu_matlab_read2geoArray(filename_guess_landuse) ;
if ~isequal(S.list2map, lu_frac.list2map)
    lu_frac = harmonize_gridlists(S, lu_frac) ;
end
lu_frac = harmonize_yearLists(S, lu_frac) ;
disp('Importing crop fractions...')
crop_frac = lpjgu_matlab_read2geoArray(filename_guess_cropfrac) ;
if ~isequal(S.list2map, crop_frac.list2map)
    crop_frac = harmonize_gridlists(S, crop_frac) ;
end
crop_frac = harmonize_yearLists(S, crop_frac) ;

% Import cell area
disp('Importing cell area...')
cell_area_YXqd = transpose(ncread(filename_staticData_quarterdeg, 'carea')) ;
cell_area_YX = aggregate_cell_area(cell_area_YXqd, xres, yres) ;
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

crop_diagnostics_options
make_timeseries(S, crop_area, crop_frac, outDir, thisFile, timeseries_opts)


%% Make maps for each crop stand type

crop_diagnostics_options
make_maps(S, crop_area, crop_frac, outDir, thisFile, map_opts)


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


% Plot area-weighted global time series for each crop stand type
function make_timeseries(S, crop_area, crop_frac, outDir, thisFile, opts)

% Process units
if opts.perarea
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
    opts.perarea_xvy = S.garr_xvy(:,IB,:) ;

    % Get crop areas
    [C, ~, IB] = intersect(theseVars, crop_area.varNames, 'stable') ;
    if ~isequal(C, theseVars)
        error('Some crop(s) in LPJ-GUESS output data not found in crop area input')
    end
    area_xvy = crop_area.garr_xvy(:,IB,:) ;

    % Get totals
    totals_xvy = opts.perarea_xvy .* area_xvy ;
    totals_xvy(:,end+1,:) = sum(totals_xvy, 2) ; %#ok<AGROW> 
    theseVars{end+1} = 'Combined' ; %#ok<AGROW>

    % Get global sums
    data_vy = squeeze(sum(totals_xvy, 1)) ;
    
    % Convert back to per-area value, if needed
    if opts.perarea
        area_xvy(:,end+1,:) = sum(area_xvy, 2) ; %#ok<AGROW> 
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
    figure('Color', 'w', 'Position', opts.figure_window_position)
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
            'lineWidth', opts.lineWidth)
        hold off
    end
    legend(thisLegend, 'Location', 'Best')
    xlabel('Year')
    ylabel(units)
    set(gca, 'fontSize', opts.fontSize)
    thisTitle = sprintf('%s global %s', thisCrop_rf, titleName) ;
    title(thisTitle)
    
    % Save figure
    outFile = fullfile(outDir, ['Timeseries ' thisTitle, '.pdf']) ;
    export_fig(outFile, '-r300')
    close
end
disp('Done.')


end


function make_maps(S, crop_area, crop_frac, outDir, thisFile, opts)

% Process units
if opts.perarea
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

% Get included years
if ~isfield(S, 'yearList')
    error('Need to code handling for S with no year dimension')
end
[C, ~, Iyears] = intersect(opts.incl_years, S.yearList, 'stable') ;
if ~isequal(C(:), opts.incl_years(:))
    error('Year list mismatch')
end
y1 = min(opts.incl_years) ;
yN = max(opts.incl_years) ;

fprintf('Mapping %d-%d:\n', y1, yN)
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
    opts.perarea_xvy = S.garr_xvy(:,IB,:) ;

    % Get crop areas
    [C, ~, IB] = intersect(theseVars, crop_area.varNames, 'stable') ;
    if ~isequal(C, theseVars)
        error('Some crop(s) in LPJ-GUESS output data not found in crop area input')
    end
    area_xvy = crop_area.garr_xvy(:,IB,:) ;

    % Get totals
    totals_xvy = opts.perarea_xvy .* area_xvy ;
    totals_xvy(:,end+1,:) = sum(totals_xvy, 2) ; %#ok<AGROW> 
    theseVars{end+1} = 'Combined' ; %#ok<AGROW>

    % Get mean over selected years
    data_xv = mean(totals_xvy(:,:,Iyears), 3) ;

    % Get areas
    area_xvy(:,end+1,:) = sum(area_xvy, 2) ; %#ok<AGROW> 
    area_xv = mean(area_xvy(:,:,Iyears), 3) ;
    
    % Convert back to per-area value, if needed
    if opts.perarea
        data_xv = data_xv ./ area_xv ;
        data_xv(area_xv == 0) = NaN ;
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
    figure('Color', 'w', 'Position', opts.figure_window_position)
    data_xv = data_xv * conversion_factor ;
    new_caxis = [nanmin(data_xv(:)) nanmax(data_xv(:))] ;
    for v = 1:length(theseVars)
        thisVar = thisLegend{v} ;
        switch thisVar
            case 'Rainfed'
                nx = 2 ;
                p = 1 ;
            case 'Irrigated'
                nx = 2 ;
                p = 2 ;
            case 'Combined' 
                nx = 1 ;
                p = 2 ;
            otherwise
                error('%s not recognized for lineStyle', thisVar)
        end
        map_YX = lpjgu_vector2map(data_xv(:,v), opts.mapSize, S.list2map) ;
        area_YX = lpjgu_vector2map(area_xv(:,v), opts.mapSize, S.list2map) ;
        map_YX(area_YX == 0) = NaN ;
        subplot_tight(2, nx, p, opts.subplot_spacing)

        worldmap('World')
        load coastlines
        plotm(coastlat,coastlon, 'k')
        pcolorm(opts.lats, -179.75:0.5:179.75, map_YX)
        framem('FLineWidth', opts.lineWidth)
        mlabel off; plabel off

        if ~isequal(new_caxis, [0 0]) && ~isequaln(new_caxis, [NaN NaN])
            caxis(new_caxis)
        end
        hcb = colorbar ;
        xlabel(hcb, units)
        set(gca, 'FontSize', opts.fontSize)
        title(thisLegend{v})
    end
    thisTitle = sprintf('%s global %s %d-%d', thisCrop_rf, titleName, y1, yN) ;
    hsgt = sgtitle(thisTitle, 'FontSize', opts.fontSize*1.25, 'FontWeight', 'bold') ;
    
    % Save figure
    outFile = fullfile(outDir, ['Map ' thisTitle '.png']) ;
    export_fig(outFile, '-r150')
    close
end
disp('Done.')


end
