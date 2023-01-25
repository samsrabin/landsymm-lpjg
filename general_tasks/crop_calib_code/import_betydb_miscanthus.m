function miscanthus_x_os = import_betydb_miscanthus(dir_data, yield_lpj, year1, yearN)
% Get an X by 2 array:
%    - Cols: observed (from betyDB), simulated (not actually Miscanthus, but some proxy)
%    - Rows: each one gridcell, rainfed or irrigated
% This breaks with the usual calibration factor method of aggregating to country level.
% This also breaks with the usual method of each year being its own point---means across
% all of year1–yearN are calculated here.

fprintf('Getting Miscanthus yields for %d-%d...\n', year1, yearN)

% Get years of interest
I = find(yield_lpj.yearList >= year1 & yield_lpj.yearList <= yearN) ;
Nyears = yearN - year1 + 1 ;
if length(I) ~= Nyears
    error('Expected %d years in %d-%d; found %d', Nyears, year1, yearN, length(I))
end
yield_lpj.maps_YXvy = yield_lpj.maps_YXvy(:,:,:,I) ;
yield_lpj.yearList = yield_lpj.yearList(I) ;

% Get Miscanthus-"equivalent" LPJ-GUESS crops
miscanthus_equivalents = {'CerealsC4', 'TeCo'} ;
miscanthus_equiv_ri = {} ;
for q = 1:length(miscanthus_equivalents)
    this_equiv = miscanthus_equivalents{q} ;
    is_this_equiv = contains(yield_lpj.varNames, this_equiv) ;
    if ~any(is_this_equiv)
        continue
    end
    if sum(is_this_equiv) ~= 2
        error('Expected 0 or 2 yield_lpj.varNames containing %s; got %d', sum(is_this_equiv))
    end
    these_equivs = yield_lpj.varNames(is_this_equiv) ;
    
    % Ensure these_equivs are {rainfed, irrigated}
    these_equivs_irrCodes = strrep(these_equivs, this_equiv, '') ;
    isIrr = contains(these_equivs_irrCodes, 'i') ; % Fragile but works for now
    if sum(isIrr) ~= 1
        error('Expected 1 irrigated this_equiv; got %d', sum(isIrr))
    end
    [isIrr, I] = sort(isIrr) ;
    miscanthus_equiv_ri = these_equivs(I) ;

    warning('Using %s and %s as proxies for rainfed and irrigated Miscathus, respectively.', ...
        miscanthus_equiv_ri{1}, miscanthus_equiv_ri{2})
    break
end
if isempty(miscanthus_equiv_ri)
    miscanthus_equivalents_char = sprintf('%s ', miscanthus_equivalents{:}) ;
    error('No "equivalents" for Miscanthus found (checked: %s)', miscanthus_equivalents_char)
end
[~, ~, IB] = intersect(miscanthus_equiv_ri, yield_lpj.varNames, 'stable') ;
yield_lpj.maps_YXvy = yield_lpj.maps_YXvy(:,:,IB,:) ;
yield_lpj.varNames = yield_lpj.varNames(IB) ;

% Set up lat/lon maps
if (isfield(yield_lpj, 'lat_extent') && ~isequal(yield_lpj.lat_extent, [-90 90])) ...
        || (isfield(yield_lpj, 'lat_orient') && ~strcmp(yield_lpj.lat_orient, 'center')) ...
        || ~isequal(size(yield_lpj.maps_YXvy, 1:2), [360 720])
    warning('Edit import_betydb_miscanthus() to handle non-standard grids')
    yield_map_YXi = [] ;
    return
end
res = 0.5 ;
lats = (-(90-res/2):res:(90-res/2))' ;
lats_map = repmat(lats,[1 360/res]) ;
lons = (-(180-res/2):res:(180-res/2)) ;
lons_map = repmat(lons,[180/res 1]) ;
if ~isequal(size(lats_map), size(lons_map))
    error('Map size mismatch for lon/lat maps')
end
Nlats = length(lats) ;
Nlons = length(lons) ;
if ~isequal(size(lats_map), size(yield_lpj.maps_YXvy, 1:2))
    error('Map size mismatch between LPJ-GUESS and betyDB')
end

% Import betyDB table
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
betydb_data = readtable(fullfile(dir_data, 'other', 'betydb_grass_yields2.csv')) ;

% Include only rows with known irrigation in period of interest
betydb_data = betydb_data( ...
    betydb_data.year >= year1 ...
    & betydb_data.year <= yearN ...
    & ~strcmp(betydb_data.irrig, 'NA'),:) ;

% Extract some info
betydb_yields = betydb_data.mean_yield ;
betydb_irrBool = strcmp(betydb_data.irrig, 'TRUE') ;

% Aggregate Miscanthus yields to LPJ-GUESS gridcell means
yield_map_RF = nan(Nlats,Nlons) ;
yield_map_IR = nan(Nlats,Nlons) ;
for ii = 1:Nlons
    thisLon = lons(ii) ;
    inThisLon = (betydb_data.lon>=thisLon-res/2 & betydb_data.lon<thisLon+res/2) ;
    if any(inThisLon)
        for jj = 1:Nlats
            thisLat = lats(jj) ;
            inThisLat = (betydb_data.lat>=thisLat-res/2 & betydb_data.lat<thisLat+res/2) ;
            if any(inThisLon & inThisLat)
                if any(inThisLon & inThisLat & ~betydb_irrBool)
                    yield_map_RF(jj,ii) = mean(betydb_yields(inThisLon & inThisLat & ~betydb_irrBool)) ;
                end
                if any(inThisLon & inThisLat & betydb_irrBool)
                    yield_map_IR(jj,ii) = mean(betydb_yields(inThisLon & inThisLat & betydb_irrBool)) ;
                end
            end
        end
    end
end

rf_x_os = get_x_os(yield_lpj, isIrr==0, yield_map_RF) ;
ir_x_os = get_x_os(yield_lpj, isIrr==1, yield_map_IR) ;
miscanthus_x_os = cat(1, rf_x_os, ir_x_os) ;

disp('Done.')

end


function x_os = get_x_os(yield_lpj, inclIrr, obs_map)

sim_YX = mean(yield_lpj.maps_YXvy(:,:,inclIrr,:), 4) ;
isok = ~isnan(obs_map) & ~isnan(sim_YX) ;
x_os = [obs_map(isok) sim_YX(isok)] ;

end