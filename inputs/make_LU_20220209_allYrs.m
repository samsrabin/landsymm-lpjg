%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate land use files for yield table runs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use historical land use (with tiny areas of everything) until y0_past, at
% which point everything is converted to pasture (again, aside from tiny
% refugia of other LUs). After one year, the gridcell's LU is split up
% evenly among all stands.

list_lu = {'NATURAL', 'CROPLAND', 'PASTURE', 'FOREST'} ;
list_forest = {'forC'} ;

% When is the year for which we want the first outputs from the experimental stands?
y1_expt_list = 1850:20:2090 ;

% Which land use versions are we using?
remapVer = '10_g2p' ; % CerealsC3, CerealsC4, Rice, OilNfix, OilOther, 
                      % StarchyRoots, Pulses, Sugarbeet, Sugarcane, 
                      % FruitAndVeg

% No historical land use?
use_historical_lu = true ;

dir_out = sprintf('/Users/Shared/LandSyMM/inputs/LU_20220209/%s', remapVer) ;

% Output behaviors/settings
outPrec = 6 ;
overwrite = true ;
outWidth = 1 ;
delimiter = ' ' ;
fancy = false ;


%% Setup

% Create output directory, if needed
if ~exist(dir_out, 'dir')
    mkdir(dir_out)
end


%% Get ISIMIP3 historical land use
disp('Reading input files...')

verbose = false ;

% Determine files to use
if contains(remapVer, '_g2p')
    LUdir = sprintf('/Users/Shared/G2P/inputs/LU/remaps_v%s', ...
        remapVer) ;
else
    LUdir = sprintf('/Users/Shared/PLUM/input/remaps_v%s', ...
        remapVer) ;
end
file_lu = sprintf('%s/LU.remapv%s.txt', ...
    LUdir, remapVer);
file_cf = sprintf('%s/cropfracs.remapv%s.txt', ...
    LUdir, remapVer);

% Read input files
lu_in = lpjgu_matlab_read2geoArray(file_lu) ;
Ncells = length(lu_in.list2map) ;
Nyears = length(lu_in.yearList) ;
if use_historical_lu
    cf_in = lpjgu_matlab_read2geoArray(file_cf) ;
    
    % Check
    if ~isequal(cf_in.lonlats, lu_in.lonlats)
        error('Need to align gridlists of land use and cropfrac files')
    end
    if isfield(cf_in, 'yearList')
        error('This code relies on crop fractions not varying over time')
    end
else
    lu_in = unneeded_LU_to_ntrl(lu_in, list_lu, NaN) ;
    lu_in = add_zeros_for_missing_LUs(list_lu, lu_in) ;
    % Make Ncells*Nlu*0 (empty) array that will be added to in next step
    lu_in.garr_xvy = nan(Ncells, length(list_lu), 0) ;
    lu_in.yearList = [] ;
    Nyears = 0 ;
end

disp('Done.')


%% Set up output structures
disp('Setting up output structures...')

if use_historical_lu
    lu_out = lu_in ;
    cf_out = cf_in ;

    % Ensure sum of LU fractions to 1
    lu_sums_x1y = sum(lu_out.garr_xvy,2) ;
    if any(any(abs(lu_sums_x1y - 1) > 10^-outPrec))
        lu_out.garr_xvy = lu_out.garr_xvy ...
            ./ repmat(lu_sums_x1y, [1 size(lu_out.garr_xvy, 2) 1]) ;
        sum_to_1_test = sum(lu_out.garr_xvy,2) - 1 ;
        if any(any(abs(sum_to_1_test) > 10^-outPrec))
            error('Land use fractions do not sum to 1 (max abs. deviation %g) even after trying to fix', max(max(abs(sum_to_1_test))))
        end
    end
    
    % Move any unneeded LU types into NATURAL
    lu_out = unneeded_LU_to_ntrl(lu_out, list_lu, outPrec) ;

    % Trim any unneeded crop types
    is_unneeded = nansum(cf_out.garr_xv, 1) == 0 ;
    if any(is_unneeded)
        cf_out.garr_xv(:,is_unneeded) = [] ;
        cf_out.varNames(is_unneeded) = [] ;
    end
    
    % Ensure that every cell has at least some CerealsC3
    minCrop = 10^-outPrec ;
    noCrop_x1y = lu_out.garr_xvy(:,strcmp(lu_out.varNames, 'CROPLAND'),:) < minCrop ;
    never_cropped = all(noCrop_x1y, 3) ;
    cf_out.garr_xv(never_cropped,:) = 0 ;
    i_CerealsC3 = find(strcmp(cf_out.varNames, 'CerealsC3')) ;
    cf_out.garr_xv(never_cropped, i_CerealsC3) = 1 ;
    if any(noCrop_x1y(:))
        if any(any(abs(sum_to_1_test) > 10^-outPrec))
            error('Land use fractions do not sum to 1 (max abs. deviation %g) before ensuring that every cell has at least some crop', max(max(abs(sum_to_1_test))))
        end
        lu_out = donate_land_to_crop(lu_out, noCrop_x1y, 'PASTURE', minCrop, verbose) ;
        noCrop_x1y = lu_out.garr_xvy(:,strcmp(lu_out.varNames, 'CROPLAND'),:) < minCrop ;
        if any(noCrop_x1y(:))
            lu_out = donate_land_to_crop(lu_out, noCrop_x1y, 'NATURAL', minCrop, verbose) ;
            noCrop_x1y = lu_out.garr_xvy(:,strcmp(lu_out.varNames, 'CROPLAND'),:) < minCrop ;
            if any(noCrop_x1y(:))
                error('There are still %d cell-years with no cropland', length(find(noCrop_x1y)))
            end
        end
    end
    sum_to_1_test = sum(lu_out.garr_xvy,2) - 1 ;
    if any(any(abs(sum_to_1_test) > 10^-outPrec))
        error('Land use fractions do not sum to 1 (max abs. deviation %g)', max(max(abs(sum_to_1_test))))
    end
    needsCerealsC3 = cf_out.garr_xv(:, strcmp(cf_out.varNames, 'CerealsC3')) == 0 ;
    if any(needsCerealsC3)
        maxCropArea = max(cf_out.garr_xv, [], 2) ;
        for c = 1:length(cf_out.varNames)
            thisCrop = cf_out.varNames{c} ;
            if strcmp(thisCrop, 'CerealsC3')
                continue
            end
            thisIsMax = cf_out.garr_xv(:,c) == maxCropArea ;
            if ~any(thisIsMax & needsCerealsC3)
                continue
            end
            cf_out.garr_xv(thisIsMax & needsCerealsC3, i_CerealsC3) = minCrop ;
            cf_out.garr_xv(thisIsMax & needsCerealsC3, c) = ...
                cf_out.garr_xv(thisIsMax & needsCerealsC3, c) - minCrop ;

            needsCerealsC3 = cf_out.garr_xv(:, strcmp(cf_out.varNames, 'CerealsC3')) == 0 ;
            if ~any(needsCerealsC3)
                break
            end
        end
    end
    if any(needsCerealsC3)
        error('Some cell(s) still missing CerealsC3?')
    end

    
    % Add zeros for missing LUs
    lu_out = add_zeros_for_missing_LUs(list_lu, lu_out) ;

end

disp('Done.')


%% Process

y1_expt_list_desc = sort(y1_expt_list,'descend') ;
for y1_expt = y1_expt_list_desc

    y0_expt = y1_expt - 1 ; % Because establishment happens at END of first year
    y1_past = y0_expt - 1 ;

    fprintf('Processing %dpast... ', y1_past)

    if use_historical_lu
        % Trim any unneeded years
        lu_out.garr_xvy(:,:,lu_out.yearList >= y1_past) = [] ;
        lu_out.yearList(lu_out.yearList >= y1_past) = [] ;
        Nyears = length(lu_out.yearList) ;
        if length(lu_out.yearList) > length(unique(lu_out.yearList))
            error('Non-unique value(s) in lu_out.yearList')
        end
    else
        lu_out = lu_in ;
    end

    % Add extra year at beginning, if needed
    if Nyears == 0
        if ~isempty(lu_out.garr_xvy)
            error('Nyears==0 but lu_out.garr_xvy not empty somehow')
        elseif ~isempty(lu_out.yearList)
            error('Nyears==0 but lu_out.yearList not empty somehow')
        end
        tmp_xv = zeros(size(lu_out.garr_xvy, 1:2)) ;
        tmp_xv(:,strcmp(lu_out.varNames, 'NATURAL')) = 1 ;
        lu_out.garr_xvy = cat(3, lu_out.garr_xvy, tmp_xv) ;
        clear tmp_xv
        lu_out.yearList = y1_past - 1 ;
        Nyears = length(lu_out.yearList) ;
        if length(lu_out.yearList) > length(unique(lu_out.yearList))
            error('Non-unique value(s) in lu_out.yearList after adding to beginning')
        end
    end

    % Add extra year, if needed
    if max(lu_out.yearList) < y1_past-1
        lu_out.garr_xvy = cat(3, lu_out.garr_xvy, lu_out.garr_xvy(:,:,end)) ;
        lu_out.yearList(Nyears+1) = y1_past - 1 ;
        Nyears = Nyears + 1 ;
        if length(lu_out.yearList) ~= size(lu_out.garr_xvy, 3); error('Nyear mismatch'); end
        if length(lu_out.yearList) > length(unique(lu_out.yearList))
            error('Non-unique value(s) in lu_out.yearList after adding to end')
        end
    end
    
    % Add 1 year of pasture, then all LUs
    if Nyears > 0
        lu_out.yearList(end+1) = y1_past ;
    end
    Nyears = length(lu_out) ;
    lu_out.garr_xvy(:,:,end+1) = 0 ;
    lu_out.garr_xvy(:,strcmp(lu_out.varNames,'PASTURE'),end) = 1 ;
    Nlu_out = length(list_lu) ;
    Nforest = length(list_forest) ;
    thisDenom = Nlu_out - 1 + Nforest ;
    lu_out.yearList(end+1) = y0_expt ;
    lu_out.garr_xvy(:,:,end+1) = 1 / thisDenom ;
    lu_out.garr_xvy(:,strcmp(lu_out.varNames,'FOREST'),end) = ...
        Nforest / thisDenom ;

    % Sanity checks
    if ~isequal(sort(list_lu), sort(lu_out.varNames))
        error('Mismatch between list_lu and lu_out.varNames')
    end
    sum_to_1_test = sum(lu_out.garr_xvy,2) - 1 ;
    if any(any(abs(sum_to_1_test) > 1e-6))
        error('Land use fractions do not sum to 1 (max abs. deviation %g)', max(max(abs(sum_to_1_test))))
    end
    if length(lu_out.yearList) > length(unique(lu_out.yearList))
        error('Non-unique value(s) in lu_out.yearList')
    end

    % Save cropfrac if this is the first time through the loop
    if exist('cf_out', 'var') && y1_expt == y1_expt_list_desc(1)

        cf_out_header = [{'Lon', 'Lat'}, cf_out.varNames] ;
        file_cf_out = sprintf('%s/remap%s_someCrop_cropfrac.txt', ...
            dir_out, strrep(remapVer, '_', '')) ;
        
        fprintf('Saving cropfrac... ')
        lpjgu_matlab_saveTable(cf_out_header, cf_out, file_cf_out, ...
            'outPrec', outPrec, ...
            'outWidth', outWidth, ...
            'delimiter', delimiter, ...
            'overwrite', overwrite, ...
            'fancy', fancy, ...
            'verbose', verbose) ;

    end % if first year in loop (saving cropfrac)
    
    % Save land use
    if use_historical_lu
        basename_prefix = sprintf('remap%s_someCrop_', strrep(remapVer, '_', '')) ;
    else
        basename_prefix = '../' ;
    end
    file_lu_out = sprintf('%s/%s%dpast_%dall_LU.txt', ...
        dir_out, basename_prefix, y1_past, y0_expt) ;
    lu_out_header = [{'Lon', 'Lat', 'Year'}, lu_out.varNames] ;
    fprintf('Saving LU... ')
    lpjgu_matlab_saveTable(lu_out_header, lu_out, file_lu_out, ...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'verbose', verbose) ;
    disp('Done.')

end % loop through years

disp('All done!')


%% FUNCTIONS

function S = unneeded_LU_to_ntrl(S, list_lu, outPrec)

[~, I] = setdiff(S.varNames, list_lu) ;
if ~isempty(I)
    sum_to_1_test = sum(S.garr_xvy,2) - 1 ;
    if ~isnan(outPrec) && any(any(abs(sum_to_1_test) > 10^-outPrec))
        error('Land use fractions do not sum to 1 (max abs. deviation %g) before moving unneeded LU types into NATURAL', max(max(abs(sum_to_1_test))))
    end
    S.garr_xvy(:,strcmp(S.varNames, 'NATURAL'),:) = ...
        S.garr_xvy(:,strcmp(S.varNames, 'NATURAL'),:) ...
        + sum(S.garr_xvy(:,I,:), 2) ;
    S.garr_xvy(:,I,:) = [] ;
    S.varNames(I) = [] ;
end

end


function S = add_zeros_for_missing_LUs(list_lu, S)

Ncells = length(S.list2map) ;
Nyears = length(S.yearList) ;

missing_LUs = setdiff(list_lu, S.varNames) ;
if ~isempty(missing_LUs)
    if ~isempty(S.garr_xvy)
        S.garr_xvy = cat(2, S.garr_xvy, ...
            zeros(Ncells, length(missing_LUs), Nyears)) ;
    end
    S.varNames = [S.varNames missing_LUs] ;
end

end

function lu = donate_land_to_crop(lu, noCrop_x1y, donor, minCrop, verbose)

% Sanity checks
if any(any(any(lu.garr_xvy < 0)))
    error('Negative value(s) in lu.garr_xvy (worst: %g)', min(min(min(lu.garr_xvy))))
elseif any(any(any(isnan(lu.garr_xvy))))
    error('NaN(s) in lu.garr_xvy')
end

I_crop = find(strcmp(lu.varNames, 'CROPLAND')) ;
if length(I_crop) ~= 1
    error('Expected 1 I_crop; found %d', length(I_crop))
end

I_donor = find(strcmp(lu.varNames, donor)) ;
if length(I_donor) ~= 1
    error('Expected 1 I_donor; found %d', length(I_donor))
end

lu_xCy = lu.garr_xvy(:,I_crop,:) ;
lu_xDy = lu.garr_xvy(:,I_donor,:) ;

okDonor_x1y = lu_xDy >= 2*minCrop ;
noCrop_okDonor_x1y = noCrop_x1y & okDonor_x1y ;

Nbad_before = length(find(noCrop_x1y)) ;

if any(noCrop_okDonor_x1y(:))

    lu_xDy(noCrop_okDonor_x1y) = ...
        lu_xDy(noCrop_okDonor_x1y) - minCrop ;
    lu_xCy(noCrop_okDonor_x1y) = ...
        lu_xCy(noCrop_okDonor_x1y) + minCrop ;
    
    Nbad_after = length(find(lu_xCy < minCrop)) ;
    
    msg = sprintf('# of no-crop cell-years: %d -> %d', ...
        Nbad_before, ...
        Nbad_after) ;
    
    if Nbad_after > Nbad_before
        error(msg) %#ok<SPERR>
    end
    
    if verbose
        fprintf('%s\n', msg)
    end
    
    lu.garr_xvy(:,I_donor,:) = lu_xDy ;
    lu.garr_xvy(:,I_crop,:) = lu_xCy ;

end

% Sanity checks
if any(any(any(lu.garr_xvy < 0)))
    error('Negative value(s) in lu.garr_xvy (worst: %g)', min(min(min(lu.garr_xvy))))
elseif any(any(any(isnan(lu.garr_xvy))))
    error('NaN(s) in lu.garr_xvy')
end

end


function S = import_lu(y0_past, filename)

% Import
S = lpjgu_matlab_read2geoArray(filename) ;

% Restrict to needed years
if isfield(S, 'yearList')
    S.garr_xvy(:,:,S.yearList >= y0_past) = [] ;
    S.yearList(S.yearList >= y0_past) = [] ;
end

end

