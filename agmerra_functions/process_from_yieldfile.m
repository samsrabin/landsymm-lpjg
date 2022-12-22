function process_from_yieldfile(file_yield, runInfo, dir_out)

% Setup
irrList_in = {'', 'i'} ;
irrList_out = {'_noirr', '_firr'} ;
Nirr = length(irrList_in) ;
latgrid = 89.75:-0.5:-89.75;
longrid = -179.75:0.5:179.75;

if ~exist(dir_out, 'dir')
    mkdir(dir_out)
end

% Input checks
if length(irrList_out) ~= Nirr
    error('Length mismatch between input vs. output irrigation lists')
end

% Import yield
disp('Importing yield...')
yield = lpjgu_matlab_read2geoArray(file_yield, ...
    'verboseIfNoMat', false) ;

% Convert to t/ha
yield.garr_xvy = 10 * yield.garr_xvy ;
theseUnits = 't ha-1 yr-1' ;

% Get gridlist
gridlist = struct() ;
gridlist.list_to_map = yield.list2map ;
gridlist.lonlats = yield.lonlats ;
gridlist.mask_YX = false(360, 720) ;
gridlist.mask_YX(gridlist.list_to_map) = true ;

% Trim to yearList_out
[C, ~, IB] = intersect(runInfo.yearList_out, yield.yearList, 'stable') ;
% Append blank year, if needed
if isequal(runInfo.yearList_out, 1980:2010) ...
        && isequal(shiftdim(C), shiftdim(1980:2009))
    yield.yearList = 1980:2010 ;
    yield.garr_xvy(:,:,end+1) = zeros(size(yield.garr_xvy(:,:,end))) ;
elseif ~isequal(shiftdim(C), shiftdim(runInfo.yearList_out))
    error('Not all years in runInfo.yearList_out present in yield file (%d-%d)', ...
        min(yield.yearList), max(yield.yearList))
else
    yield.garr_xvy = yield.garr_xvy(:,:,IB) ;
end

% Get max wheats, if needed
if length(intersect(yield.varNames, ...
        strcat({'CerealsC3s', 'CerealsC3w'}, irrList_in{1}))) == 2
    for ii = 1:Nirr
        thisIrr_in = irrList_in{ii} ;
        thisCropIrr_new = ['CerealsC3' thisIrr_in] ;
        inCrops = strcat({'CerealsC3s', 'CerealsC3w'}, thisIrr_in) ;
        [~, IA] = intersect(yield.varNames, inCrops) ;
        if length(IA) ~= 2
            error('??')
        end
        maxYield_g1y = max(yield.garr_xvy(:,IA,:), [], 2) ;
        yield.garr_xvy(:,IA,:) = [] ;
        yield.varNames(IA) = [] ;
        yield.garr_xvy(:,end+1,:) = maxYield_g1y ;
        yield.varNames{end+1} = thisCropIrr_new ;
    end
end

% Get crop lists
cropIrrList_in = yield.varNames ;
cropIrrList_in(strcmp(cropIrrList_in, 'ExtraCrop')) = [] ;
cropList_in = cropIrrList_in ;
for ii = 1:Nirr
    thisIrr_in = irrList_in{ii} ;
    if isempty(thisIrr_in)
        continue
    end
    % Delete this irrigation indicator
    pattern = [thisIrr_in '$'] ;
    isThisIrr = ~cellfun(@isempty, regexp(cropList_in, pattern)) ;
    cropList_in(isThisIrr) = regexprep(cropList_in(isThisIrr), pattern, '') ;
end
cropList_in = unique(cropList_in, 'stable') ;
Ncrops = length(cropList_in) ;

for c = 1:Ncrops
    thisCrop_in = cropList_in{c} ;
    
    % Get output crop
    switch thisCrop_in
        case {'CerealsC3'}
            thisCrop_out = 'whe' ;
        case {'CerealsC4'}
            thisCrop_out = 'mai' ;
        case {'Oilcrops'}
            thisCrop_out = 'soy' ;
        case {'Rice'}
            thisCrop_out = 'ric' ;
        otherwise
            error('thisCrop_in %s not recognized', thisCrop_in)
    end
    
    % Get output variable
    thisVar_out = sprintf('yield_%s', thisCrop_out) ;
    
    for ii = 1:Nirr
        thisIrr_in = irrList_in{ii} ;
        thisCropIrr_in = sprintf('%s%s', ...
            thisCrop_in, thisIrr_in) ;
        fprintf('%s...\n', thisCropIrr_in)
        
        % Get yield of this crop-irr
        v = find(strcmp(yield.varNames, thisCropIrr_in)) ;
        if length(v) ~= 1
            error('Expected to find 1 match of %s in yield.varNames; found %d', ...
                thisCropIrr_in, length(v))
        end
        data_out_xy = squeeze(yield.garr_xvy(:,v,:)) ;
        
        % Get output filename
        thisIrr_out = irrList_out{ii} ;
        file_out = sprintf('%s/%s_%s_hist_fullharm_%s_%s_%s_annual_%d_%d.nc', ...
            dir_out, ...
            lower(runInfo.modelname), runInfo.climate_forcing, ...
            strrep(thisIrr_out, '_', ''), 'yield', thisCrop_out, ...
            min(runInfo.yearList_out), max(runInfo.yearList_out)) ;
        
        % Save
        save_file_by10s(file_out, runInfo, latgrid, longrid, ...
            'single', thisVar_out, theseUnits, ...
            data_out_xy, gridlist, true)
        
    end
end

disp('Done')

end