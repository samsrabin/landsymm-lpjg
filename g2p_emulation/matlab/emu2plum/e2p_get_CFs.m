function T = e2p_get_CFs(cropList, ggcm, cfDir, combineCrops, verbose)

T = table() ;

failure_is_not_an_option = ~isempty(combineCrops) ;

% Make sure the directory exists
if ~exist(cfDir, 'dir')
    if failure_is_not_an_option
        error('Calibration factor directory (%s) not found.', ...
            cfDir)
    else
        warning('Calibration factor directory (%s) not found; skipping figure.', ...
            cfDir)
        return
    end
end

% Find the latest file for this GGCM
thePattern = sprintf('%s/*%s*.csv', cfDir, ggcm) ;
cf_files = dir(thePattern) ;
if isempty(cf_files)
    if failure_is_not_an_option
        error('No calibration factor CSV found matching %s.', ...
            thePattern)
    else
        warning('No calibration factor CSV found matching %s; skipping figure.', ...
            thePattern)
        return
    end
end
[~, I] = sort([cf_files.datenum]) ;
cf_file = sprintf('%s/%s', ...
    cf_files(I(end)).folder, cf_files(I(end)).name) ;

% Read and sort
T = readtable(cf_file) ;
[C, ~, IB] = intersect(cropList, table2cell(T(:,1)), 'stable') ;
if ~isequal(shiftdim(C), shiftdim(cropList))
    if failure_is_not_an_option
        error('Calibration factor CSV %s does not contain all values in cropList.', ...
            cf_file)
    else
        warning('Calibration factor CSV %s does not contain all values in cropList. Skipping figure.', ...
            cf_file)
        return
    end
end
T = T(IB,:) ;

if verbose
    fprintf('%s calibration factors from\n', ggcm)
    disp(cf_file)
    disp(T)
end

    
end