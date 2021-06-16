function cf = e2p_get_CFs(cropList, ggcm, cfDir, verbose)

cf = [] ;

% Make sure the directory exists
if ~exist(cfDir, 'dir')
    warning('Calibration factor directory (%s) not found; skipping figure.', ...
        cfDir)
    return
end

% Find the latest file for this GGCM
thePattern = sprintf('%s/*%s*.csv', cfDir, ggcm) ;
cf_files = dir(thePattern) ;
if isempty(cf_files)
    warning('No calibration factor CSV found matching %s; skipping figure.', ...
        thePattern)
    return
end
[~, I] = sort([cf_files.datenum]) ;
cf_file = sprintf('%s/%s', ...
    cf_files(I(end)).folder, cf_files(I(end)).name) ;

% Read and sort
T = readtable(cf_file) ;
[C, ~, IB] = intersect(cropList, table2cell(T(:,1)), 'stable') ;
if ~isequal(shiftdim(C), shiftdim(cropList))
    warning('Calibration factor CSV %s does not contain all values in cropList. Skipping figure.', ...
        cf_file)
    return
end
T = T(IB,:) ;

if verbose
    fprintf('%s calibration factors from\n', ggcm)
    disp(cf_file)
    disp(T)
end

cf = table2array(T(:,2)) ;


    
end