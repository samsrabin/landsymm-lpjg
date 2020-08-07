function data_bl_lpj = e2p_import_bl_lpj( ...
    baseline_y1, baseline_yN, topDir_lpj, which_file, get_unneeded, ...
    gridlist_target)

ii = 0 ;
y1 = baseline_y1 ;
yN = 0 ;
while yN < baseline_yN
    ii = ii + 1 ;
    
    % Get this directory
    [thisDir, yN] = e2p_get_thisDir(topDir_lpj, y1, baseline_yN) ;
    
    thisDir_path = sprintf('%s/%s', thisDir.folder, thisDir.name) ;
    subdirs = dir([thisDir_path '/out*']) ;
    if ~isempty(subdirs)
        thisDir = sprintf('%s/%s', subdirs(end).folder, subdirs(end).name) ;
    else
        thisDir = sprintf('%s/%s', thisDir.folder, thisDir.name) ;
    end
    clear subdirs
    
    % Import
    filename = sprintf('%s/%s.out', thisDir, which_file) ;
    if ii==1
        data_bl_lpj = lpjgu_matlab_read2geoArray(filename, 'verboseIfNoMat', false, ...
            'target', gridlist_target) ;
    else
        tmp = lpjgu_matlab_read2geoArray(filename, 'verboseIfNoMat', false, ...
            'target', gridlist_target) ;
        data_bl_lpj.garr_xv = data_bl_lpj.garr_xv * (ii-1)/ii + tmp.garr_xv * 1/ii ;
        clear tmp
    end
    
    % Prepare for next iteration
    y1 = yN + 1 ;
end

% Trim unneeded variables
data_bl_lpj = e2p_trim_unneeded(data_bl_lpj, get_unneeded) ;


end