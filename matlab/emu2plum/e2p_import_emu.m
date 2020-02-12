function [data_bl_emu, data_fu_emu, Ntpers] = e2p_import_emu( ...
    topDir_emu_bl, topDir_emu_fu, gcm, ggcm, rcp, which_file, ...
    baseline_yN, future_yN)

% Baseline
dir_emu = sprintf('%s/%s', topDir_emu_bl, ggcm) ;
if ~exist(dir_emu,'dir')
    error('Emulator directory not found (%s)', dir_emu) ;
end
filename = sprintf('%s/%s.out', dir_emu, which_file) ;
data_bl_emu = lpjgu_matlab_read2geoArray(filename, 'verboseIfNoMat', false) ;

% Future: Get number of time periods
dir_emu = sprintf('%s/%s/%s/%s', topDir_emu_fu, gcm, ggcm, rcp) ;
if ~exist(dir_emu,'dir')
    error('Emulator directory not found (%s)', dir_emu) ;
end
y1 = baseline_yN + 1 ;
yN = 0 ;
Ntpers = 0 ;
while yN < future_yN
    Ntpers = Ntpers + 1 ;
    [~, yN] = e2p_get_thisDir(dir_emu, y1, future_yN) ;
    y1 = yN + 1 ;
end

% Future: Import
y1 = baseline_yN + 1 ;
yN = 0 ;
ii = 0 ;
while yN < future_yN
    ii = ii + 1 ;
    [thisDir, yN] = e2p_get_thisDir(dir_emu, y1, future_yN) ;
    
    filename = sprintf('%s/%s/%s.out', thisDir.folder, thisDir.name, which_file) ;
    if ii==1
        data_fu_emu = lpjgu_matlab_read2geoArray(filename, 'verboseIfNoMat', true) ;
        data_fu_emu.y1s = nan(Ntpers,1) ;
        data_fu_emu.yNs = nan(Ntpers,1) ;
        if ~isequal(data_bl_emu.varNames,data_fu_emu.varNames)
            error('Add code to deal with mismatch of baseline vs. future varNames')
        end
        data_fu_emu.garr_xvt = nan([size(data_fu_emu.garr_xv) Ntpers]) ;
        data_fu_emu.garr_xvt(:,:,1) = data_fu_emu.garr_xv ;
        data_fu_emu = rmfield(data_fu_emu, 'garr_xv') ;
    else
        tmp = lpjgu_matlab_read2geoArray(filename, 'verboseIfNoMat', false) ;
        if ~isequal(tmp.varNames, data_fu_emu.varNames)
            error('Mismatch in variable names between time periods in future emulator outputs')
        end
        data_fu_emu.garr_xvt(:,:,ii) = tmp.garr_xv ;
        clear tmp
    end
    data_fu_emu.y1s(ii) = y1 ;
    data_fu_emu.yNs(ii) = yN ;
    
    y1 = yN + 1 ;
end


end