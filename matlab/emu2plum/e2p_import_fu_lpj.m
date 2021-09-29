function data_fu_lpj = e2p_import_fu_lpj( ...
    baseline_yN, future_ts, future_yN, topDir_lpj, which_file, ...
    varargin)

gridlist_target = [] ;
if ~isempty(varargin)
    gridlist_target = varargin{1} ;
    if length(varargin) > 1
        error('e2p_import_fu_lpj() takes at most 1 optional argument (gridlist_target)')
    end
end

% Future: Set up
y1 = baseline_yN + 1 ;
ii = 0 ;
test_y1s = (baseline_yN+1):future_ts:future_yN ;
test_yNs = (baseline_yN+future_ts):future_ts:future_yN ;
test_Nyears = zeros(size(test_y1s)) ;
Ntpers = length(test_y1s) ;

% Future: Import
while y1 < future_yN
    ii = ii + 1 ;
    
    % Get this directory
    [thisDir, yN] = e2p_get_thisDir(topDir_lpj, y1, future_yN, which_file) ;
%     fprintf('%d-%d\n', y1, yN) ;
    thisDir_path = sprintf('%s/%s', thisDir.folder, thisDir.name) ;
    subdirs = dir([thisDir_path '/out*']) ;
    if ~isempty(subdirs)
        thisDir = sprintf('%s/%s', subdirs(end).folder, subdirs(end).name) ;
    else
        thisDir = sprintf('%s/%s', thisDir.folder, thisDir.name) ;
    end
    clear subdirs
    
    % Which time period are we dealing with?
    tper_i = 0 ;
    for t = 1:Ntpers
        if y1 >= test_y1s(t) && yN <= test_yNs(t)
            tper_i = t ;
            break
        end
    end
    if tper_i == 0
        error('Error finding tper_i')
    end
%     fprintf('tper_i %d\n', tper_i)
    
    filename = sprintf('%s/%s.out', thisDir, which_file) ;
    if isempty(gridlist_target)
        tmp = lpjgu_matlab_read2geoArray(filename, 'verboseIfNoMat', false) ;
    else
        tmp = lpjgu_matlab_read2geoArray(filename, 'verboseIfNoMat', false, ...
            'target', gridlist_target) ;
    end
    if ii==1
        data_fu_lpj.garr_xvt = zeros([size(tmp.garr_xv) Ntpers]) ;
        data_fu_lpj.lonlats = tmp.lonlats ;
        data_fu_lpj.list2map = tmp.list2map ;
        data_fu_lpj.varNames = tmp.varNames ;
    else
        if ~isequal(tmp.varNames, data_fu_lpj.varNames)
            error('Mismatch in variable names between time periods in future LPJ-GUESS outputs')
        end
    end
    
    % Weights for weighted average
    Nyears_existing = test_Nyears(tper_i) ;
    Nyears_new = yN - y1 + 1 ;
    Nyears_total = Nyears_existing + Nyears_new ;
    weight_existing = Nyears_existing / Nyears_total ;
    weight_new = Nyears_new / Nyears_total ;
    
    % Save
    data_fu_lpj.garr_xvt(:,:,tper_i) = ...
        data_fu_lpj.garr_xvt(:,:,tper_i) * weight_existing + ...
        tmp.garr_xv * weight_new ;
    
    clear tmp
    test_Nyears(tper_i) = test_Nyears(tper_i) + Nyears_new ;
    y1 = yN + 1 ;
end

% Check that each time step included 
if any(test_Nyears ~= future_ts)
    error('At least one time step did not include %d years'' worth of data', future_ts)
end

% Sanity check
if max(max(max(data_fu_lpj.garr_xvt))) == 0
    error('LPJ-GUESS %s data is all zeros on import', ...
        which_file)
end

% Trim unneeded variables
if any(strcmp({'yield', 'gsirrigation'}, which_file))
    data_fu_lpj = e2p_trim_unneeded(data_fu_lpj) ;
else
    if strcmp(which_file, 'anpp')
        cols2keep = {'C3G', 'C4G'} ;
    elseif strcmp(which_file, 'tot_runoff')
        cols2keep = {'Surf'} ;
    else
        error('which_file %s not recognized', which_file)
    end
    [~, IA] = intersect(data_fu_lpj.varNames, cols2keep) ;
    if length(IA) ~= length(cols2keep)
        error('Expected to find %d matches in cols2keep; found %d', ...
            length(cols2keep), length(IA))
    end
    data_fu_lpj.garr_xvt = data_fu_lpj.garr_xvt(:,IA,:) ;
    data_fu_lpj.varNames = data_fu_lpj.varNames(IA) ;
end

% Sanity check
if max(max(max(data_fu_lpj.garr_xvt))) == 0
    error('LPJ-GUESS %s data is all zeros after trimming unneeded variables', ...
        which_file)
end


end
