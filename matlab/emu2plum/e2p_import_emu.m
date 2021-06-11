function [data_bl_emu, data_fu_emu] = e2p_import_emu( ...
    topDir_emu, thisGCM, thisEmu, thisSSP, which_file, ...
    cropList_in, gridlist, ts1_list, tsN_list, Nlist, ...
    baseline_yN, future_yN_emu, irrList_in, irrList_out, ...
	emuVer, adaptation)

% Translate which_file to Christoph's tokens
if strcmp(which_file, 'yield')
    which_file = 'yields' ;
elseif strcmp(which_file, 'gsirrigation')
    which_file = 'iwd' ;
else
    error('which_file %s not recognized', which_file)
end

Ncells = length(gridlist.list_to_map) ;
Ntpers = length(ts1_list) ;
NN = length(Nlist) ;
Ncrops_in = length(cropList_in) ;
yearList_in = (baseline_yN+1):future_yN_emu ;
Nyears = length(yearList_in) ;

% Get short crop name tokens
cropList_in_short = cropList_in ;
for c = 1:Ncrops_in
    thisCrop = cropList_in{c} ;
    thisCrop_short = e2p_get_thisCrop_short(thisCrop) ;
    cropList_in_short(strcmp(cropList_in_short, thisCrop)) = {thisCrop_short} ;
end

% Which crops does this emulator have for this GCM-SSP?
cropList_in_this_yield = get_cropList_in_this('yield', ...
    topDir_emu, emuVer, thisSSP, thisGCM, thisEmu, adaptation, ...
    cropList_in)  ;
cropList_in_this_iwd = get_cropList_in_this('iwd', ...
    topDir_emu, emuVer, thisSSP, thisGCM, thisEmu, adaptation, ...
    cropList_in)  ;
if ~isequal(sort(cropList_in_this_yield), sort(cropList_in_this_iwd))
    msg = 'yield and iwd files have different crop lists:' ;
    msg = sprintf('%s\n          yield   iwd', msg) ;
    cropList_in_this = unique([cropList_in_this_yield ; cropList_in_this_iwd]) ;
    for c = 1:length(cropList_in_this)
        thisCrop = cropList_in_this{c} ;
        thisCrop_short = e2p_get_thisCrop_short(thisCrop) ;
        msg = sprintf('%s\n    %s       %d     %d', ...
            msg, thisCrop_short, ...
            any(strcmp(cropList_in_this_yield, thisCrop)), ...
            any(strcmp(cropList_in_this_iwd, thisCrop))) ;
    end
    error(msg)
end
cropList_in_this = cropList_in_this_iwd ;
Ncrops_in_this = length(cropList_in_this) ;

% Set up input arrays
data_in_bl_xcni = nan(Ncells, Ncrops_in_this, NN, 2) ;
data_in_fu_xtcni = nan(Ncells, Ntpers, Ncrops_in_this, NN, 2) ;
if strcmp(which_file, 'iwd')
    data_in_bl_xcni(:,:,:,strcmp(irrList_in,'rf')) = 0 ;
    data_in_fu_xtcni(:,:,:,:,strcmp(irrList_in,'rf')) = 0 ;
elseif ~strcmp(which_file, 'yields')
    error('which_file %s not recognized', which_file)
end

% Import
% fprintf('Importing %s %s %s (%d of %d)...\n', ...
%     thisGCM, thisSSP, thisEmu, Ngse_worked, Ngse)
for c = 1:Ncrops_in_this
    
    % Get this crop's long and short names
    thisCrop = cropList_in_this{c} ;
    thisCrop_short = cropList_in_short(strcmp(cropList_in, thisCrop)) ;
    if length(thisCrop_short) ~= 1
        error('Error parsing short version of %s', thisCrop)
    end
    thisCrop_short = thisCrop_short{1} ;
    
    % Loop through N levels
    for n = 1:NN
        thisN = Nlist(n) ;
        
        % Get directories and file search patterns
        thisDir_bl = sprintf('%s/%s/CMIP6/A%d_N%d_%s/baseline/%s', ...
            topDir_emu, which_file, adaptation, thisN, emuVer, thisEmu) ;
        pattern_in_BL = sprintf('%s/*_%s_%s_*baseline*', ...
            thisDir_bl, thisCrop, thisEmu) ;
        thisDir_fu = sprintf('%s/%s/CMIP6/A%d_N%d_%s/%s/%s', ...
            topDir_emu, which_file, adaptation, thisN, emuVer, thisSSP, thisEmu) ;
        if ~exist(thisDir_fu, 'dir')
            error('thisDir_fu does not exist: %s', thisDir_fu)
        end
        if ~exist(thisDir_bl, 'dir')
            warning('baseline/ directory not found; looking in thisDir_fu (%s)', thisDir_fu)
            thisDir_bl = thisDir_fu ;
            pattern_in_BL = sprintf('%s/*_%s_*_%s_%s_*baseline*', ...
                thisDir_bl, thisGCM, thisCrop, thisEmu) ;
        end
        
        % Get file with data for baseline
        files = dir(pattern_in_BL) ;
        if length(files) ~= 1
            error('thisDir_bl ok, but error finding thisFile_BL: %d found', length(files))
        end
        thisFile_BL = sprintf('%s/%s', files.folder, files.name) ;
        
        % Get file with data for future
        pattern_in = sprintf('%s/*_%s_*_%s_%s_*movingwindow*', ...
            thisDir_fu, thisGCM, thisCrop, thisEmu) ;
        files = dir(pattern_in) ;
        if length(files) ~= 1
            error('thisDir_fu ok, but error finding thisFile: %d found', length(files))
        end
        thisFile = sprintf('%s/%s', files.folder, files.name) ;
        
%         % Print status
%         nread = nread + 1 ;
%         fprintf('    %d of %d (%s %s %s %s N%d)...\n', ...
%             nread, Nfiles/3, ...
%             thisGCM, thisSSP, thisEmu, thisCrop, thisN)
        
        % Read yield files
        if strcmp(which_file, 'yields')
            for ii = 1:2
                thisIrr = irrList_in{ii} ;
                thisVar = sprintf('yield_%s_%s', thisIrr, thisCrop_short) ;
                yieldBL_in_YX = flipud(transpose(ncread(thisFile_BL, thisVar))) ;
                data_in_bl_xcni(:,c,n,ii) = yieldBL_in_YX(gridlist.list_to_map) ;
                yield_in_YXy = flip(permute(ncread(thisFile, thisVar), [2 1 3]), 1) ;
                if size(yield_in_YXy,3) ~= Nyears
                    error('Expected %d years in yield file but found %d', ...
                        Nyears, size(yield_in_YXy,3))
                end
                data_in_fu_xtcni(:,:,c,n,ii) = ...
                    process_timesteps(ts1_list, tsN_list, ...
                    yearList_in, yield_in_YXy, gridlist);
                clear yieldBL_in_YX yield_in_YXy
            end
        elseif strcmp(which_file, 'iwd')
            ii = find(strcmp(irrList_in,'ir')) ;
            if length(ii) ~= 1
                error('Error finding ii: %d found', length(ii))
            end
            
            % Read IWD file
            thisVar = sprintf('iwd_%s', thisCrop_short) ;
            iwdBL_in_YX = flipud(transpose(ncread(thisFile_BL, thisVar))) ;
            data_in_bl_xcni(:,c,n,ii) = iwdBL_in_YX(gridlist.list_to_map) ;
            iwd_in_YXy = flip(permute(ncread(thisFile, thisVar), [2 1 3]), 1) ;
            if size(iwd_in_YXy,3) ~= Nyears
                error('Expected %d years in IWD file but found %d', ...
                    Nyears, size(iwd_in_YXy,3))
            end
            data_in_fu_xtcni(:,:,c,n,ii) = ...
                process_timesteps(ts1_list, tsN_list, ...
                yearList_in, iwd_in_YXy, gridlist);
            clear iwdBL_in_YX iwd_in_YXy
        else
            error('which_file %s not recognized', which_file)
        end
    end % Loop through N values
end % Loop through crops

% Get full names of variables
cropIrrList_in = [shiftdim(cropList_in_this) ; shiftdim(strcat(cropList_in_this, 'i'))] ;
cropIrrNlist_in = {} ;
for c = 1:length(cropIrrList_in)
    thisCropIrr = cropIrrList_in{c} ;
    for n = 1:NN
        thisN = pad(num2str(Nlist(n)), 3, 'left', '0') ;
        thisCropIrrN = [thisCropIrr thisN] ;
        cropIrrNlist_in{end+1} = thisCropIrrN ; %#ok<AGROW>
    end
end

% Combine crop*Nfert*irrig
clear thisCrop*
data_out_bl_xv = nan(Ncells, length(cropIrrNlist_in)) ;
data_out_fu_xtv = nan(Ncells, Ntpers, length(cropIrrNlist_in)) ;
for c = 1:Ncrops_in_this
    thisCrop_in = cropList_in_this{c} ;
    for ii = 1:2
        thisCropIrr = [thisCrop_in irrList_out{ii}] ;
        for n = 1:NN
            thisN = pad(num2str(Nlist(n)), 3, 'left', '0') ;
            thisCropIrrN = [thisCropIrr thisN] ;
            thisIndex = find(strcmp(cropIrrNlist_in, thisCropIrrN)) ;
            if length(thisIndex) ~= 1
                error('Error finding index of %s in cropIrrNlist_in', thisCropIrrN)
            end
            data_out_bl_xv(:,thisIndex) = data_in_bl_xcni(:,c,n,ii) ;
            data_out_fu_xtv(:,:,thisIndex) = data_in_fu_xtcni(:,:,c,n,ii) ;
        end
    end
end
            
% Create required output structs, converting tons/ha to kg/m2
tpha_to_kgpm2 = 0.1 ;

data_bl_emu.varNames = cropIrrNlist_in ;
data_bl_emu.list2map = gridlist.list_to_map ;
data_bl_emu.incl_years = e2p_get_incl_years(thisFile_BL) ;
data_bl_emu.garr_xv = data_out_bl_xv * tpha_to_kgpm2;

data_fu_emu.varNames = cropIrrNlist_in ;
data_fu_emu.list2map = gridlist.list_to_map ;
data_fu_emu.y1s = ts1_list ;
data_fu_emu.yNs = tsN_list ;
data_fu_emu.incl_years = e2p_get_incl_years(thisFile) ;
data_fu_emu.garr_xvt = permute(data_out_fu_xtv, [1 3 2]) * tpha_to_kgpm2;

% Fill any missing "future" years by weighting with baseline
t = 1 ;
while min(data_fu_emu.incl_years) > ts1_list(t) && t <= Ntpers
    yearList_pd = ts1_list(t):tsN_list(t) ;
    Nyears_pd = length(yearList_pd) ;
    missing_years = setdiff(yearList_pd, data_fu_emu.incl_years) ;
    if isempty(missing_years)
        error('How are no years missing here?')
    end
    Nmissing = length(missing_years) ;
    Npresent = Nyears_pd - Nmissing ;
    if min(data_bl_emu.incl_years) > missing_years(1) || max(data_bl_emu.incl_years) < missing_years(end)
        error('Unable to fill missing "future" years %d-%d with baseline because baseline range is %d-%d', ...
            missing_years(1), missing_years(end), data_bl_emu.incl_years(1), data_bl_emu.incl_years(end))
    end
    wt_fu = Npresent / Nyears_pd ;
    wt_bl = 1 - wt_fu ;
    warning('Filling missing "future" years %d-%d (%d/%d in period %d-%d) with emulated baseline mean', ...
        missing_years(1), missing_years(end), ...
        Nmissing, Nyears_pd, ts1_list(t), tsN_list(t))
    data_fu_emu.garr_xvt(:,:,t) = ...
        data_fu_emu.garr_xvt(:,:,t) * wt_fu ...
        + data_bl_emu.garr_xv * wt_bl ;
    t = t+1 ;
end

% Make sure emulated future range covers the latest period
if max(data_fu_emu.incl_years) < tsN_list(end)
    error('Last period ends %d but future emulation only goes to %d', ...
        tsN_list(end), max(data_fu_emu.incl_years))
end


end


function out_xt = process_timesteps(ts1_list, tsN_list, yearList_in, ...
    in_YXy, gridlist)

Ntpers = length(ts1_list) ;
out_xt = nan(length(gridlist.list_to_map), Ntpers) ;

for t = 1:Ntpers
    thisY1 = ts1_list(t) ;
    thisYN = tsN_list(t) ;
    okyrs = yearList_in>=thisY1 & yearList_in<=thisYN ;
    tmp_YX = mean(in_YXy(:,:,okyrs),3) ;
    out_xt(:,t) = tmp_YX(gridlist.list_to_map) ;
end

end


function cropList_in_this = get_cropList_in_this(iwd_or_yield, ...
    topDir_emu, emuVer, thisSSP, thisGCM, thisEmu, adaptation, ...
    cropList_in)

pattern = sprintf('%s/**/*%s_%s*%s_A%d*%s_baseline*_%s.nc4', ...
    topDir_emu, thisSSP, thisGCM, thisEmu, adaptation, iwd_or_yield, emuVer) ;
cropList_in_this = dir(pattern) ;
if isempty(cropList_in_this)
    msgStruct.messsage = ['No files found matching ' pattern] ;
    msgStruct.identifier = 'e2p:e2p_import_emu:noFilesFound' ;
%     error(msgStruct)
    error(msgStruct.identifier, msgStruct.messsage)
    stop
end
cropList_in_this = {cropList_in_this.name}' ;
cropList_in_this = regexprep(cropList_in_this, ...
    ['cmip6_' thisSSP '_' thisGCM '_r\di\dp\d_'], '') ;
expression = sprintf('_%s_A%d_(N%s_)?emulated_%s_baseline_1980_2010_average_%s.nc4', ...
    thisEmu, adaptation, '\d+', iwd_or_yield, emuVer) ;
cropList_in_this = unique(regexprep(cropList_in_this, ...
    expression, '')) ;
if any(~startsWith(cropList_in_this, cropList_in) | ~endsWith(cropList_in_this, cropList_in))
    cropList_in_this %#ok<NOPRT>
    error('Error parsing crops present in %s for %s %s', ...
        thisEmu, thisGCM, thisSSP)
end


end
