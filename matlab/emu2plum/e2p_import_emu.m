function [data_bl_emu, data_fu_emu] = e2p_import_emu( ...
    topDir_emu, thisGCM, thisEmu, thisSSP, which_file, ...
    cropList_in, gridlist, ts1_list, tsN_list, Nlist, ...
    baseline_yN, future_yN_emu, irrList_in, irrList_out)

% Translate which_file to Christoph's tokens
if strcmp(which_file, 'yield')
    which_file = 'yields' ;
elseif strcmp(which_file, 'gsirrigation')
    which_file = 'iwd' ;
else
    error('which_file %s not recognized', which_file)
end

Ncells = length(gridlist.list_to_map) ;
Nts = length(ts1_list) ;
NN = length(Nlist) ;
Ncrops_in = length(cropList_in) ;
yearList_in = (baseline_yN+1):future_yN_emu ;
Nyears = length(yearList_in) ;

% Get short crop name tokens
cropList_in_short = cropList_in ;
for c = 1:Ncrops_in
    thisCrop = cropList_in{c} ;
    switch thisCrop
        case 'spring_wheat'
            thisCrop_short = 'swh' ;
        case 'winter_wheat'
            thisCrop_short = 'wwh' ;
        case 'maize'
            thisCrop_short = 'mai' ;
        case 'soy'
            thisCrop_short = 'soy' ;
        case 'rice'
            thisCrop_short = 'ric' ;
        otherwise
            error('Can''t parse %s into shortname', thisCrop)
    end
    cropList_in_short(strcmp(cropList_in_short, thisCrop)) = {thisCrop_short} ;
end

% Which crops does this emulator have for this GCM-SSP?
cropList_in_this = dir(sprintf('%s/**/*%s_%s*%s*iwd_baseline*.nc4', ...
    topDir_emu, thisSSP, thisGCM, thisEmu)) ;
if isempty(cropList_in_this)
    msgStruct.messsage = sprintf( ...
        'No files found matching %s/**/*%s_%s*%s*iwd_baseline*.nc4', ...
        topDir_emu, thisSSP, thisGCM, thisEmu) ;
    msgStruct.identifier = 'e2p:e2p_import_emu:noFilesFound' ;
%     error(msgStruct)
    error(msgStruct.identifier, msgStruct.messsage)
    stop
end
cropList_in_this = {cropList_in_this.name}' ;
cropList_in_this = regexprep(cropList_in_this, ...
    ['cmip6_' thisSSP '_' thisGCM '_r\di\dp\d_'], '') ;
cropList_in_this = unique(regexprep(cropList_in_this, ...
    ['_' thisEmu '_A0_N\d+_emulated_iwd_baseline_1980_2010_average_v2.0.nc4'],'')) ;
if any(~startsWith(cropList_in_this, cropList_in) | ~endsWith(cropList_in_this, cropList_in))
    cropList_in_this %#ok<NOPRT>
    error('Error parsing crops present in %s for %s %s', ...
        thisEmu, thisGCM, thisSSP)
end
Ncrops_in_this = length(cropList_in_this) ;

% Set up input arrays
data_in_bl_xcni = nan(Ncells, Ncrops_in_this, NN, 2) ;
data_in_fu_xtcni = nan(Ncells, Nts, Ncrops_in_this, NN, 2) ;
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
        
        % Get file with data for baseline
        thisDir = sprintf('%s/%s/CMIP6/A0_N%d/%s/%s', ...
            topDir_emu, which_file, thisN, thisSSP, thisEmu) ;
        if ~exist(thisDir, 'dir')
            error('thisDir does not exist: %s', thisDir)
        end
        pattern_in_BL = sprintf('%s/*_%s_*_%s_%s_*baseline*', ...
            thisDir, thisGCM, thisCrop, thisEmu) ;
        files = dir(pattern_in_BL) ;
        if length(files) > 1
            error('thisDir_yield ok, but error finding thisFile_BL: %d found', length(files))
        end
        thisFile_BL = sprintf('%s/%s', files.folder, files.name) ;
        
        % Get file with data for future
        pattern_in = sprintf('%s/*_%s_*_%s_%s_*movingwindow*', ...
            thisDir, thisGCM, thisCrop, thisEmu) ;
        files = dir(pattern_in) ;
        if length(files) ~= 1
            error('thisDir ok, but error finding thisFile: %d found', length(files))
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
            ii = strcmp(irrList_in,'ir') ;
            
            % Read IWD file
            thisVar = sprintf('iwd_%s', thisCrop_short) ;
            iwdBL_in_YX = flipud(transpose(ncread(thisFile_BL, thisVar))) ;
            data_in_bl_xcni(:,c,n,ii) = iwdBL_in_YX(gridlist.list_to_map) ;
            iwd_in_YXy = flip(permute(ncread(thisFile, thisVar), [2 1 3]), 1) ;
            if size(iwd_in_YXy,3) ~= Nyears
                error('Expected %d years in IWD file but found %d', ...
                    Nyears, size(iwd_in_YXy,3))
            end
            if length(ii) ~= 1
                error('Error finding ii: %d found', length(ii))
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
data_out_fu_xtv = nan(Ncells, Nts, length(cropIrrNlist_in)) ;
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
            
% Create required output arrays, converting tons/ha to kg/m2
tpha_to_kgpm2 = 0.1 ;

data_bl_emu.varNames = cropIrrNlist_in ;
data_bl_emu.list2map = gridlist.list_to_map ;
data_bl_emu.y1s = ts1_list ;
data_bl_emu.yNs = tsN_list ;
data_bl_emu.garr_xv = data_out_bl_xv * tpha_to_kgpm2;

data_fu_emu.varNames = cropIrrNlist_in ;
data_fu_emu.list2map = gridlist.list_to_map ;
data_fu_emu.y1s = ts1_list ;
data_fu_emu.yNs = tsN_list ;
data_fu_emu.garr_xvt = permute(data_out_fu_xtv, [1 3 2]) * tpha_to_kgpm2;


end


function out_xt = process_timesteps(ts1_list, tsN_list, yearList_in, ...
    in_YXy, gridlist)

Nts = length(ts1_list) ;
out_xt = nan(length(gridlist.list_to_map), Nts) ;

for t = 1:Nts
    thisY1 = ts1_list(t) ;
    thisYN = tsN_list(t) ;
    okyrs = yearList_in>=thisY1 & yearList_in<=thisYN ;
    tmp_YX = mean(in_YXy(:,:,okyrs),3) ;
    out_xt(:,t) = tmp_YX(gridlist.list_to_map) ;
end

end