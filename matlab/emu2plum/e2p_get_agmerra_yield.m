function S_out = e2p_get_agmerra_yield(...
    varNames_emu, topDir_phase2, ggcm, list2map, lonlats, getN, adaptation, ...
    which_file, required)

warning('on','all')

% Get info
Nvars_emu = length(varNames_emu) ;
fileVar = get_fileVar(which_file, ggcm, adaptation) ;

% Set up output structure
S_out.varNames = varNames_emu ;
S_out.list2map = list2map ;
S_out.lonlats = lonlats ;
S_out.garr_xv = nan(length(list2map),Nvars_emu) ;

% Set up conversion factor
if strcmp(which_file, 'yield')
    conv_fact = 0.1 ; % tons/ha to kg/m2
elseif strcmp(which_file, 'gsirrigation')
    conv_fact = 1 ;
else
    error('which_file %s not recognized', which_file)
end

% Get AgMERRA-forced yield or irrigation water applied
for v = 1:Nvars_emu
    
    % Get info about this crop
    thisCropIrr = getbasenamei(varNames_emu{v}) ;
    thisN = str2double(getN(varNames_emu{v})) ;
    isirrig = strcmp(thisCropIrr(end),'i') ;
    thisCrop = thisCropIrr ;
    if isirrig
        thisCrop = thisCrop(1:end-1) ;
    end
    
    % No irrigation water ever applied for non-irrigated crops
    if ~isirrig && strcmp(which_file, 'gsirrigation')
        S_out.garr_xv(:,v) = 0 ;
        continue
    end
    
    % Get filename
    switch thisCrop
        case 'maize'; thisCrop_short = 'mai' ;
        case 'soy'; thisCrop_short = 'soy' ;
        case 'rice'; thisCrop_short = 'ric' ;
        case 'spring_wheat'; thisCrop_short = 'swh' ;
        case 'winter_wheat'; thisCrop_short = 'wwh' ;
        otherwise; error('thisCrop (%s) not recognized', thisCrop)
    end
    fileVar = get_fileVar(which_file, ggcm, adaptation) ;
    thisFile = sprintf('%s/%s/A%d/%s/%s_agmerra_fullharm_%s_%s_global_annual_1980_2010_C360_T0_W0_N%d_A%d.nc4', ...
        topDir_phase2, thisCrop, adaptation, fileVar, lower(ggcm), ...
        fileVar, thisCrop_short, thisN, ...
        adaptation) ;
    if isirrig
        thisFile = strrep(thisFile, 'W0', 'Winf') ;
    end
    nofilefound = false ;
    if ~exist(thisFile, 'file')
        thisAN = sprintf('A%d_N%d', adaptation, thisN) ;
        msg = sprintf('Phase 2 AgMERRA %s %s not found, nor were any substitutes', ...
            thisCropIrr, thisAN) ;
        msg = sprintf('%s (A%d_N%d, A%d_NNA, or A%d_NNA)', ...
            msg, ~adaptation, thisN, adaptation, ~adaptation) ;
        % Try A[other]
        tmp_thisAdapt = ~adaptation ;
        tmp_thisN = thisN ;
        fileVar = get_fileVar(which_file, ggcm, tmp_thisAdapt) ;
        thisFile = sprintf('%s/%s/A%d/%s/%s_agmerra_fullharm_%s_%s_global_annual_1980_2010_C360_T0_W0_N%d_A%d.nc4', ...
            topDir_phase2, thisCrop, tmp_thisAdapt, fileVar, lower(ggcm), ...
            fileVar, thisCrop_short, thisN, ...
            tmp_thisAdapt) ;
        if ~exist(thisFile, 'file')
            % Try NNA
            tmp_thisAdapt = adaptation ;
            tmp_thisN = NaN ;
            fileVar = get_fileVar(which_file, ggcm, tmp_thisAdapt) ;
            thisFile = sprintf('%s/%s/A%d/%s/%s_agmerra_fullharm_%s_%s_global_annual_1980_2010_C360_T0_W0_NNA_A%d.nc4', ...
                topDir_phase2, thisCrop, tmp_thisAdapt, fileVar, lower(ggcm), ...
                fileVar, thisCrop_short, ...
                tmp_thisAdapt) ;
            if ~exist(thisFile, 'file')
                % Try A[other] NNA
                tmp_thisAdapt = ~adaptation ;
                fileVar = get_fileVar(which_file, ggcm, tmp_thisAdapt) ;
                thisFile = sprintf('%s/%s/A%d/%s/%s_agmerra_fullharm_%s_%s_global_annual_1980_2010_C360_T0_W0_NNA_A%d.nc4', ...
                    topDir_phase2, thisCrop, tmp_thisAdapt, fileVar, lower(ggcm), ...
                    fileVar, thisCrop_short, ...
                    tmp_thisAdapt) ;
                if ~exist(thisFile, 'file')
                    if required
                        error('Could not find baseline AgMERRA %s for %s', ...
                            which_file, thisCrop)
                    else
                        nofilefound = true ;
                        warning('%s\n  Will not check for missing or low baseline AgMERRA %s!', ...
                            msg, which_file) ;
                    end
                end
            end
        end
        
        % Warn if using a substitute file
        if ~nofilefound
            if isnan(tmp_thisN)
                replacementAN = sprintf('A%d_NNA', tmp_thisAdapt) ;
            else
                replacementAN = sprintf('A%d_N%d', tmp_thisAdapt, thisN) ;
            end
            warning('Phase 2 AgMERRA %s %s not found; using %s instead', ...
                thisCropIrr, thisAN, replacementAN) ;
        end
    end
    
    % Import (exclude last timestep to avoid incomplete final seasons)
    if nofilefound
        S_out.garr_xv(:,v) = Inf ;
    else
        tmp_XYt = ncread(thisFile,sprintf('%s_%s', fileVar, thisCrop_short)) ;
        tmp_YX = flipud(transpose(nanmean(tmp_XYt(:,:,1:end-1), 3))) ;
        tmp = tmp_YX(list2map) ;
        if ~any(~isnan(tmp))
            warning('%s is all NaN', thisCrop)
        end
        S_out.garr_xv(:,v) = tmp * conv_fact ;
    end
    
end

end


function fileVar = get_fileVar(which_file, ggcm, adaptation)

if strcmp(which_file, 'yield')
    fileVar = 'yield' ;
elseif strcmp(which_file, 'gsirrigation')
    if strcmp(ggcm, 'EPIC-TAMU')
        % Misnamed
        fileVar = 'pirrw' ;
    else
        fileVar = 'pirrww' ;
    end
else
    error('which_file %s not recognized', which_file)
end


end
