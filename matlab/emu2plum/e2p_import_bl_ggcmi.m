function S_out = e2p_import_bl_ggcmi(...
    varNames_emu, topDir_phase2, topDir_emu, ggcm, list2map, lonlats, ...
    adaptation, which_file, force_consistent_baseline)

warning('on','all')

% Get info
Nvars_emu = length(varNames_emu) ;

% Set up output structure
S_out.varNames = varNames_emu ;
S_out.list2map = list2map ;
S_out.lonlats = lonlats ;
S_out.garr_xv = nan(length(list2map),Nvars_emu) ;
S_out.actually_emu = false(size(varNames_emu)) ;

% Set up conversion factor
if strcmp(which_file, 'yield')
    conv_fact = 0.1 ; % tons/ha to kg/m2
elseif strcmp(which_file, 'gsirrigation')
    conv_fact = 1 ;
else
    error('which_file %s not recognized', which_file)
end

% Find "best" possible files
[file_list, fileVar_list, A_used, sim_used] = get_files( ...
    varNames_emu, adaptation, which_file, topDir_phase2, topDir_emu, ...
    ggcm, force_consistent_baseline, false) ;

% If you had to resort to emulator for ANY treatment of a crop, use
% emulator for ALL treatments of that crop
if force_consistent_baseline
    varNames_emu_asCrop = getbasename(varNames_emu) ;
    cropList_emu = unique(varNames_emu_asCrop) ;
    for c = 1:length(cropList_emu)
        thisCrop = cropList_emu{c} ;
        isThisCrop = strcmp(varNames_emu_asCrop, thisCrop) ;
        if any(~sim_used(isThisCrop))
            warning('Using emulated (not simulated) baseline for %s', thisCrop)
            sim_used(isThisCrop) = 0 ;
            % Find emulator files
            [tmp_file_list, tmp_fileVar_list, tmp_A_used, tmp_sim_used] = get_files( ...
                varNames_emu(isThisCrop), adaptation, which_file, topDir_phase2, topDir_emu, ...
                ggcm, force_consistent_baseline, true) ;
            file_list(isThisCrop) = tmp_file_list ;
            fileVar_list(isThisCrop) = tmp_fileVar_list ;
            A_used(isThisCrop) = tmp_A_used ;
            sim_used(isThisCrop) = tmp_sim_used ;
            % Sanity check
            if any(sim_used(isThisCrop))
                error('???')
            end
        end
    end
end

% Warn, if you didn't already
if force_consistent_baseline
    for v = 1:Nvars_emu
        % Get info about this crop
        [thisCropIrr, thisN, thisAN] = ...
            get_crop_info(varNames_emu{v}, adaptation) ;
        this_sim_used = sim_used(v) ;
        this_A_used = A_used(v) ;
        
        if this_A_used ~= adaptation
            if this_sim_used
                sim_or_em = 'simulated' ;
            else
                sim_or_em = '*emulated*' ;
            end
            replacementAN = sprintf('A%d_N%d', this_A_used, thisN) ;
            warning('Phase 2 %s %s %s not found; using %s instead', ...
                sim_or_em, thisCropIrr, thisAN, replacementAN) ;
        end
    end
end

% Import
for v = 1:Nvars_emu
    
    % Get info about this crop
    [~, ~, ~, isirrig, thisCrop, thisCrop_short] = ...
        get_crop_info(varNames_emu{v}, adaptation) ;
    
    % No irrigation water ever applied for non-irrigated crops
    if ~isirrig && strcmp(which_file, 'gsirrigation')
        S_out.garr_xv(:,v) = 0 ;
        continue
    end
    
    % Get information about this file
    thisFile = file_list{v} ;
    fileVar = fileVar_list{v} ;
    using_emulated = ~sim_used(v) ;
    
    % Import (exclude last timestep to avoid incomplete final seasons)
    varname = sprintf('%s_%s', fileVar, thisCrop_short) ;
    if using_emulated == 1
        S_out.actually_emu(v) = true ;
        if strcmp(which_file, 'yield')
            if isirrig
                varname = strrep(varname, '_', '_ir_') ;
            else
                varname = strrep(varname, '_', '_rf_') ;
            end
        end
        tmp_YX = flipud(transpose(ncread( ...
            thisFile, ...
            varname))) ;
    else
        tmp_XYt = ncread(thisFile, varname) ;
        tmp_YX = flipud(transpose(nanmean(tmp_XYt(:,:,1:end-1), 3))) ;
    end
    tmp = tmp_YX(list2map) ;
    if ~any(~isnan(tmp))
        warning('%s is all NaN', thisCrop)
    end
    S_out.garr_xv(:,v) = tmp * conv_fact ;
    
end

S_out.incl_years = e2p_get_incl_years(thisFile) ;

end


function [thisCropIrr, thisN, thisAN, isirrig, thisCrop, thisCrop_short] = ...
    get_crop_info(thisVar, adaptation)

thisCropIrr = getbasenamei(thisVar) ;
thisN = getN_num(thisVar) ;
thisAN = sprintf('A%d_N%d', adaptation, thisN) ;
isirrig = strcmp(thisCropIrr(end),'i') ;
thisCrop = thisCropIrr ;
if isirrig
    thisCrop = thisCrop(1:end-1) ;
end
thisCrop_short = e2p_get_thisCrop_short(thisCrop) ;

end


function fileVar = get_fileVar(which_file, ggcm)

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


function [thisFile, fileVar, using_emulated] = get_file(...
    topDir_phase2, topDir_emu, thisCrop, ...
    adaptation, ggcm, thisCrop_short, thisN, isirrig, which_file, ...
    allow_emulated)

% Parse filename parts
fileVar = get_fileVar(which_file, ggcm) ;
if isnan(thisN)
    N_char = 'NNA' ;
else
    N_char = ['N' num2str(thisN)] ;
end
if isirrig
    W_char = 'Winf' ;
else
    W_char = 'W0' ;
end

% First, try actual phase2 result
% OR skip if FORCING emulated
thisFile = '' ;
if allow_emulated < 2
    thisDir_phase2 = sprintf('%s/%s/A%d/%s', ...
        topDir_phase2, thisCrop, adaptation, fileVar) ;
    thisFile = sprintf(['%s/%s_agmerra_fullharm_%s_%s_global_annual_' ...
        '1980_2010_C360_T0_%s_%s_A%d.nc4'], ...
        thisDir_phase2, lower(ggcm), ...
        fileVar, thisCrop_short, W_char, N_char, adaptation) ;
    using_emulated = 0 ;
end

% If not found, see if there's an emulated version
if isempty(thisFile) || ~exist(thisFile, 'file')
    if allow_emulated > 0
        
        if isnan(thisN)
            error('Deal with N=NaN when looking for emulated replacement for Phase 2 baseline')
        end
        
        % Different variable name in emulator outputs than in Phase2 sims
        if ~strcmp(fileVar, 'yield')
            fileVar = 'iwd' ;
        end
        
        A_N = sprintf('A%d_N%d', adaptation, thisN) ;
        
        % Should be the same for any CMIP version, GCM, and SSP
        thisPattern = sprintf(['%s/%s/CMIP6/%s/ssp126/%s/cmip6_ssp126_' ...
            'UKESM1-0-LL_r1i1p1_%s_%s_%s_emulated_%s_baseline_1980_2010_average*.nc4'], ...
            topDir_emu, fileVar, A_N, ggcm, thisCrop, ggcm, A_N, fileVar) ;
        
        % Find matches
        S = dir(thisPattern) ;
        if isempty(S)
            % If still not found, save as empty char
            thisFile = '' ;
            using_emulated = -1 ;
        elseif length(S) > 1
            error('%d matches found for %s (expected 0 or 1)', ...
                length(S), thisPattern)
        else
            using_emulated = 1 ;
            thisFile = sprintf('%s/%s', S.folder, S.name) ;
            if ~exist(thisFile, 'file')
                error('???')
            end
        end
    else
        thisFile = '' ;
        using_emulated = -1 ;
    end
end

end


function [file_list, fileVar_list, A_used, sim_used] = get_files( ...
    varNames_emu, adaptation, which_file, topDir_phase2, topDir_emu, ...
    ggcm, force_consistent_baseline, force_emulator)

Nvars_emu = length(varNames_emu) ;
file_list = cell(Nvars_emu, 1) ;
fileVar_list = cell(Nvars_emu, 1) ;
A_used = adaptation * ones(Nvars_emu, 1) ;
sim_used = ones(Nvars_emu, 1) ;
for v = 1:Nvars_emu
    
    % Get info about this crop
    [thisCropIrr, thisN, thisAN, isirrig, thisCrop, thisCrop_short] = ...
        get_crop_info(varNames_emu{v}, adaptation) ;
    
    % No irrigation water ever applied for non-irrigated crops
    if ~isirrig && strcmp(which_file, 'gsirrigation')
        file_list{v} = '' ;
        fileVar_list{v} = '' ;
        continue
    end
    
    % Look for file, not allowing emulated version
    % OR skip if forcing to use emulated
    thisFile = '' ;
    if ~force_emulator
        [thisFile, fileVar] = get_file( ...
            topDir_phase2, topDir_emu, thisCrop, ...
            adaptation, ggcm, thisCrop_short, thisN, isirrig, which_file, ...
            false) ;
    end
    if isempty(thisFile)
        
        % Try A[other], again not allowing emulated version
        % OR skip if forcing to use emulated
        tmp_thisAdapt = ~adaptation ;
        tmp_thisN = thisN ;
        if ~force_emulator
            [thisFile, fileVar, using_emulated] = get_file( ...
                topDir_phase2, topDir_emu, thisCrop, ...
                tmp_thisAdapt, ggcm, thisCrop_short, tmp_thisN, isirrig, which_file, ...
                false) ;
        end
        
        if isempty(thisFile)
            
            % Look for canonical file again, this time allowing emulated
            % version
            [thisFile, fileVar, using_emulated] = get_file( ...
                topDir_phase2, topDir_emu, thisCrop, ...
                adaptation, ggcm, thisCrop_short, thisN, isirrig, which_file, ...
                1+force_emulator) ;
            
            if isempty(thisFile)
                % Try A[other], this time allowing emulated version
                tmp_thisAdapt = ~adaptation ;
                tmp_thisN = thisN ;
                [thisFile, fileVar, using_emulated] = get_file( ...
                    topDir_phase2, topDir_emu, thisCrop, ...
                    tmp_thisAdapt, ggcm, thisCrop_short, tmp_thisN, isirrig, which_file, ...
                    1+force_emulator) ;
                
                if isempty(thisFile)
                    
                    if force_consistent_baseline
                        error('Finding NaN fertilization files not yet coded for force_consistent_baseline')
                    end
                    
                    % Try NNA
                    tmp_thisAdapt = adaptation ;
                    tmp_thisN = NaN ;
                    [thisFile, fileVar, using_emulated] = get_file( ...
                        topDir_phase2, topDir_emu, thisCrop, ...
                        tmp_thisAdapt, ggcm, thisCrop_short, tmp_thisN, isirrig, which_file, ...
                        false) ;
                    
                    if isempty(thisFile)
                        % Try A[other] NNA
                        tmp_thisAdapt = ~adaptation ;
                        [thisFile, fileVar, using_emulated] = get_file( ...
                            topDir_phase2, topDir_emu, thisCrop, ...
                            tmp_thisAdapt, ggcm, thisCrop_short, tmp_thisN, isirrig, which_file, ...
                            false) ;
                        
                        if isempty(thisFile)
                            error('Could not find baseline file %s for %s', ...
                                which_file, thisCropIrr)
                        end
                    end
                end
            end
        end
        
        sim_used(v) = ~using_emulated ;
        
        % Warn if using a substitute file
        if isnan(tmp_thisN)
            replacementAN = sprintf('A%d_NNA', tmp_thisAdapt) ;
        else
            replacementAN = sprintf('A%d_N%d', tmp_thisAdapt, thisN) ;
            A_used(v) = tmp_thisAdapt ;
        end
        if ~force_consistent_baseline
            if using_emulated == 1
                warning('Phase 2 simulated %s %s not found; using *EMULATED* %s instead', ...
                    thisCropIrr, thisAN, replacementAN) ;
            else
                warning('Phase 2 simulated %s %s not found; using %s instead', ...
                    thisCropIrr, thisAN, replacementAN) ;
            end
        end
    end
    
    % Save information about this file
    file_list{v} = thisFile ;
    fileVar_list{v} = fileVar ;
    
end % Loop through variables


end

