function S_out = e2p_get_agmerra_baseline(...
    varNames_emu, topDir_phase2, topDir_emu, ggcm, list2map, lonlats, ...
    adaptation, which_file)

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

% Get AgMERRA-forced yield or irrigation water applied
for v = 1:Nvars_emu
    
    % Get info about this crop
    thisCropIrr = getbasenamei(varNames_emu{v}) ;
    thisN = str2double(getN(varNames_emu{v})) ;
    thisAN = sprintf('A%d_N%d', adaptation, thisN) ;
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
    
    thisCrop_short = e2p_get_thisCrop_short(thisCrop) ;
    
    % Look for file, not allowing emulated version
    [thisFile, fileVar, using_emulated] = get_file( ...
        topDir_phase2, topDir_emu, thisCrop, ...
        adaptation, ggcm, thisCrop_short, thisN, isirrig, which_file, ...
        false) ;
    nofilefound = false ;
    replacementAN = '' ;
    if isempty(thisFile)
        msg = sprintf('Phase 2 AgMERRA %s %s not found, nor were any substitutes', ...
            thisCropIrr, thisAN) ;
        msg = sprintf('%s (A%d_N%d, A%d_NNA, or A%d_NNA)', ...
            msg, ~adaptation, thisN, adaptation, ~adaptation) ;
        
        % Try A[other], again not allowing emulated version
        tmp_thisAdapt = ~adaptation ;
        tmp_thisN = thisN ;
        [thisFile, fileVar, using_emulated] = get_file( ...
            topDir_phase2, topDir_emu, thisCrop, ...
            tmp_thisAdapt, ggcm, thisCrop_short, tmp_thisN, isirrig, which_file, ...
            false) ;
        
        if isempty(thisFile)
            
            % Look for canonical file again, this time allowing emulated
            % version
            [thisFile, fileVar, using_emulated] = get_file( ...
                topDir_phase2, topDir_emu, thisCrop, ...
                adaptation, ggcm, thisCrop_short, thisN, isirrig, which_file, ...
                true) ;
            
            if isempty(thisFile)
                % Try A[other], this time allowing emulated version
                tmp_thisAdapt = ~adaptation ;
                tmp_thisN = thisN ;
                [thisFile, fileVar, using_emulated] = get_file( ...
                    topDir_phase2, topDir_emu, thisCrop, ...
                    tmp_thisAdapt, ggcm, thisCrop_short, tmp_thisN, isirrig, which_file, ...
                    true) ;
                
                if isempty(thisFile)
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
                            error('Could not find baseline AgMERRA %s for %s', ...
                                which_file, thisCropIrr)
                        end
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
            if using_emulated == 1
                warning('Phase 2 simulated %s %s not found; using *EMULATED* %s instead', ...
                    thisCropIrr, thisAN, replacementAN) ;
            else
                warning('Phase 2 simulated %s %s not found; using %s instead', ...
                    thisCropIrr, thisAN, replacementAN) ;
            end
        end
    end
    
    % Import (exclude last timestep to avoid incomplete final seasons)
    if nofilefound
        S_out.garr_xv(:,v) = Inf ;
    else
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
    
end

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
thisDir_phase2 = sprintf('%s/%s/A%d/%s', ...
    topDir_phase2, thisCrop, adaptation, fileVar) ;
thisFile = sprintf(['%s/%s_agmerra_fullharm_%s_%s_global_annual_' ...
    '1980_2010_C360_T0_%s_%s_A%d.nc4'], ...
    thisDir_phase2, lower(ggcm), ...
    fileVar, thisCrop_short, W_char, N_char, adaptation) ;
using_emulated = 0 ;

% If not found, see if there's an emulated version
if allow_emulated
    if ~exist(thisFile, 'file')
        
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
        
    end
else
    thisFile = '' ;
    using_emulated = -1 ;
end



end
