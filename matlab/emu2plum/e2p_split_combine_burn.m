function [S, combineCrops_out, actually_emuBL_char, cf_table_out, ...
excl_vecs_after] = ...
    e2p_split_combine_burn(S, combineCrops_in, ...
    cropList_out, cf_table_in, varargin)
% If tied, goes to alphabetically-first source type.

excl_vecs = {} ;
if ~isempty(varargin)
    excl_vecs = varargin{1} ;
    if length(varargin) > 1
        error('Maximum 1 optional argument: excl_vecs')
    end
end

% Verbosity
% 0: Silent
% 1: E.g., "max_wheati <- max(spring_wheati, winter_wheati)"
% 2: E.g., "max_wheati0010 <- max(spring_wheati0010, winter_wheati0010)"
%          "max_wheati0060 <- max(spring_wheati0060, winter_wheati0060)"
%          etc.
verbose = 1 ;

varNames = S.varNames ;
varNames_orig = varNames ;
data_in = S.garr_xvt ;

if isempty(cf_table_in)
    cf_table_out = table() ;
else
    cropList_cf = cf_table_in.Crop ;
    cf_values = cf_table_in.calibration_factor ;
end

cropList_data = unique(getbasename(varNames)) ;

actually_emuBL_char = {} ;
if isfield(S, 'actually_emuBL_char')
    actually_emuBL_char = S.actually_emuBL_char ;
    if any(cellfun(@isempty, actually_emuBL_char))
        error('Empty member(s) in actually_emuBL_char')
    end
end

data_out = data_in ;
Ncombines = length(combineCrops_in) ;
combineCrops_out = combineCrops_in ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop through specified combinations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, ~, ~, cropListI_before] = e2p_get_names(varNames, []) ;
excl_vecs_before = excl_vecs ;

for t = 1:Ncombines
    
    % Get info about this transformation
    combineCrops_this = combineCrops_in(t) ;
    combineCrops_sources = combineCrops_this.sourceCrops_cf ;
    combineCrops_dest = combineCrops_this.destCrop ;
    combineCrops_sources_orig = combineCrops_sources ;
    combineCrops_dest_orig = combineCrops_dest ;
    
    % Do we need to burn in calibration factors?
    if isempty(cf_table_in)
        needs_burnin = false ;
    else
        sources_in_cf = intersect(cropList_cf, combineCrops_sources) ;
        if any(strcmp(cropList_cf, combineCrops_dest))
            if ~isempty(sources_in_cf)
                error(['If destination crop is in calibration factor list, ', ...
                    'then no source crop may be.'])
            end
            needs_burnin = false ;
        elseif length(sources_in_cf) ~= length(combineCrops_sources)
            error(['If destination crop is not in calibration factor list, ', ...
                'then all sources must be (%d/%d found)'], ...
                length(sources_in_cf), length(combineCrops_sources))
        else
            needs_burnin = true ;
        end
    end
    
    % Change target and destination names, if needed.
    C = intersect(cropList_data, combineCrops_sources) ;
    if length(C) ~= length(combineCrops_sources)
        % First look for LPJ-GUESS-style crop names
        try
            I_data = e2p_translate_crops_lpj2out(...
                cropList_data, combineCrops_sources) ;
            
            % Failing that, look for GGCMI-style crop names
        catch ME1
            try
                I_data = e2p_translate_crops_agm2out(...
                    cropList_data, combineCrops_sources) ;
            catch ME2
                error(['%s: Error finding combineCrops_sources in cropList_cf: ' ...
                    'Expected %d, found %d. Unable to find substitute crops, experiencing ' ...
                    'the following errors:\n\n%s\n\n%s'], ...
                    combineCrops_dest, length(combineCrops_sources), length(C), ...
                    ME1.message, ME2.message)
            end
        end
        combineCrops_sources = cropList_data(I_data) ;
        combineCrops_this.sourceCrops_cf = combineCrops_sources ;
    end
    
    if verbose == 1
        print_msg(combineCrops_dest, combineCrops_sources, needs_burnin)
    end
    
    % Get variable names of all instances of sourceA
    sourceA = combineCrops_sources{1} ;
    varNames_sourceA = varNames(contains(varNames, sourceA)) ;
    NsourceA = length(varNames_sourceA) ;
    if NsourceA == 0
        error('No variables found matching %s', sourceA)
    end
    
    % Set up array for tracking which source was max
    tmp_size = size(data_in) ;
    tmp_size(2) = NsourceA ;
    which_is_max = nan(tmp_size) ;
    max_is_zero = false(tmp_size) ;
    
    for w = 1:NsourceA
        
        % Get indices
        thisA = varNames_sourceA{w} ;
        [i_theseABetc, i_thisM, varNames] = ...
            e2p_combineCropInds(thisA, varNames, combineCrops_this) ;
        cropList_data = unique(getbasename(varNames)) ;
        
        if verbose == 2
            print_msg(varNames{i_thisM}, varNames(i_theseABetc), needs_burnin)
        end
        
        % Get maxima
        [M, I] = max(data_in(:,i_theseABetc,:), [], 2) ;
        I(isnan(M)) = NaN ;
        which_is_max(:,w,:) = I ;
        max_is_zero(:,w,:) = M==0 ;
        data_out(:,i_thisM,:) = M ;
        
        check_noeffect_preexisting( ...
                data_in, varNames_orig, ...
                data_out, varNames)
        
        % Assign character-based designator for whether it was simulated or
        % emulated
        if isfield(S, 'actually_emuBL_char')
            actually_emuBL_these = actually_emuBL_char(i_theseABetc) ;
            if length(unique(actually_emuBL_these)) == 1
                actually_emuBL_char{i_thisM} = actually_emuBL_these{1} ; %#ok<AGROW>
            else
                tmp = '' ;
                for s = 1:length(i_theseABetc)
                    tmp = [tmp actually_emuBL_char{i_theseABetc(1)}] ; %#ok<AGROW>
                end
                actually_emuBL_char{i_thisM} = tmp ; %#ok<AGROW>
            end
            if any(cellfun(@isempty, actually_emuBL_char))
                error('Empty member(s) in actually_emuBL_char')
            end
        end
        
        % Update exclusion arrays
        [~, ~, ~, cropListI_after] = e2p_get_names(varNames, [], true) ;
        if ~isempty(excl_vecs) && ~isequal(cropListI_before, cropListI_after)
            excl_vecs_after = ...
                e2p_update_excl_arrays( ...
                varNames, cropListI_before, cropListI_after, ...
                i_thisM, i_theseABetc, excl_vecs_before) ;
%             keyboard
%             table(cropListI_after', excl_vecs_after{4}')
            excl_vecs_before = excl_vecs_after ;
        end
        cropListI_before = cropListI_after ;
        
    end
    
    % Save results to combineCrops_out
    combineCrops_this.whichmax_xvt = which_is_max ;
    combineCrops_this.maxzero_xvt = max_is_zero ;
    combineCrops_this.varNames_sourceA = varNames_sourceA ;
    combineCrops_out(t) = combineCrops_this ;
    
    % "Burn in" calibration factors?
    if needs_burnin
        
        if length(cropList_cf) ~= length(cf_values)
            error('Length mismatch between cropList_cf and cf_values')
        end
            
        % Multiply outputs by calibration factors
        for w = 1:NsourceA
            
            % Get index of target variable
            thisA = varNames_sourceA{w} ;
            [~, i_thisM_burning, varNames_burning] = ...
                e2p_combineCropInds(thisA, varNames, combineCrops_this) ;
            if ~isequal(varNames, varNames_burning)
                error('varNames should not have changed just now...')
            end
        
            data_out_x1t = data_out(:,i_thisM_burning,:) ;
            for s = 1:length(i_theseABetc)
                thisSource_cf = combineCrops_sources_orig{s} ;
                thisCF = cf_values(strcmp(cropList_cf, thisSource_cf)) ;
                if length(thisCF) ~= 1
                    error('Error finding thisCF for %s: expected 1, found %d', ...
                        thisSource_cf, length(thisCF))
                end
                thisIsMax = which_is_max(:,w,:) == s ;
                data_out_x1t(thisIsMax) = data_out_x1t(thisIsMax) * thisCF ;
            end
            data_out(:,i_thisM_burning,:) = data_out_x1t ;
            
            check_noeffect_preexisting( ...
                data_in, varNames_orig, ...
                data_out, varNames)
            
        end
        
        % Combine calibration factors
        [~,IA] = intersect(cropList_cf, combineCrops_sources_orig) ;
        cf_values(IA) = [] ;
        cropList_cf(IA) = [] ;
        cf_values(end+1) = 1 ; %#ok<AGROW>
        cropList_cf{end+1} = combineCrops_dest_orig ; %#ok<AGROW>
    end
    
end

% Sort new calibration factors list
if ~isempty(cf_table_in)
    [~, ~, IB] = intersect(cropList_out, cropList_cf, 'stable') ;
    cropList_cf = cropList_cf(IB) ;
    cf_values = cf_values(IB) ;
    if ~isequal(shiftdim(cropList_cf), shiftdim(cropList_out))
        error('After all crop combining and burning-in, cropList_cf and cropList_out should be identical')
    end
    cf_table_out = table(cropList_cf, cf_values, ...
        'VariableNames', cf_table_in.Properties.VariableNames) ;
end

check_noeffect_preexisting( ...
    data_in, varNames_orig, ...
    data_out, varNames)

% Make output structure
S = rmfield(S, 'garr_xvt') ;
S.varNames = varNames ;
if isfield(S, 'actually_emuBL_char')
    S.actually_emuBL_char = actually_emuBL_char ;
end
S.garr_xvt = data_out ;

end


function print_msg(thisDest, theseSources, needs_burnin)

dispStr = sprintf('%s <- max(%s)', ...
    thisDest, strjoin(theseSources, ', ')) ;
if needs_burnin
    dispStr = [dispStr '; burning in calib. factors'] ;
end

disp(dispStr)

end


function check_noeffect_preexisting( ...
    data_in, varNames_orig, ...
    data_out, varNames)
% Make sure that no variables present in data_in differ between data_in
% and data_out.

[C, IA, IB] = intersect(varNames_orig, varNames) ;
if ~isequaln(data_in(:,IA,:), data_out(:,IB,:))
    affected_variables = {} ;
    for cc = 1:length(C)
        if ~isequaln(data_in(:,IA(cc),:), data_out(:,IB(cc),:))
            affected_variables = [affected_variables C(cc)] ;
        end
    end
    affected_variables
    error('Why were pre-existing variables affected?')
end

end



