function [Sout, excl_vecs_after] = e2p_apply_max(...
    Sin, combineCrops_file, varargin)

excl_vecs = {} ;
if ~isempty(varargin)
    excl_vecs = varargin{1} ;
    if length(varargin) > 1
        error('Maximum 1 optional argument: excl_vecs')
    end
end

Sout = Sin ;

actually_emuBL_char = {} ;
if isfield(Sin, 'actually_emuBL')
    actually_emuBL_char = cell(size(Sin.actually_emuBL)) ;
    actually_emuBL_char(~Sin.actually_emuBL) = {'sim'} ;
    actually_emuBL_char(Sin.actually_emuBL) = {'*EMU*'} ;
elseif isfield(Sin, 'actually_emuBL_char')
    actually_emuBL_char = Sin.actually_emuBL_char ;
end

% Load existing info
load(combineCrops_file, 'combineCrops') ;
if isfield(Sin,'garr_xv')
    out_xvt = Sin.garr_xv ;
elseif isfield(Sin,'garr_xvt')
    out_xvt = Sin.garr_xvt ;
else
    error('???')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop through specified combinations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, ~, ~, cropListI_before] = e2p_get_names(Sout.varNames, []) ;
excl_vecs_before = excl_vecs ;

for t = 1:length(combineCrops)
    
    % Get info about this transformation
    combineCrops_this = combineCrops(t) ;
    which_is_max = combineCrops_this.whichmax_xvt ;
    varNames_sourceA = combineCrops_this.varNames_sourceA ;
    
    NsourceA = length(varNames_sourceA) ;
    
    for w = 1:NsourceA
        
        % Get indices
        thisA = varNames_sourceA{w} ;
        [i_theseABetc, i_thisM, Sout.varNames] = ...
            e2p_combineCropInds(thisA, Sout.varNames, combineCrops_this) ;
        
        if isfield(Sin,'garr_xv')
            data_theseABetc_xvt = Sin.garr_xv(:,i_theseABetc) ;
            tmp_x1t_size = size(Sin.garr_xv) ;
        elseif isfield(Sin,'garr_xvt')
            data_theseABetc_xvt = Sin.garr_xvt(:,i_theseABetc,:) ;
            tmp_x1t_size = size(Sin.garr_xvt) ;
        else
            error('???')
        end
        
        % Assuming this function is only ever being applied for irrigation,
        % we can safely fill NaN values for rainfed crops with zeros.
        this_is_rainfed = strcmp(getbasename(thisA), getbasenamei(thisA)) ;
        if this_is_rainfed
            if any(any(any(data_theseABetc_xvt > 0)))
                error('Expected all zeros for irrigation of %s', ...
                    getbasenamei(thisA))
            else
                out_xvt(:,i_thisM,:) = 0 ;
            end
        else
            tmp_x1t_size(2) = 1 ;
            tmp_x1t = nan(tmp_x1t_size) ;
            for s = 1:length(i_theseABetc)
                thisIsMax = which_is_max(:,w,:) == s ;
                data_thisSource_x1t = data_theseABetc_xvt(:,s,:) ;
                tmp_x1t(thisIsMax) = data_thisSource_x1t(thisIsMax) ;
            end
            out_xvt(:,i_thisM,:) = tmp_x1t ;
        end
        
        % Assign character-based designator for whether it was simulated or
        % emulated
        if ~isempty(actually_emuBL_char)
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
        end
        
        % Update exclusion arrays
        [~, ~, ~, cropListI_after] = e2p_get_names(Sout.varNames, [], true) ;
        if ~isempty(excl_vecs) && ~isequal(cropListI_before, cropListI_after)
            excl_vecs_after = ...
                e2p_update_excl_arrays( ...
                Sout.varNames, cropListI_before, cropListI_after, ...
                i_thisM, i_theseABetc, excl_vecs_before) ;
%             keyboard
%             table(cropListI_after', excl_vecs_after{4}')
            excl_vecs_before = excl_vecs_after ;
        end
        cropListI_before = cropListI_after ;
        
    end
    
end

if ~isempty(actually_emuBL_char)
    Sout.actually_emuBL_char = actually_emuBL_char ;
end

if isfield(Sin,'garr_xv')
    Sout.garr_xv = out_xvt ;
elseif isfield(Sin,'garr_xvt')
    Sout.garr_xvt = out_xvt ;
else
    error('???')
end

end