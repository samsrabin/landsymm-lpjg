function Sout = e2p_apply_max_wheat(Sin, combineCrops_file)

Sout = Sin ;

actually_emuBL_char = {} ;
if isfield(Sin, 'actually_emuBL')
    actually_emuBL_char = cell(size(Sin.actually_emuBL)) ;
    actually_emuBL_char(~Sin.actually_emuBL) = {'sim'} ;
    actually_emuBL_char(Sin.actually_emuBL) = {'*EMU*'} ;
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

for t = 1:size(combineCrops, 1)
    
    % Get info about this transformation
    combineCrops_row = combineCrops(t,:) ;
%     combineCrops_sources = combineCrops_row{2} ;
%     combineCrops_dest = combineCrops_row{1} ;
    which_is_max = combineCrops_row{:,3} ;
    varNames_sourceA = combineCrops_row{:,4} ;
    
    NsourceA = length(varNames_sourceA) ;
    
    for w = 1:NsourceA
        
        % Get indices
        thisA = varNames_sourceA{w} ;
        [i_theseABetc, i_thisM, Sout.varNames] = ...
            e2p_combineCropInds(thisA, Sout.varNames, combineCrops_row) ;
        
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
        
%         % Assign character-based designator for whether it was simulated or
%         % emulated
%         if isfield(Sin, 'actually_emuBL')
%             if Sin.actually_emuBL(i_thisWW) == Sin.actually_emuBL(i_thisSW)
%                 actually_emuBL_char{i_thisMW} = actually_emuBL_char{i_thisWW} ; %#ok<AGROW>
%             else
%                 actually_emuBL_char{i_thisMW} = ...
%                     [actually_emuBL_char{i_thisWW} actually_emuBL_char{i_thisSW}] ; %#ok<AGROW>
%             end
%         end

        % Assign character-based designator for whether it was simulated or
        % emulated
        if isfield(Sin, 'actually_emuBL') && ~isempty(Sin.actually_emuBL)
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
    end
    
end

if isfield(Sin, 'actually_emuBL')
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