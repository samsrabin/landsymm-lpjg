function [data_out, varNames, combineCrops_out, actually_emuBL_char] = ...
    e2p_get_max_wheat(data_in, varNames, combineCrops_in, ...
    cropList_mid, cropList_out, cf_list, varargin)
% If tied, goes to alphabetically-first source type.

cropList_cf = cropList_mid ;
actually_emuBL = [] ;
if ~isempty(varargin)
    actually_emuBL = varargin{1} ;
    if length(varargin) > 1
        error('Maximum 1 optional argument: actually_emuBL')
    end
end

actually_emuBL_char = {} ;
if ~isempty(actually_emuBL)
    actually_emuBL_char = cell(size(actually_emuBL)) ;
    actually_emuBL_char(~actually_emuBL) = {'sim'} ;
    actually_emuBL_char(actually_emuBL) = {'*EMU*'} ;
end

data_out = data_in ;
Ncombines = size(combineCrops_in, 1) ;
combineCrops_out = [combineCrops_in cell(Ncombines, 2)] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop through specified combinations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:Ncombines
    
    % Get info about this transformation
    combineCrops_row = combineCrops_in(t,:) ;
    combineCrops_sources = combineCrops_row{2} ;
    combineCrops_dest = combineCrops_row{1} ;
    
    % Change target and destination names for wheat, if needed
    combineCrops_dest_orig = combineCrops_dest ;
    varNames_base = unique(getbasename(varNames)) ;
    C = intersect(varNames_base, combineCrops_sources) ;
    if length(C) ~= length(combineCrops_sources)
        if ~strcmp(combineCrops_dest, 'CerealsC3') || ~isempty(C)
            error('%s: Error finding combineCrops_sources in cropList_mid: Expected %d, found %d', ...
                combineCrops_dest, length(combineCrops_sources), length(C))
        else
            combineCrops_sources2 = {'spring_wheat', 'winter_wheat'} ;
            C2 = intersect(varNames_base, combineCrops_sources2) ;
            if length(C2) ~= length(combineCrops_sources2)
                st = dbstack ;
                error('%s for %s: Error finding combineCrops_sources in cropList_mid: Expected %d, found %d', ...
                    st.name, combineCrops_dest, length(combineCrops_sources), length(C))
            end
            combineCrops_sources = combineCrops_sources2 ;
            combineCrops_dest = 'max_wheat' ;
            combineCrops_row = {combineCrops_dest, combineCrops_sources} ;
        end
    end
    
    % If combination not needed, then skip
    if isempty(C) && any(strcmp(cropList_mid, combineCrops_dest))
%         fprintf('Skipping %s\n', combineCrops_dest)
        continue
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
    
    for w = 1:NsourceA
        
        % Get indices
        thisA = varNames_sourceA{w} ;
        [i_theseABetc, i_thisM, varNames] = ...
            e2p_combineCropInds(thisA, varNames, combineCrops_row) ;
        
        % Get maxima
        [M, I] = max(data_in(:,i_theseABetc,:), [], 2) ;
        I(isnan(M)) = NaN ;
        which_is_max(:,w,:) = I ;
        data_out(:,i_thisM,:) = M ;
        
        % Assign character-based designator for whether it was simulated or
        % emulated
        if ~isempty(actually_emuBL)
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
    
    % Save results to combineCrops_out
    combineCrops_row{:,3} = which_is_max ;
    combineCrops_row{:,4} = varNames_sourceA ;
    combineCrops_out(t,:) = combineCrops_row ;
    
    % "Burn in" calibration factors?
    if ~any(strcmp(cropList_mid, combineCrops_dest_orig))
        
        if length(cropList_mid) ~= length(cf_list)
            error('Length mismatch between cropList_mid and cf_list')
        end
        
        error('Need to add code for burn-in')
        
    end
    
    
end

% if ~isequal(cropList_cf, cropList_out)
%     error('After all crop combining and burning-in, cropList_cf and cropList_out should be identical')
% end

