function [S, combineCrops_out, actually_emu_char, cf_list] = ...
    e2p_split_combine_burn(S, combineCrops_in, ...
    cropList_mid, cropList_out, cf_list, varargin)
% If tied, goes to alphabetically-first source type.

% Verbosity
% 0: Silent
% 1: E.g., "max_wheati <- max(spring_wheati, winter_wheati)"
% 2: E.g., "max_wheati0010 <- max(spring_wheati0010, winter_wheati0010)"
%          "max_wheati0060 <- max(spring_wheati0060, winter_wheati0060)"
%          etc.
verbose = 0 ;

varNames = S.varNames ;
data_in = S.garr_xvt ;

cropList_cf = cropList_mid ;
cropList_data = unique(getbasename(varNames)) ;

actually_emu = [] ;
if ~isempty(varargin)
    actually_emu = varargin{1} ;
    if length(varargin) > 1
        error('Maximum 1 optional argument: actually_emu')
    end
end

actually_emu_char = {} ;
if ~isempty(actually_emu)
    actually_emu_char = cell(size(actually_emu)) ;
    actually_emu_char(~actually_emu) = {'sim'} ;
    actually_emu_char(actually_emu) = {'*EMU*'} ;
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
    combineCrops_sources_orig = combineCrops_sources ;
    combineCrops_dest_orig = combineCrops_dest ;
    
    % If combination not needed, then skip
    C = intersect(cropList_data, combineCrops_sources) ;
    
    % Change target and destination names, if needed.
    % (Initially looks for CerealsC3s and CerealsC3w, as in LPJ-GUESS,
    % but other models use spring_wheat and winter_wheat.)
    needs_burnin = false ;
    if length(C) ~= length(combineCrops_sources)
        if ~strcmp(combineCrops_dest, 'CerealsC3')
            needs_burnin = true ;
            try
                I_data = e2p_translate_crops_lpj2out(...
                    cropList_data, combineCrops_sources) ;
            catch ME1
                try
                    I_data = e2p_translate_crops_agm2out(...
                        cropList_data, combineCrops_sources) ;
                catch ME2
                    error(['%s: Error finding combineCrops_sources in cropList_mid: ' ...
                        'Expected %d, found %d. Unable to find substitute crops, experiencing ' ...
                        'the following errors:\n\n%s\n\n%s'], ...
                        combineCrops_dest, length(combineCrops_sources), length(C), ...
                        ME1.message, ME2.message)
                end
            end
            combineCrops_sources = cropList_data(I_data) ;
            combineCrops_row{:,2} = combineCrops_sources ;
        else
            if ~isempty(C)
                error('How did you find %d wheats?', length(C))
            end
            combineCrops_sources2 = {'spring_wheat', 'winter_wheat'} ;
            C2 = intersect(cropList_data, combineCrops_sources2) ;
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
    
    for w = 1:NsourceA
        
        % Get indices
        thisA = varNames_sourceA{w} ;
        [i_theseABetc, i_thisM, varNames] = ...
            e2p_wheatInds(thisA, varNames, combineCrops_row) ;
        cropList_data = unique(getbasename(varNames)) ;
        
        if verbose == 2
            print_msg(varNames{i_thisM}, varNames(i_theseABetc), needs_burnin)
        end
        
        % Get maxima
        [M, I] = max(data_in(:,i_theseABetc,:), [], 2) ;
        I(isnan(M)) = NaN ;
        which_is_max(:,w,:) = I ;
        data_out(:,i_thisM,:) = M ;
        
        % Assign character-based designator for whether it was simulated or
        % emulated
        if ~isempty(actually_emu)
            actually_emu_these = actually_emu_char(i_theseABetc) ;
            if length(unique(actually_emu_these)) == 1
                actually_emu_char{i_thisM} = actually_emu_these{1} ; %#ok<AGROW>
            else
                tmp = '' ;
                for s = 1:length(i_theseABetc)
                    tmp = [tmp actually_emu_char{i_theseABetc(1)}] ; %#ok<AGROW>
                end
                actually_emu_char{i_thisM} = tmp ; %#ok<AGROW>
            end
        end
        
    end
    
    % Save results to combineCrops_out
    combineCrops_row{:,3} = which_is_max ;
    combineCrops_row{:,4} = varNames_sourceA ;
    combineCrops_out(t,:) = combineCrops_row ;
    
    % "Burn in" calibration factors?
    if needs_burnin
        
        if length(cropList_cf) ~= length(cf_list)
            error('Length mismatch between cropList_mid and cf_list')
        end
            
        % Multiply outputs by calibration factors
        for w = 1:NsourceA
            data_out_x1t = data_out(:,w,:) ;
            for s = 1:length(i_theseABetc)
                thisSource_cf = combineCrops_sources_orig{s} ;
                thisCF = cf_list(strcmp(cropList_cf, thisSource_cf)) ;
                if length(thisCF) ~= 1
                    error('Error finding thisCF for %s: expected 1, found %d', ...
                        thisSource_cf, length(thisCF))
                end
                thisIsMax = which_is_max(:,w,:) == s ;
                data_out_x1t(thisIsMax) = data_out_x1t(thisIsMax) * thisCF ;
            end
            data_out(:,w,:) = data_out_x1t ;
        end
        
        % Combine calibration factors
        [~,IA] = intersect(cropList_cf, combineCrops_sources_orig) ;
        cf_list(IA) = [] ;
        cropList_cf(IA) = [] ;
        cf_list(end+1) = 1 ; %#ok<AGROW>
        cropList_cf{end+1} = combineCrops_dest_orig ; %#ok<AGROW>
    end
    
end

% Sort new calibration factors list
[~, ~, IB] = intersect(cropList_out, cropList_cf, 'stable') ;
cropList_cf = cropList_cf(IB) ;
cf_list = cf_list(IB) ;

if ~isequal(cropList_cf, cropList_out)
    error('After all crop combining and burning-in, cropList_cf and cropList_out should be identical')
end

% Make output structure
S = rmfield(S, 'garr_xvt') ;
S.varNames = varNames ;
S.garr_xvt = data_out ;

end



function print_msg(thisDest, theseSources, needs_burnin)

% str = [thisDest ' <- max('] ;
% for s = 1:length(theseSources)
%     str = [str theseSources{s}] ;
%     if s < length(theseSources)
%         str = [str ', '] ;
%     else
%         str = [str ')'] ;
%     end
% end

dispStr = sprintf('%s <- max(%s)', ...
    thisDest, strjoin(theseSources, ', ')) ;
if needs_burnin
    dispStr = [dispStr '; burning in calib. factors'] ;
end

disp(dispStr)

end

