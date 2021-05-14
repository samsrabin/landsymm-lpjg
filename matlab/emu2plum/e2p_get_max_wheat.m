function [data_out, varNames, is_ww_max, winter_wheats, actually_emu_char] = ...
    e2p_get_max_wheat(data_in, varNames, varargin)
% If tied, goes to spring wheat.

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

% Get variable names of all winter wheats
winter_wheats = varNames(contains(varNames,{'winter_wheat', 'CerealsC3w'})) ;
Nww = length(winter_wheats) ;
if Nww == 0
    error('No winter wheats found!')
end

% Set up array for tracking whether winter wheat was max
tmp_size = size(data_in) ;
tmp_size(2) = Nww ;
is_ww_max = false(tmp_size) ;

for w = 1:Nww
    
    % Get indices
    thisWW = winter_wheats{w} ;
    [i_thisWW, i_thisSW, i_thisMW, varNames] = ...
        e2p_wheatInds(thisWW, varNames) ;
    
    % Which wheat is max?
    is_ww_max(:,w,:) = ... 
        data_in(:,i_thisWW,:) > data_in(:,i_thisSW,:) ;
    
    % Above is always false if either winter or spring wheat is NaN. Here,
    % make sure that winter wheat is correctly registered as maximum if it
    % has a yield and spring wheat has NaN.
    tmp = is_ww_max(:,w,:) ;
    tmp(isnan(data_in(:,i_thisSW,:)) & ~isnan(data_in(:,i_thisWW,:))) = true ;
    is_ww_max(:,w,:) = tmp ;
    
    % Assign max to i_thisMW
    data_out(:,i_thisMW,:) = ...
        max(data_in(:,[i_thisWW i_thisSW],:), [], 2) ;
    
    % Assign character-based designator for whether it was simulated or
    % emulated
    if ~isempty(actually_emu)
        if actually_emu(i_thisWW) == actually_emu(i_thisSW)
            actually_emu_char{i_thisMW} = actually_emu_char{i_thisWW} ;
        else
            actually_emu_char{i_thisMW} = ...
                [actually_emu_char{i_thisWW} actually_emu_char{i_thisSW}] ;
        end
    end
    
end