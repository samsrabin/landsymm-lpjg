function [data_out, varNames, is_ww_max, winter_wheats] = e2p_get_max_wheat(data_in, varNames)
% If tied, goes to spring wheat.

data_out = data_in ;

% Get variable names of all winter wheats
winter_wheats = varNames(contains(varNames,'winter_wheat')) ;
Nww = length(winter_wheats) ;

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
    tmp(:,isnan(data_in(:,i_thisSW,:)) & ~isnan(data_in(:,i_thisSW,:)),:) = true ;
    is_ww_max(:,w,:) = tmp ;
    
    % Assign max to i_thisMW
    data_out(:,i_thisMW,:) = ...
        max(data_in(:,[i_thisWW i_thisSW],:), [], 2) ;
    
end