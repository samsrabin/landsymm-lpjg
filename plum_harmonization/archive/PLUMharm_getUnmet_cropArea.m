function [total_unmet_agri_YXv, ...
    mid_y1_2deg_agri_YXv, mid_y1_2deg_ntrl_YX] = ...
    PLUMharm_getUnmet_cropArea(...
    mid_y1_2deg_agri_YXv, vegd_2deg_y1_YX)

Nagri = size(mid_y1_2deg_agri_YXv,3) ;
vegd_2deg_y1_YXv = repmat(vegd_2deg_y1_YX,[1 1 Nagri]) ;

% Check that neither crop nor pasture exceeds nonBare area
% USING FIX TO MAKE UNMET POSITIVE
unmetagri1_YXv = (mid_y1_2deg_agri_YXv - vegd_2deg_y1_YXv) .* (mid_y1_2deg_agri_YXv > vegd_2deg_y1_YXv) ;
mid_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv - unmetagri1_YXv ;

% Check that neither crop nor pasture is negative
% (unmet is NEGATIVE)
unmetagri0_YXv = mid_y1_2deg_agri_YXv .* (mid_y1_2deg_agri_YXv < 0) ;
mid_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv - unmetagri0_YXv ;

% Check that nonBare area is conserved
% (unmet is POSITIVE)
exceedLand_YX = sum(mid_y1_2deg_agri_YXv,3) - vegd_2deg_y1_YX ;
exceedLand_YXv = repmat(exceedLand_YX,[1 1 Nagri]) ;
unmetagri_YXv = mid_y1_2deg_agri_YXv ./ repmat(sum(mid_y1_2deg_agri_YXv,3) + 1e-12,[1 1 Nagri]) ...
    .* exceedLand_YXv .* (exceedLand_YXv>0) ;
mid_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv - unmetagri_YXv ;
mid_y1_2deg_ntrl_YX = vegd_2deg_y1_YX - sum(mid_y1_2deg_agri_YXv,3) ;

% Compute the total amount of crop or pasture increase / decrease that is not able to be met within the 2 degree gridcells
total_unmet_agri_YXv = unmetagri1_YXv + unmetagri0_YXv + unmetagri_YXv ;


end