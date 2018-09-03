function [total_unmet_crop_YX, total_unmet_past_YX, ...
    mid_y1_2deg_crop_YX, mid_y1_2deg_past_YX, mid_y1_2deg_ntrl_YX] = ...
    PLUMharm_getUnmet(...
    mid_y1_2deg_crop_YX, mid_y1_2deg_past_YX, vegd_2deg_y1_YX)

    % Check that neither crop nor pasture exceeds nonBare area
    % USING FIX TO MAKE UNMET POSITIVE
    unmetcrop1_YX = (mid_y1_2deg_crop_YX - vegd_2deg_y1_YX) .* (mid_y1_2deg_crop_YX > vegd_2deg_y1_YX) ;
    unmetpast1_YX = (mid_y1_2deg_past_YX - vegd_2deg_y1_YX) .* (mid_y1_2deg_past_YX > vegd_2deg_y1_YX) ;
    mid_y1_2deg_crop_YX = mid_y1_2deg_crop_YX - unmetcrop1_YX ;
    mid_y1_2deg_past_YX = mid_y1_2deg_past_YX - unmetpast1_YX ;
    
    % Check that neither crop nor pasture is negative
    % (unmet is NEGATIVE)
    unmetcrop0_YX = mid_y1_2deg_crop_YX .* (mid_y1_2deg_crop_YX < 0) ;
    unmetpast0_YX = mid_y1_2deg_past_YX .* (mid_y1_2deg_past_YX < 0) ;
    mid_y1_2deg_crop_YX = mid_y1_2deg_crop_YX - unmetcrop0_YX ;
    mid_y1_2deg_past_YX = mid_y1_2deg_past_YX - unmetpast0_YX ;
    
    % Check that nonBare area is conserved
    % (unmet is POSITIVE)
    exceedLand_YX = mid_y1_2deg_crop_YX + mid_y1_2deg_past_YX - vegd_2deg_y1_YX ;
    unmetcrop_YX = mid_y1_2deg_crop_YX ./ (mid_y1_2deg_crop_YX + mid_y1_2deg_past_YX + 1e-12) ...
        .* exceedLand_YX .* (exceedLand_YX>0) ;
    unmetpast_YX = mid_y1_2deg_past_YX ./ (mid_y1_2deg_crop_YX + mid_y1_2deg_past_YX + 1e-12) ...
        .* exceedLand_YX .* (exceedLand_YX>0) ;
    mid_y1_2deg_crop_YX = mid_y1_2deg_crop_YX - unmetcrop_YX ;
    mid_y1_2deg_past_YX = mid_y1_2deg_past_YX - unmetpast_YX ;
    mid_y1_2deg_agri_YX = mid_y1_2deg_crop_YX + mid_y1_2deg_past_YX ;
    mid_y1_2deg_ntrl_YX = vegd_2deg_y1_YX - mid_y1_2deg_agri_YX ;
    
    % Compute the total amount of crop or pasture increase / decrease that is not able to be met within the 2 degree gridcells
    total_unmet_crop_YX = unmetcrop1_YX + unmetcrop0_YX + unmetcrop_YX ;
    total_unmet_past_YX = unmetpast1_YX + unmetpast0_YX + unmetpast_YX ;


end