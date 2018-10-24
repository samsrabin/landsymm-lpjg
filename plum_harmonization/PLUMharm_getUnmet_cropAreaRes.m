function [unmet_YXv, ...
    mid_y1_2deg_agri_YXv, mid_y1_2deg_ntrl_YX] = ...
    PLUMharm_getUnmet_cropAreaRes(...
    mid1_y1_2deg_agri_YXv, base_2deg_vegd_YX, ...
    resArea_2deg_YX, sum_agri_y0_YX, debugIJ_2deg, dbCrop)

% Properly correct for zero-denominators? FALSE is default LUH1 behavior,
% which just adds 1e-12 to denominator.
proper_zero_denoms = true ;

% Misc.
Nagri = size(mid1_y1_2deg_agri_YXv,3) ;

% Get non-reserved area
unresArea_2deg_YX = base_2deg_vegd_YX - resArea_2deg_YX ;
if ~isempty(find(unresArea_2deg_YX>base_2deg_vegd_YX,1))
    error('How is unreserved area > base vegetated area?')
end
max_agri_y1_YX = max(unresArea_2deg_YX,sum_agri_y0_YX) ;
unresArea_2deg_YXv = repmat(unresArea_2deg_YX,[1 1 Nagri]) ;

% Only use output agri from now on
mid_y1_2deg_agri_YXv = mid1_y1_2deg_agri_YXv ;

do_debug = ~isempty(debugIJ_2deg) ;
if do_debug
    dbI = debugIJ_2deg(1) ;
    dbJ = debugIJ_2deg(2) ;
    fprintf('PLUMharm_getUnmet_cropAreaRes for cell (%d,%d):\n',dbI,dbJ) ;
    disp('   Initial values:')
    fprintf('      %s\t%0.4e\n',pad('agri_y0:',20),sum_agri_y0_YX(dbI,dbJ)) ;
    fprintf('      %s\t%0.4e\n',pad('agri_y1_mid:',20),sum(mid_y1_2deg_agri_YXv(dbI,dbJ,:),3)) ;
    if ~isempty(dbCrop)
        fprintf('      %s\t%0.4e\n',pad('this_y1_mid:',20),mid_y1_2deg_agri_YXv(dbI,dbJ,dbCrop)) ;
    end
    mid1_y1_2deg_ntrl_YX = base_2deg_vegd_YX - sum(mid1_y1_2deg_agri_YXv,3) ;
    fprintf('      %s\t%0.4e\n',pad('ntrl_y1_mid:',20),mid1_y1_2deg_ntrl_YX(dbI,dbJ)) ;
    fprintf('      %s\t%0.4e\n',pad('vegd_base:',20),base_2deg_vegd_YX(dbI,dbJ)) ;
end

% Check that no agri LU exceeds non-reserved area. Unmet will need to go into
% agri LUs of neighboring cells.
% USING FIX TO MAKE UNMET POSITIVE
% unmetA_YXv = (mid_y1_2deg_agri_YXv - base_2deg_vegd_YXv) .* (mid_y1_2deg_agri_YXv > base_2deg_vegd_YXv) ;
unmetA_YXv = (mid_y1_2deg_agri_YXv - unresArea_2deg_YXv) .* (mid_y1_2deg_agri_YXv > unresArea_2deg_YXv) ;
mid_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv - unmetA_YXv ;
if do_debug
    disp('   After check that no LU exceeds non-reserved area:')
    fprintf('      %s\t%0.4e\n',pad('agri_y1_out:',20),sum(mid_y1_2deg_agri_YXv(dbI,dbJ,:),3)) ;
    fprintf('      %s\t%0.4e\n',pad('unmetA:',20),sum(unmetA_YXv(dbI,dbJ,:),3)) ;
    if ~isempty(dbCrop)
        fprintf('      %s\t%0.4e\n',pad('this_y1_out:',20),mid_y1_2deg_agri_YXv(dbI,dbJ,dbCrop)) ;
        fprintf('      %s\t%0.4e\n',pad('unmetA_this:',20),unmetA_YXv(dbI,dbJ,dbCrop)) ;
    end
end

% Check that no agri LU is negative. Unmet will need to come from agri LUs
% of neighboring cells.
% (unmet is NEGATIVE)
unmetB_YXv = mid_y1_2deg_agri_YXv .* (mid_y1_2deg_agri_YXv < 0) ;
mid_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv - unmetB_YXv ;
if do_debug
    disp('   After check that no LU is negative:')
    fprintf('      %s\t%0.4e\n',pad('agri_y1_out:',20),sum(mid_y1_2deg_agri_YXv(dbI,dbJ,:),3)) ;
    fprintf('      %s\t%0.4e\n',pad('unmetB:',20),sum(unmetB_YXv(dbI,dbJ,:),3)) ;
    if ~isempty(dbCrop)
        fprintf('      %s\t%0.4e\n',pad('this_y1_out:',20),mid_y1_2deg_agri_YXv(dbI,dbJ,dbCrop)) ;
        fprintf('      %s\t%0.4e\n',pad('unmetB_this:',20),unmetB_YXv(dbI,dbJ,dbCrop)) ;
    end
end


% Check that sum(agri_y1) does not exceed non-reserved area, and if it
% does, make sure it does not exceed sum(agri_y0). I.e., make sure that it
% doesn't encroach MORE onto reserved land. NOTE that this also takes care
% of situation where ntrl is negative, because that indicates that agri
% exceeds vegd, which means that agri *definitely* exceeds non-reserved. It
% also, of course, takes care of conservation of vegetated area.
% (unmet is POSITIVE)
sum_agri_y1_YX = sum(mid_y1_2deg_agri_YXv,3) ;
isOK_YX = sum_agri_y1_YX <= max_agri_y1_YX ;
j = 0 ;
unmetD_YXv = zeros(size(mid_y1_2deg_agri_YXv)) ;
while any(any(~isOK_YX))
    j = j + 1 ;
    if j>50
        error('Possible infinite loop in "Check that sum(agri_y1) does not exceed non-reserved area."')
    end
    exceedLandR_YX = (sum_agri_y1_YX - max_agri_y1_YX) .* (sum_agri_y1_YX > max_agri_y1_YX) ;
    exceedLandR_YXv = repmat(exceedLandR_YX,[1 1 Nagri]) ;
    if proper_zero_denoms
        weights_YXv = mid_y1_2deg_agri_YXv ./ repmat(sum_agri_y1_YX,[1 1 Nagri]) ;
        unmetDtmp_YXv = weights_YXv .* exceedLandR_YXv ;
        unmetDtmp_YXv(repmat(sum_agri_y1_YX,[1 1 Nagri])==0) = 0 ;
    else
        unmetDtmp_YXv = mid_y1_2deg_agri_YXv ./ repmat(sum_agri_y1_YX + 1e-12,[1 1 Nagri]) ...
            .* exceedLandR_YXv ;
    end
    unmetD_YXv = unmetD_YXv + unmetDtmp_YXv ;
    mid_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv - unmetDtmp_YXv ;
    sum_agri_y1_YX = sum(mid_y1_2deg_agri_YXv,3) ;
    isOK_YX = sum_agri_y1_YX <= max_agri_y1_YX ;
end
mid_y1_2deg_ntrl_YX = base_2deg_vegd_YX - sum(mid_y1_2deg_agri_YXv,3) ;

% Compute the total amount of crop or pasture increase / decrease that is not able to be met within the 2 degree gridcells
unmet_YXv = unmetA_YXv + unmetB_YXv + unmetD_YXv;

if do_debug
    disp('   After check that sum(agri_y1) is compatible with reserved area:')
    fprintf('      %s\t%0.4e\n',pad('agri_y1_out:',20),sum(mid_y1_2deg_agri_YXv(dbI,dbJ,:),3)) ;
    fprintf('      %s\t%0.4e\n',pad('unmetD:',20),sum(unmetDtmp_YXv(dbI,dbJ,:),3)) ;
    fprintf('      %s\t%0.4e\n',pad('total_unmet:',20),sum(unmet_YXv(dbI,dbJ,:),3)) ;
    if ~isempty(dbCrop)
        fprintf('      %s\t%0.4e\n',pad('this_y1_out:',20),mid_y1_2deg_agri_YXv(dbI,dbJ,dbCrop)) ;
        fprintf('      %s\t%0.4e\n',pad('unmetD_this:',20),unmetDtmp_YXv(dbI,dbJ,dbCrop)) ;
        fprintf('      %s\t%0.4e\n',pad('total_unmet_this:',20),unmet_YXv(dbI,dbJ,dbCrop)) ;
    end
    fprintf('      %s\t%0.4e\n',pad('ntrl_y1_out:',20),mid_y1_2deg_ntrl_YX(dbI,dbJ)) ;
end


end