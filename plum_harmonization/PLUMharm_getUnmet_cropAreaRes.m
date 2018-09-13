function [total_unmet_agri_YXv, ...
    out_y1_2deg_agri_YXv, out_y1_2deg_ntrl_YX] = ...
    PLUMharm_getUnmet_cropAreaRes(...
    mid_y1_2deg_agri_YXv, base_vegd_2deg_YX, ...
    resArea_2deg_YX, sum_agri_y0_YX, debugIJ_2deg)

% Properly correct for zero-denominators? FALSE is default LUH1 behavior,
% which just adds 1e-12 to denominator.
proper_zero_denoms = true ;

Nagri = size(mid_y1_2deg_agri_YXv,3) ;
base_vegd_2deg_YXv = repmat(base_vegd_2deg_YX,[1 1 Nagri]) ;
mid_y1_2deg_ntrl_YX = base_vegd_2deg_YX - sum(mid_y1_2deg_agri_YXv,3) ;

do_debug = ~isempty(debugIJ_2deg) ;
if do_debug
    i = debugIJ_2deg(1) ;
    j = debugIJ_2deg(2) ;
    disp('Initial values:')
    fprintf('   %s\t%0.4e\n',pad('agri_y0:',20),sum_agri_y0_YX(i,j)) ;
    fprintf('   %s\t%0.4e\n',pad('agri_y1_mid:',20),sum(mid_y1_2deg_agri_YXv(i,j,:),3)) ;
    fprintf('   %s\t%0.4e\n',pad('ntrl_y1_mid:',20),mid_y1_2deg_ntrl_YX(i,j)) ;
    fprintf('   %s\t%0.4e\n',pad('vegd_base:',20),base_vegd_2deg_YX(i,j)) ;
end

% Only use out_y1_2deg_*_YXv from here on (i.e., not
% mid_y1_2deg_*_YXv or base_vegd_2deg_YX)
out_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv ;
out_y1_2deg_ntrl_YX = mid_y1_2deg_ntrl_YX ;
out_y1_vegd_2deg_YX = base_vegd_2deg_YX ;

% Check that no LU exceeds nonBare area
% USING FIX TO MAKE UNMET POSITIVE
unmetagri1_YXv = (out_y1_2deg_agri_YXv - base_vegd_2deg_YXv) .* (out_y1_2deg_agri_YXv > base_vegd_2deg_YXv) ;
out_y1_2deg_agri_YXv = out_y1_2deg_agri_YXv - unmetagri1_YXv ;
out_y1_2deg_ntrl_YX = out_y1_vegd_2deg_YX - sum(out_y1_2deg_agri_YXv,3) ;
out_y1_2deg_vegd_YX = sum(out_y1_2deg_agri_YXv,3) + out_y1_2deg_ntrl_YX ;
if do_debug
    disp('After check that no LU exceeds nonBare area:')
    fprintf('   %s\t%0.4e\n',pad('agri_y1_out:',20),sum(out_y1_2deg_agri_YXv(i,j,:),3)) ;
    fprintf('   %s\t%0.4e\n',pad('ntrl_y1_out:',20),out_y1_2deg_ntrl_YX(i,j)) ;
    fprintf('   %s\t%0.4e\n',pad('vegd_y1_out:',20),out_y1_2deg_vegd_YX(i,j)) 
end

% Check that no LU is negative
% (unmet is NEGATIVE)
unmetagri0_YXv = out_y1_2deg_agri_YXv .* (out_y1_2deg_agri_YXv < 0) ;
out_y1_2deg_agri_YXv = out_y1_2deg_agri_YXv - unmetagri0_YXv ;
out_y1_2deg_ntrl_YX = out_y1_vegd_2deg_YX - sum(out_y1_2deg_agri_YXv,3) ;
out_y1_2deg_vegd_YX = sum(out_y1_2deg_agri_YXv,3) + out_y1_2deg_ntrl_YX ;
if do_debug
    disp('After check that no LU is negative:')
    fprintf('   %s\t%0.4e\n',pad('agri_y1_out:',20),sum(out_y1_2deg_agri_YXv(i,j,:),3)) ;
    fprintf('   %s\t%0.4e\n',pad('ntrl_y1_out:',20),out_y1_2deg_ntrl_YX(i,j)) ;
    fprintf('   %s\t%0.4e\n',pad('vegd_y1_out:',20),out_y1_2deg_vegd_YX(i,j)) 
end

% Check that vegetated area is conserved
% (unmet is POSITIVE)
exceedLand_YX = sum(out_y1_2deg_agri_YXv,3) - base_vegd_2deg_YX ;
exceedLand_YXv = repmat(exceedLand_YX,[1 1 Nagri]) ;
if proper_zero_denoms
    unmetagri_YXv = out_y1_2deg_agri_YXv ./ repmat(sum(out_y1_2deg_agri_YXv,3),[1 1 Nagri]) ...
        .* exceedLand_YXv .* (exceedLand_YXv>0) ;
    unmetagri_YXv(repmat(sum(out_y1_2deg_agri_YXv,3),[1 1 Nagri])==0) = 0 ;
else
    unmetagri_YXv = out_y1_2deg_agri_YXv ./ repmat(sum(out_y1_2deg_agri_YXv,3) + 1e-12,[1 1 Nagri]) ...
        .* exceedLand_YXv .* (exceedLand_YXv>0) ;
end
out_y1_2deg_agri_YXv = out_y1_2deg_agri_YXv - unmetagri_YXv ;
out_y1_2deg_ntrl_YX = out_y1_vegd_2deg_YX - sum(out_y1_2deg_agri_YXv,3) ;
out_y1_2deg_vegd_YX = sum(out_y1_2deg_agri_YXv,3) + out_y1_2deg_ntrl_YX ;
if do_debug
    disp('After check that vegetated area is conserved:')
    fprintf('   %s\t%0.4e\n',pad('agri_y1_out:',20),sum(out_y1_2deg_agri_YXv(i,j,:),3)) ;
    fprintf('   %s\t%0.4e\n',pad('ntrl_y1_out:',20),out_y1_2deg_ntrl_YX(i,j)) ;
    fprintf('   %s\t%0.4e\n',pad('vegd_y1_out:',20),out_y1_2deg_vegd_YX(i,j)) 
end

% Check that sum(agri_y1) does not exceed non-reserved area, and if it
% does, make sure it does not exceed sum(agri_y0). I.e., make sure that it
% doesn't encroach MORE onto reserved land. NOTE that this also takes care
% of situation where ntrl is negative, because that indicates that agri
% exceeds vegd, which means that agri *definitely* exceeds non-reserved.
% (unmet is POSITIVE)
unresArea_2deg_YX = base_vegd_2deg_YX - resArea_2deg_YX ;
sum_agri_y1_YX = sum(out_y1_2deg_agri_YXv,3) ;
% if do_debug
%     disp('Unaffected:')
%     fprintf('   %s\t%0.4e\n',pad('vegd_y1_out:',20),out_y1_2deg_vegd_YX(i,j)) ;
%     fprintf('   %s\t%0.4e\n',pad('reserved:',20),resArea_2deg_YX(i,j)) ;
%     fprintf('   %s\t%0.4e\n',pad('unreserved:',20),unresArea_2deg_YX(i,j)) ;
%     fprintf('   %s\t%0.4e\n',pad('agri_y0:',20),sum_agri_y0_YX(i,j)) ;
%     fprintf('   %s\t%0.4e\n',pad('ntrl_y0:',20),base_vegd_2deg_YX(i,j)-sum_agri_y0_YX(i,j)) ;
% end
% if do_debug
%     disp('Before fix:')
%     fprintf('   %s\t%0.4e\n',pad('agri_y1:',20),sum_agri_y1_YX(i,j)) ;
%     fprintf('   %s\t%0.4e\n',pad('ntrl_y1:',20),out_y1_2deg_ntrl_YX(i,j)) ;
%     fprintf('   %s\t%0.4e\n',pad('unmet:',20),sum(unmetagri1_YXv(i,j,:) + unmetagri0_YXv(i,j,:) + unmetagri_YXv(i,j,:),3)) ;
% end
isOK_YX = sum_agri_y1_YX<=sum_agri_y0_YX | sum_agri_y1_YX<=unresArea_2deg_YX ;
tmp_YX = max(unresArea_2deg_YX,sum_agri_y0_YX) ;
exceedLandR_YX = zeros(size(tmp_YX)) ;
exceedLandR_YX(~isOK_YX) = sum_agri_y1_YX(~isOK_YX) - tmp_YX(~isOK_YX) ;
exceedLandR_YXv = repmat(exceedLandR_YX,[1 1 Nagri]) ;
if proper_zero_denoms
    unmetagriR_YXv = out_y1_2deg_agri_YXv ./ repmat(sum_agri_y1_YX,[1 1 Nagri]) ...
        .* exceedLandR_YXv ;
    unmetagriR_YXv(repmat(sum_agri_y1_YX,[1 1 Nagri])==0) = 0 ;
else
    unmetagriR_YXv = out_y1_2deg_agri_YXv ./ repmat(sum_agri_y1_YX + 1e-12,[1 1 Nagri]) ...
        .* exceedLandR_YXv ;
end
% if do_debug
%     disp('During:')
%     fprintf('   %s\t%0.4e\n',pad('exceedLandR:',20),exceedLandR_YX(i,j)) ;
%     fprintf('   %s\t%0.4e\n',pad('unmetagriR:',20),sum(unmetagriR_YXv(i,j,:),3)) ;
% end
out_y1_2deg_agri_YXv = out_y1_2deg_agri_YXv - unmetagriR_YXv ;
out_y1_2deg_ntrl_YX = out_y1_vegd_2deg_YX - sum(out_y1_2deg_agri_YXv,3) ;
out_y1_2deg_vegd_YX = sum(out_y1_2deg_agri_YXv,3) + out_y1_2deg_ntrl_YX ;

% Compute the total amount of crop or pasture increase / decrease that is not able to be met within the 2 degree gridcells
total_unmet_agri_YXv = unmetagri1_YXv + unmetagri0_YXv + unmetagri_YXv + unmetagriR_YXv;

% if do_debug
%     disp('After fix:')
%     fprintf('   %s\t%0.4e\n',pad('agri_y1:',20),sum(out_y1_2deg_agri_YXv(i,j,:),3)) ;
%     fprintf('   %s\t%0.4e\n',pad('ntrl_y1:',20),out_y1_2deg_ntrl_YX(i,j)) ;
%     fprintf('   %s\t%0.4e\n',pad('unmet:',20),sum(total_unmet_agri_YXv(i,j,:),3)) ;
% end
if do_debug
    disp('After check that sum(agri_y1) is compatible with reserved area:')
    fprintf('   %s\t%0.4e\n',pad('agri_y1_out:',20),sum(out_y1_2deg_agri_YXv(i,j,:),3)) ;
    fprintf('   %s\t%0.4e\n',pad('ntrl_y1_out:',20),out_y1_2deg_ntrl_YX(i,j)) ;
    fprintf('   %s\t%0.4e\n',pad('vegd_y1_out:',20),out_y1_2deg_vegd_YX(i,j)) 
end


end