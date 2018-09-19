function [unmetMgmt_YXv, mid2_y1_2deg_mgmt_YXv] = ...
    PLUMharm_getUnmet_mgmt(...
    out_y0_2deg_mgmt_YXv, out_y0_2deg_crop_YXv, ...
    in_y0_2deg_mgmt_YXv, in_y1_2deg_mgmt_YXv, ...
    in_y0_2deg_crop_YXv, in_y1_2deg_crop_YXv, mid_y1_2deg_crop_YXv, ...
    max_mgmt, debugIJ_2deg, mgmtName, dbCrop)

% How much tolerance do we have for TOTAL management error?
mgmtTot_tol = 1e-6 ;

% Get debugging info
do_debug = ~isempty(debugIJ_2deg) ;
if do_debug
    dbi = debugIJ_2deg(1) ;
    dbj = debugIJ_2deg(2) ;
    fprintf('PLUMharm_getUnmet_mgmt %s, [%d %d] LU %d (mgmtTot values):\n',...
        mgmtName, dbi, dbj, dbCrop) ;
end

% Get totals
out_y0_2deg_mgmtTot_YXv = out_y0_2deg_mgmt_YXv .* out_y0_2deg_crop_YXv ;
in_y0_2deg_mgmtTot_YXv = in_y0_2deg_mgmt_YXv .* in_y0_2deg_crop_YXv ;
in_y1_2deg_mgmtTot_YXv = in_y1_2deg_mgmt_YXv .* in_y1_2deg_crop_YXv ;

% Get delta
mgmtTot_d_YXv = in_y1_2deg_mgmtTot_YXv - in_y0_2deg_mgmtTot_YXv ;

% Debugging
if do_debug
    fprintf('   Initial delta:    %0.4e\n', mgmtTot_d_YXv(dbi,dbj,dbCrop)) ;
end

% Do not allow delta that would drop mgmt application to <0
% (i.e., too much loss of MGMTAPP: UNMET SHOULD BE NEGATIVE)
isTooMuchLoss = out_y0_2deg_mgmtTot_YXv + mgmtTot_d_YXv < -mgmtTot_tol ;
i = 0 ;
unmetMgmtA_YXv = zeros(size(out_y0_2deg_mgmtTot_YXv)) ;
while any(isTooMuchLoss(:))
    i = i + 1 ;
    if i==50
        error('Infinite loop?')
    end
    tooMuchByThis = (out_y0_2deg_mgmtTot_YXv + mgmtTot_d_YXv) .* isTooMuchLoss ;
    unmetMgmtA_YXv = unmetMgmtA_YXv + tooMuchByThis ;
    mgmtTot_d_YXv = mgmtTot_d_YXv - tooMuchByThis ;
    isTooMuchLoss = out_y0_2deg_mgmtTot_YXv + mgmtTot_d_YXv < -mgmtTot_tol ;
end
mid1_y1_2deg_mgmtTot_YXv = out_y0_2deg_mgmtTot_YXv + mgmtTot_d_YXv ;
mid1_y1_2deg_mgmtTot_YXv(mid1_y1_2deg_mgmtTot_YXv<0) = 0 ;   % Code above should ensure these negative values are tiny

% Debugging
if do_debug
    fprintf('   After avoiding <0:   %0.4e\n', mid1_y1_2deg_mgmtTot_YXv(dbi,dbj,dbCrop)) ;
    fprintf('              UnmetA:   %0.4e\n', unmetMgmtA_YXv(dbi,dbj,dbCrop)) ;
end

% Do not allow delta that would increase mgmt application to above maximum
% per-area rate, considering what the NEW areas look like. (Unmet should be
% POSITIVE)
max_mgmt_YXv = repmat(permute(max_mgmt,[2 3 1]),[size(mid1_y1_2deg_mgmtTot_YXv,1) size(mid1_y1_2deg_mgmtTot_YXv,2) 1]) ;
max_mgmtTot_YXv = max_mgmt_YXv .* mid_y1_2deg_crop_YXv ;
if do_debug
    fprintf('   Crop area:   %0.4e\n', mid_y1_2deg_crop_YXv(dbi,dbj,dbCrop)) ;
    fprintf('   Max mgmt:    %0.4e\n', max_mgmtTot_YXv(dbi,dbj,dbCrop)) ;
end
isTooMuchNow = mid1_y1_2deg_mgmtTot_YXv > max_mgmtTot_YXv ;
tooMuchByThis = (mid1_y1_2deg_mgmtTot_YXv - max_mgmtTot_YXv) .* isTooMuchNow ;
unmetMgmtB_YXv = tooMuchByThis ;
mid2_y1_2deg_mgmtTot_YXv = mid1_y1_2deg_mgmtTot_YXv ;
mid2_y1_2deg_mgmtTot_YXv(isTooMuchNow) = max_mgmtTot_YXv(isTooMuchNow) ;
isTooMuchNow = mid2_y1_2deg_mgmtTot_YXv > max_mgmtTot_YXv ;
if any(isTooMuchNow(:))
    error('Go back to WHILE loop for upper-limiting mgmt inputs?')
end

% Get function outputs
mid2_y1_2deg_mgmt_YXv = mid2_y1_2deg_mgmtTot_YXv ./ mid_y1_2deg_crop_YXv ;
mid2_y1_2deg_mgmt_YXv(mid_y1_2deg_crop_YXv==0) = 0 ;
mid2_y1_2deg_mgmt_YXv(mid2_y1_2deg_mgmt_YXv>max_mgmt_YXv) = max_mgmt_YXv(mid2_y1_2deg_mgmt_YXv>max_mgmt_YXv) ; % Code above should make sure these excesses are tiny
unmetMgmt_YXv = unmetMgmtA_YXv + unmetMgmtB_YXv ;

% Debugging
if do_debug
    fprintf('   After avoiding >max: %0.4e\n', mid2_y1_2deg_mgmtTot_YXv(dbi,dbj,dbCrop)) ;
    fprintf('                UnmetB: %0.4e\n', unmetMgmtB_YXv(dbi,dbj,dbCrop)) ;
    fprintf('                 Unmet: %0.4e\n', unmetMgmt_YXv(dbi,dbj,dbCrop)) ;
end

% Sanity check
if any(isnan(mid2_y1_2deg_mgmt_YXv(:)))
    error('Some NaN in mid2_y1_2deg_mgmt_YXv!')
elseif any(isnan(unmetMgmt_YXv(:)))
    error('Some NaN in unmetMgmt_YXv!')
end

end