function debug_global_mgmts(...
    mapsA_mgmt_YXv, mapsB_mgmt_YXv, unmet_mgmt_YXv, ...
    mapsA_area_YXv, mapsB_area_YXv, ...
    units, convfact, ...
    msg_intro, msgA, msgB, msg_mgmt, ...
    LPJGcrops, dbCrop, thisYear)

msgWidth = 45 ;

mapsA_YXv = convfact * (mapsA_mgmt_YXv .* mapsA_area_YXv) ;
mapsB_YXv = convfact * (mapsB_mgmt_YXv .* mapsB_area_YXv + unmet_mgmt_YXv) ;

disp(' ')
disp(msg_intro)

% Crop of interest
if ~isempty(dbCrop)
    tmp0 = sum(sum(mapsA_YXv(:,:,dbCrop))) ;
    tmp1 = sum(sum(mapsB_YXv(:,:,dbCrop))) ;
    fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob %s %s (%s):', msgA, thisYear-1, LPJGcrops{dbCrop}, msg_mgmt, units),msgWidth,'left'), tmp0)
    fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob %s %s (%s):', msgB, thisYear,   LPJGcrops{dbCrop}, msg_mgmt, units),msgWidth,'left'), tmp1)
    fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)
end

% All crops 
tmp0 = sum(mapsA_YXv(:)) ;
tmp1 = sum(mapsB_YXv(:)) ;
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob total %s (%s):', msgA, thisYear-1, msg_mgmt, units),msgWidth,'left'), tmp0)
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob total %s (%s):', msgB, thisYear,   msg_mgmt, units),msgWidth,'left'), tmp1)
fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)


end