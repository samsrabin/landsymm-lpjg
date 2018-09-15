function debug_global_areas(...
    mapsA_YXv, mapsB_YXv, ...
    msg_intro, msgA, msgB, ...
    LPJGcrops, isAgri, dbCrop, thisYear)

msgWidth = 45 ;
units = 'm2' ;

disp(' ')
disp(msg_intro)
if ~isempty(dbCrop)
    tmp0 = sum(sum(mapsA_YXv(:,:,dbCrop))) ;
    tmp1 = sum(sum(mapsB_YXv(:,:,dbCrop))) ;
    fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob %s area (%s):', msgA, thisYear-1, LPJGcrops{dbCrop}, units),msgWidth,'left'), tmp0)
    fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob %s area (%s):', msgB, thisYear,   LPJGcrops{dbCrop}, units),msgWidth,'left'), tmp1)
    fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)
end

% Agricultural 
tmp0 = sum(sum(sum(mapsA_YXv(:,:,isAgri)))) ;
tmp1 = sum(sum(sum(mapsB_YXv(:,:,isAgri)))) ;
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob agri area (%s):', msgA, thisYear-1, units),msgWidth,'left'), tmp0)
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob agri area (%s):', msgB, thisYear,   units),msgWidth,'left'), tmp1)
fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)

% All land
tmp0 = sum(mapsA_YXv(:)) ;
tmp1 = sum(mapsB_YXv(:)) ;
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob land area (%s):', msgA, thisYear-1, units),msgWidth,'left'), tmp0)
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob land area (%s):', msgB, thisYear,   units),msgWidth,'left'), tmp1)
fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)

% % Vegetated
% tmp = sum(sum(out_y0_2deg_vegd_YX(:,:,1:end-1))) ;
% fprintf('%s_%d_2deg glob vegd area (%s): %0.4e\n', thisYear-1, units, tmp)
% out_y1_2deg_agri_YX = sum(out_y1_2deg_agri_YXv,3) ;
% tmp = max(out_y1_2deg_agri_YX(:) - out_y0_2deg_vegd_YX(:)) ;
% fprintf('Max agri-vegd area (km2): %0.1f\n', tmp)
% clear tmp*


end