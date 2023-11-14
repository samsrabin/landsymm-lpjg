function debug_global_areas(...
    mapsA_YXv, mapsB_YXv, ...
    msg_intro, msgA, msgB, ...
    LUnames, isCrop, isAgri, dbCrop, thisYear)

msgWidth = 45 ;
units = 'm2' ;

disp(' ')
disp(msg_intro)
if ~isempty(dbCrop)
    tmp0 = sum(sum(mapsA_YXv(:,:,dbCrop))) ;
    tmp1 = sum(sum(mapsB_YXv(:,:,dbCrop))) ;
    fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob %s area (%s):', msgA, thisYear-1, LUnames{dbCrop}, units),msgWidth,'left'), tmp0)
    fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob %s area (%s):', msgB, thisYear,   LUnames{dbCrop}, units),msgWidth,'left'), tmp1)
    fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)
end

% Cropland 
tmp0 = sum(sum(sum(mapsA_YXv(:,:,isCrop)))) ;
tmp1 = sum(sum(sum(mapsB_YXv(:,:,isCrop)))) ;
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob crop area (%s):', msgA, thisYear-1, units),msgWidth,'left'), tmp0)
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob crop area (%s):', msgB, thisYear,   units),msgWidth,'left'), tmp1)
fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)

% Pasture 
tmp0 = sum(sum(sum(mapsA_YXv(:,:,strcmp(LUnames,'PASTURE'))))) ;
tmp1 = sum(sum(sum(mapsB_YXv(:,:,strcmp(LUnames,'PASTURE'))))) ;
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob past area (%s):', msgA, thisYear-1, units),msgWidth,'left'), tmp0)
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob past area (%s):', msgB, thisYear,   units),msgWidth,'left'), tmp1)
fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)

% Agricultural 
tmp0 = sum(sum(sum(mapsA_YXv(:,:,isAgri)))) ;
tmp1 = sum(sum(sum(mapsB_YXv(:,:,isAgri)))) ;
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob agri area (%s):', msgA, thisYear-1, units),msgWidth,'left'), tmp0)
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob agri area (%s):', msgB, thisYear,   units),msgWidth,'left'), tmp1)
fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)

% Natural 
tmp0 = sum(sum(sum(mapsA_YXv(:,:,strcmp(LUnames,'NATURAL'))))) ;
tmp1 = sum(sum(sum(mapsB_YXv(:,:,strcmp(LUnames,'NATURAL'))))) ;
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob ntrl area (%s):', msgA, thisYear-1, units),msgWidth,'left'), tmp0)
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob ntrl area (%s):', msgB, thisYear,   units),msgWidth,'left'), tmp1)
fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)

% Vegetated 
tmp0 = sum(sum(sum(mapsA_YXv(:,:,~strcmp(LUnames,'BARREN'))))) ;
tmp1 = sum(sum(sum(mapsB_YXv(:,:,~strcmp(LUnames,'BARREN'))))) ;
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob vegd area (%s):', msgA, thisYear-1, units),msgWidth,'left'), tmp0)
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob vegd area (%s):', msgB, thisYear,   units),msgWidth,'left'), tmp1)
fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)

% Barren 
tmp0 = sum(sum(sum(mapsA_YXv(:,:,strcmp(LUnames,'BARREN'))))) ;
tmp1 = sum(sum(sum(mapsB_YXv(:,:,strcmp(LUnames,'BARREN'))))) ;
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob bare area (%s):', msgA, thisYear-1, units),msgWidth,'left'), tmp0)
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob bare area (%s):', msgB, thisYear,   units),msgWidth,'left'), tmp1)
fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)

% All land
tmp0 = sum(mapsA_YXv(:)) ;
tmp1 = sum(mapsB_YXv(:)) ;
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob land area (%s):', msgA, thisYear-1, units),msgWidth,'left'), tmp0)
fprintf('%s %0.4e\n', pad(sprintf('%s_%d glob land area (%s):', msgB, thisYear,   units),msgWidth,'left'), tmp1)
fprintf('%s %0.4e\n\n', pad(sprintf('diff (%s):',units),msgWidth,'left'), tmp1-tmp0)


end