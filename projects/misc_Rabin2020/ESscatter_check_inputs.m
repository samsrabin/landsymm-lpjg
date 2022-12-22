function ESscatter_check_inputs(...
    mapstructX, thisVarX, ...
    mapstructY, thisVarY, ...
    countries_YX)

if ~any(strcmp(mapstructX.varNames, thisVarX))
    error('Variable %s not found in mapstructX!', thisVarX)
elseif ~any(strcmp(mapstructY.varNames, thisVarY))
    error('Variable %s not found in mapstructY!', thisVarY)
end
if ~isfield(mapstructX, 'maps_YXvyB')
    error('maps_YXvyB not found in mapstructX!')
elseif ~isfield(mapstructX, 'maps_YXvyr')
    error('maps_YXvyr not found in mapstructX!')
elseif ~isfield(mapstructY, 'maps_YXvyB')
    error('maps_YXvyB not found in mapstructY!')
elseif ~isfield(mapstructY, 'maps_YXvyr')
    error('maps_YXvyr not found in mapstructY!')
end
size_XB = size(mapstructX.maps_YXvyB(:,:,1,:)) ;
size_Xr = size(mapstructX.maps_YXvyr(:,:,1,:,:)) ;
if ~isequal(size_XB(1:2), size(countries_YX))
    error('~isequal(size_XB(1:2), size(countries_YX))')
elseif ~isequal(size_Xr(1:2), size(countries_YX))
    error('~isequal(size_Xr(1:2), size(countries_YX))')
end
if ~isequal(size_XB, size(mapstructY.maps_YXvyB(:,:,1,:)))
    error('~isequal(size_XB, size(mapstructY.maps_YXvyB(:,:,1,:)))')
elseif ~isequal(size_Xr, size(mapstructY.maps_YXvyr(:,:,1,:,:)))
    error('~isequal(size_Xr, size(mapstructY.maps_YXvyr(:,:,1,:,:)))')
end


end