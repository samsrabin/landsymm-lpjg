function land_area_YX_out = aggregate_land_area(land_area_YX_in,xres_out,yres_out)

% Get info first
xres_in = 360 / size(land_area_YX_in,2) ;
yres_in = 180 / size(land_area_YX_in,1) ;
if xres_out < xres_in
    error(['xres_out ' num2str(xres_out) ' is too small given xres_in=' num2str(xres_in).'])
    
elseif yres_out < yres_in
    error(['yres_out ' num2str(yres_out) ' is too small given yres_in=' num2str(yres_in) '.'])
    
elseif xres_in~=xres_out || yres_in~=yres_out
    
    xratio = xres_out / xres_in ;
    yratio = yres_out / yres_in ;
    if ~isint(xratio) || ~isint(yratio)
        error(['xratio (' num2str(xratio) ') and yratio (' num2str(xratio) ') must both be integers!'])
    end
    
    % Aggregate X
    land_area_YX_tmp = zeros(size(land_area_YX_in,1),size(land_area_YX_in,2)/xratio) ;
    for i = 1:xratio
        land_area_YX_tmp = land_area_YX_tmp + land_area_YX_in(:,i:xratio:end) ;
    end
    
    % Check
    globalland_area_in = sum(sum(land_area_YX_in)) ;
    if abs(sum(sum(land_area_YX_tmp)) - globalland_area_in)>1e-6
        error('Error in X aggregation.')
    end
    
    % Aggregate Y
    land_area_YX_out = zeros(size(land_area_YX_in,1)/yratio,size(land_area_YX_tmp,2)) ;
    for j = 1:yratio
        land_area_YX_out = land_area_YX_out + land_area_YX_tmp(j:yratio:end,:) ;
    end
    
    % Check
    if abs(sum(sum(land_area_YX_out)) - globalland_area_in)>1e-6
        error('Error in Y aggregation.')
    end
    
else
    land_area_YX_out = land_area_YX_in ;
    
end


end