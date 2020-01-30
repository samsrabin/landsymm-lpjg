function [xres_out,yres_out] = lpjgu_process_resolution(xres_in,yres_in,...
    in_lons,in_lats,verboseIfNoMat,verbose)

if xres_in>0 && yres_in>0
    xres_out = xres_in ;
    yres_out = yres_in ;
elseif xres_in>0 && ~(yres_in>0)
    % If only xres provided, set yres to same value
    xres_out = xres_in ;
    yres_out = xres_in ;
elseif ~(xres_in>0) && yres_in>0
    % If only yres provided, set xres to same value
    xres_out = yres_in ;
    yres_out = yres_in ;
else
    % Determine X and Y resolution
    if verboseIfNoMat || verbose
        disp('      Determining X and Y resolution...')
    end
    if ~(xres_in>0)
        unique_lats = unique(in_lats) ;
        xres_out = min(abs(unique_lats(1:end-1)-unique_lats(2:end))) ;
    end
    if ~(yres_in>0)
        unique_lons = unique(in_lons) ;
        yres_out = min(abs(unique_lons(1:end-1)-unique_lons(2:end))) ;
    end
    if verboseIfNoMat || verbose
        disp(['      Assuming X res. = ' num2str(xres_out) ', Y res. = ' num2str(yres_out)])
    end
end



end