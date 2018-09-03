function out_YXv = PLUMharm_aggregate_mgmt(in_YXv,inArea_YXv,in_res,out_res)

res_ratio = out_res / in_res ;
if ~isint(res_ratio)
    error('out_res/in_res must be an integer!')
end

% Where there is no area, set mgmt to NaN
in_YXv(inArea_YXv==0) = NaN ;

% Convert NaNs to zeros
in_YXv(isnan(in_YXv)) = 0 ;

tmp = zeros(size(in_YXv,1), ...
            size(in_YXv,2)/res_ratio, ...
            size(in_YXv,3)) ;
tmpArea = tmp ;
for j = 1:res_ratio
    % Area in new column
    newArea = inArea_YXv(:,j:res_ratio:end,:) ;
    % The current value of mgmt input, weighted by area totaled so far as
    % fraction of so-far plus new area
    tmpWtd = tmp .* tmpArea./(tmpArea+newArea) ;
    % The value of mgmt input in the new column, weighted by area in new
    % column as fraction of so-far plus new area
    newWtd = in_YXv(:,j:res_ratio:end,:) .* newArea./(tmpArea+newArea) ;
    % Where no so-far or new area, set old and new values to zero.
    tmpWtd(tmpArea+newArea==0) = 0 ;
    newWtd(tmpArea+newArea==0) = 0 ;
    % Update aggregated management input and area
    tmp = tmpWtd + newWtd ;
    tmpArea = tmpArea + newArea ;
end

out_YXv = zeros(size(in_YXv,1)/res_ratio, ...
               size(in_YXv,2)/res_ratio, ...
               size(in_YXv,3)) ;
outArea_YX = out_YXv ;
for i = 1:res_ratio
    % Area in new row
    newArea = tmpArea(i:res_ratio:end,:,:) ;
    % The current value of mgmt input, weighted by area totaled so far as
    % fraction of so-far plus new area
    out2wtd = out_YXv .* outArea_YX./(outArea_YX+newArea) ;
    % The value of mgmt input in the new row, weighted by area in new
    % row as fraction of so-far plus new area
    new2wtd = tmp(i:res_ratio:end,:,:) .* newArea./(outArea_YX+newArea) ;
    % Where no so-far or new area, set old and new values to zero.
    out2wtd(outArea_YX+newArea==0) = 0 ;
    new2wtd(outArea_YX+newArea==0) = 0 ;
    % Update aggregated management input and area
    out_YXv = out2wtd + new2wtd;
    outArea_YX = outArea_YX + newArea ;
end

% Where there is no area, set mgmt to NaN
out_YXv(outArea_YX==0) = NaN ;


end