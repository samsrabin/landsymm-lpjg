function out_array = PLUMharm_aggregate_mgmt(in_array,inArea_array,in_res,out_res)

Ndims_in = length(find(size(in_array)>1)) ;
if Ndims_in>4
    error('Adapt code to work with >4 dimensions!')
end

res_ratio = out_res / in_res ;
if ~isint(res_ratio)
    error('out_res/in_res must be an integer!')
end

% Where there is no area, set mgmt to NaN
in_array(inArea_array==0) = NaN ;

% Convert NaNs to zeros
in_array(isnan(in_array)) = 0 ;

tmp = zeros(size(in_array,1), ...
    size(in_array,2)/res_ratio, ...
    size(in_array,3),size(in_array,4)) ;
tmpArea = tmp ;
for j = 1:res_ratio
    % Area in new column
    newArea = inArea_array(:,j:res_ratio:end,:,:) ;
    % The current value of mgmt input, weighted by area totaled so far as
    % fraction of so-far plus new area
    tmpWtd = tmp .* tmpArea./(tmpArea+newArea) ;
    % The value of mgmt input in the new column, weighted by area in new
    % column as fraction of so-far plus new area
    newWtd = in_array(:,j:res_ratio:end,:,:) .* newArea./(tmpArea+newArea) ;
    % Where no so-far or new area, set old and new values to zero.
    tmpWtd(tmpArea+newArea==0) = 0 ;
    newWtd(tmpArea+newArea==0) = 0 ;
    % Update aggregated management input and area
    tmp = tmpWtd + newWtd ;
    tmpArea = tmpArea + newArea ;
end

out_array = zeros(size(in_array,1)/res_ratio, ...
    size(in_array,2)/res_ratio, ...
    size(in_array,3),size(in_array,4)) ;
outArea_array = out_array ;
for i = 1:res_ratio
    % Area in new row
    newArea = tmpArea(i:res_ratio:end,:,:,:) ;
    % The current value of mgmt input, weighted by area totaled so far as
    % fraction of so-far plus new area
    out2wtd = out_array .* outArea_array./(outArea_array+newArea) ;
    % The value of mgmt input in the new row, weighted by area in new
    % row as fraction of so-far plus new area
    new2wtd = tmp(i:res_ratio:end,:,:,:) .* newArea./(outArea_array+newArea) ;
    % Where no so-far or new area, set old and new values to zero.
    out2wtd(outArea_array+newArea==0) = 0 ;
    new2wtd(outArea_array+newArea==0) = 0 ;
    % Update aggregated management input and area
    out_array = out2wtd + new2wtd;
    outArea_array = outArea_array + newArea ;
end


% Where there is no area, set mgmt to NaN
out_array(outArea_array==0) = NaN ;


end