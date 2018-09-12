function out_array = PLUMharm_aggregate(in_array,in_res,out_res)

Ndims_in = length(find(size(in_array)>1)) ;
if Ndims_in>4
    error('Adapt code to work with >4 dimensions!')
end

res_ratio = out_res / in_res ;
if ~isint(res_ratio)
    error('out_res/in_res must be an integer!')
end

tmp = zeros(size(in_array,1), ...
    size(in_array,2)/res_ratio, ...
    size(in_array,3),size(in_array,4)) ;
for j = 1:res_ratio
    tmp = tmp + in_array(:,j:res_ratio:end,:,:) ;
end

out_array = zeros(size(in_array,1)/res_ratio, ...
    size(in_array,2)/res_ratio, ...
    size(in_array,3),size(in_array,4)) ;
for i = 1:res_ratio
    out_array = out_array + tmp(i:res_ratio:end,:,:,:) ;
end


end