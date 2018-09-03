function out_YX = PLUMharm_aggregate(in_YX,in_res,out_res)

res_ratio = out_res / in_res ;
if ~isint(res_ratio)
    error('out_res/in_res must be an integer!')
end

tmp = zeros(size(in_YX,1), ...
            size(in_YX,2)/res_ratio, ...
            size(in_YX,3)) ;
for j = 1:res_ratio
    tmp = tmp + in_YX(:,j:res_ratio:end,:) ;
end

out_YX = zeros(size(in_YX,1)/res_ratio, ...
               size(in_YX,2)/res_ratio, ...
               size(in_YX,3)) ;
for i = 1:res_ratio
    out_YX = out_YX + tmp(i:res_ratio:end,:,:) ;
end


end