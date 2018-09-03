function out_YXv = PLUMharm_aggregate(in_YXv,in_res,out_res)

res_ratio = out_res / in_res ;
if ~isint(res_ratio)
    error('out_res/in_res must be an integer!')
end

tmp = zeros(size(in_YXv,1), ...
            size(in_YXv,2)/res_ratio, ...
            size(in_YXv,3)) ;
for j = 1:res_ratio
    tmp = tmp + in_YXv(:,j:res_ratio:end,:) ;
end

out_YXv = zeros(size(in_YXv,1)/res_ratio, ...
               size(in_YXv,2)/res_ratio, ...
               size(in_YXv,3)) ;
for i = 1:res_ratio
    out_YXv = out_YXv + tmp(i:res_ratio:end,:,:) ;
end


end