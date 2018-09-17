function out_YXv = PLUMharm_fixTinyNegs(in_YXv, landArea_YXv)
% Rounding errors can result in small negative values. Fix.

Nlu = size(in_YXv, 3) ;

tolerance = 1 ;   % m2
if min(in_YXv(:)) < -tolerance
    error('"Large" negative values of area!')
end

out_YXv = in_YXv ;
if any(out_YXv(:)<0)
    out_YXv = out_YXv ./ landArea_YXv ;
    out_YXv(out_YXv<0) = 0 ;
    tmp_YX = sum(out_YXv,3) ;
    j = 0 ;
    while max(max(abs(tmp_YX-1))) > 3*eps
        j = j+1 ;
        if j>50
            error('Possible infinite loop in fixing tiny negative areas!')
        end
        out_YXv = out_YXv ./ repmat(tmp_YX,[1 1 Nlu]) ;
        tmp_YX = sum(out_YXv,3) ;
    end
    out_YXv = out_YXv .* landArea_YXv ;
    out_YXv(landArea_YXv==0) = 0 ;
    
end



end