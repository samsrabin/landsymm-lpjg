function out_YXv = PLUMharm_fixTinyNegs(in_YXv, landArea_YXv)
% Rounding errors can result in small negative values. Fix.

Nlu = size(in_YXv, 3) ;

tolerance = 1 ;   % m2
if min(in_YXv(:)) < -tolerance
    error('"Large" negative values of area! Min %0.1e', min(in_YXv(:)))
end

out_YXv = in_YXv ;

isBad_YX = any(out_YXv<0,3) ;
if any(any(isBad_YX))
    Nbad = length(find(isBad_YX)) ;
    isBad_YXv = repmat(isBad_YX,[1 1 Nlu]) ;
    tmp = out_YXv(isBad_YXv) ;
    tmp_landArea = landArea_YXv(isBad_YXv) ;
    tmp = tmp ./ tmp_landArea ;
    if any(isnan(tmp))
        error('How do you have NaN in tmp?')
    end
    tmp_xv = reshape(tmp,[Nbad Nlu]) ;
    tmp_xv(tmp_xv<0) = 0 ;
    tmp_xSum = sum(tmp_xv,2) ;
    j = 0 ;
    while max(abs(tmp_xSum-1)) > 3*eps
        j = j+1 ;
        if j>50
            error('Possible infinite loop in fixing tiny negative areas!')
        end
        tmp_xv = tmp_xv ./ repmat(tmp_xSum,[1 Nlu]) ;
        tmp_xSum = sum(tmp_xv,2) ;
    end
    out_YXv(isBad_YXv) = tmp_xv(:) .* tmp_landArea ;
end


end