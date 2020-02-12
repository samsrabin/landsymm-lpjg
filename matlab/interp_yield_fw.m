function [out_rf_YXv, out_ir_YXv] = interp_yield_fw(...
    f, ...
    croparea_ha_YXy, Nlevels, ...
    yield_rf_YXyn, yield_ir_YXyn)

if length(Nlevels)~=3
    error('This only works with 3 N levels at the moment!')
end

% Do not allow yield_Nmid < yield_Nmin
yield_rf_YXyn(:,:,:,2) = avoid_moreIn_lessOut(yield_rf_YXyn(:,:,:,1), yield_rf_YXyn(:,:,:,2), 'yield_Nmin', 'yield_Nmid') ;
yield_ir_YXyn(:,:,:,2) = avoid_moreIn_lessOut(yield_ir_YXyn(:,:,:,1), yield_ir_YXyn(:,:,:,2), 'yield_Nmin', 'yield_Nmid') ;

% Do not allow yield_Nmax < yield_Nmid
yield_rf_YXyn(:,:,:,3) = avoid_moreIn_lessOut(yield_rf_YXyn(:,:,:,2), yield_rf_YXyn(:,:,:,3), 'yield_Nmid', 'yield_Nmax') ;
yield_ir_YXyn(:,:,:,3) = avoid_moreIn_lessOut(yield_ir_YXyn(:,:,:,2), yield_ir_YXyn(:,:,:,3), 'yield_Nmid', 'yield_Nmax') ;

% Do not allow yield_rf > yield_ir
yield_ir_YXyn = avoid_moreIn_lessOut(yield_rf_YXyn, yield_ir_YXyn, 'yield_rf', 'yield_ir') ;

% Interpolate to get yield at "50%" irrigation (assumes this gets you 80%
% of the way toward the yield at 100% irrigation)
yield_mi_YXyn = yield_rf_YXyn + 0.8*(yield_ir_YXyn-yield_rf_YXyn) ;

% Ignore where there is no crop area
yield_rf_YXyn(repmat(croparea_ha_YXy,[1 1 1 size(yield_rf_YXyn,4)])==0) = NaN ;
yield_mi_YXyn(repmat(croparea_ha_YXy,[1 1 1 size(yield_mi_YXyn,4)])==0) = NaN ;
yield_ir_YXyn(repmat(croparea_ha_YXy,[1 1 1 size(yield_ir_YXyn,4)])==0) = NaN ;

% Get curve parameters
A = yield_rf_YXyn(:,:,:,1) ;
B = yield_rf_YXyn(:,:,:,end) - A ;
C = yield_ir_YXyn(:,:,:,1) - A ;
D = yield_ir_YXyn(:,:,:,end) - C - yield_rf_YXyn(:,:,:,end);
pAlpha_YXy = -log( ...
                  1-( ...
                     yield_rf_YXyn(:,:,:,2)-yield_rf_YXyn(:,:,:,1) ...
                    ) ...
                  ./( ...
                     yield_rf_YXyn(:,:,:,end)-yield_rf_YXyn(:,:,:,1) ...
                    ) ...
                 ) * (min(Nlevels)-max(Nlevels))/(median(Nlevels)-min(Nlevels)) ;
if ~isreal(pAlpha_YXy)
    error('pAlpha has an imaginary component! Problem with yields not monotonically increasing across Nferts and/or from rainfed to irrig?')
end
pBeta_YXy = -log( ...
                 1-( ...
                    yield_mi_YXyn(:,:,:,1)-yield_rf_YXyn(:,:,:,1)...
                   ) ...
                 ./( ...
                    yield_ir_YXyn(:,:,:,1)-yield_rf_YXyn(:,:,:,1) ...
                   )...
                ) * (1-0)/(0-0.5) ;
if ~isreal(pBeta_YXy)
    error('pBeta has an imaginary component! Problem with yields not monotonically increasing across Nferts and/or from rainfed to irrig?')
end

% Calculate yield for observed N levels
w = 0 ;
out_rf_YXv = (A+B.*(1-exp(-pAlpha_YXy.*f))) ...
        + C.*(1-exp(-pBeta_YXy*w)) ...
        + D.*(1-exp(-pAlpha_YXy.*f)).*(1-exp(-pBeta_YXy*w)) ;
w = 1 ;
out_ir_YXv = (A+B.*(1-exp(-pAlpha_YXy.*f))) ...
        + C.*(1-exp(-pBeta_YXy*w)) ...
        + D.*(1-exp(-pAlpha_YXy.*f)).*(1-exp(-pBeta_YXy*w)) ;

end



function in2_YXv = avoid_moreIn_lessOut(in1_YXv, in2_YXv, in1, in2)

tmpDiff = in2_YXv - in1_YXv ;
isBad = in1_YXv > in2_YXv ;
if nnz(isBad)
    Nbad = nnz(isBad) ;
    pctBad = 100*Nbad/length(find(~isnan(tmpDiff))) ;
    tmp = tmpDiff(tmpDiff<0) ;
    tmp95 = prctile(abs(tmp),95) ;
    warning([num2str(round(pctBad,1)) ...
        '% of cells have ' in1 ' > ' in2 '! 95th prctile of diff. = ' ...
        num2str(tmp95) '. Increasing ' in2 ' to ' in1 '.'])
    in2_YXv(isBad) = in1_YXv(isBad) ;
end

end