function [out_rf_YXy, out_ir_YXy] = interp_yield_f(...
    nfert_YX, adjust_bad, ...
    croparea_ha_YXy, Nlevels, ...
    yield_rf_YXyn, yield_ir_YXyn, ...
    warn_missingVal, warn_yieldLowNgtHigh)

if length(Nlevels)~=3
    error('This only works with 3 N levels at the moment!')
end

% Equation corrections/adjustments?
pAlpha_fix = true ;
adj_f = true ;
linearize_where_needed_YXy = true ;

% Convert kgN/ha to kgN/m2
Nlevels = Nlevels * 1e-4 ;

% Calculate f
f = nfert_YX ./ max(Nlevels) ;
Nyears = size(yield_rf_YXyn,3) ;
f_YXy = repmat(f, [1 1 Nyears]) ;

% Other nfert info
% Nlevels_norm = Nlevels / max(Nlevels) ;
Nlevels_norm = (Nlevels - min(Nlevels)) / range(Nlevels) ; % Correction by Peter Alexander 2018-05-19
nfert_YXy = repmat(nfert_YX,[1 1 Nyears]) ;
nfert_YXy1 = nfert_YXy < min(Nlevels) ;
nfert_YXy2 = nfert_YXy >= min(Nlevels) & nfert_YXy < Nlevels(2) ;
nfert_YXy3 = nfert_YXy >= Nlevels(2) & nfert_YXy < Nlevels(3) ;
nfert_YXy4 = nfert_YXy >= Nlevels(3) ;

% Make sure you have at least some yields >0
if isempty(find(yield_rf_YXyn>0,1))
    error('isempty(find(yield_rf_YXyn>0,1))')
end
if isempty(find(yield_ir_YXyn>0,1))
    error('isempty(find(yield_ir_YXyn>0,1))')
end

% Make sure you have no negative yields
% % % if ~isempty(find(yield_rf_YXyn<0,1))
% % %     error('~isempty(find(yield_rf_YXyn<0,1))')
% % % end
% % % if ~isempty(find(yield_ir_YXyn<0,1))
% % %     error('~isempty(find(yield_ir_YXyn<0,1))')
% % % end
if ~isempty(find(yield_rf_YXyn<0,1))
    if length(unique(yield_rf_YXyn(yield_rf_YXyn<0))) > 1
        error('length(unique(yield_rf_YXyn(yield_rf_YXyn<0))) > 1')
    else
        Nbad = length(find(yield_rf_YXyn<0)) ;
        missingval = unique(yield_rf_YXyn(yield_rf_YXyn<0)) ;
        if warn_missingVal
            warning(['Rainfed: It looks like maybe there is a missing value (' num2str(missingval) ') in ' num2str(Nbad) '/' num2str(length(find(~isnan(yield_rf_YXyn)))) ' cells. Setting matching elements to zero.'])
        end
        yield_rf_YXyn(yield_rf_YXyn==Nbad) = 0 ;
    end
end
if ~isempty(find(yield_ir_YXyn<0,1))
    if length(unique(yield_ir_YXyn(yield_ir_YXyn<0))) > 1
        error('length(unique(yield_ir_YXyn(yield_ir_YXyn<0))) > 1')
    else
        Nbad = length(find(yield_ir_YXyn<0)) ;
        missingval = unique(yield_ir_YXyn(yield_ir_YXyn<0)) ;
        warning(['Irrigated: It looks like maybe there is a missing value (' num2str(missingval) ') in ' num2str(Nbad) '/' num2str(length(find(~isnan(yield_ir_YXyn)))) ' cells. Setting matching elements to zero.'])
        yield_ir_YXyn(yield_ir_YXyn==Nbad) = 0 ;
    end
end

% If ignoring cells with yield_Nmid < yield_Nmin or yield_Nmax < yield_Nmid, do that
if ~adjust_bad
    isBad_YXyn = repmat(yield_rf_YXyn(:,:,:,1)>yield_rf_YXyn(:,:,:,2) | yield_rf_YXyn(:,:,:,2)>yield_rf_YXyn(:,:,:,3),[1 1 1 3]) ;
    yield_rf_YXyn(isBad_YXyn) = NaN ;
    isBad_YXyn = repmat(yield_ir_YXyn(:,:,:,1)>yield_ir_YXyn(:,:,:,2) | yield_ir_YXyn(:,:,:,2)>yield_ir_YXyn(:,:,:,3),[1 1 1 3]) ;
    yield_ir_YXyn(isBad_YXyn) = NaN ;
end

% If doing so, adjust yield_Nmid to fix yield_Nmid < yield_Nmin
if adjust_bad
    yield_rf_YXyn(:,:,:,2) = avoid_moreIn_lessOut(yield_rf_YXyn(:,:,:,1), yield_rf_YXyn(:,:,:,2), 'yield_Nmin', 'yield_Nmid', adjust_bad, warn_yieldLowNgtHigh) ;
    yield_ir_YXyn(:,:,:,2) = avoid_moreIn_lessOut(yield_ir_YXyn(:,:,:,1), yield_ir_YXyn(:,:,:,2), 'yield_Nmin', 'yield_Nmid', adjust_bad, warn_yieldLowNgtHigh) ;
    if isempty(find(yield_rf_YXyn>0,1))
        error('isempty(find(yield_rf_YXyn>0,1))')
    end
    if isempty(find(yield_ir_YXyn>0,1))
        error('isempty(find(yield_ir_YXyn>0,1))')
    end
end

% If doing so, adjust yield_Nmax to fix yield_Nmax < yield_Nmid
if adjust_bad
    yield_rf_YXyn(:,:,:,3) = avoid_moreIn_lessOut(yield_rf_YXyn(:,:,:,2), yield_rf_YXyn(:,:,:,3), 'yield_Nmid', 'yield_Nmax', adjust_bad, warn_yieldLowNgtHigh) ;
    yield_ir_YXyn(:,:,:,3) = avoid_moreIn_lessOut(yield_ir_YXyn(:,:,:,2), yield_ir_YXyn(:,:,:,3), 'yield_Nmid', 'yield_Nmax', adjust_bad, warn_yieldLowNgtHigh) ;
    if isempty(find(yield_rf_YXyn>0,1))
        error('isempty(find(yield_rf_YXyn>0,1))')
    end
    if isempty(find(yield_ir_YXyn>0,1))
        error('isempty(find(yield_ir_YXyn>0,1))')
    end
end

% Check that yields are either ignored or at least not decreasing with
% increased nfert
if ~isempty(find(yield_rf_YXyn(:,:,:,1)>yield_rf_YXyn(:,:,:,2),1))
    error('Somehow still some yield_rf_min > yield_rf_mid!')
end
if ~isempty(find(yield_rf_YXyn(:,:,:,2)>yield_rf_YXyn(:,:,:,3),1))
    error('Somehow still some yield_rf_mid > yield_rf_max!')
end
if ~isempty(find(yield_ir_YXyn(:,:,:,1)>yield_ir_YXyn(:,:,:,2),1))
    error('Somehow still some yield_ir_min > yield_ir_mid!')
end
if ~isempty(find(yield_ir_YXyn(:,:,:,2)>yield_ir_YXyn(:,:,:,3),1))
    error('Somehow still some yield_ir_mid > yield_ir_max!')
end

% % Ignore where there is no crop area
% yield_rf_YXyn(repmat(croparea_ha_YXy,[1 1 1 size(yield_rf_YXyn,4)])==0) = NaN ;
% yield_ir_YXyn(repmat(croparea_ha_YXy,[1 1 1 size(yield_ir_YXyn,4)])==0) = NaN ;
% if isempty(find(yield_rf_YXyn>0,1))
%     error('isempty(find(yield_rf_YXyn>0,1))')
% end
% if isempty(find(yield_ir_YXyn>0,1))
%     error('isempty(find(yield_ir_YXyn>0,1))')
% end

% Get curve parameters
A_rf_YXy = yield_rf_YXyn(:,:,:,1) ;
B_rf_YXy = yield_rf_YXyn(:,:,:,end) - A_rf_YXy ;
theNum = yield_rf_YXyn(:,:,:,2)-yield_rf_YXyn(:,:,:,1) ;
theDenom = yield_rf_YXyn(:,:,:,end)-yield_rf_YXyn(:,:,:,1) ;
term2 = (min(Nlevels)-max(Nlevels))/(Nlevels(2)-min(Nlevels)) ;
if isempty(find(theNum>0 & theDenom>0,1))
    error('isempty(find(theNum>0 & theDenom>0,1))')
end
pAlpha_rf_YXy = -log(1-theNum./theDenom) * term2 ;
if ~isreal(pAlpha_rf_YXy)
    error('pAlpha_rf has an imaginary component! Problem with yields not monotonically increasing across Nferts and/or from rainfed to irrig?')
elseif isempty(find(~isnan(pAlpha_rf_YXy),1))
    error('pAlpha_rf is all NaN!')
end

A_ir_YXy = yield_ir_YXyn(:,:,:,1) ;
B_ir_YXy = yield_ir_YXyn(:,:,:,end) - A_ir_YXy ;
pAlpha_ir_YXy = -log(1-(yield_ir_YXyn(:,:,:,2)-yield_ir_YXyn(:,:,:,1)) ...
                ./(yield_ir_YXyn(:,:,:,end)-yield_ir_YXyn(:,:,:,1))) ...
                * (min(Nlevels)-max(Nlevels))/(median(Nlevels)-min(Nlevels)) ;
if ~isreal(pAlpha_ir_YXy)
    error('pAlpha_ir has an imaginary component! Problem with yields not monotonically increasing across Nferts and/or from rainfed to irrig?')
elseif isempty(find(~isnan(pAlpha_ir_YXy),1))
    error('pAlpha_ir is all NaN!')
end

% Make corrections to equation, if doing so
f_YXy2use = f_YXy ;
if adj_f
    f_YXy2use = f_YXy-min(Nlevels_norm) ;
end
if pAlpha_fix
    pAlpha_rf_YXy = -pAlpha_rf_YXy ;
    pAlpha_ir_YXy = -pAlpha_ir_YXy ;
end

% Calculate yield for observed N levels
out_rf_YXy = (A_rf_YXy+B_rf_YXy.*(1-exp(-pAlpha_rf_YXy.*f_YXy2use))) ;
out_ir_YXy = (A_ir_YXy+B_ir_YXy.*(1-exp(-pAlpha_ir_YXy.*f_YXy2use))) ;

% Linearize where needed_YXy, if doing so
if linearize_where_needed_YXy
    % Rainfed, y1==y2
    needed_r1_YXy = yield_rf_YXyn(:,:,:,1)==yield_rf_YXyn(:,:,:,2) & yield_rf_YXyn(:,:,:,2)~=yield_rf_YXyn(:,:,:,3);
    m_YXy = (yield_rf_YXyn(:,:,:,3) - yield_rf_YXyn(:,:,:,2)) / (Nlevels_norm(3) - Nlevels_norm(2)) ;
    b_YXy = yield_rf_YXyn(:,:,:,3) - m_YXy*Nlevels_norm(3) ;
    out_rf_YXy(needed_r1_YXy) = m_YXy(needed_r1_YXy) .* f_YXy(needed_r1_YXy) + b_YXy(needed_r1_YXy) ;
    yield_rf_YXy2 = yield_rf_YXyn(:,:,:,2) ;
    yield_rf_YXy3 = yield_rf_YXyn(:,:,:,3) ;
    out_rf_YXy(needed_r1_YXy & out_rf_YXy<yield_rf_YXy2) = yield_rf_YXy2(needed_r1_YXy & out_rf_YXy<yield_rf_YXy2) ;
    out_rf_YXy(needed_r1_YXy & out_rf_YXy>yield_rf_YXy3) = yield_rf_YXy3(needed_r1_YXy & out_rf_YXy>yield_rf_YXy3) ;
    % Rainfed, y2==y3
    needed_r2_YXy = yield_rf_YXyn(:,:,:,1)~=yield_rf_YXyn(:,:,:,2) & yield_rf_YXyn(:,:,:,2)==yield_rf_YXyn(:,:,:,3);
    m_YXy = (yield_rf_YXyn(:,:,:,2) - yield_rf_YXyn(:,:,:,1)) / (Nlevels_norm(2) - Nlevels_norm(1)) ;
    b_YXy = yield_rf_YXyn(:,:,:,2) - m_YXy*Nlevels_norm(2) ;
    out_rf_YXy(needed_r2_YXy) = m_YXy(needed_r2_YXy) .* f_YXy(needed_r2_YXy) + b_YXy(needed_r2_YXy) ;
    yield_rf_YXy3 = yield_rf_YXyn(:,:,:,3) ;
    out_rf_YXy(needed_r2_YXy & out_rf_YXy<0) = 0 ;
    out_rf_YXy(needed_r2_YXy & out_rf_YXy>yield_rf_YXy3) = yield_rf_YXy3(needed_r2_YXy & out_rf_YXy>yield_rf_YXy3) ;
    % Irrigated, y1==y2
    needed_i1_YXy = yield_ir_YXyn(:,:,:,1)==yield_ir_YXyn(:,:,:,2) & yield_ir_YXyn(:,:,:,2)~=yield_ir_YXyn(:,:,:,3);
    m_YXy = (yield_ir_YXyn(:,:,:,3) - yield_ir_YXyn(:,:,:,2)) / (Nlevels_norm(3) - Nlevels_norm(2)) ;
    b_YXy = yield_ir_YXyn(:,:,:,3) - m_YXy*Nlevels_norm(3) ;
    out_ir_YXy(needed_i1_YXy) = m_YXy(needed_i1_YXy) .* f_YXy(needed_i1_YXy) + b_YXy(needed_i1_YXy) ;
    yield_ir_YXy2 = yield_ir_YXyn(:,:,:,2) ;
    yield_ir_YXy3 = yield_ir_YXyn(:,:,:,3) ;
    out_ir_YXy(needed_i1_YXy & out_ir_YXy<yield_ir_YXy2) = yield_ir_YXy2(needed_i1_YXy & out_ir_YXy<yield_ir_YXy2) ;
    out_ir_YXy(needed_i1_YXy & out_ir_YXy>yield_ir_YXy3) = yield_ir_YXy3(needed_i1_YXy & out_ir_YXy>yield_ir_YXy3) ;
    % Irrigated, y2==y3
    needed_i2_YXy = yield_ir_YXyn(:,:,:,1)~=yield_ir_YXyn(:,:,:,2) & yield_ir_YXyn(:,:,:,2)==yield_ir_YXyn(:,:,:,3);
    m_YXy = (yield_ir_YXyn(:,:,:,2) - yield_ir_YXyn(:,:,:,1)) / (Nlevels_norm(2) - Nlevels_norm(1)) ;
    b_YXy = yield_ir_YXyn(:,:,:,2) - m_YXy*Nlevels_norm(2) ;
    out_ir_YXy(needed_i2_YXy) = m_YXy(needed_i2_YXy) .* f_YXy(needed_i2_YXy) + b_YXy(needed_i2_YXy) ;
    yield_ir_YXy3 = yield_ir_YXyn(:,:,:,3) ;
    out_ir_YXy(needed_i2_YXy & out_ir_YXy<0) = 0 ;
    out_ir_YXy(needed_i2_YXy & out_ir_YXy>yield_ir_YXy3) = yield_ir_YXy3(needed_i2_YXy & out_ir_YXy>yield_ir_YXy3) ;
end
    
% If yields are the same at all N levels, interpolation gets Inf-y. Set to
% one of the yields instead.
all3same_rf_YXy = yield_rf_YXyn(:,:,:,1)==yield_rf_YXyn(:,:,:,2) & yield_rf_YXyn(:,:,:,2)==yield_rf_YXyn(:,:,:,3) ;
yield_rf_YXy1 = yield_rf_YXyn(:,:,:,1) ;
out_rf_YXy(all3same_rf_YXy) = yield_rf_YXy1(all3same_rf_YXy) ;
all3same_ir_YXy = yield_ir_YXyn(:,:,:,1)==yield_ir_YXyn(:,:,:,2) & yield_ir_YXyn(:,:,:,2)==yield_ir_YXyn(:,:,:,3) ;
yield_ir_YXy1 = yield_ir_YXyn(:,:,:,1) ;
out_ir_YXy(all3same_ir_YXy) = yield_ir_YXy1(all3same_ir_YXy) ;

% % Avoid nonsense
out_rf_YXy(out_rf_YXy<0) = 0 ;
out_ir_YXy(out_ir_YXy<0) = 0 ;



% Check that yields are either ignored or at least not decreasing with
% increased nfert
if ~isempty(find(yield_rf_YXyn(:,:,:,1)>yield_rf_YXyn(:,:,:,2),1))
    error('Somehow still some yield_rf_min > yield_rf_mid!')
end
if ~isempty(find(yield_rf_YXyn(:,:,:,2)>yield_rf_YXyn(:,:,:,3),1))
    error('Somehow still some yield_rf_mid > yield_rf_max!')
end
if ~isempty(find(yield_ir_YXyn(:,:,:,1)>yield_ir_YXyn(:,:,:,2),1))
    error('Somehow still some yield_ir_min > yield_ir_mid!')
end
if ~isempty(find(yield_ir_YXyn(:,:,:,2)>yield_ir_YXyn(:,:,:,3),1))
    error('Somehow still some yield_ir_mid > yield_ir_max!')
end



% % % % Check that interpolation is fine
% % % test_rf_YXy = nan(size(out_rf_YXy)) ;
% % % test_rf_YXy(~isnan(out_rf_YXy)) = 0 ;
% % % % Actual nfert less than minimum Nlevel
% % % test_rf_YXy(nfert_YXy1 ...
% % %     & yield_rf_YXyn(:,:,:,1)>0 ...
% % %     & out_rf_YXy >= yield_rf_YXyn(:,:,:,1) ...
% % %     & ~needed_r1_YXy ...
% % %     & ~(needed_r2_YXy ...
% % %     & out_rf_YXy==yield_rf_YXyn(:,:,:,1))) = -1 ;
% % % % % Actual nfert between Nlevels 1 and 2
% % % % test_rf_YX(nfert_YXy2 & any(out_rf_YXy < yield_rf_YXyn(:,:,:,1) | out_rf_YXy >= yield_rf_YXyn(:,:,:,2),3)) = -2 ;
% % % % test_rf_YX(nfert_YXy2 & ~any(out_rf_YXy < yield_rf_YXyn(:,:,:,1) | out_rf_YXy >= yield_rf_YXyn(:,:,:,2),3)) = 0 ;
% % % % % Actual nfert between Nlevels 2 and 3
% % % % test_rf_YX(nfert_YXy3 & any(out_rf_YXy < yield_rf_YXyn(:,:,:,2) | out_rf_YXy >= yield_rf_YXyn(:,:,:,3),3)) = -3 ;
% % % % test_rf_YX(nfert_YXy3 & ~any(out_rf_YXy < yield_rf_YXyn(:,:,:,2) | out_rf_YXy >= yield_rf_YXyn(:,:,:,3),3)) = 0 ;
% % % % Actual nfert >= maximum Nlevel
% % % test_rf_YXy(nfert_YXy4 & out_rf_YXy < yield_rf_YXyn(:,:,:,1)) = -4 ;
% % % % y_out where y_min==y_mid==y_max==0 is taken care of above
% % % test_rf_YXy(all3same_rf_YXy) = 0 ;
% % % test_rf_YX = min(test_rf_YXy,[],3) ;
% % % shademap(-test_rf_YX); caxis([0 4])
% % % 
% % % % Check that interpolation is fine
% % % test_ir_YXy = nan(size(out_ir_YXy)) ;
% % % test_ir_YXy(~isnan(out_ir_YXy)) = 0 ;
% % % % Actual nfert less than minimum Nlevel
% % % test_ir_YXy(nfert_YXy1 ...
% % %     & yield_ir_YXyn(:,:,:,1)>0 ...
% % %     & out_ir_YXy >= yield_ir_YXyn(:,:,:,1) ...
% % %     & ~needed_r1_YXy ...
% % %     & ~(needed_r2_YXy ...
% % %     & out_ir_YXy==yield_ir_YXyn(:,:,:,1))) = -1 ;
% % % % % Actual nfert between Nlevels 1 and 2
% % % % test_ir_YX(nfert_YXy2 & any(out_ir_YXy < yield_ir_YXyn(:,:,:,1) | out_ir_YXy >= yield_ir_YXyn(:,:,:,2),3)) = -2 ;
% % % % test_ir_YX(nfert_YXy2 & ~any(out_ir_YXy < yield_ir_YXyn(:,:,:,1) | out_ir_YXy >= yield_ir_YXyn(:,:,:,2),3)) = 0 ;
% % % % % Actual nfert between Nlevels 2 and 3
% % % % test_ir_YX(nfert_YXy3 & any(out_ir_YXy < yield_ir_YXyn(:,:,:,2) | out_ir_YXy >= yield_ir_YXyn(:,:,:,3),3)) = -3 ;
% % % % test_ir_YX(nfert_YXy3 & ~any(out_ir_YXy < yield_ir_YXyn(:,:,:,2) | out_ir_YXy >= yield_ir_YXyn(:,:,:,3),3)) = 0 ;
% % % % Actual nfert >= maximum Nlevel
% % % test_ir_YXy(nfert_YXy4 & out_ir_YXy < yield_ir_YXyn(:,:,:,1)) = -4 ;
% % % % y_out where y_min==y_mid==y_max==0 is taken care of above
% % % test_ir_YXy(all3same_ir_YXy) = 0 ;
% % % test_ir_YX = min(test_ir_YXy,[],3) ;
% % % shademap(-test_ir_YX); caxis([0 4])

end



function in2_YXv = avoid_moreIn_lessOut(in1_YXv, in2_YXv, in1, in2, ...
    adjust_bad, warn_yieldLowNgtHigh)

if isequaln(in1_YXv,in2_YXv)
    error('Equal arrays provided to avoid_moreIn_lessOut!')
end

tmpDiff = in2_YXv - in1_YXv ;
isBad = in1_YXv > in2_YXv ;
if nnz(isBad)
    Nbad = nnz(isBad) ;
    pctBad = 100*Nbad/length(find(~isnan(tmpDiff))) ;
    tmp = tmpDiff(tmpDiff<0) ;
    tmp95 = prctile(abs(tmp),95) ;
    if adjust_bad
        if warn_yieldLowNgtHigh
            warning([num2str(round(pctBad,1)) ...
                '% of cells have ' in1 ' > ' in2 '! 95th prctile of diff. = ' ...
                num2str(tmp95) '. Increasing ' in2 ' to ' in1 '.'])
        end
        in2_YXv(isBad) = in1_YXv(isBad) ;
    else
        if warn_yieldLowNgtHigh
            warning([num2str(round(pctBad,1)) ...
                '% of cells have ' in1 ' > ' in2 '! 95th prctile of diff. = ' ...
                num2str(tmp95) '. Ignoring these cells.'])
        end
        in1_YXv(isBad) = NaN ;
        in2_YXv(isBad) = NaN ;
    end
end

if isequaln(in1_YXv,in2_YXv)
    error('avoid_moreIn_lessOut resulted in equal arrays!')
end

end