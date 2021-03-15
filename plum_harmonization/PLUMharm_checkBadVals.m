function PLUMharm_checkBadVals(in_lu_YXv, in_nfert_YXv, in_irrig_YXv, ...
    landArea_YX, LUnames, msg, outPrec)

% LU
isBare = strcmp(LUnames,'BARREN') ;
if ~isempty(in_lu_YXv) && any(isBare)
    vegd_YX = sum(in_lu_YXv(:,:,~isBare),3) ;
end

if ~isempty(in_lu_YXv)
    in_luFrac_YXv = in_lu_YXv./repmat(landArea_YX, [1 1 size(in_lu_YXv,3)]) ;
    if any(isnan(in_lu_YXv(:)))
        error('NaN(s) in %s LU maps!', msg)
    elseif min(in_luFrac_YXv(:))<0
        error('Negative value(s) in %s LU maps! Min %0.1e', msg, min(in_luFrac_YXv(:)))
    elseif max(in_luFrac_YXv(:))>1+10^-outPrec
        error('Value(s) >1 in %s LU maps! Max %0.1e', msg, max(in_luFrac_YXv(:))-1)
    elseif max(max(sum(in_luFrac_YXv,3)))>1+10^-outPrec
        error('Sum(s) >1 in %s LU maps! Max overage %0.1e', msg, max(max(sum(in_luFrac_YXv,3)))-1)
    elseif any(isBare) && any(any(vegd_YX==0 & landArea_YX>0))
        error('Zero vegetation in cell(s) with land! Max land area %0.1e, total %0.1e.', max(landArea_YX(vegd_YX==0)), sum(landArea_YX(vegd_YX==0)))
    end
end

% Nfert
if ~isempty(in_nfert_YXv)
    if any(isnan(in_nfert_YXv(:)))
        error('NaN(s) in %s nfert maps!', msg)
    elseif any(in_nfert_YXv(:)<0)
        error('Negative value(s) in %s nfert maps! Min %0.1e', msg, min(in_nfert_YXv(:)))
    end
end

% Irrigation
if ~isempty(in_irrig_YXv)
    if any(isnan(in_irrig_YXv(:)))
        error('NaN(s) in %s irrig maps!', msg)
    elseif any(in_irrig_YXv(:)<0)
        error('Negative value(s) in %s irrig maps! Min %0.1e', msg, min(in_irrig_YXv(:)))
    elseif max(in_irrig_YXv(:))>1+10^-outPrec
        error('Value(s) >1 in %s irrig maps! Max overage %0.1e', msg, max(in_irrig_YXv(:))-1)
    end
end


end