function out_YXv = PLUMharm_interpolateMgmt(...
    in_YXv, area_YXv, landArea_YX, LPJGcrops, inpaint_method)

% Allow negative values?
allow_neg = false ;

% Allow values higher than seen in original?
% Default=TRUE for compatibility with remap_v4.m.
allow_tooHigh = true ;

% Loop through crops and interpolate
out_YXv = nan(size(in_YXv)) ;
Ncrops_lpjg = length(LPJGcrops) ;
for c = 1:Ncrops_lpjg
    if strcmp(LPJGcrops{c},'ExtraCrop') || strcmp(LPJGcrops{c},'Miscanthus')
        out_YXv(:,:,c) = zeros(size(landArea_YX)) ;
        continue
    end
    out_YXv(:,:,c) = inpaint_nans(double(in_YXv(:,:,c)), inpaint_method) ;
end

% Set management to zero where there is no land (probably unnecessary)
out_YXv(repmat(landArea_YX==0,[1 1 Ncrops_lpjg])) = 0 ;

% Do not allow negative values
if ~allow_neg && any(out_YXv(:)<0)
%     warning('Setting negative members of out_YXv to zero.')
    out_YXv(out_YXv<0) = 0 ;
end

% Do not allow values in excess of max observed for each crop
if ~allow_tooHigh
    max_mgmt = max(max(in_YXv,[],1),[],2) ;
    max_mgmt(nansum(nansum(area_YXv,1),2)==0) = 0 ;
    max_mgmt_YXv = repmat(max_mgmt,[size(landArea_YX) 1]) ;
    if any(out_YXv(:) > max_mgmt_YXv(:))
%         warning('Limiting some members of out_YXv to maximum observed before interpolation.')
        is2much = out_YXv > max_mgmt_YXv ;
        out_YXv(is2much) = max_mgmt_YXv(is2much) ;
    end
end



end