function cropf_out = PLUMharm_combineIrrig(cropf_in, LPJGcrops)

Ncrops_lpjg = length(LPJGcrops) ;

any_irrigated = ~isequal(cropf_in.varNames,LPJGcrops) ;
if ~any_irrigated
    cropf_out = cropf_in ;
else
    
    cropf_out.varNames = LPJGcrops ;
    
    if isfield(cropf_in,'maps_YXv')
        cropf_out.maps_YXv = nan(size(cropf_in.maps_YXv,1), ...
            size(cropf_in.maps_YXv,2), ...
            Ncrops_lpjg) ;
        for c = 1:Ncrops_lpjg
            thisCrop = LPJGcrops{c} ;
            cropf_out.maps_YXv(:,:,c) = sum(cropf_in.maps_YXv(:,:,contains(cropf_in.varNames,thisCrop)),3) ;
        end
    else
        cropf_out.yearList = cropf_in.yearList ;
        cropf_out.maps_YXvy = nan(size(cropf_in.maps_YXvy,1), ...
            size(cropf_in.maps_YXvy,2), ...
            Ncrops_lpjg, size(cropf_in.maps_YXvy,4)) ;
        for c = 1:Ncrops_lpjg
            thisCrop = LPJGcrops{c} ;
            cropf_out.maps_YXvy(:,:,c,:) = sum(cropf_in.maps_YXvy(:,:,contains(cropf_in.varNames,thisCrop),:),3) ;
        end
    end
end

end