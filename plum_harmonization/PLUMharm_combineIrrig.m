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
            thisCropI = [thisCrop 'i'] ;
            is_thisCrop  = strcmp(cropf_in.varNames,thisCrop) ;
            is_thisCropI = strcmp(cropf_in.varNames,thisCropI) ;
            if any(is_thisCropI)
                cropf_out.maps_YXv(:,:,c) = cropf_in.maps_YXv(:,:,is_thisCrop) ...
                    + cropf_in.maps_YXv(:,:,is_thisCropI) ;
            else
                cropf_out.maps_YXv(:,:,c) = cropf_in.maps_YXv(:,:,is_thisCrop) ;
            end
        end
    else
        cropf_out.yearList = cropf_in.yearList ;
        cropf_out.maps_YXvy = nan(size(cropf_in.maps_YXvy,1), ...
            size(cropf_in.maps_YXvy,2), ...
            Ncrops_lpjg, size(cropf_in.maps_YXvy,4)) ;
        for c = 1:Ncrops_lpjg
            thisCrop = LPJGcrops{c} ;
            thisCropI = [thisCrop 'i'] ;
            is_thisCrop  = strcmp(cropf_in.varNames,thisCrop) ;
            is_thisCropI = strcmp(cropf_in.varNames,thisCropI) ;
            if any(is_thisCropI)
                cropf_out.maps_YXvy(:,:,c,:) = cropf_in.maps_YXvy(:,:,is_thisCrop,:) ...
                    + cropf_in.maps_YXvy(:,:,is_thisCropI) ;
            else
                cropf_out.maps_YXvy(:,:,c,:) = cropf_in.maps_YXvy(:,:,is_thisCrop,:) ;
            end
        end
    end
end


end