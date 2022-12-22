function [struct_out, outFrom] = plum7ify_struct(...
    struct_in, cropTypes_conv_plum, cropTypes_conv_lpjg, ...
    varNames_lpjg, varNames_plum, method, outFrom_in)

Ncrops_lpjg = length(varNames_lpjg) ;
Ncrops_plum = length(varNames_plum) ;

if length(cropTypes_conv_lpjg) ~= length(cropTypes_conv_plum)
    error('length(cropTypes_conv_lpjg) ~= length(cropTypes_conv_plum)')
end

if Ncrops_lpjg < Ncrops_plum
    lpjg_more = true ;
elseif Ncrops_lpjg > Ncrops_plum
    lpjg_more = false ;
else
    error('Ncrops_lpjg == Ncrops_plum')
end


Nyears = length(struct_in.yearList) ;

struct_out.list_to_map = struct_in.list_to_map ;
struct_out.varNames = varNames_plum ;
struct_out.yearList = struct_in.yearList ;
struct_out.maps_YXvy = nan(size(struct_in.maps_YXvy,1),size(struct_in.maps_YXvy,2),Ncrops_plum,Nyears) ;

if strcmp(method,'max')
    outFrom.list_to_map = struct_in.list_to_map ;
    outFrom.varNames = varNames_plum ;
    outFrom.yearList = struct_in.yearList ;
    outFrom.maps_YXvy = zeros(size(struct_in.maps_YXvy,1),size(struct_in.maps_YXvy,2),Ncrops_plum,Nyears,'uint8') ;
else
    outFrom = struct() ;
end

for c = 1:Ncrops_plum
    if lpjg_more
        % Do straight 1:1 conversion from LPJ-GUESS types to PLUM types.
        % Ignores "method."
        thisCrop_lpjg = cropTypes_conv_lpjg{c} ;
        thisCrop_plum = cropTypes_conv_plum{c} ;
% % %         thisInd_lpjg = getInds(varNames_lpjg,thisCrop_lpjg,1) ;
        thisInd_lpjg = getInds(struct_in.varNames,thisCrop_lpjg,1) ;
        thisInd_plum = getInds(varNames_plum,thisCrop_plum,1) ;
        struct_out.maps_YXvy(:,:,thisInd_plum,:) = struct_in.maps_YXvy(:,:,thisInd_lpjg,:) ;
    else
        % Combine LPJ-GUESS types into PLUM types
        thisCrop_plum = varNames_plum{c} ;
        inds_conv = getInds(cropTypes_conv_plum,thisCrop_plum,Inf) ;
        thisCrop_lpjg = cropTypes_conv_lpjg(inds_conv) ;
% % %         [~,inds_lpjg] = intersect(varNames_lpjg,thisCrop_lpjg) ;
        [~,inds_lpjg] = intersect(struct_in.varNames,thisCrop_lpjg) ;
        tmp = struct_in.maps_YXvy(:,:,inds_lpjg,:) ;
        if strcmp(method,'sum')
            struct_out.maps_YXvy(:,:,c,:) = sum(tmp,3) ;
        elseif strcmp(method,'wtd_av')
            thisSum_repd = repmat(sum(tmp,3),[1 1 size(tmp,3) 1]) ;
            weights = tmp ./ thisSum_repd ;
            weights(thisSum_repd==0) = 0 ;
            struct_out.maps_YXvy(:,:,c,:) = sum(tmp .* weights,3) ;
        elseif strcmp(method,'max')
            [struct_out.maps_YXvy(:,:,c,:),I] = max(tmp,[],3) ;
            if isempty(find(~isnan(mean(struct_out.maps_YXvy(:,:,c,:),4)),1))
                % Leave outFrom.maps_YXvy(:,:,c,:) as all zeroes
                warning(['No land found for PLUM type ' thisCrop_plum '.'])
            else
                tmpOutFrom = outFrom.maps_YXvy(:,:,c,:) ;
                for i = 1:length(inds_lpjg)
                    thisLPJGind = inds_lpjg(i) ;
                    tmpOutFrom(~isnan(struct_out.maps_YXvy(:,:,c,:)) & I==i) = thisLPJGind ;
                end
                outFrom.maps_YXvy(:,:,c,:) = tmpOutFrom ;
                clear tmpOutFrom
            end
        elseif strcmp(method,'specified') ;
            thisSpecified_YX_y = outFrom_in.maps_YXvy(:,:,c,:) ;
            unique_thisSpecified = unique(thisSpecified_YX_y(thisSpecified_YX_y>0)) ;
            Nthis = length(unique_thisSpecified) ;
            if Nthis==0
                warning(['Skipping ' thisCrop_plum '.'])
            elseif Nthis==1
                struct_out.maps_YXvy(:,:,c,:) = struct_in.maps_YXvy(:,:,unique_thisSpecified,:) ;
            else
                thisOut = struct_out.maps_YXvy(:,:,c,:) ;
                for i = 1:Nthis
                    thisI = unique_thisSpecified(i) ;
                    thisIn =  struct_in.maps_YXvy(:,:,thisI,:) ;
                    thisOut(thisSpecified_YX_y==thisI) = thisIn(thisSpecified_YX_y==thisI) ;
                    clear thisI thisIn
                end
                struct_out.maps_YXvy(:,:,c,:) = thisOut ;
                clear thisOut
            end
            clear thisSpecified_YX_y unique_thisSpecified
        else
            error(['method (' method ') not recognized!'])
        end
    end
    clear tmp
end


end



function thisInd =  getInds(varNames,thisCrop,maxFound)

thisInd = find(strcmp(varNames,thisCrop)) ;

if isempty(thisInd)
    error('isempty(thisInd)') ;
elseif length(thisInd) > maxFound
    error(['length(thisInd) > ' num2str(maxFound)]) ;
end

end