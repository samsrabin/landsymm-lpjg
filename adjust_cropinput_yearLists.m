function cropinputs = adjust_cropinput_yearLists(cropinputs, yearList)

Nyears = length(yearList) ;

if isfield(cropinputs, 'maps_YXv') || isfield(cropinputs, 'maps_YXvy')

    if ~isfield(cropinputs,'yearList') && size(cropinputs.maps_YXv,4)==1
        % One-year cropinputs file
        cropinputs.maps_YXvy = repmat(cropinputs.maps_YXv,[1 1 1 Nyears]) ;
        cropinputs.yearList = yearList ;
    end
    
    if ~isequal(cropinputs.yearList,yearList)
        if min(cropinputs.yearList) <= min(yearList) && max(cropinputs.yearList) >= max(yearList)
            cropinputs.maps_YXvy = cropinputs.maps_YXvy(:,:,:,cropinputs.yearList>=min(yearList) & cropinputs.yearList<=max(yearList)) ;
    %         cropinputs.yearList = transpose(cropinputs.yearList(1):max(yearList)) ;
            cropinputs.yearList = cropinputs.yearList(cropinputs.yearList>=min(yearList) & cropinputs.yearList<=max(yearList)) ;
        elseif min(cropinputs.yearList) > min(yearList) && max(cropinputs.yearList) == max(yearList) && isequaln(cropinputs.maps_YXvy(:,:,:,1),cropinputs.maps_YXvy(:,:,:,end))
            Nmissing = length(yearList(1):cropinputs.yearList(1)-1) ;
            cropinputs.maps_YXvy = cat(4,repmat(cropinputs.maps_YXvy(:,:,:,1),[1 1 1 Nmissing]),cropinputs.maps_YXvy(:,:,:,cropinputs.yearList>=min(yearList) & cropinputs.yearList>= max(yearList))) ;
            cropinputs.yearList = yearList ;
        elseif min(cropinputs.yearList) > min(yearList) && max(cropinputs.yearList) > max(yearList) && isequaln(cropinputs.maps_YXvy(:,:,:,1),cropinputs.maps_YXvy(:,:,:,end))
            Nmissing = length(yearList(1):cropinputs.yearList(1)-1) ;
            cropinputs.maps_YXvy = cat(4,...
                repmat(cropinputs.maps_YXvy(:,:,:,1),[1 1 1 Nmissing]),...
                cropinputs.maps_YXvy(:,:,:,...
                cropinputs.yearList>=min(yearList) ...
                & cropinputs.yearList<=max(yearList)...
                )) ;
            cropinputs.yearList = yearList ;
        end
        if min(cropinputs.yearList) == min(yearList) + 5 && max(cropinputs.yearList)==max(yearList)
            warning('Adjusting cropinputs to account for 5 years'' padding at beginning.')
            cropinputs.yearList = yearList ;
            cropinputs.maps_YXvy = cat(4,repmat(cropinputs.maps_YXvy(:,:,:,1),[1 1 1 5]),cropinputs.maps_YXvy) ;
        elseif ~isequal(cropinputs.yearList,yearList)
            error('Rework cropinputs so yearLists match!')
        end
    end

elseif isfield(cropinputs, 'garr_xv') || isfield(cropinputs, 'garr_xvy')

    if ~isfield(cropinputs,'yearList') && size(cropinputs.garr_xv,3)==1
        % One-year cropinputs file
        cropinputs.garr_xvy = repmat(cropinputs.garr_xv,[1 1 Nyears]) ;
        cropinputs.yearList = yearList ;
    end

    if ~isequal(cropinputs.yearList,yearList)
        if min(cropinputs.yearList) <= min(yearList) && max(cropinputs.yearList) >= max(yearList)
            cropinputs.garr_xvy = cropinputs.garr_xvy(:,:,cropinputs.yearList>=min(yearList) & cropinputs.yearList<=max(yearList)) ;
    %         cropinputs.yearList = transpose(cropinputs.yearList(1):max(yearList)) ;
            cropinputs.yearList = cropinputs.yearList(cropinputs.yearList>=min(yearList) & cropinputs.yearList<=max(yearList)) ;
        elseif min(cropinputs.yearList) > min(yearList) && max(cropinputs.yearList) == max(yearList) && isequaln(cropinputs.garr_xvy(:,:,1),cropinputs.garr_xvy(:,:,end))
            Nmissing = length(yearList(1):cropinputs.yearList(1)-1) ;
            cropinputs.garr_xvy = cat(3,repmat(cropinputs.garr_xvy(:,:,1),[1 1 Nmissing]),cropinputs.garr_xvy(:,:,cropinputs.yearList>=min(yearList) & cropinputs.yearList>= max(yearList))) ;
            cropinputs.yearList = yearList ;
        elseif min(cropinputs.yearList) > min(yearList) && max(cropinputs.yearList) > max(yearList) && isequaln(cropinputs.garr_xvy(:,:,1),cropinputs.garr_xvy(:,:,end))
            Nmissing = length(yearList(1):cropinputs.yearList(1)-1) ;
            cropinputs.garr_xvy = cat(3,...
                repmat(cropinputs.garr_xvy(:,:,1),[1 1 Nmissing]),...
                cropinputs.garr_xvy(:,:,...
                cropinputs.yearList>=min(yearList) ...
                & cropinputs.yearList<=max(yearList)...
                )) ;
            cropinputs.yearList = yearList ;
        end
        if min(cropinputs.yearList) == min(yearList) + 5 && max(cropinputs.yearList)==max(yearList)
            warning('Adjusting cropinputs to account for 5 years'' padding at beginning.')
            cropinputs.yearList = yearList ;
            cropinputs.garr_xvy = cat(3,repmat(cropinputs.garr_xvy(:,:,1),[1 1 5]),cropinputs.garr_xvy) ;
        elseif ~isequal(cropinputs.yearList,yearList)
            error('Rework cropinputs so yearLists match!')
        end
    end

else
    error('cropinputs does not appear to contain either maps or garrays')
end

end
