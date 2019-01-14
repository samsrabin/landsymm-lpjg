function cropinputs = adjust_cropinput_yearLists(cropinputs, yearList)

Nyears = length(yearList) ;

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


end