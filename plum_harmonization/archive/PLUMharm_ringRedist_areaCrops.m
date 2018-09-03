function out_y1_2deg_agri_YXv = ...
    PLUMharm_ringRedist_areaCrops(...
    mid_y1_2deg_agri_YXv, ...
    vegd_2deg_y1_YX, ...
    total_unmet_agri_YXv, ...
    landArea_2deg_YX, ...
    do_debug, in_y0orig_2deg, in_y1orig_2deg, out_y0_2deg_agri_YXv)
% Loop through every 2-degree gridcell. If a gridcell has unmet crop
% or pasture, look for place to put this unmet amount in neighboring
% rings, starting with gridcells that are 1 unit away, then 2, etc.
% until all unmet has been displaced to new 2 degree cells. Track the
% displaced crop and pasture and the # of "rings" needed for each
% 2-degree gridcell.

Nagri = size(mid_y1_2deg_agri_YXv,3) ;

% Create grids for tracking displaced agriculture
debug_ij = [Inf Inf] ;
if do_debug
    displaced_agri_YXv = zeros(size(mid_y1_2deg_agri_YXv)) ;
    agri_rings_YXv = zeros(size(mid_y1_2deg_agri_YXv),'uint8') ;
    thisCell_list = [] ;
end

out_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv ;

ny = size(landArea_2deg_YX,1) ;
nx = size(landArea_2deg_YX,2) ;
ks = 0:(ny-1) ;
ms = 0:(nx-1) ;

for k = ks
    for m = ms
        thisCell = sub2ind(size(landArea_2deg_YX),k+1,m+1) ;
        if do_debug
            thisCell_list = [thisCell_list thisCell] ;
        end
        
        % Do it for each
        for i = 1:Nagri
            j = 0 ;   % Note that this differs from how it was done in LUH1 code (pasture started with j of last crop iteration)
            total_unmet_this_YX = total_unmet_agri_YXv(:,:,i) ;
            if do_debug
                displaced_this_YX = displaced_agri_YXv(:,:,i) ;
                this_rings_YX = agri_rings_YXv(:,:,i) ;
            else
                displaced_this_YX = [] ;
            end
            
            out_y1_2deg_this_YX = out_y1_2deg_agri_YXv(:,:,i) ;
            while abs(total_unmet_this_YX(thisCell))>1e-8
                j = j+1 ;
                if j>100
                    error('Possible infinite loop in crop ring adjustments.')
                end
                % Calculate the total agricultural (crop+past) area
                out_y1_2deg_agri_YX = sum(out_y1_2deg_agri_YXv,3) ;
                % Set up rings and indices
                this_rings_YX(thisCell) = j;
                i_k = mod(k + (-j:j),ny) + 1 ;
                i_m = mod(m + (-j:j),nx) + 1 ;
                [I_K,I_M] = meshgrid(i_k,i_m);
                tmp = reshape(cat(2,I_K',I_M'),[],2) ;
                thisRing = sub2ind(size(total_unmet_this_YX),tmp(:,1),tmp(:,2)) ;
                [out_y1_2deg_this_YX, total_unmet_this_YX, displaced_this_YX] = ...
                    PLUMharm_doRings(out_y1_2deg_this_YX, out_y1_2deg_agri_YX, total_unmet_this_YX, displaced_this_YX, ...
                    thisCell, thisRing, vegd_2deg_y1_YX) ;
                out_y1_2deg_agri_YXv(:,:,i) = out_y1_2deg_this_YX ;
                
            end % while loop
            total_unmet_agri_YXv(:,:,i) = total_unmet_this_YX ;
            if do_debug
                displaced_agri_YXv(:,:,i) = displaced_this_YX ;
                agri_rings_YXv(:,:,i) = this_rings_YX ;
                diffH_YXv = out_y1_2deg_agri_YXv - out_y0_2deg_agri_YXv ;
                diffO_YXv = diffH_YXv ;
            end
        end % for loop
        
    end % for m
end % for k

if do_debug
%     keyboard
end


end


function [out_YX, total_unmet_YX, displaced_YX] = ...
    PLUMharm_doRings( ...
    out_YX, out_agri_YX, total_unmet_YX, displaced_YX, ...
    thisCell, thisRing, landarea_notBare_YX)

%%% Internals
% avail_space = "other" area in thisRing
% total_avail_space = avail_space summed over all gridcells in ring
%%% Outputs
% displaced_YX = net area of this land type moved from or to this gridcell.
%                Negative value indicates ???

% if thisCell needs to donate THISLU to the rest of the ring
if total_unmet_YX(thisCell)>0
    avail_space = landarea_notBare_YX(thisRing) - out_agri_YX(thisRing) ;
    total_avail_space = sum(sum(avail_space .* (avail_space>0))) ;
  % if thisCell has more THISLU than can be drawn from othr in thisRing
    if total_unmet_YX(thisCell) >= total_avail_space
        out_YX(thisRing) = out_YX(thisRing) + avail_space .* (avail_space>0);
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) + avail_space .* (avail_space>0);
        end
        total_unmet_YX(thisCell) = total_unmet_YX(thisCell) - total_avail_space ;
  % else the othr in thisRing is sufficient to absorb the extra THISLU from thisCell
    else
        out_YX(thisRing) = out_YX(thisRing) + total_unmet_YX(thisCell) * avail_space ./ total_avail_space .* (avail_space>0);
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) + total_unmet_YX(thisCell) * avail_space ./ total_avail_space .* (avail_space>0);
        end
        total_unmet_YX(thisCell) = 0;
    end
    
% elseif thisCell needs THISLU donated from rest of thisRing
elseif total_unmet_YX(thisCell)<0
    avail_space = out_YX(thisRing);
    total_avail_space = sum(sum(avail_space .* (avail_space>0))) ;
  % if the THISLU available in thisRing is not sufficient to satisfy the demand of thisCell
    if total_unmet_YX(thisCell) <= -total_avail_space
        out_YX(thisRing) = out_YX(thisRing) - avail_space .* (avail_space>0);
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) - avail_space .* (avail_space>0);
        end
        total_unmet_YX(thisCell) = total_unmet_YX(thisCell) + total_avail_space ;
  % elseif the THISLU available in the rest of thisRing is sufficient to satisfy the demand of thisCell
    elseif total_unmet_YX(thisCell) > -total_avail_space
        out_YX(thisRing) = out_YX(thisRing) + total_unmet_YX(thisCell) * avail_space ./ total_avail_space .* (avail_space>0);
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) + total_unmet_YX(thisCell) * avail_space ./ total_avail_space .* (avail_space>0);
        end
        total_unmet_YX(thisCell) = 0;
    end
end




end