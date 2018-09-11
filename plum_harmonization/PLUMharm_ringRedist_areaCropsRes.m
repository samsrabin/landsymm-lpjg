function [out_y1_2deg_agri_YXv, out_y1_2deg_nfert_YXv, out_y1_2deg_irrig_YXv] = ...
    PLUMharm_ringRedist_areaCropsRes(...
    mid_y1_2deg_agri_YXv, ...
    vegd_2deg_y1_YX, ...
    total_unmet_agri_YXv, ...
    landArea_2deg_YX, ...
    debugIJ, in_y0orig_2deg, in_y1orig_2deg, out_y0_2deg_agri_YXv, ...
    nonResNtrl_YX, LUnames_agri, ...
    mid_y1_2deg_nfert_YXv, mid_y1_2deg_irrig_YXv)
% Loop through every 2-degree gridcell. If a gridcell has unmet crop
% or pasture, look for place to put this unmet amount in neighboring
% rings, starting with gridcells that are 1 unit away, then 2, etc.
% until all unmet has been displaced to new 2 degree cells. Track the
% displaced crop and pasture and the # of "rings" needed for each
% 2-degree gridcell.

Nagri = size(mid_y1_2deg_agri_YXv,3) ;
hasMgmtInput = ~contains(LUnames_agri,{'PASTURE','ExtraCrop'}) ;
doMgmt = ~isempty(mid_y1_2deg_nfert_YXv) & ~isempty(mid_y1_2deg_irrig_YXv) ;

% Create grids for tracking displaced agriculture
% debugIJ = [999 999] ;
do_debug = ~isempty(debugIJ) ;
if do_debug
    % Just to keep track of original unmet area
    total_unmet_agri_YXv_orig = total_unmet_agri_YXv ;
    % From original code, but not sure how to interpret
    displaced_agri_YXv = zeros(size(mid_y1_2deg_agri_YXv)) ;
    % How many rings did this cell have to draw from or give to?
    agri_rings_YXv = zeros(size(mid_y1_2deg_agri_YXv),'uint8') ;
    % What is the mean distance (in # rings) that area in this cell came
    % from? NOTE that this does not track area given to a cell because it
    % had too much loss: only area given to a cell because a different cell
    % had too much gain!
    meanDist_YXv = zeros(size(mid_y1_2deg_agri_YXv)) ;
    % Keep track of what cells have been processed so far
    thisCell_list = [] ;
end

out_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv ;
if doMgmt
    out_y1_2deg_nfert_YXv = mid_y1_2deg_nfert_YXv ;
    out_y1_2deg_irrig_YXv = mid_y1_2deg_irrig_YXv ;
else
    out_y1_2deg_nfert_YXv = [] ;
    out_y1_2deg_irrig_YXv = [] ;
end

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
            j = 0 ;   % Note that this differs from how it was done in LUH1 code (pasture started with j from this cell's crop redistribution)
            total_unmet_this_YX = total_unmet_agri_YXv(:,:,i) ;
            if do_debug
                displaced_this_YX = displaced_agri_YXv(:,:,i) ;
                this_rings_YX = agri_rings_YXv(:,:,i) ;
                this_meanDist_YX = meanDist_YXv(:,:,i) ;
            else
                displaced_this_YX = [] ;
                this_meanDist_YX = [] ;
            end
            
            out_y1_2deg_this_YX = out_y1_2deg_agri_YXv(:,:,i) ;
            while abs(total_unmet_this_YX(thisCell))>1e-8
                j = j+1 ;
                if j>100
                    error('Possible infinite loop in crop ring adjustments.')
                end
                if do_debug && k+1==debugIJ(1) && m+1==debugIJ(2)
                    fprintf('%s, j = %d, total_unmet_this_YX(thisCell) = %0.4g\n',...
                        LUnames_agri{i},j,total_unmet_this_YX(thisCell)) ;
%                     if strcmp(LUnames_agri{i},'PASTURE')
%                         keyboard
%                     end
                end
                % Calculate the total agricultural (crop+past) area
                out_y1_2deg_agri_YX = sum(out_y1_2deg_agri_YXv,3) ;
                % Get current management inputs
                if doMgmt && hasMgmtInput(i)
                    out_y1_2deg_this_nfert_YX = out_y1_2deg_nfert_YXv(:,:,i) ;
                    out_y1_2deg_this_irrig_YX = out_y1_2deg_irrig_YXv(:,:,i) ;
                else
                    out_y1_2deg_this_nfert_YX = [] ;
                    out_y1_2deg_this_irrig_YX = [] ;
                end
                % Set up rings and indices
                this_rings_YX(thisCell) = j;
                i_k = mod(k + (-j:j),ny) + 1 ;
                i_m = mod(m + (-j:j),nx) + 1 ;
                [I_K,I_M] = meshgrid(i_k,i_m);
                tmp = reshape(cat(2,I_K',I_M'),[],2) ;
                if j==1
                    innerCells = thisCell ;
                else % use thisRing from previous iteration of while loop
                    innerCells = thisRing ;
                end
                thisRing = sub2ind(size(total_unmet_this_YX),tmp(:,1),tmp(:,2)) ;
                [out_y1_2deg_this_YX, total_unmet_this_YX, ...
                    out_y1_2deg_this_nfert_YX, out_y1_2deg_this_irrig_YX, ...
                    displaced_this_YX, this_meanDist_YX] = ...
                    PLUMharm_doRings_areaCropsRes(...
                    out_y1_2deg_this_YX, out_y1_2deg_agri_YX, total_unmet_this_YX, ...
                    out_y1_2deg_this_nfert_YX, out_y1_2deg_this_irrig_YX, ...
                    displaced_this_YX, ...
                    thisCell, thisRing, innerCells, nonResNtrl_YX, this_meanDist_YX) ;
                out_y1_2deg_agri_YXv(:,:,i) = out_y1_2deg_this_YX ;
                if doMgmt && hasMgmtInput(i)
                    out_y1_2deg_nfert_YXv(:,:,i) = out_y1_2deg_this_nfert_YX ;
                    out_y1_2deg_irrig_YXv(:,:,i) = out_y1_2deg_this_irrig_YX ;
                end
                if do_debug
                    meanDist_YXv(:,:,i) = this_meanDist_YX ;
                end
                
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

% figure('Color','w','Position',figurePos) ;
% h = cell(Nagri,1) ;
% for v = 1:Nagri
%     h{v} = subplot_tight(3,3,v,0.025) ;
%     tmp_YX = agri_rings_YXv(:,:,v) ;
%     tmp_YX = single(tmp_YX) ;
%     tmp_YX(landArea_2deg_YX==0) = NaN ;
%     pcolor(tmp_YX); shading flat; axis equal tight off
%     caxis([0 max(agri_rings_YXv(:))]) ;
%     colorbar ;
%     title(LUnames_agri{v})
%     set(gca,'FontSize',14) ;
% end
% figure('Color','w','Position',figurePos) ;
% h = cell(Nagri,1) ;
% for v = 1:Nagri
%     h{v} = subplot_tight(3,3,v,0.025) ;
%     tmp_YX = displaced_agri_YXv(:,:,v) ;
%     tmp_YX(landArea_2deg_YX==0) = NaN ;
%     pcolor(tmp_YX); shading flat; axis equal tight off
%     colorbar ;
%     title(LUnames_agri{v})
%     set(gca,'FontSize',14) ;
% end
% keyboard


end


function [out_YX, total_unmet_YX, ...
    out_nfert_YX, out_irrig_YX, ...
    displaced_YX, this_meanDist_YX] = ...
    PLUMharm_doRings_areaCropsRes( ...
    out_YX, out_agri_YX, total_unmet_YX, ...
    out_nfert_YX, out_irrig_YX, ...
    displaced_YX, thisCell, thisRing, innerCells, ...
    nonResNtrl_YX, this_meanDist_YX)

%%% Internals
% avail_space = "other" area in thisRing
% total_avail_space = avail_space summed over all gridcells in ring
%%% Outputs
% displaced_YX = net area of this land type moved from or to this gridcell.
%                Negative value indicates ???

thisRing_notThisCell = thisRing ;
thisRing_notThisCell(thisRing_notThisCell==thisCell) = [] ;
[~,IA,~] = intersect(thisRing,innerCells,'stable') ;
thisRing_isInnerCell = false(size(thisRing)) ;
thisRing_isInnerCell(IA) = true ;
j = (sqrt(length(thisRing))-1)/2 ;   % What ring are we on?

% if thisCell needs to give THISLU to the rest of the ring
% (i.e., thisCell had more vegetated area than allowed, or more non-NATURAL
% area than allowed)
% WARNING: This makes it so that order of land uses matters. Is there a way
% to avoid this?
if total_unmet_YX(thisCell)>0
    avail_space = nonResNtrl_YX(thisRing) - out_agri_YX(thisRing) ;
    avail_space(avail_space<0) = 0 ;
    total_avail_space = sum(sum(avail_space)) ;
    if total_avail_space>0
      % if there is not enough NATURAL_UNRESERVED in thisRing to absorb the excess THISLU in thisCell
        if total_unmet_YX(thisCell) >= total_avail_space
            to_ring = avail_space ;
            total_unmet_YX(thisCell) = total_unmet_YX(thisCell) - total_avail_space ;
      % else the NATURAL_UNRESERVED in thisRing is sufficient to absorb the excess THISLU in thisCell
        else
            to_ring = total_unmet_YX(thisCell) * avail_space ./ total_avail_space ;
            total_unmet_YX(thisCell) = 0;
        end
        out_thisRing_new = out_YX(thisRing) + to_ring ;
        to_ring_notThisCell = to_ring(thisRing_notThisCell~=thisCell) ;
        % Check that you do not have to_ring too much >0 in interior of ring.
        % THIS CAN HAPPEN BECAUSE AVAIL_SPACE CAN BE >0 IN INTERIOR. I DO
        % NOT UNDERSTAND WHY. HOPEFULLY ALWAYS TINY.
        if any(to_ring(thisRing_isInnerCell)>1e-3)
            warning('How do you have to_ring>1e-3 (%0.4f) in cell(s) not on ring perimeter? (thisCell = %d)\n', ...
                    max(to_ring(thisRing_isInnerCell)), thisCell)
        end
        out_thisRingNotThisCell_new = out_YX(thisRing_notThisCell) + to_ring_notThisCell ;
        % This cell is giving area to other cells in ring. So here, we
        % update the recipient cells' management inputs by merging in those
        % of the donor cell.
        if ~isempty(out_nfert_YX)
            now_weighted = out_nfert_YX(thisRing_notThisCell) .* out_YX(thisRing_notThisCell)./out_thisRingNotThisCell_new ;
            new_weighted = out_nfert_YX(thisCell)              * to_ring_notThisCell         ./out_thisRingNotThisCell_new ;
            now_weighted(out_thisRingNotThisCell_new==0) = 0 ;
            new_weighted(out_thisRingNotThisCell_new==0) = 0 ;
            out_nfert_YX(thisRing_notThisCell) = now_weighted + new_weighted ;
        end
        if ~isempty(out_irrig_YX)
            now_weighted = out_irrig_YX(thisRing_notThisCell) .* out_YX(thisRing_notThisCell)./out_thisRingNotThisCell_new ;
            new_weighted = out_irrig_YX(thisCell)              * to_ring_notThisCell         ./out_thisRingNotThisCell_new ;
            now_weighted(out_thisRingNotThisCell_new==0) = 0 ;
            new_weighted(out_thisRingNotThisCell_new==0) = 0 ;
            out_irrig_YX(thisRing_notThisCell) = now_weighted + new_weighted ;
        end
        % Update debugging info
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) + to_ring;
            % This cell is giving area to other cells in ring. So here, we
            % update this_meanDist_YX, which keeps track of the mean
            % distance (in # rings) that area in each cell came from.
            now_weighted = this_meanDist_YX(thisRing) .* out_YX(thisRing)./out_thisRing_new ;
            new_weighted =                          j  * to_ring./out_thisRing_new ;
            now_weighted(out_thisRing_new==0) = 0 ;
            new_weighted(out_thisRing_new==0) = 0 ;
            this_meanDist_YX(thisRing) = now_weighted + new_weighted ;
        end
        out_YX(thisRing) = out_thisRing_new ;
    end
% elseif thisCell needs THISLU donated from rest of thisRing
% (i.e., thisCell didn't have enough thisLU to satisfy PLUM-specified loss)
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