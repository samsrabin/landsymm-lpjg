function [out_y1_2deg_agri_YXv, out_ntrl_YX] = ...
    PLUMharm_ringRedist_areaCropsRes(...
    mid_y1_2deg_agri_YXv, ...
    total_unmet_agri_YXv, ...
    debugIJ, conserv_tol_pct, check_name, ...
    in_y0_area_YXv, in_y1_area_YXv, ...
    out_y0_2deg_agri_YXv, out_y0_2deg_ntrl_YX, ...
    resArea_2deg_YX, LUnames_agri)
% Loop through every 2-degree gridcell. If a gridcell has unmet crop
% or pasture, look for place to put this unmet amount in neighboring
% rings, starting with gridcells that are 1 unit away, then 2, etc.
% until all unmet has been displaced to new 2 degree cells. Track the
% displaced crop and pasture and the # of "rings" needed for each
% 2-degree gridcell.

out_ntrl_YX = out_y0_2deg_ntrl_YX ;

YXv_dims = size(mid_y1_2deg_agri_YXv) ;
YX_dims = YXv_dims(1:2) ;
Nagri = YXv_dims(3) ;

% Create grids for tracking displaced agriculture
do_debug = ~isempty(debugIJ) ;
thisCell_ofInt = Inf ;
if do_debug
    % From original code, but not sure how to interpret
    displaced_agri_YXv = zeros(YXv_dims) ;
    % How many rings did this cell have to draw from or give to?
    agri_rings_YXv = zeros(YXv_dims,'uint8') ;
    % What is the mean distance (in # rings) that area in this cell came
    % from? NOTE that this does not track area given to a cell because it
    % had too much loss: only area given to a cell because a different cell
    % had too much gain!
    meanDist_YXv = zeros(YXv_dims) ;
    % thisCell by any other name
    thisCell_ofInt = sub2ind(YX_dims,debugIJ(1),debugIJ(2)) ;
end

out_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv ;

ny = YX_dims(1) ;
nx = YX_dims(2) ;
ks = 0:(ny-1) ;
ms = 0:(nx-1) ;

is_done_YXv = false(YXv_dims) ;

Nloops = 0 ;
while any(~is_done_YXv(:))
    Nloops = Nloops + 1 ;
    if Nloops>100
        error('Possible infinite loop (top "while", %s).', check_name)
    end
    is_done_YXv_begin = is_done_YXv ;
    any_skipped_v = false(Nagri,1) ;
    
    for k = ks
        for m = ms
            thisCell = sub2ind(YX_dims,k+1,m+1) ;
            
            % It's possible that this cell DIDN'T have any area
            % available for agri expansion (and thus needed to send
            % agri expansion to a different cell), but it previously
            % donated existing agri area to a different cell, and so
            % now (some of) its own donation demand can be considered
            % satisfied. NOTE that this is different from how it was
            % done in original: Originally, a ring was started even if
            % a cell had enough room to provide for all of its "unmet."
            nonResNtrl_thisCell = max(out_ntrl_YX(thisCell) - resArea_2deg_YX(thisCell), 0) ;
            avail_thisCell = nonResNtrl_thisCell - sum(out_y1_2deg_agri_YXv(k+1,m+1,:),3) ;
            unmet_thisCell = total_unmet_agri_YXv(k+1,m+1,:) ;
            if avail_thisCell>0 && any(unmet_thisCell>0)
                unmet_thisCell_wherePos = unmet_thisCell .* (unmet_thisCell>0) ;
                unmetReduction_tot = min(avail_thisCell,sum(unmet_thisCell_wherePos)) ;
                unmetReduction = unmetReduction_tot .* unmet_thisCell_wherePos./sum(unmet_thisCell_wherePos) ;
                out_y1_2deg_agri_YXv(k+1,m+1,:) = out_y1_2deg_agri_YXv(k+1,m+1,:) + unmetReduction ;
                total_unmet_agri_YXv(k+1,m+1,:) = total_unmet_agri_YXv(k+1,m+1,:) - unmetReduction ;
                out_ntrl_YX(thisCell) = out_ntrl_YX(thisCell) - sum(unmetReduction) ;
                if out_ntrl_YX(thisCell)<0
                    error('out_ntrl_YX(thisCell)<0')
                end
            end

            % Do it for each
            for i = 1:Nagri
                
                % Skip if this crop for this cell is already done
                if is_done_YXv(k+1,m+1,i)
                    continue
                end
                
                % If this cell has no unmet, skip.
                if total_unmet_agri_YXv(k+1,m+1,i)==0
                    is_done_YXv(k+1,m+1,i) = true ;
                    continue
                end
                
                j = 0 ;   % Note that this differs from how it was done in LUH1 code (pasture started with j from this cell's crop redistribution)
                if do_debug
                    displaced_this_YX = displaced_agri_YXv(:,:,i) ;
                    this_rings_YX = agri_rings_YXv(:,:,i) ;
                    this_meanDist_YX = meanDist_YXv(:,:,i) ;
                else
                    displaced_this_YX = [] ;
                    this_meanDist_YX = [] ;
                end
                
                out_this_YX = out_y1_2deg_agri_YXv(:,:,i) ;
                total_unmet_thisCell = total_unmet_agri_YXv(k+1,m+1,i) ;
                while abs(total_unmet_thisCell)>1e-8
                    
                    j = j+1 ;
                    if j>100
                        error('Possible infinite loop in crop ring adjustments: [%d %d %d]',k+1,m+1,i)
                    end
                    if do_debug && k+1==debugIJ(1) && m+1==debugIJ(2)
                        fprintf('%s, j = %d, total_unmet_this_YX(thisCell) = %0.4g\n',...
                            LUnames_agri{i},j,total_unmet_thisCell) ;
                    end
                    
                    % Calculate the total agricultural (crop+past) area
                    out_y1_2deg_agri_YX = sum(out_y1_2deg_agri_YXv,3) ;
                    
                    % Set up rings and indices
                    this_rings_YX(thisCell) = j;
                    i_k = mod(k + (-j:j),ny) + 1 ;
                    i_m = mod(m + (-j:j),nx) + 1 ;
                    [I_K,I_M] = meshgrid(i_k,i_m);
                    tmp = reshape(cat(2,I_K',I_M'),[],2) ;
                    thisRing = sub2ind(YX_dims,tmp(:,1),tmp(:,2)) ;
                    % If ring gets big enough, you can have duplicates.
                    % Remove them.
                    if j*2+1 > min(YX_dims)
                        thisRing = unique(thisRing) ;
                    end
                    
                    % If no available, continue (i.e., expand ring and
                    % start over).
                    if total_unmet_thisCell < 0
                        avail_thisRing = out_this_YX(thisRing) ;
                        total_avail_ring = sum(avail_thisRing) ;
                    elseif total_unmet_thisCell > 0
                        nonResNtrl_thisRing = max(out_ntrl_YX(thisRing) - resArea_2deg_YX(thisRing), 0) ;
                        avail_thisRing = nonResNtrl_thisRing - out_y1_2deg_agri_YX(thisRing) ;
                        avail_thisRing(avail_thisRing<0) = 0 ;
                        total_avail_ring = sum(avail_thisRing) ;
                    else
                        error('How was this not skipped?')
                    end
                    if total_avail_ring==0
                        continue
                    end
                    
                    % Try to distribute to / take from ring
                    if total_avail_ring==0
                        error('How??')
                    end
                    [out_this_YX, out_ntrl_YX, total_unmet_thisCell, ...
                        displaced_this_YX, this_meanDist_YX] = ...
                        PLUMharm_doRings_areaCropsRes(...
                        out_this_YX, out_y1_2deg_agri_YX, total_unmet_thisCell, displaced_this_YX, ...
                        thisRing, out_ntrl_YX, resArea_2deg_YX, ...
                        this_meanDist_YX, thisCell_ofInt) ;
                    out_y1_2deg_agri_YXv(:,:,i) = out_this_YX ;
                    meanDist_YXv(:,:,i) = this_meanDist_YX ;
                    
                    % Debugging: Check that area changes are (mostly) conserved
                    if do_debug
                        out_y1_2deg_agri_YXv(:,:,i) = out_this_YX ;
                        total_unmet_agri_YXv(k+1,m+1,i) = total_unmet_thisCell ;
                        bad = PLUMharm_checkCons_area(...
                            out_y0_2deg_agri_YXv, out_y1_2deg_agri_YXv, ...
                            in_y0_area_YXv, in_y1_area_YXv, ...
                            total_unmet_agri_YXv, LUnames_agri, conserv_tol_pct, ...
                            check_name) ;
                        if bad==1
                            keyboard
                        end
                    end
                    
                end % while loop
                if do_debug
                    displaced_agri_YXv(:,:,i) = displaced_this_YX ;
                    agri_rings_YXv(:,:,i) = this_rings_YX ;
                end
                
                % Record this cell as done
                is_done_YXv(k+1,m+1,i) = true ;
                if j>0
                    out_y1_2deg_agri_YXv(:,:,i) = out_this_YX ;
                    total_unmet_agri_YXv(k+1,m+1,i) = total_unmet_thisCell ;
                end
            end % for loop
            
        end % for m
    end % for k
    
    if isequal(is_done_YXv_begin,is_done_YXv)
        error('Not enough mgmt in the world!')
    end
    
    if any(~is_done_YXv(:))
        for c = 1:Nagri
            if any_skipped_v(c)
                skipped_cells = find(~is_done_YXv(:,:,c)) ;
                thisSkipped = skipped_cells(1) ;
                [thisSkipped_I, thisSkipped_J] = ind2sub([ny nx],thisSkipped) ;
                fprintf('      %s %s: Skipped cell %d (%d,%d) and %d others.\n', ...
                    LUnames_agri{c}, check_name, thisSkipped, thisSkipped_I, thisSkipped_J, length(skipped_cells)-1) ;
            end
        end
        disp('      Trying again.')
        pause(0.1)
    end
    
end % while


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


