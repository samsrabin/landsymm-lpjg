function [out_y1_2deg_mgmt_YXv, total_unmet_mgmt_YXv, notEnough] = ...
    PLUMharm_ringRedist_mgmt(...
    mid_y1_2deg_mgmt_YXv, out_y1_2deg_cropArea_YXv, ...
    total_unmet_mgmt_YXv, max_mgmt, ...
    LPJGcrops, debugIJ, ...
    out_y0_mgmt, out_y0_area_YXv, ...
    in_y0_mgmt, in_y0_area_YXv, ...
    in_y1_mgmt, in_y1_area_YXv, ...
    conserv_tol_pct, check_name, dbCrop)
% Loop through every 2-degree gridcell. If a gridcell has unmet crop
% or pasture, look for place to put this unmet amount in neighboring
% rings, starting with gridcells that are 1 unit away, then 2, etc.
% until all unmet has been displaced to new 2 degree cells. Track the
% displaced crop and pasture and the # of "rings" needed for each
% 2-degree gridcell.

if ~isequal(size(mid_y1_2deg_mgmt_YXv), size(out_y1_2deg_cropArea_YXv))
    error('mid_y1_2deg_mgmt_YXv and out_y1_2deg_cropArea_YXv must be same size!')
end

YX_dims = [size(mid_y1_2deg_mgmt_YXv,1) size(mid_y1_2deg_mgmt_YXv,2)] ;
YXv_dims = size(mid_y1_2deg_mgmt_YXv) ;
Ncrops = YXv_dims(3) ;

% Get TOTAL management inputs
mid_y1_2deg_mgmtTot_YXv = mid_y1_2deg_mgmt_YXv .* out_y1_2deg_cropArea_YXv ;

% Create grids for tracking displaced agriculture
% debugIJ = [999 999] ;
do_debug = ~isempty(debugIJ) ;
thisCell_ofInt = Inf ;
if do_debug
    % Just to keep track of original unmet inputs
% % %     total_unmet_mgmt_YXv_orig = total_unmet_mgmt_YXv ;
    % From original code, but not sure how to interpret
    displaced_mgmt_YXv = zeros(YXv_dims) ;
    % How many rings did this cell have to draw from or give to?
    agri_rings_YXv = zeros(YXv_dims,'uint8') ;
    % What is the mean distance (in # rings) that area in this cell came
    % from? NOTE that this does not track area given to a cell because it
    % had too much loss: only area given to a cell because a different cell
    % had too much gain!
    meanDist_YXv = zeros(YXv_dims) ;
    % Get k and m of interest
    dbk = debugIJ(1) - 1 ;
    dbm = debugIJ(2) - 1 ;
    % thisCell by any other name
    thisCell_ofInt = sub2ind(size(out_y0_area_YXv(:,:,1)),debugIJ(1),debugIJ(2)) ;
end

out_y1_2deg_mgmt_YXv = mid_y1_2deg_mgmt_YXv ;
out_y1_2deg_mgmtTot_YXv = mid_y1_2deg_mgmtTot_YXv ;

ny = YX_dims(1) ;
nx = YX_dims(2) ;
ks = 0:(ny-1) ;
ms = 0:(nx-1) ;

is_done_YXv = false(YXv_dims) ;
notEnough = false(Ncrops,1) ;

Nloops = 0 ;
while any(~is_done_YXv(:))
    Nloops = Nloops + 1 ;
    if Nloops>100
        error('Possible infinite loop (top "while", %s).', check_name)
    end
    is_done_YXv_begin = is_done_YXv ;
    any_skipped_v = false(Ncrops,1) ;
    check_tooMuch = true(Ncrops,1) ;
    for k = ks
        for m = ms
            
            thisCell = sub2ind(YX_dims,k+1,m+1) ;
            if do_debug
                if k==dbk && m==dbm
%                     keyboard
                end
            end
            
            % Do it for each
            for i = 1:Ncrops
                
                % Skip if this crop for this cell is already done
                if is_done_YXv(k+1,m+1,i)
                    continue
                end
                
                if do_debug && k==dbk && m==dbm && i==dbCrop
                    x = 1 ;
                end
                
                % If there's too much desired mgmt loss for the entire
                % world to handle, set all mgmt to zero, and skip from now
                % on. Only need to do this once, for first non-skipped
                % grid cell in each "while" loop.
                
                if check_tooMuch(i)
                    check_tooMuch(i) = false ;
                    
                    % Get available maps
                    out_y1_2deg_thisTot_YX = out_y1_2deg_mgmt_YXv(:,:,i) .* out_y1_2deg_cropArea_YXv(:,:,i) ;
                    max_mgmtTot_YX = max_mgmt(i) * out_y1_2deg_cropArea_YXv(:,:,i) ;
                    total_unmet_this_YX = total_unmet_mgmt_YXv(:,:,i) ;
                    
                    % Is there too much positive or negative unmet?
                    total_pos_unmet = sum(total_unmet_this_YX(total_unmet_this_YX>0)) ;
                    total_neg_unmet = sum(total_unmet_this_YX(total_unmet_this_YX<0)) ;
                    total_avail_world_toTake = sum(sum(max_mgmtTot_YX - out_y1_2deg_thisTot_YX)) ;
                    total_avail_world_toGive = sum(sum(out_y1_2deg_thisTot_YX)) ;
                    total_avail_world_toTake = total_avail_world_toTake - total_neg_unmet ;
                    total_avail_world_toGive = total_avail_world_toGive + total_pos_unmet ;
                    tooMuchPos = total_pos_unmet > total_avail_world_toTake ;
                    tooMuchNeg = -total_neg_unmet > total_avail_world_toGive ;
                    if tooMuchPos && tooMuchNeg
                        error('Is it possible to have too much positive AND negative unmet?? Think about it. If yes, then code for it.')
                    elseif tooMuchPos
                        warning(['There is not enough mgmt headroom for %s in the world to absorb the positive unmet '...
                            '(difference = %0.4e). Setting all %s mgmt to max.'], ...
                            LPJGcrops{i}, total_pos_unmet - total_avail_world_toTake, LPJGcrops{i})
                        tmpOut_YX = max_mgmtTot_YX ./ out_y1_2deg_cropArea_YXv(:,:,i) ;
                        tmpOut_YX(out_y1_2deg_cropArea_YXv(:,:,i)==0) = 0 ;
                        out_y1_2deg_mgmt_YXv(:,:,i) = tmpOut_YX ;
                        total_unmet_mgmt_YXv(:,:,i) = 0 ; % There may be a nicer way to handle this.
                        is_done_YXv(:,:,i) = true ;
                        notEnough(i) = true ;
                        continue
                    elseif tooMuchNeg
                        warning(['There is not enough mgmt applied to %s in the world to absorb the negative unmet ' ...
                            '(difference = %0.4e). Setting all %s mgmt to zero.'], ...
                            LPJGcrops{i}, -total_neg_unmet - total_avail_world_toGive, LPJGcrops{i})
                        out_y1_2deg_mgmt_YXv(:,:,i) = 0 ;
                        total_unmet_mgmt_YXv(:,:,i) = 0 ; % There may be a nicer way to handle this.
                        is_done_YXv(:,:,i) = true ;
                        notEnough(i) = true ;
                        continue
                    end
                    clear out_y1_2deg_thisTot_YX   % Just to make sure you don't try and use it later.
                end
                
                if do_debug
                    displaced_this_YX = displaced_mgmt_YXv(:,:,i) ;
                    this_rings_YX = agri_rings_YXv(:,:,i) ;
                    this_meanDist_YX = meanDist_YXv(:,:,i) ;
                else
                    displaced_this_YX = [] ;
                    this_meanDist_YX = [] ;
                end
                
                if do_debug && k==dbk && m==dbm && i==dbCrop
                    fprintf('%s, j = %d, total_unmet_mgmt_YXv(k+1,m+1,i) =\t%0.4e\n',...
                        LPJGcrops{i},0,total_unmet_mgmt_YXv(k+1,m+1,i)) ;
                end
                
                % It's possible that this cell WAS at its max mgmt in
                % mid_y1_2deg_mgmt_YXthis (and thus needed to send mgmt to a
                % different cell), but it previously donated existing mgmt to a
                % different cell, and so now (some of) its own donation demand
                % can be considered satisfied. NOTE that this is different from
                % how it was done in original (and how it's currently done in
                % ringRedist for area): Originally, a ring was started even if
                % a cell had enough room to provide for all of its "unmet."
                max_mgmtTot_thisCell = max_mgmt(i) * out_y1_2deg_cropArea_YXv(k+1,m+1,i) ;
                out_total_thisCell = out_y1_2deg_mgmt_YXv(k+1,m+1,i) * out_y1_2deg_cropArea_YXv(k+1,m+1,i) ;
                if total_unmet_mgmt_YXv(k+1,m+1,i)>0 && out_total_thisCell < max_mgmtTot_thisCell
                    unmetReduction = min(max_mgmtTot_thisCell-out_total_thisCell,total_unmet_mgmt_YXv(k+1,m+1,i)) ;
                    out_y1_2deg_mgmt_YXv(k+1,m+1,i) = (out_total_thisCell + unmetReduction) / out_y1_2deg_cropArea_YXv(k+1,m+1,i) ;
                    total_unmet_mgmt_YXv(k+1,m+1,i) = total_unmet_mgmt_YXv(k+1,m+1,i) - unmetReduction ;
                    if do_debug && k==dbk && m==dbm && i==dbCrop
                        total_avail_world = sum(sum(out_y1_2deg_mgmt_YXv(:,:,i))) ;
                        fprintf('            now total_avail_world =\t\t%0.4e\n', total_avail_world) ;
                        fprintf('        now total_unmet_mgmt_YXv(k+1,m+1,i) =\t%0.4e\n',...
                        LPJGcrops{i},0,total_unmet_mgmt_YXv(k+1,m+1,i)) ;
                    end
                end
                
                % If this cell has no unmet, skip.
                if total_unmet_mgmt_YXv(k+1,m+1,i)==0
                    is_done_YXv(k+1,m+1,i) = true ;
                    continue
                end
                
                % If this cell needs to take mgmt from elsewhere (i.e., its
                % unmet is negative), make sure there is enough mgmt in the
                % world for that to work. If not enough, skip and come back
                % later.
                out_y1_2deg_thisTot_YX = out_y1_2deg_mgmt_YXv(:,:,i) .* out_y1_2deg_cropArea_YXv(:,:,i) ;
                if total_unmet_mgmt_YXv(k+1,m+1,i) < 0
                    total_avail_world = sum(sum(out_y1_2deg_thisTot_YX)) ;
                    if -total_unmet_mgmt_YXv(k+1,m+1,i) > total_avail_world
                        if ~any_skipped_v(i)
                            any_skipped_v(i) = true ;
%                             warning('Skipping cell %d %s (taking too much).', thisCell, LPJGcrops{i}) ;
                        end
                        if do_debug && k==dbk && m==dbm && i==dbCrop
                            fprintf('                total_avail_world =\t\t%0.4e\n', total_avail_world) ;
                        end
                        continue
                    end
                    
                % If this cell needs to give mgmt away (i.e., its unmet is
                % positive), make sure there is enough headroom for this
                % mgmt in the world for that to work. If not enough, skip
                % and come back later.
                elseif total_unmet_mgmt_YXv(k+1,m+1,i) > 0
                    max_mgmtTot_YX = max_mgmt(i) * out_y1_2deg_cropArea_YXv(:,:,i) ;
                    total_avail_world = sum(sum(max_mgmtTot_YX - out_y1_2deg_thisTot_YX)) ;
                    if total_unmet_mgmt_YXv(k+1,m+1,i) > total_avail_world
                        if ~any_skipped_v(i)
                            any_skipped_v(i) = true ;
%                             warning('Skipping cell %d %s (donating too much).', thisCell, LPJGcrops{i}) ;
                        end
                        if do_debug && k==dbk && m==dbm && i==dbCrop
                            fprintf('                total_avail_world =\t\t%0.4e\n', total_avail_world) ;
                        end
                        continue
                    end
                else
                    error('If total_unmet_mgmt_YXv(k+1,m+1,i)==0, you should have skipped!')
                end
                clear out_y1_2deg_thisTot_YX   % Just to make sure you don't try and use it later.
                if do_debug && k==dbk && m==dbm && i==dbCrop
                    fprintf('                total_avail_world =\t\t%0.4e\n', total_avail_world) ;
                end
                
%                 % It's possible that this cell WAS at its max mgmt in
%                 % mid_y1_2deg_mgmt_YXthis (and thus needed to send mgmt to a
%                 % different cell), but it previously donated existing mgmt to a
%                 % different cell, and so now (some of) its own donation demand
%                 % can be considered satisfied. NOTE that this is different from
%                 % how it was done in original (and how it's currently done in
%                 % ringRedist for area): Originally, a ring was started even if
%                 % a cell had enough room to provide for all of its "unmet."
%                 max_mgmtTot_thisCell = max_mgmt(i) * out_y1_2deg_cropArea_YXv(k+1,m+1,i) ;
%                 out_total_thisCell = out_y1_2deg_mgmt_YXv(k+1,m+1,i) * out_y1_2deg_cropArea_YXv(k+1,m+1,i) ;
%                 if total_unmet_mgmt_YXv(k+1,m+1,i)>0 && out_total_thisCell < max_mgmtTot_thisCell
%                     unmetReduction = min(max_mgmtTot_thisCell-out_total_thisCell,total_unmet_mgmt_YXv(k+1,m+1,i)) ;
%                     out_y1_2deg_mgmt_YXv(k+1,m+1,i) = (out_total_thisCell + unmetReduction) / out_y1_2deg_cropArea_YXv(k+1,m+1,i) ;
%                     total_unmet_mgmt_YXv(k+1,m+1,i) = total_unmet_mgmt_YXv(k+1,m+1,i) - unmetReduction ;
% %                     if do_debug && i==dbCrop
% %                         total_avail_world = sum(sum(out_y1_2deg_mgmt_YXv(:,:,i))) ;
% %                         fprintf('            now total_avail_world =\t\t%0.4e\n', total_avail_world) ;
% %                     end
%                 end
                
                total_unmet_thisCell = total_unmet_mgmt_YXv(k+1,m+1,i) ;
                did_ring = false ;
                j = 0 ;   % Note that this differs from how it was done in LUH1 code (pasture started with j from this cell's crop redistribution)
                while abs(total_unmet_thisCell)>1e-8
                    j = j+1 ;
                    if j==1
                        did_ring = true ;
                        out_this_YX = out_y1_2deg_mgmt_YXv(:,:,i) ;
                        out_thisArea_YX = out_y1_2deg_cropArea_YXv(:,:,i) ;
                    elseif j>100
                        error('Possible infinite loop in mgmt ring adjustments (%s).',check_name)
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
                    thisRing = sub2ind(YX_dims,tmp(:,1),tmp(:,2)) ;
                    % If ring gets big enough, you can have duplicates.
                    % Remove them.
                    if j*2+1 > min(YX_dims)
                        thisRing = unique(thisRing) ;
                    end
                    
                    % If no available, continue (i.e., expand ring and
                    % start over).
                    out_thisTot_thisRing = out_this_YX(thisRing) .* out_thisArea_YX(thisRing) ;
                    if total_unmet_thisCell < 0
                        avail_thisRing = out_thisTot_thisRing ;
                        total_avail_ring = sum(avail_thisRing) ;
                    elseif total_unmet_thisCell > 0
                        avail_thisRing = max_mgmtTot_YX(thisRing) - out_thisTot_thisRing ;
                        avail_thisRing(avail_thisRing<0) = 0 ;
                        total_avail_ring = sum(avail_thisRing) ;
                    else
                        error('How was this not skipped?')
                    end
                    % Debugging
                    if do_debug && k==dbk && m==dbm && i==dbCrop
                        if total_unmet_mgmt_YXv(k+1,m+1,i) < 0
                            total_avail_world = sum(sum(out_this_YX .* out_thisArea_YX)) ;
                        else
                            total_avail_world = sum(sum(max_mgmtTot_YX - out_this_YX.*out_thisArea_YX)) ;
                        end
                        fprintf('                                j =\t\t%d\n', j) ;
                        fprintf('                total_avail_world =\t\t%0.4e\n', total_avail_world) ;
                        if do_debug && k==dbk && m==dbm && i==dbCrop
                            fprintf('%s, j = %d, total_unmet_mgmt_YXv(k+1,m+1,i) =\t%0.4e\n',...
                                LPJGcrops{i},j,total_unmet_thisCell) ;
                            fprintf('                total_avail_ring =\t\t%0.4e\n', total_avail_ring) ;
                            if total_avail_ring~=0
                                x=1;
                            end
                        end
                    end
                    if total_avail_ring==0
                        continue
                    end
                    
                    % Try to distribute to / take from ring
                    [out_this_YX, total_unmet_thisCell, ...
                        displaced_this_YX, this_meanDist_YX] = ...
                        PLUMharm_doRings_mgmt(...
                        out_this_YX, total_unmet_thisCell, ...
                        out_thisArea_YX, max_mgmt(i), displaced_this_YX, ...
                        thisCell, thisRing, innerCells, this_meanDist_YX, ...
                        conserv_tol_pct, i, do_debug, dbCrop, thisCell_ofInt) ;
                    
                    % Debugging: Update mean distance
                    if do_debug
                        meanDist_YXv(:,:,i) = this_meanDist_YX ;
                    end
                    
                    % Debugging: Check that management changes are (mostly) conserved
                    if do_debug
                        out_y1_2deg_mgmt_YXv(:,:,i) = out_this_YX ;
                        total_unmet_mgmt_YXv(k+1,m+1,i) = total_unmet_thisCell ;
                        bad = PLUMharm_checkCons_mgmt(...
                            out_y0_mgmt.maps_YXv, out_y0_area_YXv, ...
                            out_y1_2deg_mgmt_YXv, out_y1_2deg_cropArea_YXv, ...
                            in_y0_mgmt.maps_YXv, in_y0_area_YXv, ...
                            in_y1_mgmt.maps_YXv, in_y1_area_YXv, ...
                            total_unmet_mgmt_YXv, LPJGcrops, conserv_tol_pct, notEnough, ...
                            check_name, false) ;
                        if bad==1
                            keyboard
                        end
                    end

                end % while loop
                if do_debug
                    displaced_mgmt_YXv(:,:,i) = displaced_this_YX ;
                    agri_rings_YXv(:,:,i) = this_rings_YX ;
                end
                
                % Record this cell as done
                is_done_YXv(k+1,m+1,i) = true ;
                if did_ring
                    out_y1_2deg_mgmt_YXv(:,:,i) = out_this_YX ;
                    total_unmet_mgmt_YXv(k+1,m+1,i) = total_unmet_thisCell ;
                end
                
                % If doing this, need to add the following arguments to
                % function call:
% % %                             out_y0_mgmt, out_y0_area_YXv, ...
% % %                             in_y0_mgmt, in_y0_area_YXv, ...
% % %                             in_y1_mgmt, in_y1_area_YXv, ...
% % %                             conserv_tol_pct, check_name
                if do_debug
                    bad = PLUMharm_checkCons_mgmt(...
                        out_y0_mgmt.maps_YXv, out_y0_area_YXv, ...
                        out_y1_2deg_mgmt_YXv, out_y1_2deg_cropArea_YXv, ...
                        in_y0_mgmt.maps_YXv, in_y0_area_YXv, ...
                        in_y1_mgmt.maps_YXv, in_y1_area_YXv, ...
                        total_unmet_mgmt_YXv, LPJGcrops, conserv_tol_pct, notEnough, ...
                        check_name, do_warn) ;
                    if bad==1
                        keyboard
                    end
                end
                
                
            end % for loop
            
        end % for m
    end % for k
    
    % This shouldn't ever happen now.
    if isequal(is_done_YXv_begin,is_done_YXv)
        error('Not enough mgmt in the world!')
    end
    
    if any(~is_done_YXv(:))
        for c = 1:Ncrops
            if any_skipped_v(c)
                skipped_cells = find(~is_done_YXv(:,:,c)) ;
                thisSkipped = skipped_cells(1) ;
                [thisSkipped_I, thisSkipped_J] = ind2sub([ny nx],thisSkipped) ;
                fprintf('      %s %s: Skipped cell %d (%d,%d) and %d others.\n', ...
                    LPJGcrops{c}, check_name, thisSkipped, thisSkipped_I, thisSkipped_J, length(skipped_cells)-1) ;
            end
        end
        disp('      Trying again.')
        pause(0.1)
    end
    
    
    
end % while

if do_debug
%     keyboard
end

% figure('Color','w','Position',figurePos) ;
% h = cell(Ncrops,1) ;
% for v = 1:Ncrops
%     h{v} = subplot_tight(3,3,v,0.025) ;
%     tmp_YX = agri_rings_YXv(:,:,v) ;
%     tmp_YX = single(tmp_YX) ;
%     tmp_YX(landArea_2deg_YX==0) = NaN ;
%     pcolor(tmp_YX); shading flat; axis equal tight off
%     caxis([0 max(agri_rings_YXv(:))]) ;
%     colorbar ;
%     title(LPJGcrops{v})
%     set(gca,'FontSize',14) ;
% end
% figure('Color','w','Position',figurePos) ;
% h = cell(Ncrops,1) ;
% for v = 1:Ncrops
%     h{v} = subplot_tight(3,3,v,0.025) ;
%     tmp_YX = displaced_mgmt_YXv(:,:,v) ;
%     tmp_YX(landArea_2deg_YX==0) = NaN ;
%     pcolor(tmp_YX); shading flat; axis equal tight off
%     colorbar ;
%     title(LPJGcrops{v})
%     set(gca,'FontSize',14) ;
% end
% keyboard


end
