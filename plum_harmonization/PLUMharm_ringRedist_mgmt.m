function [out_y1_2deg_mgmt_YXv, total_unmet_mgmt_YXv] = ...
    PLUMharm_ringRedist_mgmt(...
    mid_y1_2deg_mgmt_YXv, out_y1_2deg_cropArea_YXv, ...
    total_unmet_mgmt_YXv, max_mgmt, ...
    LPJGcrops, debugIJ, ...
    out_y0_mgmt, out_y0_area_YXv, ...
    in_y0_mgmt, in_y0_area_YXv, ...
    in_y1_mgmt, in_y1_area_YXv, ...
    conserv_tol_pct, check_name)
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
    % Keep track of what cells have been processed so far
    thisCell_list = [] ;
    % Get k and m of interest
    dbk = debugIJ(1) - 1 ;
    dbm = debugIJ(2) - 1 ;
end
i_ofInterest = 1 ;

out_y1_2deg_mgmt_YXv = mid_y1_2deg_mgmt_YXv ;
out_y1_2deg_mgmtTot_YXv = mid_y1_2deg_mgmtTot_YXv ;

ny = YX_dims(1) ;
nx = YX_dims(2) ;
ks = 0:(ny-1) ;
ms = 0:(nx-1) ;

is_done_YXv = false(YXv_dims) ;

while any(~is_done_YXv)
    is_done_YXv_begin = is_done_YXv ;
    for k = ks
        for m = ms
            
            thisCell = sub2ind(YX_dims,k+1,m+1) ;
            if do_debug
                %             thisCell_list = [thisCell_list thisCell] ;
                %             if k==dbk && m==dbm
                %                 keyboard
                %             end
            end
            
            % Do it for each
            for i = 1:Ncrops
                
%                 % Skip if this crop is done worldwide
%                 if ~any(any(~is_done_YXv(:,:,i)))
%                     continue
%                 end
                
                % Skip if this crop for this cell is already done
                if is_done_YXv(k+1,m+1,i)
                    continue
                end
                
                % If there's too much desired mgmt loss for the entire
                % world to handle, set all mgmt to zero, and skip from now
                % on.
                if k==ks(1) && m==ms(1)
                    testSum = sum(sum(total_unmet_mgmt_YXv(:,:,i) + out_y1_2deg_mgmt_YXv(:,:,i).*out_y1_2deg_cropArea_YXv(:,:,i))) ;
                    if testSum < 0
                        warning(['There is not enough mgmt applied to %s in the world to absorb the negative unmet.\n' ...
                            'Difference = %0.4e.\n' ...
                            'Setting all %s mgmt to zero.'], ...
                            LPJGcrops{i}, -testSum, LPJGcrops{i})
                        out_y1_2deg_mgmt_YXv(:,:,i) = 0 ;
                        total_unmet_mgmt_YXv(:,:,i) = 0 ; % There may be a nicer way to handle this.
                        is_done_YXv(:,:,i) = true ;
                        continue
                    end
                end
                
                if do_debug
                    displaced_this_YX = displaced_mgmt_YXv(:,:,i) ;
                    this_rings_YX = agri_rings_YXv(:,:,i) ;
                    this_meanDist_YX = meanDist_YXv(:,:,i) ;
                else
                    displaced_this_YX = [] ;
                    this_meanDist_YX = [] ;
                end
                
                
                if do_debug && k==dbk && m==dbm && i==i_ofInterest
                    %                 keyboard
                end
                
                % If this cell has no unmet, skip.
                if total_unmet_mgmt_YXv(k+1,m+1,i)==0
                    is_done_YXv(k+1,m+1,i) = true ;
                    continue
                end
                
                % If this cell needs to take mgmt from elsewhere (i.e., its
                % unmet is negative), make sure there is enough mgmt in the
                % world for that to work. If not enough, skip and come back
                % later. If this cell is instead donating mgmt, calculate
                % total_avail_world using appropriate method.
                if total_unmet_mgmt_YXv(k+1,m+1,i) < 0
                    max_mgmtTot_YX = max_mgmt(i) * out_y1_2deg_cropArea_YXv(:,:,i) ;
                    total_avail_world = sum(sum(out_y1_2deg_mgmt_YXv(:,:,i) .* out_y1_2deg_cropArea_YXv(:,:,i))) ;
                    if -total_unmet_mgmt_YXv(k+1,m+1,i) > total_avail_world
                        warning('Skipping %d, to return later.', thisCell) ;
                        continue
                    end
                elseif do_debug && total_unmet_mgmt_YXv(k+1,m+1,i) > 0
                    max_mgmtTot_YX = max_mgmt(i) * out_y1_2deg_cropArea_YXv(:,:,i) ;
                    total_avail_world = sum(sum(max_mgmtTot_YX - out_y1_2deg_mgmt_YXv(:,:,i))) ;
                elseif do_debug
                    error('If total_unmet_mgmt_YXv(k+1,m+1,i)==0, you should have skipped!')
                end
                if do_debug && i==i_ofInterest
                    fprintf('                total_avail_world =\t\t%0.8f\n', total_avail_world) ;
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
                out_total_thisCell = out_y1_2deg_mgmt_YXv(k+1,m+1,i) ;
                if total_unmet_mgmt_YXv(k+1,m+1,i)>0 && out_total_thisCell < max_mgmtTot_thisCell
                    unmetReduction = min(max_mgmtTot_thisCell-out_total_thisCell,total_unmet_mgmt_YXv(k+1,m+1,i)) ;
                    out_y1_2deg_mgmt_YXv(k+1,m+1,i) = (out_total_thisCell + unmetReduction) / out_y1_2deg_cropArea_YXv(k+1,m+1,i) ;
                    total_unmet_mgmt_YXv(k+1,m+1,i) = total_unmet_mgmt_YXv(k+1,m+1,i) - unmetReduction ;
                    if do_debug && i==i_ofInterest
                        total_avail_world = sum(sum(out_y1_2deg_mgmt_YXv(:,:,i))) ;
                        fprintf('            now total_avail_world =\t\t%0.8f\n', total_avail_world) ;
                    end
                end
                
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
                        error('Possible infinite loop in mgmt ring adjustments.')
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
                    
                    % Debugging
                    if do_debug && i==i_ofInterest
                        fprintf('                total_avail_world =\t\t%0.8f\n', total_avail_world) ;
                        if do_debug && k==dbk && m==dbm && i==i_ofInterest
                            out_y1_2deg_thisTot_YX = out_this_YX .* out_thisArea_YX ;
                            fprintf('%s, j = %d, total_unmet_mgmt_YXv(k+1,m+1,i) =\t%0.8f\n',...
                                LPJGcrops{i},j,total_unmet_thisCell) ;
                            if total_unmet_thisCell < 0
                                total_avail_ring = sum(out_y1_2deg_thisTot_YX(thisRing)) ;
                            else
                                total_avail_ring = sum(max_mgmtTot_YX(thisRing) - out_y1_2deg_thisTot_YX(thisRing)) ;
                            end
                            fprintf('                total_avail_ring =\t\t%0.8f\n', total_avail_ring) ;
                            if total_avail_ring~=0
                                keyboard
                            end
                        end
                    end
                    
                    % Try to distribute to / take from ring
                    [out_this_YX, total_unmet_thisCell, ...
                        displaced_this_YX, this_meanDist_YX] = ...
                        PLUMharm_doRings_mgmt(...
                        out_this_YX, total_unmet_thisCell, ...
                        out_thisArea_YX, max_mgmt(i), displaced_this_YX, ...
                        thisCell, thisRing, innerCells, this_meanDist_YX, i, do_debug, i_ofInterest) ;
                    
                    if do_debug
                        meanDist_YXv(:,:,i) = this_meanDist_YX ;
                    end
                    
                    % Debugging: Check that management changes are (mostly) conserved
                    if do_debug
                        out_y1_2deg_mgmt_YXv(:,:,i) = out_this_YX ;
                        total_unmet_mgmt_YXv(k+1,m+1,i) = total_unmet_thisCell ;
                        bad = PLUMharm_check_conservation(...
                            out_y0_mgmt.maps_YXv, out_y0_area_YXv, ...
                            out_y1_2deg_mgmt_YXv, out_y1_2deg_cropArea_YXv, ...
                            in_y0_mgmt.maps_YXv, in_y0_area_YXv, ...
                            in_y1_mgmt.maps_YXv, in_y1_area_YXv, ...
                            total_unmet_mgmt_YXv, LPJGcrops, conserv_tol_pct, ...
                            check_name) ;
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
                    bad = PLUMharm_check_conservation(...
                        out_y0_mgmt.maps_YXv, out_y0_area_YXv, ...
                        out_y1_2deg_mgmt_YXv, out_y1_2deg_cropArea_YXv, ...
                        in_y0_mgmt.maps_YXv, in_y0_area_YXv, ...
                        in_y1_mgmt.maps_YXv, in_y1_area_YXv, ...
                        total_unmet_mgmt_YXv, LPJGcrops, conserv_tol_pct, ...
                        check_name) ;
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
    
end % while

if do_debug
    keyboard
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


function [out_YX, total_unmet_thisCell, ...
    displaced_YX, this_meanDist_YX] = ...
    PLUMharm_doRings_mgmt( ...
    in_YX, total_unmet_thisCell, ...
    thisArea_YX, max_mgmt, ...
    displaced_YX, thisCell, thisRing, innerCells, ...
    this_meanDist_YX, i, do_debug, i_ofInterest)

%%% Internals
% avail_mgmt = "other" area in thisRing
% total_avail_mgmt = avail_mgmt summed over all gridcells in ring
%%% Outputs
% displaced_YX = net area of this land type moved from or to this gridcell.
%                Negative value indicates ???

thisArea_thisRing = thisArea_YX(thisRing) ;
in_total_thisRing = in_YX(thisRing) .* thisArea_thisRing ;

% Set up outs
out_YX = in_YX ;

if do_debug
    [~,IA,~] = intersect(thisRing,innerCells,'stable') ;
    thisRing_isInnerCell = false(size(thisRing)) ;
    thisRing_isInnerCell(IA) = true ;
end
j = (sqrt(length(thisRing))-1)/2 ;   % What ring are we on?

% Debugging
% thisCell_ofInt = 13205 ;

% if thisCell needs to give MANAGEMENT to the rest of the ring
% (i.e., thisCell had more MANAGEMENT than allowed)
if total_unmet_thisCell>0
    
    % Get available management
    max_mgmt_thisRing = max_mgmt * thisArea_thisRing ;
    avail_mgmt = max_mgmt_thisRing - in_total_thisRing ;
    avail_mgmt(avail_mgmt<0) = 0 ;
    total_avail_mgmt = sum(avail_mgmt) ;
    
    % Sanity checks
    if max(in_total_thisRing-max_mgmt_thisRing) > 1e-3
        error('How do you have out_total > max_mgmt_thisRing (%0.4f) in this ring? (thisCell = %d, j = %d)', ...
                    max(in_total_thisRing-max_mgmt_thisRing), thisCell, j)
    end
    if do_debug && max(avail_mgmt(thisRing_isInnerCell)) > 1e-3
        error('How do you have available mgmt (%0.4f) in cell(s) not on ring perimeter? (thisCell = %d, j = %d)', ...
                    max(avail_mgmt(thisRing_isInnerCell)), thisCell, j)
    end
%     if i==1 && any(thisRing==thisCell_ofInt) && thisCell~=thisCell_ofInt% && avail_mgmt(thisRing==thisCell_ofInt)>0
%         keyboard
%     end
    if any(isnan(avail_mgmt))
        error('NaN in avail_mgmt!')
    end
    
    % If there's some available management, distribute thisCell's excess
    % out to thisRing.
    if total_avail_mgmt>0
        
        % EITHER there is not enough headroom in thisRing to absorb the excess mgmt in thisCell
        if total_unmet_thisCell >= total_avail_mgmt
            to_ring = avail_mgmt ;
            total_unmet_thisCell = total_unmet_thisCell - total_avail_mgmt ;
        
        % OR ELSE the headroom in thisRing is sufficient to absorb the excess mgmt in thisCell
        else
            to_ring = total_unmet_thisCell * avail_mgmt / total_avail_mgmt ;
            total_unmet_thisCell = 0;
        end
        
        % Save changes
        out_mgmtTotal_thisRing_new = in_total_thisRing + to_ring ;
        out_thisRing = out_mgmtTotal_thisRing_new ./ thisArea_thisRing ;
        out_thisRing(thisArea_thisRing==0) = 0 ;
        out_YX(thisRing) = out_thisRing ;

        % Update debugging info
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) + to_ring;
            % This cell is giving area to other cells in ring. So here, we
            % update this_meanDist_YX, which keeps track of the mean
            % distance (in # rings) that area in each cell came from.
            now_weighted = this_meanDist_YX(thisRing) .* in_total_thisRing./out_mgmtTotal_thisRing_new ;
            new_weighted =                          j  * to_ring./out_mgmtTotal_thisRing_new ;
            now_weighted(out_mgmtTotal_thisRing_new==0) = 0 ;
            new_weighted(out_mgmtTotal_thisRing_new==0) = 0 ;
            this_meanDist_YX(thisRing) = now_weighted + new_weighted ;
        end
%         if i==6 && any(thisRing==thisCell_ofInt) && out_mgmtTotal_thisRing_new(thisRing==thisCell_ofInt) > out_total_thisRing(thisRing==thisCell_ofInt)
%             warning('Cell %d is receiving from %d.', thisCell_ofInt, thisCell) ;
%             fprintf('out_total before: %0.4e\n', out_total_thisRing(thisRing==thisCell_ofInt)) ;
%             fprintf('out_total after: %0.4e\n', out_mgmtTotal_thisRing_new(thisRing==thisCell_ofInt)) ;
%         end
    end
    
    
% elseif thisCell needs MANAGEMENT donated from rest of thisRing
% (i.e., thisCell didn't have enough MANAGEMENT to satisfy PLUM-specified loss)
elseif total_unmet_thisCell<0
    avail_mgmt = in_total_thisRing ;
%     if i==6 && any(thisRing==thisCell_ofInt) && thisCell~=thisCell_ofInt && avail_mgmt(thisRing==thisCell_ofInt)>0
%         keyboard
%         warning('Cell %d is donating to %d.', thisCell_ofInt, thisCell) ;
%         fprintf('out_total before: %0.4e\n', out_total_thisRing(thisRing==thisCell_ofInt)) ;
%     end
    if any(isnan(avail_mgmt))
        error('NaN in avail_mgmt!')
    end
    total_avail_mgmt = sum(sum(avail_mgmt .* (avail_mgmt>0))) ;
  % if the THISLU available in thisRing is not sufficient to satisfy the demand of thisCell
    if total_unmet_thisCell <= -total_avail_mgmt
        out_total_thisRing = in_total_thisRing - avail_mgmt .* (avail_mgmt>0) ;
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) - avail_mgmt .* (avail_mgmt>0);
        end
        total_unmet_thisCell = total_unmet_thisCell + total_avail_mgmt ;
  % elseif the THISLU available in the rest of thisRing is sufficient to satisfy the demand of thisCell
    elseif total_unmet_thisCell > -total_avail_mgmt
        out_total_thisRing = in_total_thisRing + total_unmet_thisCell * avail_mgmt ./ total_avail_mgmt .* (avail_mgmt>0);
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) + total_unmet_thisCell * avail_mgmt ./ total_avail_mgmt .* (avail_mgmt>0);
        end
        total_unmet_thisCell = 0;
    end
    out_thisRing = out_total_thisRing ./ thisArea_thisRing ;
    out_thisRing(thisArea_thisRing==0) = 0 ;
    out_YX(thisRing) = out_thisRing ;
%     if i==6 && any(thisRing==thisCell_ofInt) && thisCell~=thisCell_ofInt && avail_mgmt(thisRing==thisCell_ofInt)>0
%         out_total_thisRing = out_total_YX(thisRing) ;
%         fprintf('out_total after: %0.4e\n', out_total_thisRing(thisRing==thisCell_ofInt)) ;
%     end
end



end