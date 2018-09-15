function [out_YX, total_unmet_thisCell_out, ...
    displaced_YX, this_meanDist_YX] = ...
    PLUMharm_doRings_mgmt( ...
    in_YX, total_unmet_thisCell_in, ...
    thisArea_YX, max_mgmt, ...
    displaced_YX, thisCell, thisRing, innerCells, ...
    this_meanDist_YX, conserv_tol_pct, ...
    i, do_debug, i_ofInterest, thisCell_ofInt)

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

% if thisCell needs to give MANAGEMENT to the rest of the ring
% (i.e., thisCell had more MANAGEMENT than allowed)
if total_unmet_thisCell_in>0
    
    % Get available management
    max_mgmt_thisRing = max_mgmt * thisArea_thisRing ;
    avail_mgmt = max_mgmt_thisRing - in_total_thisRing ;
    avail_mgmt(avail_mgmt<0) = 0 ;
    total_avail_mgmt = sum(avail_mgmt) ;
    
    % Sanity checks
    if max(100*(in_total_thisRing-max_mgmt_thisRing)./max_mgmt_thisRing) > conserv_tol_pct
        error('How do you have out_total > max_mgmt_thisRing (%0.4f) in this ring? (thisCell = %d, j = %d)', ...
                    max(in_total_thisRing-max_mgmt_thisRing), thisCell, j)
    end
    if do_debug && max(avail_mgmt(thisRing_isInnerCell)) > 1e3
        error('How do you have available mgmt (%0.4f) in cell(s) not on ring perimeter? (thisCell = %d, j = %d)', ...
                    max(avail_mgmt(thisRing_isInnerCell)), thisCell, j)
    end
    if any(thisRing==thisCell_ofInt) && total_avail_mgmt>0
%         keyboard
    end
    if any(isnan(avail_mgmt))
        error('NaN in avail_mgmt!')
    end
    
    % If there's some available management, distribute thisCell's excess
    % out to thisRing.
    if total_avail_mgmt>0
        
        % EITHER there is not enough headroom in thisRing to absorb the excess mgmt in thisCell
        if total_unmet_thisCell_in >= total_avail_mgmt
            to_ring = avail_mgmt ;
            total_unmet_thisCell_out = total_unmet_thisCell_in - total_avail_mgmt ;
        
        % OR ELSE the headroom in thisRing is sufficient to absorb the excess mgmt in thisCell
        else
            to_ring = total_unmet_thisCell_in * avail_mgmt / total_avail_mgmt ;
            total_unmet_thisCell_out = 0;
        end
        
        % Debugging
        if any(thisRing==thisCell_ofInt)
%             keyboard
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
    end
    
    
% elseif thisCell needs MANAGEMENT donated from rest of thisRing
% (i.e., thisCell didn't have enough MANAGEMENT to satisfy PLUM-specified loss)
elseif total_unmet_thisCell_in<0
    avail_mgmt = in_total_thisRing ;
    if any(thisRing==thisCell_ofInt)
%         keyboard
    end
    if any(isnan(avail_mgmt))
        error('NaN in avail_mgmt!')
    end
    total_avail_mgmt = sum(sum(avail_mgmt .* (avail_mgmt>0))) ;
  % if the THISLU available in thisRing is not sufficient to satisfy the demand of thisCell
    if total_unmet_thisCell_in <= -total_avail_mgmt
        out_total_thisRing = in_total_thisRing - avail_mgmt .* (avail_mgmt>0) ;
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) - avail_mgmt .* (avail_mgmt>0);
        end
        total_unmet_thisCell_out = total_unmet_thisCell_in + total_avail_mgmt ;
  % elseif the THISLU available in the rest of thisRing is sufficient to satisfy the demand of thisCell
    elseif total_unmet_thisCell_in > -total_avail_mgmt
        out_total_thisRing = in_total_thisRing + total_unmet_thisCell_in * avail_mgmt ./ total_avail_mgmt .* (avail_mgmt>0);
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) + total_unmet_thisCell_in * avail_mgmt ./ total_avail_mgmt .* (avail_mgmt>0);
        end
        total_unmet_thisCell_out = 0;
    end
    if any(thisRing==thisCell_ofInt)
%         keyboard
    end
    out_thisRing = out_total_thisRing ./ thisArea_thisRing ;
    out_thisRing(thisArea_thisRing==0) = 0 ;
    out_YX(thisRing) = out_thisRing ;
%     if i==6 && any(thisRing==thisCell_ofInt) && thisCell~=thisCell_ofInt && avail_mgmt(thisRing==thisCell_ofInt)>0
%         out_total_thisRing = out_total_YX(thisRing) ;
%         fprintf('out_total after: %0.4e\n', out_total_thisRing(thisRing==thisCell_ofInt)) ;
%     end
end

if any(thisRing==thisCell_ofInt)
%     keyboard
end



end