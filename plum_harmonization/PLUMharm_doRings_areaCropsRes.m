function [out_YX, out_ntrl_YX, total_unmet_thisCell_out, ...
    displaced_YX, this_meanDist_YX] = ...
    PLUMharm_doRings_areaCropsRes( ...
    in_YX, out_agri_YX, total_unmet_thisCell_in, displaced_YX, ...
    thisRing, in_ntrl_YX, resArea_YX, ...
    this_meanDist_YX, thisCell_ofInt)

%%% Internals
% avail_space = "other" area in thisRing
% total_avail_space = avail_space summed over all gridcells in ring
%%% Outputs
% displaced_YX = net area of this land type moved from or to this gridcell.
%                Negative value indicates ???

% Get land available for conversion to agriculture,
% as max(NATURAL-RESERVED,0)
nonResNtrl_YX = max(in_ntrl_YX - resArea_YX, 0) ;

% Set this up
out_YX = in_YX ;

% % Debugging
% if any(thisRing==thisCell_ofInt)
%     keyboard
% end

% if thisCell needs to give THISLU to the rest of the ring
% (i.e., thisCell had more vegetated area than allowed, or more non-NATURAL
% area than allowed)
% WARNING: This makes it so that order of land uses matters. Is there a way
% to avoid this?
if total_unmet_thisCell_in>0
    avail_space = nonResNtrl_YX(thisRing) - out_agri_YX(thisRing) ;
    avail_space(avail_space<0) = 0 ;
    total_avail_space = sum(sum(avail_space)) ;
    if total_avail_space>0
      % if there is not enough NATURAL_UNRESERVED in thisRing to absorb the excess THISLU in thisCell
        if total_unmet_thisCell_in >= total_avail_space
%             % Debugging
%             if any(thisRing==thisCell_ofInt)
%                 keyboard
%             end
            to_ring = avail_space ;
            total_unmet_thisCell_out = total_unmet_thisCell_in - total_avail_space ;
      % else the NATURAL_UNRESERVED in thisRing is sufficient to absorb the excess THISLU in thisCell
        else
%             % Debugging
%             if any(thisRing==thisCell_ofInt)
%                 keyboard
%             end
            to_ring = total_unmet_thisCell_in * avail_space ./ total_avail_space ;
            total_unmet_thisCell_out = 0;
        end
%         % Debugging
%         if any(thisRing==thisCell_ofInt)
%             keyboard
%         end
        out_thisRing_new = out_YX(thisRing) + to_ring ;
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) + to_ring;
            % This cell is giving area to other cells in ring. So here, we
            % update this_meanDist_YX, which keeps track of the mean
            % distance (in # rings) that area in each cell came from.
            j = (sqrt(length(thisRing))-1)/2 ;   % What ring are we on?
            now_weighted = this_meanDist_YX(thisRing) .* out_YX(thisRing)./out_thisRing_new ;
            new_weighted =                          j  * to_ring./out_thisRing_new ;
            now_weighted(out_thisRing_new==0) = 0 ;
            new_weighted(out_thisRing_new==0) = 0 ;
            this_meanDist_YX(thisRing) = now_weighted + new_weighted ;
        end
        out_YX(thisRing) = out_thisRing_new ;
    else
        keyboard
    end
% elseif thisCell needs THISLU donated from rest of thisRing
% (i.e., thisCell didn't have enough thisLU to satisfy PLUM-specified loss)
elseif total_unmet_thisCell_in<0
    avail_space = out_YX(thisRing);
    total_avail_space = sum(sum(avail_space .* (avail_space>0))) ;
  % if the THISLU available in thisRing is not sufficient to satisfy the demand of thisCell
    if total_unmet_thisCell_in <= -total_avail_space
%         % Debugging
%         if any(thisRing==thisCell_ofInt)
%             keyboard
%         end
        out_YX(thisRing) = out_YX(thisRing) - avail_space .* (avail_space>0);
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) - avail_space .* (avail_space>0);
        end
        total_unmet_thisCell_out = total_unmet_thisCell_in + total_avail_space ;
  % elseif the THISLU available in the rest of thisRing is sufficient to satisfy the demand of thisCell
    elseif total_unmet_thisCell_in > -total_avail_space
%         % Debugging
%         if any(thisRing==thisCell_ofInt)
%             keyboard
%         end
        out_YX(thisRing) = out_YX(thisRing) + total_unmet_thisCell_in * avail_space ./ total_avail_space .* (avail_space>0);
        if ~isempty(displaced_YX)
            displaced_YX(thisRing) = displaced_YX(thisRing) + total_unmet_thisCell_in * avail_space ./ total_avail_space .* (avail_space>0);
        end
        total_unmet_thisCell_out = 0;
    end
end

out_ntrl_YX = in_ntrl_YX - (out_YX - in_YX) ;

% Debugging
if any(thisRing==thisCell_ofInt)
    fprintf('resArea=%0.4e, avail=%0.4e, in_ntrl=%0.4e, out_ntrl=%0.4e, in_this=%0.4e, out_this=%0.4e\n',...
    resArea_YX(thisCell_ofInt), avail_space(thisRing==thisCell_ofInt), ...
    in_ntrl_YX(thisCell_ofInt), out_ntrl_YX(thisCell_ofInt), ...
    in_YX(thisCell_ofInt), out_YX(thisCell_ofInt)) ;
%     keyboard
end

if ~exist('total_unmet_thisCell_out','var')
    total_unmet_thisCell_out = total_unmet_thisCell_in ;
%     keyboard
end

end