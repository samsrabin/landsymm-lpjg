function [out_y1_agri_YXv, already_done] = ...
    loop_thru_agri(already_done, k1, debug_ijk, ...
    update_avail_land, proper_zero_denoms, conserv_tol_area, ...
    i, j, iy, ix, theseCells, out_y0_vegd_theseCells, agri_d, ...
    out_y0_agri_YXv, out_y0_2deg_agri_YXv, out_y1_agri_YXv)

Nagri = size(out_y0_agri_YXv,3) ;
if k1 > Nagri
    error('k1 > Nagri')
end

for k = k1:Nagri
    
    if debug_ijk(1)==i && debug_ijk(2)==j && debug_ijk(3)==k
%         [out_y0_vegd_theseCells out_y0_vegd_YX(theseCells)]
        keyboard
    end
    
    if update_avail_land
        now_agri_YX = sum(out_y1_agri_YXv(iy,ix,:),3) ;
        avail_land = out_y0_vegd_theseCells - now_agri_YX(:) ;
    end
    
    % Skip if done previously
    if already_done(k)
        out_y1_this_theseCells = out_y1_agri_YXv(iy,ix,k) ;
        out_y1_this_theseCells = out_y1_this_theseCells(:) ;
%         if debug_ijk(1)==i && debug_ijk(2)==j
%             keyboard
%         end
    else
        % Get delta for this land use. If 0, skip. Otherwise, set
        % up for rest of logic.
        this_d = agri_d(k) ;
        if this_d==0
            out_y1_agri_YXv(iy,ix,k) = out_y0_agri_YXv(iy,ix,k) ;
            already_done(k) = true ;
            if debug_ijk(1)==i && debug_ijk(2)==j
                keyboard
            end
            continue
        end
        out_y0_this_YX = out_y0_agri_YXv(iy,ix,k) ;
        out_y0_2deg_this = out_y0_2deg_agri_YXv(i,j,k) ;
        out_y1_this_theseCells = [] ;
        
        % Update this LU
        if this_d < 0
            if proper_zero_denoms
                this_p = this_d / out_y0_2deg_this ;
                this_p(out_y0_2deg_this==0) = 0 ;
            else
                this_p = this_d / (out_y0_2deg_this + 1e-12) ;
            end
            tmp = out_y0_this_YX * (1 + this_p) ;
            out_y1_this_theseCells = reshape(tmp,size(out_y1_agri_YXv(iy,ix,k))) ;
%             if debug_ijk(1)==i && debug_ijk(2)==j
%                 keyboard
%             end
        elseif this_d > 0
            if sum(avail_land) < this_d ...
                    ...&& k < Nagri && sum(agri_d((k+1):end)) < 0
                    && k < Nagri
                
%                 if sum(agri_d((k+1):end)) >= 0
%                     error('How is sum(agri_d((k+1):end)) >= 0 ??')
%                 end

                if debug_ijk(1)==i && debug_ijk(2)==j
                    disp('GOING D E  E   P    E     R')
                    keyboard
                end
                
                
                % Go ahead and do next land uses until there is enough
                % available land. RECURSION
                [out_y1_agri_YXv, already_done] = ...
                    loop_thru_agri(already_done, k+1, debug_ijk, ...
                    update_avail_land, proper_zero_denoms, conserv_tol_area, ...
                    i, j, iy, ix, theseCells, out_y0_vegd_theseCells, agri_d, ...
                    out_y0_agri_YXv, out_y0_2deg_agri_YXv, out_y1_agri_YXv) ;
                
                
                
                % How much NATURAL land is available now we have
                % sufficently decreased enough of the next LUs?
                now_agri_YX_this = sum(out_y1_agri_YXv,3) ;
                avail_land_this = out_y0_vegd_theseCells - now_agri_YX_this(theseCells) ;
                % Distribute thisLU demand to half-degree gridcells
                % based on how much NATURAL land each has
                if proper_zero_denoms
                    if sum(avail_land_this)==0
                        error('How is there no available land here??')
                    end
                    tmp = out_y0_this_YX(:) + avail_land_this / sum(avail_land_this) * this_d ;
                else
                    tmp = out_y0_this_YX(:) + avail_land_this / (sum(avail_land_this) + 1e-12) * this_d ;
                end
                out_y1_this_theseCells = reshape(tmp,size(out_y1_agri_YXv(iy,ix,k))) ;
            else % thisLU increases and (nextLU increases OR available land is enough to satisfy thisLU demand)
                % Distribute thisLU demand to half-degree gridcells based
                % on how much NATURAL land each has
                if debug_ijk(1)==i && debug_ijk(2)==j
                    keyboard
                end
                if proper_zero_denoms
                    if sum(avail_land)==0
                        error('How is there no available land here??')
                    end
                    tmp = out_y0_this_YX(:) + avail_land / sum(avail_land) * this_d ;
                else
                    tmp = out_y0_this_YX(:) + avail_land / (sum(avail_land) + 1e-12) * this_d ;
                end
                out_y1_this_theseCells = reshape(tmp,size(out_y1_agri_YXv(iy,ix,k))) ;
            end
        else
            error('How is this_d not < or > 0??')
        end
        
        % Sanity check
        if isempty(out_y1_this_theseCells)
            error('How did out_y1_this_theseCells not get set??')
        end
    end
    
    % Reshape from matrix to vector
    tmp = out_y1_this_theseCells(:) ;
    
    % Check for invalid cell areas
    if any(tmp < -conserv_tol_area)
        error('Negative members of half-deg out_y1_this_YXv!')
    elseif any(tmp > conserv_tol_area+out_y0_vegd_theseCells)
        error('Members >vegd_area in half-deg out_y1_this_YX!')
    end
    
    % NOTE: This bit of code is where I should ensure that the total
    % agricultural area does not exceed the vegetated area. I think it
    % needs to happen AFTER loop_thru_agri in top function. At the moment,
    % though, it doesn't get called, so I'm leaving it where it is.
    this_minus_vegd = tmp - out_y0_vegd_theseCells ;
    while any(this_minus_vegd > conserv_tol_area)
        error('You have not generalized this bit of code yet.')
        % Indices of where, in theseCells, there is too much cropland
        is_too_much_crop = find(this_minus_vegd > conserv_tol_area) ;
        % What is the total amount of excess cropland?
        sum_too_much_crop = sum(this_minus_vegd(is_too_much_crop)) ;
        % Cap cropland at vegetated area
        tmp( is_too_much_crop) = out_y0_vegd_theseCells(is_too_much_crop) ;
        % Recalculate available land for all theseCells
        avail_land = out_y0_vegd_theseCells - out_y1_crop_YX(theseCells) - out_y0_past_YX(theseCells) ;
        % Distribute excess cropland to cells that didn't have excess cropland, 
        % proportionally according to how much available land there is in each cell.
        tmp(~is_too_much_crop) = out_y0_crop_YX(~is_too_much_crop) + avail_land(~is_too_much_crop) / (sum(sum(avail_land)) + 1e-12) * sum_too_much_crop ;
        % Update this_minus_vegd to see if WHILE condition is still satisfied
        this_minus_vegd = tmp - out_y0_vegd_theseCells ;
    end
   
    % Save this (if previously done: redundant but safe)
    out_y1_this_theseCells = reshape(tmp,size(out_y1_agri_YXv(iy,ix,k))) ;
    out_y1_agri_YXv(iy,ix,k) = out_y1_this_theseCells ;
    already_done(k) = true ;
    
    % Another check
    if k==Nagri
        agri_YX = sum(out_y1_agri_YXv(iy,ix,:),3) ;
        if debug_ijk(1)==i && debug_ijk(2)==j
            disp(['k = ' num2str(k)])
            [agri_YX(:) out_y0_vegd_theseCells agri_YX(:)-out_y0_vegd_theseCells]
        end
        if any(agri_YX(:) > conserv_tol_area+out_y0_vegd_theseCells)
            error('Members >vegd_area in half-deg out_y1_(crop+past)_YX!')
        end
    end
    
    
end % for k = 1:Nagri


end % loop_thru_agri()