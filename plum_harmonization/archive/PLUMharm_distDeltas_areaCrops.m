function out_y1_agri_YXv = ...
    PLUMharm_distDeltas_areaCrops( ...
    landArea_YX, landArea_2deg_YX, ...
    out_y0_2deg_agri_YXv, out_y1_2deg_agri_YXv, out_y0_agri_YXv, ...
    in_y1_vegd_YX, luh2_vegd_YX, conserv_tol_pct)
% Loop through all 2-degree gridcells and distribute the new LU changes to
% the half-degree gridcells within.
%
% If thisLU change is a decrease, apply the same
% PERCENTAGE decrease to all half-degree thisLU gridcells
% within the 2-degree cell.
%
% If the thisLU change is an increase, apply this increase to
% all half-degree gridcells within the 2-degree cell proportionally
% according to available land area. Make sure that total area within a
% half-degree gridcell does not exceed 1 or 0.


% Update avail_land after every land use? FALSE is default LUH1 behavior,
% where avail_land is only calculated at the beginning.
update_avail_land = true ;

debug_ijk = [Inf Inf Inf] ;
% debug_ijk = [37   154     2] ;

Nagri = size(out_y0_agri_YXv,3) ;

out_y1_agri_YXv = out_y0_agri_YXv ;
for i = 1:size(landArea_2deg_YX,1)
    for j = 1:size(landArea_2deg_YX,2)
        %             disp(['i = ' num2str(i) ', j = ' num2str(j)])
        iy = 4*i-(0:3) ;
        ix = 4*j-3:4*j ;
        [IY,IX] = meshgrid(iy,ix);
        tmp = reshape(cat(2,IY',IX'),[],2) ;
        theseCells = sub2ind(size(landArea_YX),tmp(:,1),tmp(:,2)) ;
        clear tmp
        
        agri_d = shiftdim(out_y1_2deg_agri_YXv(i,j,:) - out_y0_2deg_agri_YXv(i,j,:)) ;
        if isequal(agri_d,zeros(Nagri,1))
            out_y1_agri_YXv(iy,ix,:) = out_y0_agri_YXv(iy,ix,:) ;
            continue
        end
        luh2_vegd_theseCells = luh2_vegd_YX(theseCells) ;
        
        
        if ~update_avail_land
            now_agri_YX = sum(out_y0_agri_YXv(iy,ix,:),3) ;
%             avail_land = in_y1_vegd_YX(theseCells) - now_agri_YX(:) ;
            avail_land = luh2_vegd_theseCells - now_agri_YX(:) ;
        end
        
        prev_done = false ;
        for k = 1:Nagri
            
            if debug_ijk(1)==i && debug_ijk(2)==j && debug_ijk(3)==k
                keyboard
            end
            
            if update_avail_land
                now_agri_YX = sum(out_y1_agri_YXv,3) ;
%                 avail_land = in_y1_vegd_YX(theseCells) - now_agri_YX(theseCells) ;
                avail_land = luh2_vegd_theseCells - now_agri_YX(theseCells) ;
            end
            
            % Skip if done previously
            if prev_done
                prev_done = false ;
                out_y1_this_theseCells = out_y1_agri_YXv(iy,ix,k) ;
                out_y1_this_theseCells = out_y1_this_theseCells(:) ;
            else
                % Get delta for this land use. If 0, skip. Otherwise, set
                % up for rest of logic.
                this_d = agri_d(k) ;
                if this_d==0
                    out_y1_agri_YXv(iy,ix,k) = out_y0_agri_YXv(iy,ix,k) ;
                    continue
                end
                out_y0_this_YX = out_y0_agri_YXv(iy,ix,k) ;
                out_y0_2deg_this = out_y0_2deg_agri_YXv(i,j,k) ;
                out_y1_this_theseCells = [] ;
                
                % Update this LU
                if this_d < 0
                    this_p = this_d / (out_y0_2deg_this + 1e-12) ;
                    tmp = out_y0_this_YX * (1 + this_p) ;
                    out_y1_this_theseCells = reshape(tmp,size(out_y1_agri_YXv(iy,ix,k))) ;
                elseif this_d > 0
                    if sum(avail_land) < this_d ...
                            && k < Nagri ...
                            && agri_d(k+1) < 0
                        % Go ahead and do next land use
                        prev_done = true ;
                        out_y0_next_YX = out_y0_agri_YXv(iy,ix,k+1) ;
                        out_y0_2deg_next_YX = out_y0_2deg_agri_YXv(i,j,k+1) ;
                        next_d = agri_d(k+1) ;
                        next_p = next_d / (out_y0_2deg_next_YX + 1e-12) ;
                        tmp = out_y0_next_YX * (1 + next_p) ;
                        out_y1_agri_YXv(iy,ix,k+1) = reshape(tmp,size(out_y1_agri_YXv(iy,ix,k))) ;
                        % How much NATURAL land is available now that
                        % nextLU has decreased?
% % %                         now_agri_YX = sum(out_y1_agri_YXv,3) ;
% % %                         avail_land_this = luh2_vegd_YX(theseCells) - now_agri_YX(theseCells) ;
                        now_agri_YX_this = sum(out_y1_agri_YXv,3) ;
                        avail_land_this = luh2_vegd_theseCells - now_agri_YX_this(theseCells) ;
                        % Distribute thisLU demand to half-degree gridcells
                        % based on how much NATURAL land each has
% % %                         tmp = out_y0_this_YX(:) + avail_land_this / (sum(avail_land_this) + 1e-12) * this_d ;
                        tmp = out_y0_this_YX(:) + avail_land_this / (sum(avail_land_this) + 1e-12) * this_d ;
                        out_y1_this_theseCells = reshape(tmp,size(out_y1_agri_YXv(iy,ix,k))) ;
                    else % thisLU increases and (nextLU increases OR available land is enough to satisfy thisLU demand)
                        % Distribute thisLU demand to half-degree gridcells based
                        % on how much NATURAL land each has
                        tmp = out_y0_this_YX(:) + avail_land / (sum(avail_land) + 1e-12) * this_d ;
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
            
            tmp = out_y1_this_theseCells(:) ;
            this_minus_vegd = tmp - luh2_vegd_theseCells ;
            while any(this_minus_vegd > 1e-3)
                stop
                is_too_much_crop = find(this_minus_vegd > 1e-3) ;
                sum_too_much_crop = sum(this_minus_vegd(is_too_much_crop)) ;
                tmp( is_too_much_crop) = luh2_vegd_theseCells(is_too_much_crop) ;
                avail_land = in_y1_vegd_YX(theseCells) - out_y1_crop_YX(theseCells) - out_y0_past_YX(theseCells) ;
                tmp(~is_too_much_crop) = out_y0_crop_YX(~is_too_much_crop) + avail_land(~is_too_much_crop) / (sum(sum(avail_land)) + 1e-12) * sum_too_much_crop ;
                this_minus_vegd = tmp - luh2_vegd_theseCells ;
            end
            
            % Check for invalid cell areas
            if any(tmp < -1e-3)
                error('Negative members of half-deg out_y1_this_YXv!')
            elseif any(tmp > 1e-3+luh2_vegd_theseCells)
                error('Members >vegd_area in half-deg out_y1_this_YX!')
            end
            
            % Save this (if previously done: redundant but safe)
            out_y1_this_theseCells = reshape(tmp,size(out_y1_agri_YXv(iy,ix,k))) ;
            out_y1_agri_YXv(iy,ix,k) = out_y1_this_theseCells ;
            
            
        end % for k = 1:Nagri
        
        % Check for invalid cell areas
        agri_YX = sum(out_y1_agri_YXv(iy,ix,:),3) ;
        if any(agri_YX(:) > 1e-3+luh2_vegd_YX(theseCells))
            error('Members >vegd_area in half-deg out_y1_(crop+past)_YX!')
        end
        
        
        % Debugging step 4
        for k = 1:Nagri
            this_d_theseCells_halfDeg_4 = sum(sum(out_y1_agri_YXv(iy,ix,k) - out_y0_agri_YXv(iy,ix,k))) ;
            if abs((this_d_theseCells_halfDeg_4-agri_d(k))/agri_d(k)*100) > conserv_tol_pct ...
                    && agri_d(k) > 1 % km2
                error(['Global ' LUnames{i} ' area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 4)'])
            end
        end
    end
end


end