function [out_y1_crop_YX, out_y1_past_YX] = ...
    PLUMharm_distDeltas( ...
    landArea_YX, landArea_2deg_YX, ...
    out_y0_2deg_crop_YX, out_y1_2deg_crop_YX, out_y0_crop_YX, ...
    out_y0_2deg_past_YX, out_y1_2deg_past_YX, out_y0_past_YX, ...
    in_y1_vegd_YX, luh2_vegd_YX, conserv_tol_pct)

out_y1_crop_YX = zeros(size(landArea_YX)) ;
out_y1_past_YX = zeros(size(landArea_YX)) ;
for i = 1:size(landArea_2deg_YX,1)
    for j = 1:size(landArea_2deg_YX,2)
        %             disp(['i = ' num2str(i) ', j = ' num2str(j)])
        iy = 4*i-(0:3) ;
        ix = 4*j-3:4*j ;
        [IY,IX] = meshgrid(iy,ix);
        tmp = reshape(cat(2,IY',IX'),[],2) ;
        theseCells = sub2ind(size(landArea_YX),tmp(:,1),tmp(:,2)) ;
        clear tmp
        
        crop_d = out_y1_2deg_crop_YX(i,j) - out_y0_2deg_crop_YX(i,j);
        past_d = out_y1_2deg_past_YX(i,j) - out_y0_2deg_past_YX(i,j);
        avail_land = in_y1_vegd_YX(theseCells) - out_y0_crop_YX(theseCells) - out_y0_past_YX(theseCells);
        
        %%%%%%%%%%%%
        %%% Crop %%%
        %%%%%%%%%%%%
        % if crop decreases
        if crop_d <= 0
            crop_p = crop_d / (out_y0_2deg_crop_YX(i,j) + 1e-12);
            out_y1_crop_YX(theseCells) = out_y0_crop_YX(theseCells) * (1 + crop_p);
        else % crop_d >0
            % if crop increases and pasture decreases and available land is less than crop demand
            if sum(sum(avail_land))<sum(sum(crop_d)) && past_d<0
                % Go ahead and do the pasture decrease
                past_p = past_d / (out_y0_2deg_past_YX(i,j) + 1e-12);
                out_y1_past_YX(theseCells) = out_y0_past_YX(theseCells) * (1 + past_p);
                if any(isnan(out_y1_past_YX(:)))
                    error('How did you get a NaN in out_y1_past_YX?')
                end
                % How much NATURAL land is available now that pasture
                % has decreased?
                avail_land_crop = luh2_vegd_YX(theseCells) - out_y0_crop_YX(theseCells) - out_y1_past_YX(theseCells);
                % Distribute crop demand to half-degree gridcells based
                % on how much NATURAL land each has
                out_y1_crop_YX(theseCells) = out_y0_crop_YX(theseCells) + avail_land_crop / (sum(sum(avail_land_crop)) + 1e-12) * crop_d ;
            else % crop increases and pasture increases or available land is enough to satisfy crop demand
                % Distribute crop demand to half-degree gridcells based
                % on how much NATURAL land each has
                out_y1_crop_YX(theseCells) = out_y0_crop_YX(theseCells) + avail_land / (sum(sum(avail_land)) + 1e-12) * crop_d ;
            end
            
            
            % Where is there more cropland than vegetated land?
            out_y1_crop_theseCells = out_y1_crop_YX(theseCells) ;
            luh2_vegd_theseCells = luh2_vegd_YX(theseCells) ;
            crop_minus_vegd = out_y1_crop_theseCells - luh2_vegd_theseCells ;
            while any(crop_minus_vegd > 1e-3)
                is_too_much_crop = find(crop_minus_vegd > 1e-3) ;
                % What is the total amount of excess cropland?
                sum_too_much_crop = sum(crop_minus_vegd(is_too_much_crop)) ;
                % Cap cropland at vegetated area
                out_y1_crop_theseCells( is_too_much_crop) = luh2_vegd_theseCells(is_too_much_crop) ;
                % Recalculate available land
                avail_land = in_y1_vegd_YX(theseCells) - out_y1_crop_YX(theseCells) - out_y0_past_YX(theseCells) ;
                % Distribute excess cropland to cells that didn't have
                % excess cropland, proportionally according to how much
                % available land there is in each cell.
                out_y1_crop_theseCells(~is_too_much_crop) = out_y0_crop_YX(~is_too_much_crop) + avail_land(~is_too_much_crop) / (sum(sum(avail_land)) + 1e-12) * sum_too_much_crop ;
                crop_minus_vegd = out_y1_crop_theseCells - luh2_vegd_theseCells ;
            end
            out_y1_crop_YX(theseCells) = out_y1_crop_theseCells ;
            
            
            
        end
        
        %%%%%%%%%%%%
        %%% Past %%%
        %%%%%%%%%%%%
        % if pasture decreases and we didn't have the situation where
        % crop was supposed to decrease but there wasn't enough NATURAL
        % land to supply it. If pasture decreases but we DID have that
        % situation, then we've already (above) updated pasture.
        if past_d<0 && ((crop_d<=0) || (sum(sum(avail_land)) >= sum(sum(crop_d))))
            past_p = past_d / (out_y0_2deg_past_YX(i,j) + 1e-12);
            out_y1_past_YX(theseCells) = out_y0_past_YX(theseCells) * (1 + past_p);
            if any(isnan(out_y1_past_YX(:)))
                error('How did you get a NaN in out_y1_past_YX?')
            end
            % elseif pasture increases
        elseif past_d >= 0
            avail_land = luh2_vegd_YX(theseCells) - out_y1_crop_YX(theseCells) - out_y0_past_YX(theseCells);
            out_y1_past_YX(theseCells) = out_y0_past_YX(theseCells) + avail_land / (sum(sum(avail_land)) + 1e-12) * past_d ;
            if any(isnan(out_y1_past_YX(:)))
                error('How did you get a NaN in out_y1_past_YX?')
            end
            out_y1_past_theseCells = out_y1_past_YX(theseCells) ;
            luh2_vegd_theseCells = luh2_vegd_YX(theseCells) ;
            past_minus_vegd = out_y1_past_theseCells - luh2_vegd_theseCells ;
            while any(past_minus_vegd > 1e-3)
                is_too_much_past = find(past_minus_vegd > 1e-3) ;
                sum_too_much_past = sum(past_minus_vegd(is_too_much_past)) ;
                out_y1_past_theseCells( is_too_much_past) = luh2_vegd_theseCells(is_too_much_past) ;
                avail_land = in_y1_vegd_YX(theseCells) - out_y1_past_YX(theseCells) - out_y0_past_YX(theseCells) ;
                out_y1_past_theseCells(~is_too_much_past) = out_y0_past_YX(~is_too_much_past) + avail_land(~is_too_much_past) / (sum(sum(avail_land)) + 1e-12) * sum_too_much_past ;
                past_minus_vegd = out_y1_past_theseCells - luh2_vegd_theseCells ;
            end
            out_y1_past_YX(theseCells) = out_y1_past_theseCells ;
            if any(isnan(out_y1_past_YX(:)))
                error('How did you get a NaN in out_y1_past_YX?')
            end
        elseif ~(past_d<0 && sum(sum(avail_land))<sum(sum(crop_d)))
            error('How did this happen?')
        end
        
        % Check for invalid cell areas
        if any(out_y1_crop_YX(theseCells) < -1e-3)
            error('Negative members of half-deg out_y1_crop_YX!')
        elseif any(out_y1_past_YX(theseCells) < -1e-3)
            error('Negative members of half-deg out_y1_past_YX!')
        end
        if any(out_y1_crop_YX(theseCells) > 1e-3+luh2_vegd_YX(theseCells))
            error('Members >vegd_area in half-deg out_y1_crop_YX!')
        elseif any(out_y1_past_YX(theseCells) > 1e-3+luh2_vegd_YX(theseCells))
            error('Members >vegd_area in half-deg out_y1_past_YX!')
        end
        if any(out_y1_crop_YX(theseCells)+out_y1_past_YX(theseCells) > 1e-3+luh2_vegd_YX(theseCells))
            error('Members >vegd_area in half-deg out_y1_(crop+past)_YX!')
        end
        
        % Debugging step 4
        crop_d_theseCells_halfDeg_4 = sum(out_y1_crop_YX(theseCells)-out_y0_crop_YX(theseCells)) ;
        if abs((crop_d_theseCells_halfDeg_4-crop_d)/crop_d*100) > conserv_tol_pct ...
                && crop_d > 1 % km2
            error(['Global crop area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 4)'])
        end
        past_d_theseCells_halfDeg_4 = sum(out_y1_past_YX(theseCells)-out_y0_past_YX(theseCells)) ;
        if abs((past_d_theseCells_halfDeg_4-past_d)/past_d*100) > conserv_tol_pct ...
                && past_d > 1 % km2
            error(['Global past area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 4)'])
        end
    end
end



end