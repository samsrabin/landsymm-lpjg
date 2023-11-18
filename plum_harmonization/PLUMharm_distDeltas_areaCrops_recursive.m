function out_y1_agri_YXv = ...
    PLUMharm_distDeltas_areaCrops_recursive( ...
    landArea_YX, landArea_2deg_YX, ...
    out_y0_2deg_agri_YXv, out_y1_2deg_agri_YXv, out_y0_agri_YXv, ...
    out_y0_vegd_YX, conserv_tol_pct, conserv_tol_area, LUnames, ...
    debug_ijk)
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

% Properly correct for zero-denominators? FALSE is default LUH1 behavior,
% which just adds 1e-12 to denominator.
proper_zero_denoms = true ;


do_debug = ~isempty(debug_ijk) ;
if do_debug
    last_IJd = Inf ;
    this_IJd = -Inf ;
    dbI = debug_ijk(1) ;
    dbJ = debug_ijk(2) ;
    disp('PLUMharm_distDeltas_areaCrops_recursive:')
else
    debug_ijk = [Inf Inf Inf] ;
end

Nagri = size(out_y0_agri_YXv,3) ;
max_diff = 0 ;

out_y1_agri_YXv = out_y0_agri_YXv ;
for i = 1:size(landArea_2deg_YX,1)
    for j = 1:size(landArea_2deg_YX,2)
        iy = 4*i-(0:3) ;
        ix = 4*j-3:4*j ;
        [IY,IX] = meshgrid(iy,ix);
        tmp = reshape(cat(2,IY',IX'),[],2) ;
        theseCells = sub2ind(size(landArea_YX),tmp(:,1),tmp(:,2)) ;
        clear tmp
        
        
        % debugging
        if do_debug
            if ~isequal(last_IJd,this_IJd)
                this_IJd = shiftdim(out_y1_2deg_agri_YXv(dbI,dbJ,:) - out_y0_2deg_agri_YXv(dbI,dbJ,:)) ;
                if isinf(last_IJd)
                    last_IJd = this_IJd ;
                end
                fprintf(['   agri_d(%d,%d): ' repmat('%0.3e ',[1 Nagri]) '\n'],dbI,dbJ,this_IJd) ;
            end
%             if i==debug_ijk(1) && j==debug_ijk(2)
%                 keyboard
%             end
        end
        
        agri_d = shiftdim(out_y1_2deg_agri_YXv(i,j,:) - out_y0_2deg_agri_YXv(i,j,:)) ;
        if isequal(agri_d,zeros(Nagri,1))
            out_y1_agri_YXv(iy,ix,:) = out_y0_agri_YXv(iy,ix,:) ;
            continue
        end
        out_y0_vegd_theseCells = out_y0_vegd_YX(theseCells) ;        
        
        if update_avail_land
            now_agri_YX = sum(out_y0_agri_YXv(iy,ix,:),3) ;
            avail_land = out_y0_vegd_theseCells - now_agri_YX(:) ;

            % Ensure that there is SOME available land if there's a net positive change in
            % agricultural area
            if sum(avail_land) == 0 && sum(agri_d) > 0
                msg = sprintf('(%d,%d): Positive sum(agri_d) but no avail_land', ...
                    i, j) ;
                % If we're talking about very small changes in area, just ignore
                if sum(abs(agri_d)) < conserv_tol_area
                    out_y1_agri_YXv(iy,ix,:) = out_y0_agri_YXv(iy,ix,:) ;
                    warning([msg '. All sum(abs(agri_d)) < conserv_tol_area %g, so setting all agri_d to 0.'], ...
                        conserv_tol_area)
                    continue
                % Otherwise, throw an error
                else
                    warning([msg ' and some abs(agri_d) >= conserv_tol_area %g; ' ...
                        'this might result in loop_thru_agri() being unable to ' ...
                        'allocate deltas\n'], ...
                        conserv_tol_area)
                end
            end
        else
            avail_land = [] ;
        end

        % Loop through all land uses (may be recursive)
        already_done = false(1,Nagri) ;
        [out_y1_agri_YXv, ~] = ...
            loop_thru_agri(already_done, 1, debug_ijk, ...
            update_avail_land, proper_zero_denoms, conserv_tol_area, ...
            i, j, iy, ix, theseCells, out_y0_vegd_theseCells, agri_d, ...
            out_y0_agri_YXv, out_y0_2deg_agri_YXv, out_y1_agri_YXv, avail_land, ...
            LUnames(1:end-2)) ;
        
        % Check for invalid cell areas
        agri_YX = sum(out_y1_agri_YXv(iy,ix,:),3) ;
        if any(agri_YX(:) > conserv_tol_area+out_y0_vegd_YX(theseCells))
            error('Members >vegd_area in half-deg out_y1_(crop+past)_YX!')
        end
        max_diff = max(max_diff,max(agri_YX(:)-out_y0_vegd_YX(theseCells))) ;
        
        % Debugging step 4
        for k = 1:Nagri
            this_d_theseCells_halfDeg_4 = sum(sum(out_y1_agri_YXv(iy,ix,k) - out_y0_agri_YXv(iy,ix,k))) ;
            if abs((this_d_theseCells_halfDeg_4-agri_d(k))/agri_d(k)*100) > conserv_tol_pct ...
                    && agri_d(k) > 1e6 % m2
                error(['Global ' LUnames{i} ' area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 4)'])
            end
        end; clear k
    end
end

% disp(['   max_diff = ' num2str(max_diff)])


end


