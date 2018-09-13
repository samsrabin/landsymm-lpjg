function [out_y1_2deg_agri_YXv, out_y1_2deg_ntrl_YX] = ...
    PLUMharm_ringRedist_areaCropsRes(...
    mid_y1_2deg_agri_YXv, ...
    total_unmet_agri_YXv, ...
    landArea_2deg_YX, ...
    debugIJ, in_y0orig_2deg, in_y1orig_2deg, out_y0_2deg_agri_YXv, ...
    out_y0_2deg_ntrl_YX, resArea_2deg_YX, LUnames_agri)
% Loop through every 2-degree gridcell. If a gridcell has unmet crop
% or pasture, look for place to put this unmet amount in neighboring
% rings, starting with gridcells that are 1 unit away, then 2, etc.
% until all unmet has been displaced to new 2 degree cells. Track the
% displaced crop and pasture and the # of "rings" needed for each
% 2-degree gridcell.

Nagri = size(mid_y1_2deg_agri_YXv,3) ;
out_y1_2deg_ntrl_YX = out_y0_2deg_ntrl_YX ;

% Create grids for tracking displaced agriculture
% debugIJ = [999 999] ;
do_debug = ~isempty(debugIJ) ;
thisCell_ofInt = Inf ;
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
    % thisCell by any other name
    thisCell_ofInt = sub2ind(size(landArea_2deg_YX),debugIJ(1),debugIJ(2)) ;
end

out_y1_2deg_agri_YXv = mid_y1_2deg_agri_YXv ;

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
                end
                % Calculate the total agricultural (crop+past) area
                out_y1_2deg_agri_YX = sum(out_y1_2deg_agri_YXv,3) ;
                % Set up rings and indices
                this_rings_YX(thisCell) = j;
                i_k = mod(k + (-j:j),ny) + 1 ;
                i_m = mod(m + (-j:j),nx) + 1 ;
                [I_K,I_M] = meshgrid(i_k,i_m);
                tmp = reshape(cat(2,I_K',I_M'),[],2) ;
                thisRing = sub2ind(size(total_unmet_this_YX),tmp(:,1),tmp(:,2)) ;
                [out_y1_2deg_this_YX, out_y1_2deg_ntrl_YX, total_unmet_this_YX, ...
                    displaced_this_YX, this_meanDist_YX] = ...
                    PLUMharm_doRings_areaCropsRes(...
                    out_y1_2deg_this_YX, out_y1_2deg_agri_YX, total_unmet_this_YX, displaced_this_YX, ...
                    thisCell, thisRing, out_y1_2deg_ntrl_YX, resArea_2deg_YX, ...
                    this_meanDist_YX, thisCell_ofInt) ;
                out_y1_2deg_agri_YXv(:,:,i) = out_y1_2deg_this_YX ;
                meanDist_YXv(:,:,i) = this_meanDist_YX ;
                
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
    keyboard
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


