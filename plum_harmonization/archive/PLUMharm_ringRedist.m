function [out_y1_2deg_crop_YX, out_y1_2deg_past_YX] = ...
    PLUMharm_ringRedist(...
    mid_y1_2deg_crop_YX, mid_y1_2deg_past_YX, ...
    vegd_2deg_y1_YX, ...
    total_unmet_crop_YX, total_unmet_past_YX, ...
    landArea_2deg_YX, do_debug)

% Loop through every 2-degree gridcell. If a gridcell has unmet crop
% or pasture, look for place to put this unmet amount in neighboring
% rings, starting with gridcells that are 1 unit away, then 2, etc.
% until all unmet has been displaced to new 2 degree cells. Track the
% displaced crop and pasture and the # of "rings" needed for each
% 2-degree gridcell.

% Create grids for tracking displaced crop and pasture
if do_debug
    displaced_crop_YX = zeros(size(landArea_2deg_YX)) ;
    displaced_past_YX = zeros(size(landArea_2deg_YX)) ;
    crop_rings_YX = zeros(size(landArea_2deg_YX)) ;
    past_rings_YX = zeros(size(landArea_2deg_YX)) ;
else
    displaced_crop_YX = [] ;
    displaced_past_YX = [] ;
end

out_y1_2deg_crop_YX = mid_y1_2deg_crop_YX ;
out_y1_2deg_past_YX = mid_y1_2deg_past_YX ;

ny = size(landArea_2deg_YX,1) ;
nx = size(landArea_2deg_YX,2) ;
ks = 0:(ny-1) ;
ms = 0:(nx-1) ;
ks = fliplr(ks) ;
ms = fliplr(ms) ;

for k = ks
    for m = ms
        thisCell = sub2ind(size(landArea_2deg_YX),k+1,m+1) ;
        
        % Do it for crop
        j = 0 ;
        while abs(total_unmet_crop_YX(thisCell))>1e-8
            j = j+1 ;
            if j>100
                error('Possible infinite loop in crop ring adjustments.')
            end
            % Calculate the total agricultural (crop+past) area
            out_y1_2deg_agri_YX = out_y1_2deg_crop_YX + out_y1_2deg_past_YX ;
            % Set up rings and indices
            if do_debug
                crop_rings_YX(thisCell) = j;
            end
            i_k = mod(k + (-j:j),ny) + 1 ;
            i_m = mod(m + (-j:j),nx) + 1 ;
            [I_K,I_M] = meshgrid(i_k,i_m);
            tmp = reshape(cat(2,I_K',I_M'),[],2) ;
            thisRing = sub2ind(size(total_unmet_crop_YX),tmp(:,1),tmp(:,2)) ;
            [out_y1_2deg_crop_YX, total_unmet_crop_YX, displaced_crop_YX] = ...
                PLUMharm_doRings(out_y1_2deg_crop_YX, out_y1_2deg_agri_YX, total_unmet_crop_YX, displaced_crop_YX, ...
                thisCell, thisRing, vegd_2deg_y1_YX) ;
        end % while loop: crop
        
        % Do it for past
        j = 0 ;   % Note that this differs from how it was done in LUH1 code (pasture started with j of last crop iteration)
        while abs(total_unmet_past_YX(thisCell))>1e-8
            j = j+1 ;
            if j>100
                error('Possible infinite loop in past ring adjustments.')
            end
            % Calculate the total agricultural (past+past) area
            out_y1_2deg_agri_YX = out_y1_2deg_crop_YX + out_y1_2deg_past_YX ;
            % Set up rings and indices
            if do_debug
                past_rings_YX(thisCell) = j;
            end
            i_k = mod(k + (-j:j),ny) + 1 ;
            i_m = mod(m + (-j:j),nx) + 1 ;
            [I_K,I_M] = meshgrid(i_k,i_m);
            tmp = reshape(cat(2,I_K',I_M'),[],2) ;
            thisRing = sub2ind(size(total_unmet_crop_YX),tmp(:,1),tmp(:,2)) ;
            [out_y1_2deg_past_YX, total_unmet_past_YX, displaced_past_YX] = ...
                PLUMharm_doRings(out_y1_2deg_past_YX, out_y1_2deg_agri_YX, total_unmet_past_YX, displaced_past_YX, ...
                thisCell, thisRing, vegd_2deg_y1_YX) ;
        end % while loop: past
        
%         % Debugging redistribution
%         if do_debug
%             out_y1_2deg_agri_YX = out_y1_2deg_crop_YX + out_y1_2deg_past_YX ;
%             out_y1_2deg_ntrl_YXtmp = vegd_2deg_y1_YX - out_y1_2deg_agri_YX ;
%             debugCell = [42 98] ; debugLU = 'CROPLAND' ;
%             if k==debugCell(1)-1 && m==debugCell(2)-1
%                 diffH = out_y1_2deg_past_YX(debugCell(1),debugCell(2)) - out_y0_2deg_past_YX(debugCell(1),debugCell(2)) ;
%                 diffO = diffO_YXv(debugCell(1),debugCell(2),strcmp(LUnames,'PASTURE')) ;
%                 diffsPAST = [diffO ; diffH] ;
%                 diffH = out_y1_2deg_crop_YX(debugCell(1),debugCell(2)) - out_y0_2deg_crop_YX(debugCell(1),debugCell(2)) ;
%                 diffO = diffO_YXv(debugCell(1),debugCell(2),strcmp(LUnames,'CROPLAND')) ;
%                 diffsCROP = [diffO ; diffH] ;
%                 
%                 diffH = out_y1_2deg_ntrl_YXtmp(debugCell(1),debugCell(2)) - in_y0_2deg.maps_YXv(debugCell(1),debugCell(2),strcmp(LUnames,'NATURAL')) ;
%                 diffO = diffO_YXv(debugCell(1),debugCell(2),strcmp(LUnames,'NATURAL')) ;
%                 diffsNTRL = [diffO ; diffH] ;
%                 
%                 tmp = table({'orig';'harm'},diffsPAST,diffsCROP,diffsNTRL) ;
%                 tmp.Properties.VariableNames = {'which','PASTURE','CROPLAND','NATURAL'} ;
%                 tmp
%                 keyboard
%             end
%         end
        
    end % for m
end % for k

if do_debug
%     keyboard
end


end