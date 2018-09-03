%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LUH1-style harmonization for PLUM outputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_year = 2010 ;
lastYear = 2012 ;

do_debug = false ;

PLUM_in_top = removeslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1') ;

outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;

conserv_tol_pct = 0.1 ;


%% Setup

PLUM_out_top = addslashifneeded([PLUM_in_top '.harm/']) ;
unix(['mkdir -p ' PLUM_out_top]) ;
PLUM_base_in = addslashifneeded([PLUM_in_top '/' num2str(base_year)]) ;
PLUM_in_top = addslashifneeded(PLUM_in_top) ;

addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;

% Make lower-left lat/lon map (for compat. with PLUM style)
lons_map_2deg = repmat(-180:2:178,[90 1]) ;
lats_map_2deg = repmat((-90:2:88)',[1 180]) ;
lons_map = repmat(-180:0.5:179.5,[360 1]) ;
lats_map = repmat((-90:0.5:89.5)',[1 720]) ;


%% Import and aggregate reference LU map

disp('Importing reference LU map...')

% Get LUH2 land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
landArea_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear *_YXqd

% Import LUH2 base_year
luh2_file = '/Users/Shared/PLUM/input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
luh2 = lpjgu_matlab_readTable_then2map(luh2_file) ;
if ~isempty(find(luh2.maps_YXvy(:,:,contains(luh2.varNames,{'URBAN','PEATLAND'}),:)>0,1))
    error('This code is not designed to handle LUH2 inputs with any URBAN or PEATLAND area!')
end
luh2_base.maps_YXv = luh2.maps_YXvy(:,:,~contains(luh2.varNames,{'URBAN','PEATLAND'}),luh2.yearList==base_year) ;
luh2_base.varNames = luh2.varNames(~contains(luh2.varNames,{'URBAN','PEATLAND'})) ;
luh2_base.list_to_map = luh2.list_to_map ;
LUnames = luh2_base.varNames ;
Nlu = length(LUnames) ;
notBare = ~strcmp(LUnames,'BARREN') ;
clear luh2

% Harmonize LUH2 mask and PLUM mask
file_in = [PLUM_base_in 'LandCoverFract.txt'] ;
S = lpjgu_matlab_readTable_then2map(file_in,'verboseIfNoMat',false,'force_mat_nosave',true) ;
mask_YX = isnan(S.maps_YXv(:,:,1)) ...
    | sum(S.maps_YXv(:,:,contains(S.varNames,{'CROPLAND','PASTURE','NATURAL'})),3)==0 ...
    | landArea_YX==0 ;
landArea_YX(mask_YX) = 0 ;
clear S

% Get repmat 0.5º land area
landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;

% Convert from "fraction of land" to "land area (km2)"
% (This also masks where needed due to harmonization of LUH2+PLUM masks)
luh2_base.maps_YXv = luh2_base.maps_YXv .* landArea_YXv ;

% Convert NaNs to zeros for addition
luh2_base.maps_YXv(isnan(luh2_base.maps_YXv)) = 0 ;

% Get LUH2 fraction that is vegetated, barren
luh2_vegd_YX = sum(luh2_base.maps_YXv(:,:,notBare),3) ;
luh2_vegdFrac_YX = luh2_vegd_YX ./ landArea_YX ;
luh2_bare_YX = luh2_base.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
luh2_bareFrac_YX = luh2_bare_YX ./ landArea_YX ;

% Aggregate from 0.5º to 2º
tmp = luh2_base.maps_YXv(:,1:4:end,:) ...
    + luh2_base.maps_YXv(:,2:4:end,:) ...
    + luh2_base.maps_YXv(:,3:4:end,:) ...
    + luh2_base.maps_YXv(:,4:4:end,:) ;
luh2_base_2deg.maps_YXv = tmp(1:4:end,:,:) ...
                        + tmp(2:4:end,:,:) ...
                        + tmp(3:4:end,:,:) ...
                        + tmp(4:4:end,:,:) ;
clear tmp
luh2_base_2deg.varNames = luh2_base.varNames ;
tmp = landArea_YX(:,1:4:end) ...
    + landArea_YX(:,2:4:end) ...
    + landArea_YX(:,3:4:end) ...
    + landArea_YX(:,4:4:end) ;
landArea_2deg_YX = tmp(1:4:end,:) ...
                 + tmp(2:4:end,:) ...
                 + tmp(3:4:end,:) ...
                 + tmp(4:4:end,:) ;
clear tmp

% Get LUH2 2-deg fraction that is vegetated, barren
luh2_vegd_2deg_YX = sum(luh2_base_2deg.maps_YXv(:,:,notBare),3) ;
luh2_vegdFrac_2deg_YX = luh2_vegd_2deg_YX ./ landArea_2deg_YX ;
luh2_bare_2deg_YX = luh2_base_2deg.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
luh2_bareFrac_2deg_YX = luh2_bare_2deg_YX ./ landArea_2deg_YX ;
             
% Get lats/lons
list2map_2deg = find(landArea_2deg_YX>0) ;
lons_2deg = lons_map_2deg(list2map_2deg) ;
lats_2deg = lats_map_2deg(list2map_2deg) ;
list2map = find(landArea_YX>0) ;
lons = lons_map(list2map) ;
lats = lats_map(list2map) ;

disp('Done importing reference LU map.')




%% Process PLUM outputs

% The years we want to produce PLUM outputs for (will begin with
% transitions from years(1)-1 to years(1)
years = base_year+(1:(lastYear-base_year)) ;

out_y0_crop_YX = luh2_base.maps_YXv(:,:,strcmp(LUnames,'CROPLAND')) ;
out_y0_past_YX = luh2_base.maps_YXv(:,:,strcmp(LUnames,'PASTURE')) ;
out_y0_2deg_crop_YX = luh2_base_2deg.maps_YXv(:,:,strcmp(LUnames,'CROPLAND')) ;
out_y0_2deg_past_YX = luh2_base_2deg.maps_YXv(:,:,strcmp(LUnames,'PASTURE')) ;
out_y0_2deg_agri_YXv = cat(3,out_y0_2deg_crop_YX,out_y0_2deg_past_YX) ;
vegdFrac_y0_YX = luh2_vegdFrac_YX ;
bareFrac_y0_YX = luh2_bareFrac_YX ;
for y = 1:length(years)
    
    thisYear = years(y) ;
    disp(num2str(thisYear)) ;
    
    % Import and process previous year, if needed
    if y==1
        file_in = [PLUM_base_in '/LandCoverFract.txt'] ;
        [in_y0, in_y0_2deg] = ...
            PLUMharm_processPLUMin(file_in, landArea_YX, LUnames, bareFrac_y0_YX) ;
        bareFrac_y0_YX = in_y0.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ./ landArea_YX ;
        in_y0_crop_YX = in_y0.maps_YXv(:,:,strcmp(in_y0.varNames,'CROPLAND')) ;
        in_y0_past_YX = in_y0.maps_YXv(:,:,strcmp(in_y0.varNames,'PASTURE')) ;
    end
    
    % Import this year and convert to area
    file_in = [PLUM_in_top num2str(thisYear) '/LandCoverFract.txt'] ;
    [in_y1, in_y1_2deg] = ...
        PLUMharm_processPLUMin(file_in, landArea_YX, LUnames, bareFrac_y0_YX) ;
    in_y1_crop_YX = in_y1.maps_YXv(:,:,strcmp(in_y1.varNames,'CROPLAND')) ;
    in_y1_past_YX = in_y1.maps_YXv(:,:,strcmp(in_y1.varNames,'PASTURE')) ;
    in_y1_ntrl_YX = in_y1.maps_YXv(:,:,strcmp(in_y1.varNames,'NATURAL')) ;
    in_y1_bare_YX = in_y1.maps_YXv(:,:,strcmp(in_y1.varNames,'BARREN')) ;
    in_y1_vegd_YX = sum(in_y1.maps_YXv(:,:,notBare),3) ;
    in_y1_vegdFrac_YX = in_y1_vegd_YX ./ landArea_YX ;
    in_y1_bareFrac_YX = in_y1_bare_YX ./ landArea_YX ;
    in_y1_2deg_vegd_YX = sum(in_y1_2deg.maps_YXv(:,:,notBare),3) ;
        
    % Import for debugging redistribution
    if do_debug
        if y==1
            file_in = [PLUM_base_in '/LandCoverFract.txt'] ;
            [~, in_y0orig_2deg] = ...
                PLUMharm_processPLUMin(file_in, landArea_YX, LUnames, []) ;
        end
        in_y1orig_2deg = in_y1_2deg ;
        diffO_YXv = in_y1orig_2deg.maps_YXv - in_y0orig_2deg.maps_YXv ;
    end
    
    % Calculate changes in PLUM crop and pasture grids at 2 degrees
    crop_d_YX = in_y1_2deg.maps_YXv(:,:,strcmp(LUnames,'CROPLAND')) ...
              - in_y0_2deg.maps_YXv(:,:,strcmp(LUnames,'CROPLAND')) ;
    past_d_YX = in_y1_2deg.maps_YXv(:,:,strcmp(LUnames,'PASTURE')) ...
              - in_y0_2deg.maps_YXv(:,:,strcmp(LUnames,'PASTURE')) ;
          
    % Apply changes to previous grid (@2deg)
    % (if y==1, out_y0_2deg_crop_YX was specified above loop)
    mid_y1_2deg_crop_YX = out_y0_2deg_crop_YX + crop_d_YX ;
    mid_y1_2deg_past_YX = out_y0_2deg_past_YX + past_d_YX ;
    
    % Check that neither crop nor pasture exceeds nonBare area;
    % check that neither crop nor pasture is negative;
    % check that nonBare area is conserved;
    % compute the total amount of crop or pasture increase / decrease that
    % is not able to be met within the 2 degree gridcells
% %     [total_unmet_crop_YX, total_unmet_past_YX, ...
% %         mid_y1_2deg_crop_YX, mid_y1_2deg_past_YX, mid_y1_2deg_ntrl_YX] = ...
% %         PLUMharm_getUnmet(...
% %         mid_y1_2deg_crop_YX, mid_y1_2deg_past_YX, ...
% %         ...in_y1_2deg_vegd_YX) ;
% %         luh2_vegd_2deg_YX) ;
    mid_y1_2deg_agri_YXv = cat(3,mid_y1_2deg_crop_YX,mid_y1_2deg_past_YX) ;
    [total_unmet_agri_YXv, ...
        mid_y1_2deg_agri_YXv, mid_y1_2deg_ntrl_YX] = ...
        PLUMharm_getUnmet_cropArea(...
        mid_y1_2deg_agri_YXv, ...
        ...in_y1_2deg_vegd_YX) ;
        luh2_vegd_2deg_YX) ;
    total_unmet_crop_YX = total_unmet_agri_YXv(:,:,1) ;
    total_unmet_past_YX = total_unmet_agri_YXv(:,:,2) ;
    clear total_unmet_agri_YXv
    mid_y1_2deg_crop_YX = mid_y1_2deg_agri_YXv(:,:,1) ;
    mid_y1_2deg_past_YX = mid_y1_2deg_agri_YXv(:,:,2) ;
    clear mid_y1_2deg_agri_YXv
    
    % Debugging step 2: Check that global area changes are (mostly) conserved
    crop_d2_YX = mid_y1_2deg_crop_YX - out_y0_2deg_crop_YX + total_unmet_crop_YX ;
    crop_d_glob_1 = sum(crop_d_YX(:)) ;
    crop_d_glob_2 = sum(crop_d2_YX(:)) ;
    if abs((crop_d_glob_2-crop_d_glob_1)/crop_d_glob_1*100) > conserv_tol_pct
        error(['Global crop area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 2)'])
    end
    past_d2_YX = mid_y1_2deg_past_YX - out_y0_2deg_past_YX + total_unmet_past_YX ;
    past_d_glob_1 = sum(past_d_YX(:)) ;
    past_d_glob_2 = sum(past_d2_YX(:)) ;
    if abs((past_d_glob_2-past_d_glob_1)/past_d_glob_1*100) > conserv_tol_pct
        error(['Global past area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 2)'])
    end
        
    % Loop through every 2-degree gridcell. If a gridcell has unmet crop
    % or pasture, look for place to put this unmet amount in neighboring
    % rings, starting with gridcells that are 1 unit away, then 2, etc.
    % until all unmet has been displaced to new 2 degree cells. Track the
    % displaced crop and pasture and the # of "rings" needed for each
    % 2-degree gridcell.
% % %     [out_y1_2deg_crop_YX, out_y1_2deg_past_YX] = PLUMharm_ringRedist(...
% % %         mid_y1_2deg_crop_YX, mid_y1_2deg_past_YX, ...
% % %         ...in_y1_2deg_vegd_YX, ...
% % %         luh2_vegd_2deg_YX, ...
% % %         total_unmet_crop_YX, total_unmet_past_YX, ...
% % %         landArea_2deg_YX, do_debug) ;
    mid_y1_2deg_agri_YXv = cat(3,mid_y1_2deg_crop_YX, mid_y1_2deg_past_YX) ;
    total_unmet_agri_YXv = cat(3,total_unmet_crop_YX, total_unmet_past_YX) ;
    out_y1_2deg_agri_YXv = ...
        PLUMharm_ringRedist_areaCrops(...
        mid_y1_2deg_agri_YXv, ...
        ...in_y1_2deg_vegd_YX, ...
        luh2_vegd_2deg_YX, ...
        total_unmet_agri_YXv, ...
        landArea_2deg_YX, do_debug) ;
    out_y1_2deg_crop_YX = out_y1_2deg_agri_YXv(:,:,1) ;
    out_y1_2deg_past_YX = out_y1_2deg_agri_YXv(:,:,2) ;
    diffH_YXv = out_y1_2deg_agri_YXv - out_y0_2deg_agri_YXv ;
    clear out_y1_2deg_agri_YXv
    
    % Debugging step 3: Check that global area changes are (mostly) conserved
    crop_d3_YX = out_y1_2deg_crop_YX - out_y0_2deg_crop_YX  ;
    crop_d_glob_3 = sum(crop_d3_YX(:)) ;
    if abs((crop_d_glob_3-crop_d_glob_1)/crop_d_glob_1*100) > conserv_tol_pct
        error(['Global crop area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 3)'])
    end
    past_d3_YX = out_y1_2deg_past_YX - out_y0_2deg_past_YX  ;
    past_d_glob_3 = sum(past_d3_YX(:)) ;
    if abs((past_d_glob_3-past_d_glob_1)/past_d_glob_1*100) > conserv_tol_pct
        error(['Global past area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 3)'])
    end
    
    
    % Loop through all 2-degree gridcells and distribute the new crop and
    % pasture changes to the half-degree gridcells within.
    %
    % If the crop or pasture change is a decrease, apply the same
    % PERCENTAGE decrease to all half-degree crop or pasture gridcells 
    % within the 2-degree cell.
    %
    % If the crop or pasture change is an increase, apply this increase to 
    % all half-degree gridcells within the 2-degree cell proportionally 
    % according to available land area. Make sure that total area within a 
    % half-degree gridcell does not exceed 1 or 0.
    %
% % %     [out_y1_crop_YX, out_y1_past_YX] = ...
% % %         PLUMharm_distDeltas( ...
% % %         landArea_YX, landArea_2deg_YX, ...
% % %         out_y0_2deg_crop_YX, out_y1_2deg_crop_YX, out_y0_crop_YX, ...
% % %         out_y0_2deg_past_YX, out_y1_2deg_past_YX, out_y0_past_YX, ...
% % %         in_y1_vegd_YX, luh2_vegd_YX, conserv_tol_pct) ;
    
    out_y0_2deg_agri_YXv = cat(3,out_y0_2deg_crop_YX,out_y0_2deg_past_YX) ;
    out_y1_2deg_agri_YXv = cat(3,out_y1_2deg_crop_YX,out_y1_2deg_past_YX) ;
    out_y0_agri_YXv = cat(3,out_y0_crop_YX,out_y0_past_YX) ;
% % %     out_y1_agri_YXv = ...
% % %         PLUMharm_distDeltas_areaCrops( ...
% % %         landArea_YX, landArea_2deg_YX, ...
% % %         out_y0_2deg_agri_YXv, out_y1_2deg_agri_YXv, out_y0_agri_YXv, ...
% % %         in_y1_vegd_YX, luh2_vegd_YX, conserv_tol_pct) ;
    out_y1_agri_YXv = ...
        PLUMharm_distDeltas_areaCrops_recursive( ...
        landArea_YX, landArea_2deg_YX, ...
        out_y0_2deg_agri_YXv, out_y1_2deg_agri_YXv, out_y0_agri_YXv, ...
        in_y1_vegd_YX, luh2_vegd_YX, conserv_tol_pct) ;
    out_y1_crop_YX = out_y1_agri_YXv(:,:,1) ;
    out_y1_past_YX = out_y1_agri_YXv(:,:,2) ;
    
    % Ensure cells in range [0,in_y1_vegd_YX] (extreme deviations were checked
    % in loop)
    out_y1_crop_YX(out_y1_crop_YX<0) = 0 ;
    out_y1_past_YX(out_y1_past_YX<0) = 0 ;
    if any(isnan(out_y1_past_YX(:)))
        error('How did you get a NaN in out_y1_past_YX?')
    end
    out_y1_crop_YX(out_y1_crop_YX>luh2_vegd_YX) = luh2_vegd_YX(out_y1_crop_YX>luh2_vegd_YX) ;
    out_y1_past_YX(out_y1_past_YX>luh2_vegd_YX) = luh2_vegd_YX(out_y1_past_YX>luh2_vegd_YX) ;
    if any(isnan(out_y1_past_YX(:)))
        error('How did you get a NaN in out_y1_past_YX?')
    end
    
    % Debugging step 5: Check that global area changes are (mostly) conserved
    crop_d1_halfDeg_YX = in_y1_crop_YX - in_y0_crop_YX ;
    crop_d5_halfDeg_YX = out_y1_crop_YX - out_y0_crop_YX  ;
    crop_d_glob_halfDeg_1 = sum(crop_d1_halfDeg_YX(:)) ;
    crop_d_glob_halfDeg_5 = sum(crop_d5_halfDeg_YX(:)) ;
    if abs((crop_d_glob_halfDeg_5-crop_d_glob_halfDeg_1)/crop_d_glob_halfDeg_1*100) > conserv_tol_pct
        error(['Global crop area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 5)'])
    end
    past_d1_halfDeg_YX = in_y1_past_YX - in_y0_past_YX ;
    past_d5_halfDeg_YX = out_y1_past_YX - out_y0_past_YX  ;
    past_d_glob_halfDeg_1 = sum(past_d1_halfDeg_YX(:)) ;
    past_d_glob_halfDeg_5 = sum(past_d5_halfDeg_YX(:)) ;
    if abs((past_d_glob_halfDeg_5-past_d_glob_halfDeg_1)/past_d_glob_halfDeg_1*100) > conserv_tol_pct
        error(['Global past area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 5)'])
    end
    
    % Get barren and natural areas
    out_y1_bare_YX = in_y1.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
    out_y1_ntrl_YX = landArea_YX ...
        - (out_y1_crop_YX + out_y1_past_YX + out_y1_bare_YX) ;
        
    % Save new file (0.5-degree)
    out_y1.maps_YXv(:,:,strcmp(LUnames,'CROPLAND')) = out_y1_crop_YX ;
    out_y1.maps_YXv(:,:,strcmp(LUnames,'PASTURE')) = out_y1_past_YX ;
    out_y1.maps_YXv(:,:,strcmp(LUnames,'NATURAL')) = out_y1_ntrl_YX ;
    out_y1.maps_YXv(:,:,strcmp(LUnames,'BARREN')) = out_y1_bare_YX ;
    out_y1.maps_YXv = out_y1.maps_YXv ./ repmat(landArea_YX,[1 1 Nlu]) ;
    out_y1.varNames = LUnames ;
    [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
    unix(['mkdir -p ' PLUM_out_top num2str(thisYear)]) ;
    file_out = [PLUM_out_top num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.txt'] ;
    lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'progress_step_pct', 20, ...
        'verbose', false) ;
    
% %     % Save new file (2-degree)
%     out_y1_2deg_ntrl_YX = in_y1_2deg_vegd_YX - (out_y1_2deg_crop_YX + out_y1_2deg_past_YX) ;
%     out_y1_2deg.maps_YXv(:,:,strcmp(LUnames,'CROPLAND')) = out_y1_2deg_crop_YX ;
%     out_y1_2deg.maps_YXv(:,:,strcmp(LUnames,'PASTURE')) = out_y1_2deg_past_YX ;
%     out_y1_2deg.maps_YXv(:,:,strcmp(LUnames,'NATURAL')) = out_y1_2deg_ntrl_YX ;
%     out_y1_2deg.maps_YXv(:,:,strcmp(LUnames,'BARREN')) = in_y1_2deg.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
%     out_y1_2deg.maps_YXv = out_y1_2deg.maps_YXv ./ repmat(landArea_2deg_YX,[1 1 Nlu]) ;
%     out_y1_2deg.varNames = LUnames ;
% %     [out_y1_2deg_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1_2deg,list2map_2deg) ;
% %     unix(['mkdir -p ' PLUM_out_top num2str(thisYear)]) ;
% %     file_out = [PLUM_out_top num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.2deg.txt'] ;
% %     lpjgu_matlab_saveTable(out_header_cell, out_y1_2deg_array, file_out,...
% %         'outPrec', outPrec, ...
% %         'outWidth', outWidth, ...
% %         'delimiter', delimiter, ...
% %         'overwrite', overwrite, ...
% %         'fancy', fancy, ...
% %         'progress_step_pct', 20, ...
% %         'verbose', false) ;
    
    % Prepare for next iteration
    if y < length(years)
        in_y0 = in_y1 ;
        in_y0_2deg = in_y1_2deg ;
        in_y0_crop_YX = in_y1_crop_YX ;
        in_y0_past_YX = in_y1_past_YX ;
        out_y0_crop_YX = out_y1_crop_YX ;
        out_y0_past_YX = out_y1_past_YX ;
        out_y0_ntrl_YX = out_y1_ntrl_YX ;
        out_y0_2deg_crop_YX = out_y1_2deg_crop_YX ;
        out_y0_2deg_past_YX = out_y1_2deg_past_YX ;
        bareFrac_y0_YX = in_y1_bareFrac_YX ;
        if do_debug
            in_y0orig_2deg = in_y1orig_2deg ;
        end
    end
    
    
end % years loop

disp('Done')






