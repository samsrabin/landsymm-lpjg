function [albedo_jan_xy, albedo_jul_xy] = ...
    get_albedo_garr(fpc, snowdepth, LU, albedo_bs_vec, pftNames)

% Check inputs
check_containSame_ignoreOrder(pftNames,fpc.varNames)
luNames = {'CROPLAND','PASTURE','NATURAL','BARREN'} ;
check_containSame_ignoreOrder(luNames,LU.varNames)

% Setup
Nyears = size(fpc.garr_xvy,3) ;

% Establish NH and snow albedos (values and comments from Andy's script)
%%% cr=crop, gr=grass, ev=evergreen (extratropical), de=deciduous, 
%%% te = tropical evergreen
albedo_cr_sumr = 0.178 ; % values from Boisier et al. 2013 based on MODIS observations for s, w and snow-covered
albedo_cr_wntr = 0.141 ;
albedo_cr_snow = 0.546 ;
albedo_gr_sumr = 0.176 ;
albedo_gr_wntr = 0.161 ;
albedo_gr_snow = 0.568 ;
albedo_ev_sumr = 0.104 ;
albedo_ev_wntr = 0.094 ;
albedo_ev_snow = 0.205 ;
albedo_de_sumr = 0.153 ;
albedo_de_wntr = 0.117 ;
albedo_de_snow = 0.244 ;
albedo_te      = 0.140 ; % very different from albedo of temperate ev forest, value chosen based on the map in Boisier et al.
albedo_te_snow = 0.205 ;
albedo_bs_snow = 0.535 ;  

% Make maps (flip winter/summer in S. Hemisphere, as in Andy's script)
is_sh = false(360,720) ;
is_sh(1:180,:) = true ;
albedo_cr_jul_vec = albedo_cr_sumr*ones(size(fpc.garr_xvy, 1), 1) ;
albedo_cr_jan_vec = albedo_cr_wntr*ones(size(fpc.garr_xvy, 1), 1) ;
albedo_gr_jul_vec = albedo_gr_sumr*ones(size(fpc.garr_xvy, 1), 1) ;
albedo_gr_jan_vec = albedo_gr_wntr*ones(size(fpc.garr_xvy, 1), 1) ;
albedo_ev_jul_vec = albedo_ev_sumr*ones(size(fpc.garr_xvy, 1), 1) ;
albedo_ev_jan_vec = albedo_ev_wntr*ones(size(fpc.garr_xvy, 1), 1) ;
albedo_de_jul_vec = albedo_de_sumr*ones(size(fpc.garr_xvy, 1), 1) ;
albedo_de_jan_vec = albedo_de_wntr*ones(size(fpc.garr_xvy, 1), 1) ;
albedo_cr_jul_vec(is_sh) = albedo_cr_wntr ;
albedo_cr_jan_vec(is_sh) = albedo_cr_sumr ;
albedo_gr_jul_vec(is_sh) = albedo_gr_wntr ;
albedo_gr_jan_vec(is_sh) = albedo_gr_sumr ;
albedo_ev_jul_vec(is_sh) = albedo_ev_wntr ;
albedo_ev_jan_vec(is_sh) = albedo_ev_sumr ;
albedo_de_jul_vec(is_sh) = albedo_de_wntr ;
albedo_de_jan_vec(is_sh) = albedo_de_sumr ;
% albedo_te_vec = albedo_te*ones(size(fpc.garr_xvy, 1), 1) ;
% albedo_cr_snow_vec = albedo_cr_snow*ones(size(fpc.garr_xvy, 1), 1) ;
% albedo_gr_snow_vec = albedo_gr_snow*ones(size(fpc.garr_xvy, 1), 1) ;
% albedo_ev_snow_vec = albedo_ev_snow*ones(size(fpc.garr_xvy, 1), 1) ;
% albedo_de_snow_vec = albedo_de_snow*ones(size(fpc.garr_xvy, 1), 1) ;
% albedo_te_snow_vec = albedo_te_snow*ones(size(fpc.garr_xvy, 1), 1) ;

% LU_crop_xvy = LU.garr_xvy(:,strcmp(LU.varNames,'CROPLAND'),:) ;
% LU_past_xvy = LU.garr_xvy(:,strcmp(LU.varNames,'PASTURE'),:) ;
% LU_ntrl_xvy = LU.garr_xvy(:,strcmp(LU.varNames,'NATURAL'),:) ;
LU_vegd_xvy = sum(LU.garr_xvy(:,contains(LU.varNames,{'CROPLAND','PASTURE','NATURAL'}),:),2) ;
LU_crop_xvy = LU.garr_xvy(:,strcmp(LU.varNames,'CROPLAND'),:) ./ LU_vegd_xvy ;
LU_past_xvy = LU.garr_xvy(:,strcmp(LU.varNames,'PASTURE'),:) ./ LU_vegd_xvy ;
LU_ntrl_xvy = LU.garr_xvy(:,strcmp(LU.varNames,'NATURAL'),:) ./ LU_vegd_xvy ;

% Classify
is_cr = strcmp(fpc.varNames,'Crop_sum') ;
is_gr_ntrl = find_members(fpc.varNames,{'C3G','C4G'}) ;
is_gr_past = find_members(fpc.varNames,{'PC3G','PC4G'}) ;
% is_gr_past = contains(fpc.varNames,{'PC3G','PC4G','C3G_pas','C4G_pas'}) ;
is_ev = find_members(fpc.varNames,{'BNE','BINE','TeNE','TeBE'}) ;
is_de = find_members(fpc.varNames,{'BNS','TeBS','IBS','TrBR'}) ;
is_te = find_members(fpc.varNames,{'TrBE','TrIBE'}) ;

% Get cover maps
fpc_map_cr_xy = squeeze(sum(fpc.garr_xvy(:,is_cr,:),2) .* LU_crop_xvy) ;
fpc_map_gr_xy = squeeze(sum(fpc.garr_xvy(:,is_gr_ntrl,:),2) .* LU_ntrl_xvy ...
                       + sum(fpc.garr_xvy(:,is_gr_past,:),2) .* LU_past_xvy ...
                       ) ;
fpc_map_ev_xy = squeeze(sum(fpc.garr_xvy(:,is_ev,:),2) .* LU_ntrl_xvy) ;
fpc_map_de_xy = squeeze(sum(fpc.garr_xvy(:,is_de,:),2) .* LU_ntrl_xvy) ;
fpc_map_te_xy = squeeze(sum(fpc.garr_xvy(:,is_te,:),2) .* LU_ntrl_xvy) ;
% fpc_map_bs_xy = 1 - (fpc_map_cr_xy + fpc_map_gr_xy + fpc_map_ev_xy +...
%                       fpc_map_de_xy + fpc_map_te_xy) ;
fpc_map_bs_xy = squeeze( ...
    (1 - min(1,fpc.garr_xvy(:,strcmp(fpc.varNames,'Natural_sum'),:))) .* LU_ntrl_xvy ...
  + (1 - min(1,fpc.garr_xvy(:,strcmp(fpc.varNames,'Crop_sum'),:)))    .* LU_crop_xvy ...
  + (1 - min(1,fpc.garr_xvy(:,strcmp(fpc.varNames,'Pasture_sum'),:))) .* LU_past_xvy ...
  ) ;
                  
% Get snow fraction from snow depth (after Wang & Zeng, 2010, J. Appl.
% Meteor. Clim., Eq. 17)
snowdepth_jan_xy = squeeze(snowdepth.garr_xvy(:,1,:)) ;
snowdepth_jul_xy = squeeze(snowdepth.garr_xvy(:,7,:)) ;
snowfrac_jan_xy = snowdepth_jan_xy ./ (0.01+snowdepth_jan_xy) ;
snowfrac_jul_xy = snowdepth_jul_xy ./ (0.01+snowdepth_jul_xy) ;

% Get albedo maps
% albedo_bs_snow = 0.9 ;
% warning(['Assuming albedo of snow on bare soil is ' num2str(albedo_bs_snow) '!'])
albedo_jan_xy = (1-snowfrac_jan_xy) .* (...
                      fpc_map_cr_xy .* repmat(albedo_cr_jan_vec,[1 Nyears]) ...
                    + fpc_map_gr_xy .* repmat(albedo_gr_jan_vec,[1 Nyears]) ...
                    + fpc_map_ev_xy .* repmat(albedo_ev_jan_vec,[1 Nyears]) ...
                    + fpc_map_de_xy .* repmat(albedo_de_jan_vec,[1 Nyears]) ...
                    + fpc_map_te_xy * albedo_te ...
                    + fpc_map_bs_xy .* repmat(albedo_bs_vec,[1 Nyears]) ...
                 ) + ...
                 snowfrac_jan_xy .* (...
                      fpc_map_cr_xy * albedo_cr_snow ...
                    + fpc_map_gr_xy * albedo_gr_snow ...
                    + fpc_map_ev_xy * albedo_ev_snow ...
                    + fpc_map_de_xy * albedo_de_snow ...
                    + fpc_map_te_xy * albedo_te_snow ...
                    + fpc_map_bs_xy * albedo_bs_snow ...
                 );
albedo_jul_xy = (1-snowfrac_jul_xy) .* (...
                      fpc_map_cr_xy .* repmat(albedo_cr_jul_vec,[1 Nyears]) ...
                    + fpc_map_gr_xy .* repmat(albedo_gr_jul_vec,[1 Nyears]) ...
                    + fpc_map_ev_xy .* repmat(albedo_ev_jul_vec,[1 Nyears]) ...
                    + fpc_map_de_xy .* repmat(albedo_de_jul_vec,[1 Nyears]) ...
                    + fpc_map_te_xy * albedo_te ...
                    + fpc_map_bs_xy .* repmat(albedo_bs_vec,[1 Nyears]) ...
                 ) + ...
                 snowfrac_jul_xy .* (...
                      fpc_map_cr_xy * albedo_cr_snow ...
                    + fpc_map_gr_xy * albedo_gr_snow ...
                    + fpc_map_ev_xy * albedo_ev_snow ...
                    + fpc_map_de_xy * albedo_de_snow ...
                    + fpc_map_te_xy * albedo_te_snow ...
                    + fpc_map_bs_xy * albedo_bs_snow ...
                 );

    function check_containSame_ignoreOrder(A,B)
        % Makes sure that A and B contain the same elements, with no
        % repeats. Order does not matter.
        for i = 1:length(B)
            n = length(find(strcmp(A,B{i}))) ;
            if n ~= 1
                error([B{i} ' not found exactly once in A!'])
            end
        end
        for i = 1:length(A)
            n = length(find(strcmp(B,A{i}))) ;
            if n ~= 1
                error([A{i} ' not found exactly once in B!'])
            end
        end
    end

    function tf_out = find_members(A,B)
        % Finds indices of elements in A that match elements in B
        tf_out = false(size(A)) ;
        for i = 1:length(B)
            is_match = find(strcmp(A,B{i})) ;
            if isempty(is_match)
                error([B{i} ' not found in A!'])
            elseif length(is_match)>1
                error([B{i} ' repeated in A!'])
            end
            tf_out(is_match) = true ;
        end
        
    end


end