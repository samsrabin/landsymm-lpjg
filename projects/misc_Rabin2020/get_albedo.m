function [albedo_jan_YXy,albedo_jul_YXy] = ...
    get_albedo(fpc, snowdepth, LU, albedo_bs_map, pftNames)

% Check inputs
check_containSame_ignoreOrder(pftNames,fpc.varNames)
luNames = {'CROPLAND','PASTURE','NATURAL','BARREN'} ;
check_containSame_ignoreOrder(luNames,LU.varNames)

% Setup
Nyears = size(fpc.maps_YXvy,4) ;

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
albedo_cr_jul_map = albedo_cr_sumr*ones(360,720) ;
albedo_cr_jan_map = albedo_cr_wntr*ones(360,720) ;
albedo_gr_jul_map = albedo_gr_sumr*ones(360,720) ;
albedo_gr_jan_map = albedo_gr_wntr*ones(360,720) ;
albedo_ev_jul_map = albedo_ev_sumr*ones(360,720) ;
albedo_ev_jan_map = albedo_ev_wntr*ones(360,720) ;
albedo_de_jul_map = albedo_de_sumr*ones(360,720) ;
albedo_de_jan_map = albedo_de_wntr*ones(360,720) ;
albedo_cr_jul_map(is_sh) = albedo_cr_wntr ;
albedo_cr_jan_map(is_sh) = albedo_cr_sumr ;
albedo_gr_jul_map(is_sh) = albedo_gr_wntr ;
albedo_gr_jan_map(is_sh) = albedo_gr_sumr ;
albedo_ev_jul_map(is_sh) = albedo_ev_wntr ;
albedo_ev_jan_map(is_sh) = albedo_ev_sumr ;
albedo_de_jul_map(is_sh) = albedo_de_wntr ;
albedo_de_jan_map(is_sh) = albedo_de_sumr ;
% albedo_te_map = albedo_te*ones(360,720) ;
% albedo_cr_snow_map = albedo_cr_snow*ones(360,720) ;
% albedo_gr_snow_map = albedo_gr_snow*ones(360,720) ;
% albedo_ev_snow_map = albedo_ev_snow*ones(360,720) ;
% albedo_de_snow_map = albedo_de_snow*ones(360,720) ;
% albedo_te_snow_map = albedo_te_snow*ones(360,720) ;

% LU_crop_YXvy = LU.maps_YXvy(:,:,strcmp(LU.varNames,'CROPLAND'),:) ;
% LU_past_YXvy = LU.maps_YXvy(:,:,strcmp(LU.varNames,'PASTURE'),:) ;
% LU_ntrl_YXvy = LU.maps_YXvy(:,:,strcmp(LU.varNames,'NATURAL'),:) ;
LU_vegd_YXvy = sum(LU.maps_YXvy(:,:,contains(LU.varNames,{'CROPLAND','PASTURE','NATURAL'}),:),3) ;
LU_crop_YXvy = LU.maps_YXvy(:,:,strcmp(LU.varNames,'CROPLAND'),:) ./ LU_vegd_YXvy ;
LU_past_YXvy = LU.maps_YXvy(:,:,strcmp(LU.varNames,'PASTURE'),:) ./ LU_vegd_YXvy ;
LU_ntrl_YXvy = LU.maps_YXvy(:,:,strcmp(LU.varNames,'NATURAL'),:) ./ LU_vegd_YXvy ;

% Classify
is_cr = strcmp(fpc.varNames,'Crop_sum') ;
is_gr_ntrl = find_members(fpc.varNames,{'C3G','C4G'}) ;
is_gr_past = find_members(fpc.varNames,{'PC3G','PC4G'}) ;
% is_gr_past = contains(fpc.varNames,{'PC3G','PC4G','C3G_pas','C4G_pas'}) ;
is_ev = find_members(fpc.varNames,{'BNE','BINE','TeNE','TeBE'}) ;
is_de = find_members(fpc.varNames,{'BNS','TeBS','IBS','TrBR'}) ;
is_te = find_members(fpc.varNames,{'TrBE','TrIBE'}) ;

% Get cover maps
fpc_map_cr_YXy = squeeze(sum(fpc.maps_YXvy(:,:,is_cr,:),3) .* LU_crop_YXvy) ;
fpc_map_gr_YXy = squeeze(sum(fpc.maps_YXvy(:,:,is_gr_ntrl,:),3) .* LU_ntrl_YXvy ...
                       + sum(fpc.maps_YXvy(:,:,is_gr_past,:),3) .* LU_past_YXvy ...
                       ) ;
fpc_map_ev_YXy = squeeze(sum(fpc.maps_YXvy(:,:,is_ev,:),3) .* LU_ntrl_YXvy) ;
fpc_map_de_YXy = squeeze(sum(fpc.maps_YXvy(:,:,is_de,:),3) .* LU_ntrl_YXvy) ;
fpc_map_te_YXy = squeeze(sum(fpc.maps_YXvy(:,:,is_te,:),3) .* LU_ntrl_YXvy) ;
% fpc_map_bs_YXy = 1 - (fpc_map_cr_YXy + fpc_map_gr_YXy + fpc_map_ev_YXy +...
%                       fpc_map_de_YXy + fpc_map_te_YXy) ;
fpc_map_bs_YXy = squeeze( ...
    (1 - min(1,fpc.maps_YXvy(:,:,strcmp(fpc.varNames,'Natural_sum'),:))) .* LU_ntrl_YXvy ...
  + (1 - min(1,fpc.maps_YXvy(:,:,strcmp(fpc.varNames,'Crop_sum'),:)))    .* LU_crop_YXvy ...
  + (1 - min(1,fpc.maps_YXvy(:,:,strcmp(fpc.varNames,'Pasture_sum'),:))) .* LU_past_YXvy ...
  ) ;
                  
% Get snow fraction from snow depth (after Wang & Zeng, 2010, J. Appl.
% Meteor. Clim., Eq. 17)
snowdepth_jan_YXy = squeeze(snowdepth.maps_YXvy(:,:,1,:)) ;
snowdepth_jul_YXy = squeeze(snowdepth.maps_YXvy(:,:,7,:)) ;
snowfrac_jan_YXy = snowdepth_jan_YXy ./ (0.01+snowdepth_jan_YXy) ;
snowfrac_jul_YXy = snowdepth_jul_YXy ./ (0.01+snowdepth_jul_YXy) ;

% Get albedo maps
% albedo_bs_snow = 0.9 ;
% warning(['Assuming albedo of snow on bare soil is ' num2str(albedo_bs_snow) '!'])
albedo_jan_YXy = (1-snowfrac_jan_YXy) .* (...
                      fpc_map_cr_YXy .* repmat(albedo_cr_jan_map,[1 1 Nyears]) ...
                    + fpc_map_gr_YXy .* repmat(albedo_gr_jan_map,[1 1 Nyears]) ...
                    + fpc_map_ev_YXy .* repmat(albedo_ev_jan_map,[1 1 Nyears]) ...
                    + fpc_map_de_YXy .* repmat(albedo_de_jan_map,[1 1 Nyears]) ...
                    + fpc_map_te_YXy * albedo_te ...
                    + fpc_map_bs_YXy .* repmat(albedo_bs_map,[1 1 Nyears]) ...
                 ) + ...
                 snowfrac_jan_YXy .* (...
                      fpc_map_cr_YXy * albedo_cr_snow ...
                    + fpc_map_gr_YXy * albedo_gr_snow ...
                    + fpc_map_ev_YXy * albedo_ev_snow ...
                    + fpc_map_de_YXy * albedo_de_snow ...
                    + fpc_map_te_YXy * albedo_te_snow ...
                    + fpc_map_bs_YXy * albedo_bs_snow ...
                 );
albedo_jul_YXy = (1-snowfrac_jul_YXy) .* (...
                      fpc_map_cr_YXy .* repmat(albedo_cr_jul_map,[1 1 Nyears]) ...
                    + fpc_map_gr_YXy .* repmat(albedo_gr_jul_map,[1 1 Nyears]) ...
                    + fpc_map_ev_YXy .* repmat(albedo_ev_jul_map,[1 1 Nyears]) ...
                    + fpc_map_de_YXy .* repmat(albedo_de_jul_map,[1 1 Nyears]) ...
                    + fpc_map_te_YXy * albedo_te ...
                    + fpc_map_bs_YXy .* repmat(albedo_bs_map,[1 1 Nyears]) ...
                 ) + ...
                 snowfrac_jul_YXy .* (...
                      fpc_map_cr_YXy * albedo_cr_snow ...
                    + fpc_map_gr_YXy * albedo_gr_snow ...
                    + fpc_map_ev_YXy * albedo_ev_snow ...
                    + fpc_map_de_YXy * albedo_de_snow ...
                    + fpc_map_te_YXy * albedo_te_snow ...
                    + fpc_map_bs_YXy * albedo_bs_snow ...
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