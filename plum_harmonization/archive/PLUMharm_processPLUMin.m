function [S, S_2deg] = PLUMharm_processPLUMin(...
    file_in, landArea_YX, LUnames, bareFrac_y0_YX)

Nlu = length(LUnames) ;
landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;
notBare = ~strcmp(LUnames,'BARREN') ;

% Import
S = lpjgu_matlab_readTable_then2map(file_in,'verboseIfNoMat',false,'force_mat_nosave',true) ;

% Move URBAN into BARREN; remove URBAN
S.maps_YXv(:,:,strcmp(S.varNames,'BARREN')) = ...
    sum(S.maps_YXv(:,:,contains(S.varNames,{'BARREN','URBAN'})),3) ;
S.maps_YXv(:,:,strcmp(S.varNames,'URBAN')) = [] ;
S.varNames(strcmp(S.varNames,'URBAN')) = [] ;

% % Move URBAN into NATURAL; remove URBAN
% S.maps_YXv(:,:,strcmp(S.varNames,'NATURAL')) = ...
%     sum(S.maps_YXv(:,:,contains(S.varNames,{'NATURAL','URBAN'})),3) ;
% S.maps_YXv(:,:,strcmp(S.varNames,'URBAN')) = [] ;
% S.varNames(strcmp(S.varNames,'URBAN')) = [] ;

% Ensure same ordering of variables in LUH and PLUM
[~,~,IB] = intersect(LUnames,S.varNames,'stable') ;
S.maps_YXv = S.maps_YXv(:,:,IB) ;
S.varNames = S.varNames(IB) ;

% Match mask with overall mask
S.maps_YXv(landArea_YXv==0) = NaN ;

% Harmonize vegetated+barren fractions
% Will not work if any vegdFrac_y1_YX==0!
if ~isempty(bareFrac_y0_YX)
    vegdFrac_y1_YX = sum(S.maps_YXv(:,:,notBare),3) ;
    if min(vegdFrac_y1_YX(~isnan(vegdFrac_y1_YX))) == 0
        error('"Harmonize vegetated+barren fractions" will not work if any vegdFrac_y1_YX==0!')
    end
    bareFrac_y1_YX = S.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
    bare_diff_YX = bareFrac_y1_YX - bareFrac_y0_YX ;
    
    nonBare_weights_YXv = S.maps_YXv(:,:,notBare) ./ repmat(vegdFrac_y1_YX, [1 1 length(find(notBare))]) ;
    nonBare_weights_YXv(repmat(vegdFrac_y1_YX,[1 1 length(find(notBare))])==0) = 0 ;
    
    S.maps_YXv(:,:,notBare) = S.maps_YXv(:,:,notBare) ...
        + (repmat(bare_diff_YX,[1 1 length(find(notBare))]) ...
        .* nonBare_weights_YXv) ;
    S.maps_YXv(:,:,strcmp(LUnames,'BARREN')) = ...
        S.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ...
        - bare_diff_YX ;
end

% Convert map to km2 (with zeros instead of NaNs)
S.maps_YXv(isnan(S.maps_YXv)) = 0 ;
S.maps_YXv = S.maps_YXv .* landArea_YXv ;

% Get map at 2deg
if isequal(size(S.maps_YXv),[90 180 Nlu])
    % Already at 2-degree
    warning('Already at 2-degree resolution! S.maps_YXv will be all NaN.')
    S_2deg = S ;
    S.maps_YXv = nan(360,720,Nlu) ;
else
    tmp = S.maps_YXv(:,1:4:end,:) ...
        + S.maps_YXv(:,2:4:end,:) ...
        + S.maps_YXv(:,3:4:end,:) ...
        + S.maps_YXv(:,4:4:end,:) ;
    S_2deg.maps_YXv = tmp(1:4:end,:,:) ...
        + tmp(2:4:end,:,:) ...
        + tmp(3:4:end,:,:) ...
        + tmp(4:4:end,:,:) ;
    S_2deg.varNames = S.varNames ;
end


end