function [S, S_2deg] = PLUMharm_processPLUMin_areaCrops(...
    file_in_lcf, landArea_YX, LUnames, bareFrac_y0_YX, ...
    PLUMtoLPJG, LPJGcrops,norm2extra)

PUTURBANHERE = 'BARREN' ;

Nlu = length(LUnames) ;
landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;
notBare = ~strcmp(LUnames,'BARREN') ;

% Import LandCoverFract.txt
S_lcf = lpjgu_matlab_readTable_then2map(file_in_lcf,'verboseIfNoMat',false,'force_mat_nosave',true) ;

% Move URBAN into PUTURBANHERE; remove URBAN
if any(strcmp(S_lcf.varNames,'URBAN'))
    S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,PUTURBANHERE)) = ...
        sum(S_lcf.maps_YXv(:,:,contains(S_lcf.varNames,{PUTURBANHERE,'URBAN'})),3) ;
    S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,'URBAN')) = [] ;
    S_lcf.varNames(strcmp(S_lcf.varNames,'URBAN')) = [] ;
end

% Import detailed LandUse.txt
file_in_dtl = strrep(file_in_lcf,'LandCoverFract','LandUse') ;
if exist(file_in_dtl,'file')
    S_dtl = lpjgu_matlab_readTable_then2map(file_in_dtl,'verboseIfNoMat',false,'force_mat_nosave',true) ;
    PLUMcrops = S_dtl.varNames ;
    is_actual_PLUMcrops = contains(PLUMcrops,'_A') ...
        & ~contains(PLUMcrops,'ruminants') ...
        & ~contains(PLUMcrops,'monogastrics') ...
        & ~contains(PLUMcrops,'pasture') ;
    PLUMcrops = PLUMcrops(is_actual_PLUMcrops) ;
    PLUMcrops = strrep(PLUMcrops,'_A','') ;
    S_cropa.maps_YXv = S_dtl.maps_YXv(:,:,is_actual_PLUMcrops) .* S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,'CROPLAND')) ;
else
    file_in_cropfrac = strrep(file_in_lcf,'LandCoverFract','CropFract') ;
    S_cropf = lpjgu_matlab_readTable_then2map(file_in_cropfrac,'verboseIfNoMat',false,'force_mat_nosave',true) ;
    PLUMcrops = S_cropf.varNames ;
    S_cropa.maps_YXv = S_cropf.maps_YXv .* S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,'CROPLAND')) ;
end
% Translate PLUM crop names to LPJ-GUESS crop names, if necessary
[~,IA] = intersect(LPJGcrops,PLUMcrops,'stable') ;
if ~isequal(IA,shiftdim(1:length(LPJGcrops)))
    for c = 1:length(PLUMcrops)
        thisCrop = PLUMcrops{c} ;
        PLUMcrops{c} = LPJGcrops{find(strcmp(PLUMtoLPJG,thisCrop))} ;
    end
end

% Move some % of normal crops into ExtraCrop
% (Set norm2extra = 0 if reading already-harmonized files!)
S_cropa.maps_YXv(:,:,strcmp(PLUMcrops,'ExtraCrop')) = S_cropa.maps_YXv(:,:,strcmp(PLUMcrops,'ExtraCrop')) ...
    + norm2extra*sum(S_cropa.maps_YXv(:,:,~strcmp(PLUMcrops,'ExtraCrop')),3) ;
S_cropa.maps_YXv(:,:,~strcmp(PLUMcrops,'ExtraCrop')) = S_cropa.maps_YXv(:,:,~strcmp(PLUMcrops,'ExtraCrop')) ...
    - norm2extra*S_cropa.maps_YXv(:,:,~strcmp(PLUMcrops,'ExtraCrop')) ;

% Concatenate crops and other LUs
S.varNames = [PLUMcrops S_lcf.varNames(~strcmp(S_lcf.varNames,'CROPLAND'))] ;
S.maps_YXv = cat(3, S_cropa.maps_YXv, S_lcf.maps_YXv(:,:,~strcmp(S_lcf.varNames,'CROPLAND'))) ;

% Ensure same ordering of variables in LUH and PLUM
if ~isequal(LUnames,S.varNames)
    [~,~,IB] = intersect(LUnames,S.varNames,'stable') ;
%     [C,IA,IB] = intersect(LUnames,S.varNames,'stable')
    if length(IB)~=length(LUnames)
        error('length(IB)~=length(LUnames)')
    end
    S.maps_YXv = S.maps_YXv(:,:,IB) ;
    S.varNames = S.varNames(IB) ;
end

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