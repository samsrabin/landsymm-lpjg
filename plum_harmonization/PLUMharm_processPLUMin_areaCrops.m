function [S, S_nfert, S_irrig, ...
          S_2deg, S_nfert_2deg, S_irrig_2deg, ...
          latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
          max_orig_nfert] = ...
    PLUMharm_processPLUMin_areaCrops(...
        file_in_lcf, landArea_YX, landArea_2deg_YX, LUnames, bareFrac_y0_YX, ...
        latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
        PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method)

PUTURBANHERE = 'BARREN' ;
cf_kgNha_kgNm2 = 1e-4 ;
useLatestPLUMmgmt = ~isempty(latestPLUMin_nfert_2deg_YXv) ;
doInterp = ~isempty(inpaint_method) ;

Nlu = length(LUnames) ;
landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;
notBare = ~strcmp(LUnames,'BARREN') ;
isCrop = ~strcmp(LUnames,'NATURAL') & ~strcmp(LUnames,'PASTURE') & notBare ;

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
    is_actual_PLUMcrops = contains(S_dtl.varNames,'_A') ...
        & ~contains(S_dtl.varNames,'ruminants') ...
        & ~contains(S_dtl.varNames,'monogastrics') ...
        & ~contains(S_dtl.varNames,'pasture') ;
    inds_PLUMcrops_nfert = ...
        contains(S_dtl.varNames,'_FQ') ...
     & ~contains(S_dtl.varNames,'ruminants') ...
     & ~contains(S_dtl.varNames,'monogastrics') ...
     & ~contains(S_dtl.varNames,'pasture')  ;
    inds_PLUMcrops_irrig = ...
        contains(S_dtl.varNames,'_II') ...
     & ~contains(S_dtl.varNames,'ruminants') ...
     & ~contains(S_dtl.varNames,'monogastrics') ...
     & ~contains(S_dtl.varNames,'pasture') ;
    PLUMcrops = S_dtl.varNames(is_actual_PLUMcrops) ;
    PLUMcrops = strrep(PLUMcrops,'_A','') ;
    S_cropa.maps_YXv = S_dtl.maps_YXv(:,:,is_actual_PLUMcrops) .* S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,'CROPLAND')) ;
    S_nfert.maps_YXv = S_dtl.maps_YXv(:,:,inds_PLUMcrops_nfert) ;
    S_irrig.maps_YXv = S_dtl.maps_YXv(:,:,inds_PLUMcrops_irrig) ;
else
    file_in_cropfrac = strrep(file_in_lcf,'LandCoverFract','CropFract') ;
    S_cropf = lpjgu_matlab_readTable_then2map(file_in_cropfrac,'verboseIfNoMat',false,'force_mat_nosave',true) ;
    PLUMcrops = S_cropf.varNames ;
    S_cropa.maps_YXv = S_cropf.maps_YXv .* S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,'CROPLAND')) ;
    file_in_nfert = strrep(file_in_lcf,'LandCoverFract','Fert') ;
    S_nfert = lpjgu_matlab_readTable_then2map(file_in_nfert,'verboseIfNoMat',false,'force_mat_nosave',true) ;
    file_in_irrig = strrep(file_in_lcf,'LandCoverFract','Irrig') ;
    S_irrig = lpjgu_matlab_readTable_then2map(file_in_irrig,'verboseIfNoMat',false,'force_mat_nosave',true) ;
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
if ~isequal(LPJGcrops,PLUMcrops)
    [~,~,IB] = intersect(LPJGcrops,PLUMcrops,'stable') ;
    if length(IB)~=length(LPJGcrops)
        error('length(IB)~=length(LPJGcrops)')
    end
    S_nfert.maps_YXv = S_nfert.maps_YXv(:,:,IB) ;
    S_irrig.maps_YXv = S_irrig.maps_YXv(:,:,IB) ;
end
S_nfert.varNames = LPJGcrops ;
S_irrig.varNames = LPJGcrops ;

% Ensure management NaN where area == 0
S_cropArea.maps_YXv = S.maps_YXv(:,:,isCrop) ;
S_nfert.maps_YXv(S_cropArea.maps_YXv==0) = NaN ;
S_irrig.maps_YXv(S_cropArea.maps_YXv==0) = NaN ;

% Convert kgN/ha to kgN/m2
S_nfert.maps_YXv = S_nfert.maps_YXv * cf_kgNha_kgNm2 ;

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

% Get maps at 2deg
if isequal(size(S.maps_YXv),[90 180 Nlu])
    % Already at 2-degree
    warning('Already at 2-degree resolution! S.maps_YXv will be all NaN.')
    S_2deg = S ;
    S.maps_YXv = nan(360,720,Nlu) ;
    S_nfert_2deg = S_nfert ;
    S_irrig_2deg = S_irrig ;
else
    S_2deg.maps_YXv = PLUMharm_aggregate(S.maps_YXv,0.5,2) ;
    S_2deg.varNames = S.varNames ;
    S_nfert_2deg.maps_YXv = PLUMharm_aggregate_mgmt(S_nfert.maps_YXv,S.maps_YXv(:,:,isCrop),0.5,2) ;
    S_nfert_2deg.varNames = S_nfert.varNames ;
    S_irrig_2deg.maps_YXv = PLUMharm_aggregate_mgmt(S_irrig.maps_YXv,S.maps_YXv(:,:,isCrop),0.5,2) ;
    S_irrig_2deg.varNames = S_irrig.varNames ;
end

% Update "most recent PLUMin management."
% Returns -1 if latest is -1 (i.e., thisCell has never had thisCrop in PLUM
% record)
if useLatestPLUMmgmt
    latestPLUMin_nfert_2deg_YXv = nanmax(latestPLUMin_nfert_2deg_YXv,S_nfert_2deg.maps_YXv) ;
    latestPLUMin_irrig_2deg_YXv = nanmax(latestPLUMin_irrig_2deg_YXv,S_irrig_2deg.maps_YXv) ;
end

% Interpolate management inputs
max_orig_nfert = squeeze(max(max(S_nfert_2deg.maps_YXv,[],1),[],2)) ;
if doInterp
    S_nfert_2deg2.maps_YXv = PLUMharm_interpolateMgmt(...
        S_nfert_2deg.maps_YXv, S_2deg.maps_YXv(:,:,isCrop), landArea_2deg_YX,...
        LPJGcrops, inpaint_method) ;
    S_irrig_2deg2.maps_YXv = PLUMharm_interpolateMgmt(...
        S_irrig_2deg.maps_YXv, S_2deg.maps_YXv(:,:,isCrop), landArea_2deg_YX,...
        LPJGcrops, inpaint_method) ;
elseif useLatestPLUMmgmt
    S_nfert_2deg2 = S_nfert_2deg ;
    S_irrig_2deg2 = S_irrig_2deg ;
end

% In cells that previously had thisCrop but now do not, replace
% interpolated value with latestPLUMin value
if useLatestPLUMmgmt
    repl_YXv = S_2deg.maps_YXv(:,:,isCrop)==0 & latestPLUMin_nfert_2deg_YXv>=0 ;
    S_nfert_2deg2.maps_YXv(repl_YXv) = latestPLUMin_nfert_2deg_YXv(repl_YXv) ;
    repl_YXv = S_2deg.maps_YXv(:,:,isCrop)==0 & latestPLUMin_irrig_2deg_YXv>=0 ;
    S_irrig_2deg2.maps_YXv(repl_YXv) = latestPLUMin_irrig_2deg_YXv(repl_YXv) ;
end

% Save managmeent
if doInterp || useLatestPLUMmgmt
    S_nfert_2deg.maps_YXv = S_nfert_2deg2.maps_YXv ;
    S_irrig_2deg.maps_YXv = S_irrig_2deg2.maps_YXv ;
end

% Make sure you have no NaNs in half-degree mgmt arrays
S_nfert.maps_YXv(S.maps_YXv(:,:,isCrop)==0) = 0 ;
S_irrig.maps_YXv(S.maps_YXv(:,:,isCrop)==0) = 0 ;
if any(isnan(S_nfert.maps_YXv(:)))
    error('NaN in S_nfert.maps_YXv!')
elseif any(isnan(S_irrig.maps_YXv(:)))
    error('NaN in S_irrig.maps_YXv!')
end


end