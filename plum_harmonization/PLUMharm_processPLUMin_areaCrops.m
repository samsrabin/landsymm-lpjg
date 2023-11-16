function [S, S_nfert, S_irrig, ...
          S_2deg, S_nfert_2deg, S_irrig_2deg, ...
          latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
          max_orig_nfert] = ...
    PLUMharm_processPLUMin_areaCrops(...
        file_in_lcf, landArea_YX, landArea_2deg_YX, LUnames, bareFrac_y0_YX, ...
        latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
        PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method, ...
        fruitveg_sugar_2oil, varargin)

outPrec = 6 ;
if ~isempty(varargin)
    outPrec = varargin{1} ;
    if length(varargin) > 1
        error('PLUMharm_processPLUMin_areaCrops() takes at most one optional argument (outPrec).')
    end
end

PUTURBANHERE = 'BARREN' ;
cf_kgNha_kgNm2 = 1e-4 ;
useLatestPLUMmgmt = ~isempty(latestPLUMin_nfert_2deg_YXv) ;
doInterp = ~isempty(inpaint_method) ;
combineCrops = isempty(PLUMtoLPJG) ;
if combineCrops
    S_nfert = [] ;
    S_irrig = [] ;
    S_nfert_2deg = [] ;
    S_irrig_2deg = [] ;
    latestPLUMin_nfert_2deg_YXv = [] ;
    latestPLUMin_irrig_2deg_YXv = [] ;
    max_orig_nfert = [] ;
end

Nlu = length(LUnames) ;
landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;
notBare = ~strcmp(LUnames,'BARREN') ;
isCrop = ~strcmp(LUnames,'NATURAL') & ~strcmp(LUnames,'PASTURE') & notBare ;

% Import LandCoverFract.txt
S_lcf = lpjgu_matlab_readTable_then2map(file_in_lcf,'verboseIfNoMat',false,'force_mat_save',true) ;

% Move URBAN into PUTURBANHERE; remove URBAN
if any(strcmp(S_lcf.varNames,'URBAN'))
    S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,PUTURBANHERE)) = ...
        sum(S_lcf.maps_YXv(:,:,contains(S_lcf.varNames,{PUTURBANHERE,'URBAN'})),3) ;
    S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,'URBAN')) = [] ;
    S_lcf.varNames(strcmp(S_lcf.varNames,'URBAN')) = [] ;
end

% Import detailed LandUse.txt
file_in_dtl = strrep(file_in_lcf,'LandCoverFract','LandUse') ;
if exist(file_in_dtl,'file') || exist([file_in_dtl '.gz'],'file')
    S_dtl = lpjgu_matlab_readTable_then2map(file_in_dtl,'verboseIfNoMat',false,'force_mat_save',true) ;
%     PLUMcrops = S_dtl.varNames ;
    is_actual_PLUMcrops = contains(S_dtl.varNames,'_A') ...
        & ~contains(S_dtl.varNames,'ruminants') ...
        & ~contains(S_dtl.varNames,'monogastrics') ...
        & ~contains(S_dtl.varNames,'pasture') ;
    if ~combineCrops
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
    end
    PLUMcrops = S_dtl.varNames(is_actual_PLUMcrops) ;
    PLUMcrops = strrep(PLUMcrops,'_A','') ;
    S_cropa.maps_YXv = S_dtl.maps_YXv(:,:,is_actual_PLUMcrops) .* S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,'CROPLAND')) ;
    if ~combineCrops
        S_nfert.maps_YXv = S_dtl.maps_YXv(:,:,inds_PLUMcrops_nfert) ;
        S_irrig.maps_YXv = S_dtl.maps_YXv(:,:,inds_PLUMcrops_irrig) ;
        % If there is irrigation or fertilization present on SetAside, remove
        % it!
        if ~any(strcmp(PLUMcrops,'setaside'))
            error('\nsetaside not present %s', file_in_dtl)
        end
        tmpN = S_nfert.maps_YXv(:,:,strcmp(PLUMcrops,'setaside')) ;
        tmpI = S_irrig.maps_YXv(:,:,strcmp(PLUMcrops,'setaside')) ;
        if any(any(tmpN>0)) || any(any(tmpI>0))
            tmpN(tmpN>0) = 0 ;
            tmpI(tmpI>0) = 0 ;
            S_nfert.maps_YXv(:,:,strcmp(PLUMcrops,'setaside')) = tmpN ;
            S_irrig.maps_YXv(:,:,strcmp(PLUMcrops,'setaside')) = tmpI ;
        end
        clear tmpN tmpI
    end
else
    file_in_cropfrac = strrep(file_in_lcf,'LandCoverFract','CropFract') ;
    S_cropf = lpjgu_matlab_readTable_then2map(file_in_cropfrac,'verboseIfNoMat',false,'force_mat_save',true) ;
    PLUMcrops = S_cropf.varNames ;
    S_cropa.maps_YXv = S_cropf.maps_YXv .* S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,'CROPLAND')) ;
    clear S_cropf
    if ~combineCrops
        file_in_nfert = strrep(file_in_lcf,'LandCoverFract','Fert') ;
        S_nfert = lpjgu_matlab_readTable_then2map(file_in_nfert,'verboseIfNoMat',false,'force_mat_save',true) ;
        file_in_irrig = strrep(file_in_lcf,'LandCoverFract','Irrig') ;
        S_irrig = lpjgu_matlab_readTable_then2map(file_in_irrig,'verboseIfNoMat',false,'force_mat_save',true) ;
        if ~any(strcmp(PLUMcrops,'setaside'))
            error('\nsetaside not present, probably because you''re not using LandUse.txt (which wasn''t found): %s', file_in_cropfrac)
        end
    end
end

if combineCrops
    S_cropa.maps_YXv = sum(S_cropa.maps_YXv, 3) ;
    PLUMcrops = LPJGcrops ;
else
    
    % Move fruitveg and/or sugar into their proxy crops, if needed
    if fruitveg_sugar_2oil
        fakelist_plum = {'fruitveg', 'sugar'} ;
        fakelist_lpjg = {'Oilcrops', 'Oilcrops'} ;
        for c = 1:length(fakelist_lpjg)
            thisCrop_plum = fakelist_plum{c} ;
            if ~any(strcmp(PLUMcrops, thisCrop_plum))
                continue
            end
            thisCrop_lpjg = fakelist_lpjg{c} ;
            fprintf('Moving %s into %s...\n', thisCrop_plum, thisCrop_lpjg)
            % Do area first so you can get weighting
            data_real_YX = S_cropa.maps_YXv(:,:,strcmp(PLUMcrops, PLUMtoLPJG{strcmp(LPJGcrops, thisCrop_lpjg)})) ;
            data_fake_YX = S_cropa.maps_YXv(:,:,strcmp(PLUMcrops, thisCrop_plum)) ;
            wt_real_YX = data_real_YX ./ (data_real_YX + data_fake_YX) ;
            wt_fake_YX = data_fake_YX ./ (data_real_YX + data_fake_YX) ;
            S_cropa.maps_YXv(:,:,strcmp(PLUMcrops, PLUMtoLPJG{strcmp(LPJGcrops, thisCrop_lpjg)})) ...
                = data_real_YX + data_fake_YX ;
            S_cropa.maps_YXv(:,:,strcmp(PLUMcrops, thisCrop_plum)) = [] ;
            % Nfert
            data_real_YX = S_nfert.maps_YXv(:,:,strcmp(PLUMcrops, PLUMtoLPJG{strcmp(LPJGcrops, thisCrop_lpjg)})) ;
            data_fake_YX = S_nfert.maps_YXv(:,:,strcmp(PLUMcrops, thisCrop_plum)) ;
            S_nfert.maps_YXv(:,:,strcmp(PLUMcrops, PLUMtoLPJG{strcmp(LPJGcrops, thisCrop_lpjg)})) ...
                = data_real_YX.*wt_real_YX + data_fake_YX.*wt_fake_YX ;
            S_nfert.maps_YXv(:,:,strcmp(PLUMcrops, thisCrop_plum)) = [] ;
            % Irrigation
            data_real_YX = S_irrig.maps_YXv(:,:,strcmp(PLUMcrops, PLUMtoLPJG{strcmp(LPJGcrops, thisCrop_lpjg)})) ;
            data_fake_YX = S_irrig.maps_YXv(:,:,strcmp(PLUMcrops, thisCrop_plum)) ;
            S_irrig.maps_YXv(:,:,strcmp(PLUMcrops, PLUMtoLPJG{strcmp(LPJGcrops, thisCrop_lpjg)})) ...
                = data_real_YX.*wt_real_YX + data_fake_YX.*wt_fake_YX ;
            S_irrig.maps_YXv(:,:,strcmp(PLUMcrops, thisCrop_plum)) = [] ;
            clear data_real_YX data_fake_YX
            % Remove fake from PLUMcrops
            PLUMcrops(strcmp(PLUMcrops, thisCrop_plum)) = [] ;
        end
    end
    
    % If there is a mismatch in crop lists, throw an error
    if ~isequal(PLUMcrops, intersect(PLUMcrops, PLUMtoLPJG, 'stable'))
        error('Mismatch between crops in PLUMcrops and PLUMtoLPJG')
    end
    
    % If there is irrigation or fertilization present on SetAside, throw an
    % error!
    tmpN = S_nfert.maps_YXv(:,:,strcmp(PLUMcrops,'setaside')) ;
    tmpI = S_irrig.maps_YXv(:,:,strcmp(PLUMcrops,'setaside')) ;
    isBad = any(any(tmpN>0)) || any(any(tmpI>0)) ;
    if isBad
        error('Fertilization and/or irrigation on SetAside!')
    end
    clear tmpN tmpI
    
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
end

% Concatenate crops and other LUs
S.varNames = [PLUMcrops S_lcf.varNames(~strcmp(S_lcf.varNames,'CROPLAND'))] ;
S.maps_YXv = cat(3, S_cropa.maps_YXv, S_lcf.maps_YXv(:,:,~strcmp(S_lcf.varNames,'CROPLAND'))) ;

% Ensure same ordering of variables in baseline LU and PLUM
if ~isequal(LUnames,S.varNames)
    missing_types = setdiff(LUnames,S.varNames) ;
    if ~isempty(missing_types)
        if fruitveg_sugar_2oil && isequal(shiftdim(missing_types), shiftdim({'FruitAndVeg', 'Sugar'}))
            error('FruitAndVeg and Sugar are missing from PLUM types; do you need to set fruitveg_sugar_2oil to false?')
        else
            error('Types in baseline LU but missing from PLUM: %s', sprintf('%s ', missing_types))
        end
    end
    [~,~,IB] = intersect(LUnames,S.varNames,'stable') ;
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
    if ~combineCrops
        S_nfert.maps_YXv = S_nfert.maps_YXv(:,:,IB) ;
        S_irrig.maps_YXv = S_irrig.maps_YXv(:,:,IB) ;
    end
end
if ~combineCrops
    S_nfert.varNames = LPJGcrops ;
    S_irrig.varNames = LPJGcrops ;
end

if ~combineCrops
    % Ensure management NaN where area == 0
    S_cropArea.maps_YXv = S.maps_YXv(:,:,isCrop) ;
    S_nfert.maps_YXv(S_cropArea.maps_YXv==0) = NaN ;
    S_irrig.maps_YXv(S_cropArea.maps_YXv==0) = NaN ;
    
    % Convert kgN/ha to kgN/m2
    S_nfert.maps_YXv = S_nfert.maps_YXv * cf_kgNha_kgNm2 ;
end

% Match mask with overall mask
S.maps_YXv(landArea_YXv==0) = NaN ;

% Harmonize vegetated+barren fractions
if ~isempty(bareFrac_y0_YX)
    bareFrac_y1_YX = sum(S.maps_YXv(:,:,~notBare), 3) ;
    if ~isequaln(bareFrac_y0_YX, bareFrac_y1_YX)
        diff_YX = (bareFrac_y1_YX - bareFrac_y0_YX) .* landArea_YX*1e-6 ;
        total_bare_y0 = nansum(nansum(bareFrac_y0_YX .* landArea_YX*1e-6)) ;
        netDiff = nansum(nansum(diff_YX)) ;
        grossDiff = nansum(nansum(abs(diff_YX))) ;
        netDiff_pct = 100 * netDiff / total_bare_y0 ;
        grossDiff_pct = 100 * grossDiff / total_bare_y0 ;
        disp('    Harmonizing vegd+bare fractions')
        fprintf('        Sum bare  diff  = %g km2 (%0.2f%% of y0)\n', netDiff, netDiff_pct)
        fprintf('        Sum bare |diff| = %g km2 (%0.2f%% of y0)\n', grossDiff, grossDiff_pct)
        vegdFrac_y1_YX = sum(S.maps_YXv(:,:,notBare),3) ;
        vegdFrac_y1_YXrep = repmat(vegdFrac_y1_YX, [1 1 sum(notBare)]) ;
        vegdFrac_y1_YXv = S.maps_YXv(:,:,notBare) ./ vegdFrac_y1_YXrep ;
        vegdFrac_y1_YXv(vegdFrac_y1_YXrep==0) = 0 ;
        S.maps_YXv(:,:,notBare) = vegdFrac_y1_YXv .* (1 - bareFrac_y0_YX) ;
        S.maps_YXv(:,:,~notBare) = bareFrac_y0_YX ;
        tol = 1e-12 ;
        maxdiff = nanmax(nanmax(abs(sum(S.maps_YXv, 3) - 1))) ;
        if maxdiff > tol
            error('Land use fractions don''t sum to 1 within tolerance %g; max abs. difference %g', ...
                tol, maxdiff)
        end
    end
end

% Convert map to km2 (with zeros instead of NaNs)
S.maps_YXv(isnan(S.maps_YXv)) = 0 ;
S.maps_YXv = S.maps_YXv .* landArea_YXv ;

% Check for bad values
if ~combineCrops
    tmp_nfert_YXv = S_nfert.maps_YXv ;
    tmp_nfert_YXv(isnan(tmp_nfert_YXv)) = 0 ;
    tmp_irrig_YXv = S_irrig.maps_YXv ;
    tmp_irrig_YXv(isnan(tmp_irrig_YXv)) = 0 ;
    PLUMharm_checkBadVals([], tmp_nfert_YXv, tmp_irrig_YXv, ...
        [], LUnames, 'import_mgmt', outPrec) ;
    clear tmp_*
end

% Get maps at 2deg
if isequal(size(S.maps_YXv),[90 180 Nlu])
    % Already at 2-degree
    warning('Already at 2-degree resolution! S.maps_YXv will be all NaN.')
    S_2deg = S ;
    S.maps_YXv = nan(360,720,Nlu) ;
    if ~combineCrops
        S_nfert_2deg = S_nfert ;
        S_irrig_2deg = S_irrig ;
    end
else
    S_2deg.maps_YXv = PLUMharm_aggregate(S.maps_YXv,0.5,2) ;
    S_2deg.varNames = S.varNames ;
    if ~combineCrops
        S_nfert_2deg.maps_YXv = PLUMharm_aggregate_mgmt(S_nfert.maps_YXv,S.maps_YXv(:,:,isCrop),0.5,2) ;
        S_nfert_2deg.varNames = S_nfert.varNames ;
        S_irrig_2deg.maps_YXv = PLUMharm_aggregate_mgmt(S_irrig.maps_YXv,S.maps_YXv(:,:,isCrop),0.5,2) ;
        S_irrig_2deg.varNames = S_irrig.varNames ;
    end
end

if ~combineCrops
    % Update "most recent PLUMin management."
    % Returns -1 if latest is -1 (i.e., thisCell has never had thisCrop in PLUM
    % record)
    if useLatestPLUMmgmt
        latestPLUMin_nfert_2deg_YXv = nanmax(latestPLUMin_nfert_2deg_YXv,S_nfert_2deg.maps_YXv) ;
        latestPLUMin_irrig_2deg_YXv = nanmax(latestPLUMin_irrig_2deg_YXv,S_irrig_2deg.maps_YXv) ;
    end
    
    % Check for bad values
    tmp_nfert_YXv = S_nfert_2deg.maps_YXv ;
    tmp_nfert_YXv(isnan(tmp_nfert_YXv)) = 0 ;
    tmp_irrig_YXv = S_irrig_2deg.maps_YXv ;
    tmp_irrig_YXv(isnan(tmp_irrig_YXv)) = 0 ;
    PLUMharm_checkBadVals([], tmp_nfert_YXv, tmp_irrig_YXv, [], ...
        LUnames, 'import_mgmt', outPrec) ;
    clear tmp_*
    
    % Interpolate management inputs
    max_orig_nfert = squeeze(max(max(S_nfert_2deg.maps_YXv,[],1),[],2)) ;
    if doInterp
        S_nfert_2deg2.maps_YXv = PLUMharm_interpolateMgmt(...
            S_nfert_2deg.maps_YXv, S_2deg.maps_YXv(:,:,isCrop), landArea_2deg_YX,...
            LPJGcrops, inpaint_method) ;
        S_irrig_2deg2.maps_YXv = PLUMharm_interpolateMgmt(...
            S_irrig_2deg.maps_YXv, S_2deg.maps_YXv(:,:,isCrop), landArea_2deg_YX,...
            LPJGcrops, inpaint_method) ;
        % Irrigation needs to be maximum 1.
        S_irrig_2deg2.maps_YXv(S_irrig_2deg2.maps_YXv>1) = 1 ;
        % Check for bad values
        PLUMharm_checkBadVals([], S_nfert_2deg2.maps_YXv, S_irrig_2deg2.maps_YXv, ...
            [], LUnames, 'import_mgmt_2deg2', outPrec) ;
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
        % Check for bad values
        PLUMharm_checkBadVals([], S_nfert_2deg2.maps_YXv, S_irrig_2deg2.maps_YXv, ...
            [], LUnames, 'import_mgmt_2deg2.2', outPrec) ;
    end
    
    % Save management
    if doInterp || useLatestPLUMmgmt
        S_nfert_2deg.maps_YXv = S_nfert_2deg2.maps_YXv ;
        S_irrig_2deg.maps_YXv = S_irrig_2deg2.maps_YXv ;
    end
end


end
