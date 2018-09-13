%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LUH1-style harmonization for PLUM outputs, at cropType level %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_year = 2010 ;
year1 = 2046 ;
yearN = 2100 ;

% Directory for PLUM outputs
PLUM_in_toptop = {...
                  '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1' ;
%                   '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP3.v10.s1' ;
%                   '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP4.v10.s1' ;
%                   '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v10.s1' ;
                  } ;
              
% Save?
save_halfDeg_mat = true ;
save_2deg_mat = false ;
save_halfDeg_txt = false ;
save_2deg_txt = false ;

% Coordinates of 2-degree cell to debug (leave empty for no debug)
debugIJ_2deg = [] ;
% debugIJ_2deg = [29 105] ; % i=6
% debugIJ_2deg = [6 1] ;
% debugIJ_2deg = [73 98] ;
% debugIJ_2deg = [48 90] ;

% Print verbose messages?
verbose = false ;


%% Setup

addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;

% Method for inpaint_nans()
inpaint_method = 0 ;

% Use "latest PLUM management" in cells that don't have thisCrop thisYear
% but did in a previous year? FALSE = rely solely on interpolation of
% thisYear's thisMgmt.
useLatestPLUMmgmt = false ;

% Save details
save_halfDeg_any = save_halfDeg_mat || save_halfDeg_txt ;
save_2deg_any = save_2deg_mat || save_2deg_txt ;
save_any = save_halfDeg_any || save_2deg_any ;

% Output file details
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;

% How much divergence from PLUM-original transition is acceptable? (%)
conserv_tol_pct = 0.2 ;
% How much divergence from PLUM-original transition is acceptable? (m2)
% Also used for checking whether agricultural area exceeds vegetated area
% in each gridcell.
conserv_tol_area = 1e3 ;

% The fraction of real crops (by area) not included in PLUM's 6
% non-energyCrops commodities. Multiply PLUM area by (1-norm2extra) to
% get area of the crop itself without the norm2extra fraction included.
norm2extra = 0.177 ;

% Debug?
do_debug = ~isempty(debugIJ_2deg) ;
clear db2i db2j
if do_debug
    db2i = debugIJ_2deg(1) ;
    db2j = debugIJ_2deg(2) ;
end
debug_header = ':\t\tArea\t\tFrac\n' ;


%% Import reference data

doHarm = true ;
run('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/plum_harmonization/PLUMharm_importRefData.m') ;


%% Do it

% The years we want to produce PLUM outputs for (will begin with
% transitions from years(1)-1 to years(1)
years = year1:yearN ;
if year1 <= base_year
    error('year1 (%d) must be > base_year (%d)!\n', year1, base_year)
elseif year1 > yearN
    error('year1 (%d) must be <= yearN (%d)!\n', year1, yearN)
end

for d = 1:length(PLUM_in_toptop)

    %%%%%%%%%%%%%
    %%% Setup %%%
    %%%%%%%%%%%%%

    % Get directories
    PLUM_in_top = removeslashifneeded(PLUM_in_toptop{d}) ;
    disp(PLUM_in_top)
    PLUM_out_top = addslashifneeded([PLUM_in_top '.harm/']) ;
    unix(['mkdir -p ' PLUM_out_top]) ;
    PLUM_in_top = addslashifneeded(PLUM_in_top) ;

    % Read fraction of VEGETATED protected by...
    %%% PLUM's minimum natural fraction ("rate")
    [~,r] = unix(['grep "MIN_NATURAL_RATE" ' PLUM_in_top 'config.properties | sed "s/MIN_NATURAL_RATE=//"  | tr -d ''\n''']) ;
    min_natural_rate = str2double(r) ;
    resFrac_minN_YX = min_natural_rate * ones(size(resFrac_prot_YX)) ;
    %%% Total
    resFrac_YX = 1 - (1-resFrac_terr_YX).*(1-resFrac_prot_YX)*(1-min_natural_rate) ;
    % Convert to area
    resArea_YX = resFrac_YX .* base_vegd_YX ;
    resArea_2deg_YX = PLUMharm_aggregate(resArea_YX,0.5,2) ;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Process PLUM outputs %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Debugging output
    if do_debug
        fprintf('Center Lon %.2f, Lat %.2f, land area %0.4g\n\n',1+lons_map_2deg(db2i,db2j),1+lats_map_2deg(db2i,db2j),round(landArea_2deg_YX(db2i,db2j),4)) ;
        PLUMharm_debugOut('LUH2','areas',base_2deg.maps_YXv,landArea_2deg_YX,debugIJ_2deg,LUnames)
        PLUMharm_debugOut('LUH2','Nfert',base_nfertTot_2deg.maps_YXv,base_2deg.maps_YXv(:,:,isCrop),debugIJ_2deg,LPJGcrops)
    end

    bareFrac_y0_YX = base_bareFrac_YX ;
    if useLatestPLUMmgmt
        latestPLUMin_nfert_2deg_YXv = -1*ones([size(landArea_2deg_YX) Ncrops_lpjg]) ;
        latestPLUMin_irrig_2deg_YXv = -1*ones([size(landArea_2deg_YX) Ncrops_lpjg]) ;
    else
        latestPLUMin_nfert_2deg_YXv = [] ;
        latestPLUMin_irrig_2deg_YXv = [] ;
    end
    clear out_* in_*
    for y = 1:length(years)
        
        thisYear = years(y) ;
        disp(num2str(thisYear)) ;
        tic ;
                
        % Import previous harmonized year, if needed
        if thisYear-1 == base_year
            out_y0 = base ;
%             disp('out_y0 from luh2')
            out_y0_2deg = base_2deg ;
            out_y0_agri_YXv = base.maps_YXv(:,:,isAgri) ;
            out_y0_2deg_agri_YXv = base_2deg.maps_YXv(:,:,isAgri) ;
            out_y0_nfert_YXv = base_nfert.maps_YXv ;
            out_y0_nfert_2deg = base_nfert_2deg ;
            out_y0_irrig_YXv = base_irrig.maps_YXv ;
            out_y0_irrig_2deg = base_irrig_2deg ;
        elseif y==1
            file_in = [removeslashifneeded(PLUM_in_top) '.harm/' num2str(thisYear-1) 'post.base' num2str(base_year) '.mat'] ;
%             disp(['*y0* from ' file_in])
            load(file_in) ;
        end
        
        % Import for debugging redistribution
        if do_debug
            if y==1
                file_in = [PLUM_base_in '/LandCoverFract.txt'] ;
                [~, ~, ~, in_y0orig_2deg, in_y0orig_nfert_2deg, in_y0orig_irrig_2deg, ...
                    ~, ~] = ...
                    PLUMharm_processPLUMin_areaCrops(file_in, landArea_YX, landArea_2deg_YX, ...
                    LUnames, bareFrac_y0_YX, ...
                    latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
                    PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method) ;
            end
        else
            diffO_YXv = [] ;
            in_y0orig_2deg = [] ;
            in_y1orig_2deg = [] ;
        end

        % Import and process previous year, if needed
        if ~exist('in_y0','var')
            file_in = [PLUM_in_top num2str(thisYear-1) '/LandCoverFract.txt'] ;
%             disp(['in_y0 from ' file_in])
            [in_y0, ~, ~, in_y0_2deg, in_y0_nfert_2deg, in_y0_irrig_2deg, ...
                latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
                max_orig_nfert_y0] = ...
                PLUMharm_processPLUMin_areaCrops(file_in, landArea_YX, landArea_2deg_YX, ...
                LUnames, bareFrac_y0_YX, latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
                PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method) ;
            bareFrac_y0_YX = in_y0.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ./ landArea_YX ;
            in_y0_agri_YXv = in_y0.maps_YXv(:,:,isAgri) ;
            in_y0_2deg_agri_YXv = in_y0_2deg.maps_YXv(:,:,isAgri) ;
            in_y0_2deg_vegd_YX = sum(in_y0_2deg.maps_YXv(:,:,notBare),3) ;
        end

        % Import this year and convert to area
        file_in = [PLUM_in_top num2str(thisYear) '/LandCoverFract.txt'] ;
%         disp(['in_y1 from ' file_in])
        [in_y1, ~, ~, in_y1_2deg, in_y1_nfert_2deg, in_y1_irrig_2deg, ...
            latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
            max_orig_nfert_y1] = ...
            PLUMharm_processPLUMin_areaCrops(file_in, landArea_YX, landArea_2deg_YX, ...
            LUnames, bareFrac_y0_YX, latestPLUMin_nfert_2deg_YXv, latestPLUMin_irrig_2deg_YXv, ...
            PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method) ;
        in_y1_agri_YXv = in_y1.maps_YXv(:,:,isAgri) ;
        in_y1_2deg_agri_YXv = in_y1_2deg.maps_YXv(:,:,isAgri) ;
        in_y1_ntrl_YX = in_y1.maps_YXv(:,:,strcmp(in_y1.varNames,'NATURAL')) ;
        in_y1_bare_YX = in_y1.maps_YXv(:,:,strcmp(in_y1.varNames,'BARREN')) ;
        in_y1_vegd_YX = sum(in_y1.maps_YXv(:,:,notBare),3) ;
        in_y1_vegdFrac_YX = in_y1_vegd_YX ./ landArea_YX ;
        in_y1_bareFrac_YX = in_y1_bare_YX ./ landArea_YX ;
        in_y1_2deg_vegd_YX = sum(in_y1_2deg.maps_YXv(:,:,notBare),3) ;
        
        % Debugging
        if verbose
            disp(' ')
            disp('Initial import')
            tmp0 = 1e-6*1e-6*sum(sum(sum(in_y0_2deg.maps_YXv(:,:,isAgri)))) ;
            tmp1 = 1e-6*1e-6*sum(sum(sum(in_y1_2deg.maps_YXv(:,:,isAgri)))) ;
            fprintf('in_%d glob agri area (million km2): %0.1f\n', thisYear-1, tmp0)
            fprintf('in_%d glob agri area (million km2): %0.1f\n', thisYear,   tmp1)
            fprintf('                  diff (million km2): %0.1f\n\n', tmp1-tmp0)
            tmp0 = 1e-6*1e-6*sum(in_y0_2deg.maps_YXv(:)) ;
            tmp1 = 1e-6*1e-6*sum(in_y1_2deg.maps_YXv(:)) ;
            fprintf('in_%d glob land area (million km2): %0.1f\n', thisYear-1, tmp0)
            fprintf('in_%d glob land area (million km2): %0.1f\n', thisYear,   tmp1)
            fprintf('                  diff (million km2): %0.1f\n\n', tmp1-tmp0)
            tmp = 1e-6*1e-6*sum(base_vegd_2deg_YX(:)) ;
            fprintf('base_2deg glob vegd area (million km2): %0.1f\n', tmp)
            clear tmp*
        end
        
%         % Debugging
%         if verbose
%             disp(' ')
%             disp('Initial import')
%             tmp0 = 1e-6*1e-3*sum(sum(in_y0_nfert_2deg.maps_YXv(:,:,2) .* in_y0_2deg.maps_YXv(:,:,2))) ;
%             tmp1 = 1e-6*1e-3*sum(sum(in_y1_nfert_2deg.maps_YXv(:,:,2) .* in_y1_2deg.maps_YXv(:,:,2))) ;
%             fprintf('in_%d glob Nfert (Mt): %0.3f\n', thisYear-1, tmp0)
%             fprintf('in_%d glob Nfert (Mt): %0.3f\n', thisYear,   tmp1)
%             fprintf('              diff (Mt): %0.3f\n\n', tmp1-tmp0)
%         end
        
        % Import for debugging redistribution
        if do_debug
            in_y1orig_2deg = in_y1_2deg ;
            diffO_YXv = in_y1orig_2deg.maps_YXv - in_y0orig_2deg.maps_YXv ;
        end
        
        % Make sure that total global loss of a land use does not exceed
        % the total global area of that land use.
        % NOTE: v10.SSP4 has Miscanthus in 2010, and some gridcells lose
        % Miscanthus from 2010 to 2011. This cannot be satisfied using my
        % 2010 maps because I have no Miscanthus. However, it shouldn't be
        % a problem as long as the grid cells PLUM specifies as losing
        % Miscanthus lose ALL of their Miscanthus.
        for i = 1:length(LUnames)
            thisLU = LUnames{i} ;
            out0_thisLU_YX = out_y0.maps_YXv(:,:,i) ;
            in0_thisLU_YX = in_y0.maps_YXv(:,:,i) ;
            in1_thisLU_YX = in_y1.maps_YXv(:,:,i) ;
            inDiff_thisLU_YX = in1_thisLU_YX - in0_thisLU_YX ;
            globLoss_thisLU = sum(inDiff_thisLU_YX(inDiff_thisLU_YX<0)) ;
            globArea_thisLU = sum(out0_thisLU_YX(:)) ;
            NETglobChg_thisLU = sum(inDiff_thisLU_YX(:)) ;
            if -globLoss_thisLU > globArea_thisLU
                if do_debug
                    disp(LUnames{i})
                    fprintf('globArea_y0_thisLu = %0.4e\n',sum(in0_thisLU_YX(:))) ;
                    fprintf('globArea_y1_thisLu = %0.4e\n',sum(in1_thisLU_YX(:))) ;
                    fprintf('NETglobChg_thisLU  = %0.4e\n',NETglobChg_thisLU) ;
                    fprintf('globLoss_thisLU    = %0.4e\n',globLoss_thisLU) ;
                    fprintf('globArea_thisLU    = %0.4e\n',globArea_thisLU) ;
                end
                if strcmp(thisLU,'Miscanthus') ...
                && thisYear-1==base_year
                    if sum(in1_thisLU_YX(in0_thisLU_YX>0)) == 0
                        warning('PLUM has Miscanthus in 2010, and some gridcells lose Miscanthus from 2010 to 2011. This cannot be satisfied using my 2010 maps because I have no Miscanthus. However, it shouldn''t be a problem because the grid cells PLUM specifies as losing Miscanthus lose ALL it.')
                    else
                        warning('PLUM has Miscanthus in 2010, and some gridcells lose Miscanthus from 2010 to 2011. This cannot be satisfied using my 2010 maps because I have no Miscanthus. But not all grid cells PLUM specifies as losing Miscanthus lose ALL of it! Difference of %0.1g m2. Ignoring.',sum(in1_thisLU_YX(in0_thisLU_YX>0)))
                    end
                    tmp = in_y0.maps_YXv(:,:,i) ;
                    tmp(inDiff_thisLU_YX<0) = in1_thisLU_YX(inDiff_thisLU_YX<0) ;
                    in_y0.maps_YXv(:,:,i) = tmp ;
                    in_y0_2deg.maps_YXv = PLUMharm_aggregate(in_y0.maps_YXv,0.5,2) ;
                    in_y0_agri_YXv = in_y0.maps_YXv(:,:,isAgri) ;
                    in_y0_2deg_agri_YXv = in_y0_2deg.maps_YXv(:,:,isAgri) ;
                elseif 0 <= NETglobChg_thisLU && NETglobChg_thisLU < conserv_tol_area
                    warning('Total global loss of %s (%0.4e m2) exceeds its y0 area (%0.4e m2). LIKELY WILL cause infinite loop in ringRedist.\n',thisLU,-globLoss_thisLU,globArea_thisLU) ;
                else
                    warning('Total global loss of %s (%0.4e m2) exceeds its y0 area (%0.4e m2). (((May))) cause infinite loop in ringRedist.\n',thisLU,-globLoss_thisLU,globArea_thisLU) ;
                end
                pause(0.1)
            end
        end; clear *_thisLU*
        
        % Get maximum allowed N fertilization (maximum seen for this crop
        % in any PLUM output thus far, or in LUH2, or in any harmonized
        % output thus far [although no harmonized output should exceed max
        % anyway)
        max_harm_nfert_y0 = squeeze(max(max(out_y0_nfert_2deg.maps_YXv,[],1),[],2)) ;
        if ~exist('max_nfert_y0','var')
            max_nfert_y1 = max([max_orig_nfert_y0 max_orig_nfert_y1 max_harm_nfert_y0],[],2) ;
        else
            max_nfert_y1 = max([max_orig_nfert_y0 max_orig_nfert_y1 max_harm_nfert_y0 max_nfert_y0],[],2) ;
        end
        max_nfert_y1_YXv = repmat(permute(max_nfert_y1,[2 3 1]),[size(landArea_2deg_YX) 1]) ;
        tooMuchNfert_YXv = in_y0_nfert_2deg.maps_YXv > max_nfert_y1_YXv ;
                
        % Do not allow invalid management inputs. This will only affect
        % cells that did not have any of this crop in the original PLUM
        % outputs.
        in_y0_nfert_2deg.maps_YXv(in_y0_nfert_2deg.maps_YXv<0) = 0 ;
        in_y1_nfert_2deg.maps_YXv(in_y1_nfert_2deg.maps_YXv<0) = 0 ;
        in_y0_nfert_2deg.maps_YXv(tooMuchNfert_YXv) = max_nfert_y1_YXv(tooMuchNfert_YXv) ;
        in_y1_nfert_2deg.maps_YXv(tooMuchNfert_YXv) = max_nfert_y1_YXv(tooMuchNfert_YXv) ;
        in_y0_irrig_2deg.maps_YXv(in_y0_irrig_2deg.maps_YXv<0) = 0 ;
        in_y1_irrig_2deg.maps_YXv(in_y1_irrig_2deg.maps_YXv<0) = 0 ;
        in_y0_irrig_2deg.maps_YXv(in_y0_irrig_2deg.maps_YXv>1) = 1 ;
        in_y1_irrig_2deg.maps_YXv(in_y1_irrig_2deg.maps_YXv>1) = 1 ;
        
%         disp('After limiting to [0, max]')
%         tmp0 = 1e-6*1e-3*sum(sum(in_y0_nfert_2deg.maps_YXv(:,:,2) .* in_y0_2deg.maps_YXv(:,:,2))) ;
%         tmp1 = 1e-6*1e-3*sum(sum(in_y1_nfert_2deg.maps_YXv(:,:,2) .* in_y1_2deg.maps_YXv(:,:,2))) ;
%         fprintf('in_%d glob Nfert (Mt): %0.3f\n', thisYear-1, tmp0)
%         fprintf('in_%d glob Nfert (Mt): %0.3f\n', thisYear,   tmp1)
%         fprintf('              diff (Mt): %0.3f\n\n', tmp1-tmp0)
        
        % Debugging output
        if do_debug
            PLUMharm_debugOut_deltas('iny0_to_iny1','areas',in_y0_2deg.maps_YXv,in_y1_2deg.maps_YXv,debugIJ_2deg,LUnames)
            PLUMharm_debugOut_deltas('iny0_to_iny1','Nfert', ...
                in_y0_nfert_2deg.maps_YXv .* in_y0_2deg.maps_YXv(:,:,isCrop), ...
                in_y1_nfert_2deg.maps_YXv .* in_y1_2deg.maps_YXv(:,:,isCrop), ...
                debugIJ_2deg, LPJGcrops)
        end
        
        % Calculate changes in PLUM agri grids at 2 degrees
        %%% Negative indicates LOSS of thisAgri area or mgmt
        agri_d_YXv = in_y1_2deg_agri_YXv - in_y0_2deg_agri_YXv ;

        % Apply changes to previous grid (@2deg)
        mid1_y1_2deg_agri_YXv = out_y0_2deg_agri_YXv + agri_d_YXv ;
        
        % Debugging output
        if do_debug
            mid1_y1_2deg_ntrl_YX = out_y0_2deg.maps_YXv(:,:,strcmp(out_y0_2deg.varNames,'NATURAL')) ...
                                 - sum(agri_d_YXv,3) ;
            PLUMharm_debugOut_deltas('outy0_to_mid1y1', 'areas',...
                out_y0_2deg.maps_YXv, ...
                cat(3,mid1_y1_2deg_agri_YXv,mid1_y1_2deg_ntrl_YX,base_bare_2deg_YX), ...
                debugIJ_2deg,LUnames) ;
        end

        % Get land available for conversion to agriculture,
        % as max(NATURAL-RESERVED,0)
        out_y0_2deg_ntrl_YX = out_y0_2deg.maps_YXv(:,:,strcmp(out_y0_2deg.varNames,'NATURAL')) ;
        nonResNtrl_YX = max(out_y0_2deg_ntrl_YX - resArea_2deg_YX, 0) ;

        % Check that neither crop nor pasture exceeds nonBare area;
        % check that neither crop nor pasture is negative;
        % check that nonBare area is conserved;
        % compute the total amount of crop or pasture increase / decrease that
        % is not able to be met within the 2 degree gridcells
        %%% Negative total_unmet indicates TOO MUCH LOSS of thisAgri area
        [total_unmet_agri_YXv, ...
            mid_y1_2deg_agri_YXv, mid_y1_2deg_ntrl_YX] = ...
            PLUMharm_getUnmet_cropAreaRes(...
            mid1_y1_2deg_agri_YXv, base_vegd_2deg_YX, ...
            resArea_2deg_YX, sum(out_y0_2deg_agri_YXv,3), debugIJ_2deg) ;
                
        % Check 2: Check that global area changes are (mostly) conserved
        for i = 1:Nagri
            agri_d2_YX = mid_y1_2deg_agri_YXv(:,:,i) - out_y0_2deg_agri_YXv(:,:,i) + total_unmet_agri_YXv(:,:,i) ;
            agri_d_YX = agri_d_YXv(:,:,i) ;
            agri_d_glob_1 = sum(agri_d_YX(:)) ;
            agri_d_glob_2 = sum(agri_d2_YX(:)) ;
            if abs((agri_d_glob_2-agri_d_glob_1)/agri_d_glob_1*100) > conserv_tol_pct
                error(['Global ' LUnames{i} ' area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (step 2)'])
            end
        end
        
        % Check 2: Check that cells' vegetated area does not change much
        mid_y1_2deg_vegd_YX = sum(mid_y1_2deg_agri_YXv,3) + out_y0_2deg_ntrl_YX - sum(agri_d_YXv,3) + sum(total_unmet_agri_YXv,3) ;
        nonCons_pct_YX = 100*(mid_y1_2deg_vegd_YX - base_vegd_2deg_YX) ./ base_vegd_2deg_YX ;
        if max(abs(nonCons_pct_YX(:))) > conserv_tol_pct
            error('By-cell vegetated area not conserved! Worst errors %0.1e and %0.1e.\n', ...
                min(nonCons_pct_YX(:)), max(nonCons_pct_YX(:))) ;
        end
        
        % Debugging
        if verbose
            disp(' ')
            disp('After addition of deltas')
            tmp0 = 1e-6*1e-6*sum(sum(sum(out_y0_2deg.maps_YXv(:,:,isAgri)))) ;
            tmp1 = 1e-6*1e-6*sum(mid_y1_2deg_agri_YXv(:)) ;
            fprintf('out_%d glob agri area (million km2): %0.1f\n', thisYear-1, tmp0)
            fprintf('mid_%d glob agri area (million km2): %0.1f\n', thisYear,   tmp1)
            fprintf('                    diff (million km2): %0.1f\n\n', tmp1-tmp0)

            tmp0 = 1e-6*1e-6*sum(out_y0_2deg.maps_YXv(:)) ;
            tmp1 = 1e-6*1e-6*sum(sum(sum(mid_y1_2deg_agri_YXv,3) ...
                                   + sum(total_unmet_agri_YXv,3) ...
                                   + mid_y1_2deg_ntrl_YX ...
                                   + in_y1_2deg.maps_YXv(:,:,~notBare))) ;
            fprintf('out_%d glob land area (million km2): %0.1f\n', thisYear-1, tmp0)
            fprintf('mid_%d glob land area (million km2): %0.1f\n', thisYear,   tmp1)
            fprintf('                    diff (million km2): %0.1f\n\n', tmp1-tmp0)
            tmp = 1e-6*1e-6*sum(base_vegd_2deg_YX(:)) ;
            fprintf('base_2deg glob vegd area (million km2): %0.1f\n', tmp)
            clear tmp*
        end

        % Debugging output
        if do_debug
            PLUMharm_debugOut_deltas('outy0_to_midy1', 'areas', ...
                out_y0_2deg.maps_YXv, ...
                cat(3,mid_y1_2deg_agri_YXv,mid_y1_2deg_ntrl_YX,base_bare_2deg_YX), ...
                debugIJ_2deg,LUnames) ;
        end
        
        % Loop through every 2-degree gridcell. If a gridcell has unmet crop
        % or pasture, look for place to put this unmet amount in neighboring
        % rings, starting with gridcells that are 1 unit away, then 2, etc.
        % until all unmet has been displaced to new 2 degree cells. Track the
        % displaced crop and pasture and the # of "rings" needed for each
        % 2-degree gridcell.
        [out_y1_2deg_agri_YXv, out_y1_2deg_ntrl_YX] = ...
            PLUMharm_ringRedist_areaCropsRes(...
            mid_y1_2deg_agri_YXv, ...
            total_unmet_agri_YXv, ...
            landArea_2deg_YX, ...
            [], in_y0orig_2deg, in_y1orig_2deg, out_y0_2deg_agri_YXv, ...
            mid_y1_2deg_ntrl_YX, resArea_2deg_YX, LUnames_agri) ;
        out_y1_2deg_vegd_YX = sum(cat(3,out_y1_2deg_agri_YXv,out_y1_2deg_ntrl_YX),3) ;
        
        % Debugging
        if verbose
            disp(' ')
            disp('After ringRedist')
            tmp0 = 1e-6*1e-6*sum(sum(sum(out_y0_2deg.maps_YXv(:,:,isAgri)))) ;
            tmp1 = 1e-6*1e-6*sum(sum(sum(out_y1_2deg_agri_YXv))) ;
            fprintf('out_%d glob agri area (million km2): %0.1f\n', thisYear-1, tmp0)
            fprintf('out_%d glob agri area (million km2): %0.1f\n', thisYear,   tmp1)
            fprintf('                    diff (million km2): %0.1f\n\n', tmp1-tmp0)
            tmp0 = 1e-6*1e-6*sum(out_y0_2deg.maps_YXv(:)) ;
            tmp1 = 1e-6*1e-6*sum(sum(out_y1_2deg_vegd_YX ...
                                   + in_y1_2deg.maps_YXv(:,:,~notBare))) ;
            fprintf('out_%d glob land area (million km2): %0.1f\n', thisYear-1, tmp0)
            fprintf('out_%d glob land area (million km2): %0.1f\n', thisYear,   tmp1)
            fprintf('                    diff (million km2): %0.1f\n\n', tmp1-tmp0)
            tmp = 1e-6*1e-6*sum(base_vegd_2deg_YX(:)) ;
            fprintf('base_2deg glob vegd area (million km2): %0.1f\n', tmp)
            out_y1_2deg_agri_YX = sum(out_y1_2deg_agri_YXv,3) ;
            tmp = max(out_y1_2deg_agri_YX(:) - base_vegd_2deg_YX(:)) ;
            fprintf('Max agri-vegd area (km2): %0.1f\n', tmp)
            clear tmp*
        end
        
        % Debugging output
        if do_debug
            PLUMharm_debugOut_deltas('outy0_to_outy1', 'areas', ...
                out_y0_2deg.maps_YXv, ...
                cat(3,out_y1_2deg_agri_YXv, out_y1_2deg_ntrl_YX, base_bare_2deg_YX), ...
                debugIJ_2deg,LUnames) ;
        end
        
        % Check 3: Check that agricultural area does not exceed vegetated
        % area
        out_y1_2deg_agri_YX = sum(out_y1_2deg_agri_YXv,3) ;
        if any(out_y1_2deg_agri_YX(:) - base_vegd_2deg_YX(:) > conserv_tol_area)
            error('More agricultural than vegetated area! (Check 3)')
        end

        % Check 3: Check that global area changes are (mostly) conserved
        diffH_YXv = out_y1_2deg_agri_YXv - out_y0_2deg_agri_YXv ;
        for i = 1:Nagri
            agri_d_YX = agri_d_YXv(:,:,i) ;
            agri_d_glob_1 = sum(agri_d_YX(:)) ;
            agri_d3_YX = out_y1_2deg_agri_YXv(:,:,i) - out_y0_2deg_agri_YXv(:,:,i) ;
            agri_d_glob_3 = sum(agri_d3_YX(:)) ;
            if abs((agri_d_glob_3-agri_d_glob_1)/agri_d_glob_1*100) > conserv_tol_pct
                warning(['Global ' LUnames{i} ' area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (' num2str((agri_d_glob_3-agri_d_glob_1)/agri_d_glob_1*100) '%) (step 3)'])
            end
        end
        
        % Make sure that managements are within acceptable bounds,
        % considering (a) how much each cell had available to lose and (b)
        % how much room each cell has to gain before exceeding max rate
        [total_unmet_nfert_YXv, mid1_y1_2deg_nfert_YXv] = ...
            PLUMharm_getUnmet_mgmt(...
                out_y0_nfert_2deg.maps_YXv, out_y0_2deg.maps_YXv(:,:,isCrop), ...
                in_y0_nfert_2deg.maps_YXv, in_y1_nfert_2deg.maps_YXv, ...
                in_y0_2deg.maps_YXv(:,:,isCrop), in_y1_2deg.maps_YXv(:,:,isCrop), ...
                out_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), max_nfert_y1) ;
        [total_unmet_irrig_YXv, mid1_y1_2deg_irrig_YXv] = ...
            PLUMharm_getUnmet_mgmt(...
                out_y0_irrig_2deg.maps_YXv, out_y0_2deg.maps_YXv(:,:,isCrop), ...
                in_y0_irrig_2deg.maps_YXv, in_y1_irrig_2deg.maps_YXv, ...
                in_y0_2deg.maps_YXv(:,:,isCrop), in_y1_2deg.maps_YXv(:,:,isCrop), ...
                out_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ones(size(max_nfert_y1))) ;
            
        % Do not allow invalid management inputs.
        if any(mid1_y1_2deg_nfert_YXv(:)<-1e-6)
            error('Negative value(s) in mid1_y1_2deg_nfert_YXv!')
        elseif any(mid1_y1_2deg_irrig_YXv(:)<0)
            error('Negative value(s) in mid1_y1_2deg_irrig_YXv!')
        elseif any(mid1_y1_2deg_irrig_YXv(:)>1)
            error('Value(s) >1 in mid1_y1_2deg_irrig_YXv!')
        end
            
        % Check 2b: Check that management changes are (mostly) conserved
%         disp('Check 2b')
%         tmp0 = 1e-6*1e-3*sum(sum(in_y0_nfert_2deg.maps_YXv(:,:,2) .* in_y0_2deg.maps_YXv(:,:,2))) ;
%         tmp1 = 1e-6*1e-3*sum(sum(in_y1_nfert_2deg.maps_YXv(:,:,2) .* in_y1_2deg.maps_YXv(:,:,2))) ;
%         fprintf('in_%d glob Nfert (Mt): %0.3f\n', thisYear-1, tmp0)
%         fprintf('in_%d glob Nfert (Mt): %0.3f\n', thisYear,   tmp1)
%         fprintf('            diff (Mt): %0.3f\n\n', tmp1-tmp0)
%         tmp0 = 1e-6*1e-3*sum(sum(out_y0_2deg_agri_YXv(:,:,2) .* out_y0_nfert_2deg.maps_YXv(:,:,2))) ;
%         tmp1 = 1e-6*1e-3*sum(sum(out_y1_2deg_agri_YXv(:,:,2) .* mid1_y1_2deg_nfert_YXv(:,:,2) + total_unmet_nfert_YXv(:,:,2))) ;
%         fprintf(' out_%d glob Nfert (Mt): %0.3f\n', thisYear-1, tmp0)
%         fprintf('mid1_%d glob Nfert (Mt): %0.3f\n', thisYear,   tmp1)
%         fprintf('              diff (Mt): %0.3f\n\n', tmp1-tmp0)
        PLUMharm_check_conservation(...
            out_y0_nfert_2deg.maps_YXv, out_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            mid1_y1_2deg_nfert_YXv, out_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y0_nfert_2deg.maps_YXv, in_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y1_nfert_2deg.maps_YXv, in_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            total_unmet_nfert_YXv, LPJGcrops, conserv_tol_pct, ...
            '2b nfert') ;
        PLUMharm_check_conservation(...
            out_y0_irrig_2deg.maps_YXv, out_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            mid1_y1_2deg_irrig_YXv, out_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y0_irrig_2deg.maps_YXv, in_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y1_irrig_2deg.maps_YXv, in_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            total_unmet_irrig_YXv, LPJGcrops, conserv_tol_pct, ...
            '2b irrig') ;
        
% % %         if do_debug
% % %             PLUMharm_debugOut_deltas('outy0_to_midy1', 'Nfert', ...
% % %                 out_y0_2deg_nfert_YXv.*out_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
% % %                 mid1_y1_2deg_nfert_YXv.*out_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
% % %                 debugIJ_2deg,LPJGcrops) ;
% % %         end
                
        % Do ring redistribution for management inputs
        [out_y1_2deg_nfert_YXv, total_unmet2_nfert_YXv] = PLUMharm_ringRedist_mgmt(...
            mid1_y1_2deg_nfert_YXv, out_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            total_unmet_nfert_YXv, max_nfert_y1, ...
            LPJGcrops, [], ...
            out_y0_nfert_2deg, out_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y0_nfert_2deg, in_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y1_nfert_2deg, in_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            conserv_tol_pct, ...
            'ringRedist nfert') ;
%         tmp0 = 1e-6*1e-3*sum(sum(out_y0_2deg.maps_YXv(:,:,2) .* out_y0_nfert_2deg.maps_YXv(:,:,2))) ;
%         tmp1 = 1e-6*1e-3*sum(sum(out_y1_2deg_agri_YXv(:,:,2) .* out_y1_2deg_nfert_YXv(:,:,2) + total_unmet2_nfert_YXv(:,:,2))) ;
%         fprintf(' out_%d glob Nfert (Mt): %0.3f\n', thisYear-1, tmp0)
%         fprintf('mid1_%d glob Nfert (Mt): %0.3f\n', thisYear,   tmp1)
%         fprintf('              diff (Mt): %0.3f\n\n', tmp1-tmp0)
        
        [out_y1_2deg_irrig_YXv, total_unmet2_irrig_YXv] = PLUMharm_ringRedist_mgmt(...
            mid1_y1_2deg_irrig_YXv, out_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            total_unmet_irrig_YXv, ones(size(max_nfert_y1)), ...
            LPJGcrops, [], ...
            out_y0_irrig_2deg, out_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y0_irrig_2deg, in_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y1_irrig_2deg, in_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            conserv_tol_pct, ...
            'ringRedist irrig') ;
        
        % Check 3b: Check that management changes are (mostly) conserved
        PLUMharm_check_conservation(...
            out_y0_nfert_2deg.maps_YXv, out_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            out_y1_2deg_nfert_YXv, out_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y0_nfert_2deg.maps_YXv, in_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y1_nfert_2deg.maps_YXv, in_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            total_unmet2_nfert_YXv, LPJGcrops, conserv_tol_pct, ...
            '3b nfert') ;
        PLUMharm_check_conservation(...
            out_y0_irrig_2deg.maps_YXv, out_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            out_y1_2deg_irrig_YXv, out_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y0_irrig_2deg.maps_YXv, in_y0_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            in_y1_irrig_2deg.maps_YXv, in_y1_2deg_agri_YXv(:,:,~strcmp(LUnames_agri,'PASTURE')), ...
            total_unmet2_irrig_YXv, LPJGcrops, conserv_tol_pct, ...
            '3b irrig') ;
        
        % Do not allow invalid management inputs.
        if any(out_y1_2deg_nfert_YXv(:)<0)
            error('Negative value(s) in out_y1_2deg_nfert_YXv!')
        elseif any(out_y1_2deg_irrig_YXv(:)<0)
            error('Negative value(s) in out_y1_2deg_irrig_YXv!')
        elseif any(out_y1_2deg_irrig_YXv(:)>1+1e-6)
            error('Value(s) >1 in out_y1_2deg_irrig_YXv!')
        end

        % Loop through all 2-degree gridcells and distribute the new crop and
        % pasture area changes to the half-degree gridcells within.
        %
        % If the crop or pasture change is a decrease, apply the same
        % PERCENTAGE decrease to all half-degree crop or pasture gridcells 
        % within the 2-degree cell.
        %
        % If the crop or pasture change is an increase, apply this increase to 
        % all half-degree gridcells within the 2-degree cell proportionally 
        % according to available land area. Make sure that total area within a 
        % half-degree gridcell does not exceed 1 or 0.
        out_y1_agri_YXv = ...
            PLUMharm_distDeltas_areaCrops_recursive( ...
            landArea_YX, landArea_2deg_YX, ...
            out_y0_2deg_agri_YXv, out_y1_2deg_agri_YXv, out_y0_agri_YXv, ...
            base_vegd_YX, base_vegd_YX, conserv_tol_pct, conserv_tol_area, LUnames) ;
        
        % Debugging
        if verbose
            disp(' ')
            disp('Now at half-degree')
            tmp0 = 1e-6*1e-6*sum(sum(sum(out_y0.maps_YXv(:,:,isAgri)))) ;
            tmp1 = 1e-6*1e-6*sum(out_y1_agri_YXv(:)) ;
            fprintf('out_%d glob agri area (million km2): %0.1f\n', thisYear-1, tmp0)
            fprintf('out_%d glob agri area (million km2): %0.1f\n', thisYear,   tmp1)
            fprintf('                    diff (million km2): %0.1f\n\n', tmp1-tmp0)
            out_y1_ntrl_YX_tmp = landArea_YX ...
                - (sum(out_y1_agri_YXv,3) + in_y1.maps_YXv(:,:,~notBare)) ;
            tmp0 = 1e-6*1e-6*sum(out_y0.maps_YXv(:)) ;
            tmp1 = 1e-6*1e-6*sum(sum(sum(out_y1_agri_YXv,3) ...
                                   + out_y1_ntrl_YX_tmp ...
                                   + in_y1.maps_YXv(:,:,~notBare))) ;
            fprintf('out_%d glob land area (million km2): %0.1f\n', thisYear-1, tmp0)
            fprintf('out_%d glob land area (million km2): %0.1f\n', thisYear,   tmp1)
            fprintf('                    diff (million km2): %0.1f\n\n', tmp1-tmp0)
            tmp = 1e-6*1e-6*sum(base_vegd_2deg_YX(:)) ;
            fprintf('base_2deg glob vegd area (million km2): %0.1f\n', tmp)
            clear tmp*
        end
        
        % Distribute management inputs from 2-degree to half-degree cells.
        out_y1_nfert_YXv = PLUMharm_distMgmt(out_y1_2deg_nfert_YXv,2,0.5) ;
        out_y1_irrig_YXv = PLUMharm_distMgmt(out_y1_2deg_irrig_YXv,2,0.5) ;

        % Ensure cells in range [0,in_y1_vegd_YX] (extreme deviations were checked
        % in loop)
        if any(isnan(out_y1_agri_YXv(:)))
            error('How did you get a NaN in out_y1_agri_YXv?')
        end
        out_y1_agri_YXv(out_y1_agri_YXv<0) = 0 ;
        base_vegd_YXv = repmat(base_vegd_YX,[1 1 Nagri]) ;
        out_y1_agri_YXv(out_y1_agri_YXv>base_vegd_YXv) = out_y1_agri_YXv(out_y1_agri_YXv>base_vegd_YXv) ;
        clear base_vegd_YXv

        % Check 5: Check that global area changes are (mostly) conserved
        for i = 1:Nagri
            this_d1_halfDeg_YX = in_y1_agri_YXv(:,:,i) - in_y0_agri_YXv(:,:,i) ;
            this_d5_halfDeg_YX = out_y1_agri_YXv(:,:,i) - out_y0_agri_YXv(:,:,i) ;
            this_d_glob_halfDeg_1 = sum(this_d1_halfDeg_YX(:)) ;
            this_d_glob_halfDeg_5 = sum(this_d5_halfDeg_YX(:)) ;
            if abs((this_d_glob_halfDeg_5-this_d_glob_halfDeg_1)/this_d_glob_halfDeg_1*100) > conserv_tol_pct
                error(['Global ' LUnames{i} ' area changes are not conserved to within ' num2str(conserv_tol_pct) '%! (' num2str((this_d_glob_halfDeg_5-this_d_glob_halfDeg_1)/this_d_glob_halfDeg_1*100) '%) (step 5)'])
            end
        end
        
        % Check 5: Check that global mgmt changes are (mostly) conserved
        for i = 1:Ncrops_lpjg
            this_d1_halfDeg_YX = in_y1_nfert_2deg.maps_YXv(:,:,i) .* in_y1_2deg.maps_YXv(:,:,i) ...
                               - in_y0_nfert_2deg.maps_YXv(:,:,i) .* in_y0_2deg.maps_YXv(:,:,i) ;
            this_d5_halfDeg_YX = out_y1_nfert_YXv(:,:,i).*out_y1_agri_YXv(:,:,i) ...
                               - out_y0_nfert_YXv(:,:,i).*out_y0_agri_YXv(:,:,i) ;
            this_d_glob_halfDeg_1 = sum(this_d1_halfDeg_YX(:)) ;
            this_d_glob_halfDeg_5 = sum(this_d5_halfDeg_YX(:)) ;
            if abs((this_d_glob_halfDeg_5-this_d_glob_halfDeg_1)/this_d_glob_halfDeg_1*100) > conserv_tol_pct
                error(['Global ' LUnames{i} ' nfert changes are not conserved to within ' num2str(conserv_tol_pct) '%! (' num2str((this_d_glob_halfDeg_5-this_d_glob_halfDeg_1)/this_d_glob_halfDeg_1*100) '%) (step 5)'])
            end
        end
        for i = 1:Ncrops_lpjg
            this_d1_halfDeg_YX = in_y1_irrig_2deg.maps_YXv(:,:,i) .* in_y1_2deg.maps_YXv(:,:,i) ...
                               - in_y0_irrig_2deg.maps_YXv(:,:,i) .* in_y0_2deg.maps_YXv(:,:,i) ;
            this_d5_halfDeg_YX = out_y1_irrig_YXv(:,:,i).*out_y1_agri_YXv(:,:,i) ...
                               - out_y0_irrig_YXv(:,:,i).*out_y0_agri_YXv(:,:,i) ;
            this_d_glob_halfDeg_1 = sum(this_d1_halfDeg_YX(:)) ;
            this_d_glob_halfDeg_5 = sum(this_d5_halfDeg_YX(:)) ;
            if abs((this_d_glob_halfDeg_5-this_d_glob_halfDeg_1)/this_d_glob_halfDeg_1*100) > conserv_tol_pct
                error(['Global ' LUnames{i} ' irrig changes are not conserved to within ' num2str(conserv_tol_pct) '%! (' num2str((this_d_glob_halfDeg_5-this_d_glob_halfDeg_1)/this_d_glob_halfDeg_1*100) '%) (step 5)'])
            end
        end

        % Get land use areas
        out_y1_bare_YX = in_y1.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
        out_y1_ntrl_YX = landArea_YX ...
            - (sum(out_y1_agri_YXv,3) + out_y1_bare_YX) ;
        out_y1_past_YX = out_y1_agri_YXv(:,:,isAgri_isPast) ;
        out_y1_crop_YX = sum(out_y1_agri_YXv(:,:,~isAgri_isPast),3) ;
        out_y1_2deg_bare_YX = in_y1_2deg.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
        out_y1_2deg_past_YX = out_y1_2deg_agri_YXv(:,:,isAgri_isPast) ;
        out_y1_2deg_crop_YX = sum(out_y1_2deg_agri_YXv(:,:,~isAgri_isPast),3) ;
        if ~exist('out_y1_2deg','var')
            out_y1_2deg.varNames = out_y0_2deg.varNames ;
        end
        out_y1_2deg.maps_YXv = cat(3, out_y1_2deg_agri_YXv, out_y1_2deg_ntrl_YX, out_y1_2deg_bare_YX) ;

%         tmp0 = 1e-6*1e-3*sum(sum(out_y0_2deg.maps_YXv(:,:,2) .* out_y0_nfert_2deg.maps_YXv(:,:,2))) ;
%         tmp1 = 1e-6*1e-3*sum(sum(out_y1_2deg.maps_YXv(:,:,2) .* out_y1_2deg_nfert_YXv(:,:,2))) ;
%         fprintf('out_%d glob Nfert (Mt): %0.3f\n', thisYear-1, tmp0)
%         fprintf('out_%d glob Nfert (Mt): %0.3f\n', thisYear,   tmp1)
%         fprintf('             diff (Mt): %0.3f\n\n', tmp1-tmp0)
%         keyboard
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Write output files %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if save_any
            disp(['  Done processing (' toc_hms(toc) '). Now writing.'])
        end
        
        if save_halfDeg_any
            % Save new LandCoverFract.txt (0.5-degree)
            unix(['mkdir -p ' PLUM_out_top num2str(thisYear)]) ;
            out_y1.varNames = {'PASTURE','CROPLAND','NATURAL','BARREN'} ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'CROPLAND')) = out_y1_crop_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'PASTURE')) = out_y1_past_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'NATURAL')) = out_y1_ntrl_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'BARREN')) = out_y1_bare_YX ;
            out_y1.maps_YXv = out_y1.maps_YXv ./ repmat(landArea_YX,[1 1 length(out_y1.varNames)]) ;
            if save_halfDeg_mat
                file_out = [PLUM_out_top num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.mat'] ;
                save(file_out,'out_y1') ;
            end
            if save_halfDeg_txt
                [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
                file_out = [PLUM_out_top num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.txt'] ;
                lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                    'outPrec', outPrec, ...
                    'outWidth', outWidth, ...
                    'delimiter', delimiter, ...
                    'overwrite', overwrite, ...
                    'fancy', fancy, ...
                    'progress_step_pct', 20, ...
                    'verbose', false) ;
            end
            clear out_y1
            
            % Save new CropFract.txt (0.5-degree)
            out_y1.maps_YXv = out_y1_agri_YXv(:,:,~isAgri_isPast) ./ repmat(out_y1_crop_YX,[1 1 length(LPJGcrops)]) ;
            out_y1.maps_YXv(repmat(out_y1_crop_YX,[1 1 length(LPJGcrops)])==0) = 0 ;
            out_y1.varNames = LPJGcrops ;
            if save_halfDeg_mat
                file_out = [PLUM_out_top num2str(thisYear) '/CropFract.base' num2str(base_year) '.mat'] ;
                save(file_out,'out_y1') ;
            end
            if save_halfDeg_txt
                [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
                file_out = [PLUM_out_top num2str(thisYear) '/CropFract.base' num2str(base_year) '.txt'] ;
                lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                    'outPrec', outPrec, ...
                    'outWidth', outWidth, ...
                    'delimiter', delimiter, ...
                    'overwrite', overwrite, ...
                    'fancy', fancy, ...
                    'progress_step_pct', 20, ...
                    'verbose', false) ;
            end
            clear out_y1
            
            % Save new Fert.txt (0.5-degree)
            % Convert from kg/m2 to kg/ha for compatibility with original
            % PLUM style.
            out_y1.maps_YXv = 1e4*out_y1_nfert_YXv ;
            out_y1.varNames = LPJGcrops ;
            if save_halfDeg_mat
                file_out = [PLUM_out_top num2str(thisYear) '/Fert.base' num2str(base_year) '.mat'] ;
                save(file_out,'out_y1') ;
            end
            if save_halfDeg_txt
                [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
                file_out = [PLUM_out_top num2str(thisYear) '/Fert.base' num2str(base_year) '.txt'] ;
                lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                    'outPrec', outPrec, ...
                    'outWidth', outWidth, ...
                    'delimiter', delimiter, ...
                    'overwrite', overwrite, ...
                    'fancy', fancy, ...
                    'progress_step_pct', 20, ...
                    'verbose', false) ;
            end
            clear out_y1
            
            % Save new Irrig.txt (0.5-degree)
            out_y1.maps_YXv = out_y1_irrig_YXv ;
            out_y1.varNames = LPJGcrops ;
            if save_halfDeg_mat
                file_out = [PLUM_out_top num2str(thisYear) '/Irrig.base' num2str(base_year) '.mat'] ;
                save(file_out,'out_y1') ;
            end
            if save_halfDeg_txt
                [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
                file_out = [PLUM_out_top num2str(thisYear) '/Irrig.base' num2str(base_year) '.txt'] ;
                lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                    'outPrec', outPrec, ...
                    'outWidth', outWidth, ...
                    'delimiter', delimiter, ...
                    'overwrite', overwrite, ...
                    'fancy', fancy, ...
                    'progress_step_pct', 20, ...
                    'verbose', false) ;
            end
            clear out_y1
            
        end

        if save_2deg_any
            % Save new LandCoverFract.txt (2-degree)
            unix(['mkdir -p ' PLUM_out_top num2str(thisYear)]) ;
            out_y1.varNames = {'PASTURE','CROPLAND','NATURAL','BARREN'} ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'CROPLAND')) = out_y1_2deg_crop_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'PASTURE')) = out_y1_2deg_past_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'NATURAL')) = out_y1_2deg_ntrl_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'BARREN')) = out_y1_2deg_bare_YX ;
            out_y1.maps_YXv = out_y1.maps_YXv ./ repmat(landArea_2deg_YX,[1 1 length(out_y1.varNames)]) ;
            if save_2deg_mat
                file_out = [PLUM_out_top num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.2deg.txt'] ;
                save(file_out,'out_y1') ;
            end
            if save_2deg_txt
                [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
                file_out = [PLUM_out_top num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.2deg.txt'] ;
                lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                    'outPrec', outPrec, ...
                    'outWidth', outWidth, ...
                    'delimiter', delimiter, ...
                    'overwrite', overwrite, ...
                    'fancy', fancy, ...
                    'progress_step_pct', 20, ...
                    'verbose', false) ;
            end
            clear out_y1

            % Save new CropFract.txt (2-degree)
            out_y1.maps_YXv = out_y1_2deg_agri_YXv(:,:,~isAgri_isPast) ./ repmat(out_y1_2deg_crop_YX,[1 1 length(LPJGcrops)]) ;
            out_y1.maps_YXv(repmat(out_y1_2deg_crop_YX,[1 1 length(LPJGcrops)])==0) = 0 ;
            out_y1.varNames = LPJGcrops ;
            if save_2deg_mat
                file_out = [PLUM_out_top num2str(thisYear) '/CropFract.base' num2str(base_year) '.2deg.txt'] ;
                save(file_out,'out_y1') ;
            end
            if save_2deg_txt
                [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
                file_out = [PLUM_out_top num2str(thisYear) '/CropFract.base' num2str(base_year) '.2deg.txt'] ;
                lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                    'outPrec', outPrec, ...
                    'outWidth', outWidth, ...
                    'delimiter', delimiter, ...
                    'overwrite', overwrite, ...
                    'fancy', fancy, ...
                    'progress_step_pct', 20, ...
                    'verbose', false) ;
            end
            clear out_y1
            
            % Save new Fert.txt (2-degree)
            % Convert from kg/m2 to kg/ha for compatibility with original
            % PLUM style.
            out_y1.maps_YXv = 1e4*out_y1_2deg_nfert_YXv ;
            out_y1.varNames = LPJGcrops ;
            if save_2deg_mat
                file_out = [PLUM_out_top num2str(thisYear) '/Fert.base' num2str(base_year) '.2deg.txt'] ;
                save(file_out,'out_y1') ;
            end
            if save_2deg_txt
                [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
                file_out = [PLUM_out_top num2str(thisYear) '/Fert.base' num2str(base_year) '.2deg.txt'] ;
                lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                    'outPrec', outPrec, ...
                    'outWidth', outWidth, ...
                    'delimiter', delimiter, ...
                    'overwrite', overwrite, ...
                    'fancy', fancy, ...
                    'progress_step_pct', 20, ...
                    'verbose', false) ;
            end
            clear out_y1
            
            % Save new Irrig.txt (2-degree)
            out_y1.maps_YXv = out_y1_2deg_irrig_YXv ;
            out_y1.maps_YXv(repmat(out_y1_2deg_crop_YX,[1 1 length(LPJGcrops)])==0) = 0 ;
            out_y1.varNames = LPJGcrops ;
            if save_2deg_mat
                file_out = [PLUM_out_top num2str(thisYear) '/Irrig.base' num2str(base_year) '.2deg.txt'] ;
                save(file_out,'out_y1') ;
            end
            if save_2deg_txt
                [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
                file_out = [PLUM_out_top num2str(thisYear) '/Irrig.base' num2str(base_year) '.2deg.txt'] ;
                lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                    'outPrec', outPrec, ...
                    'outWidth', outWidth, ...
                    'delimiter', delimiter, ...
                    'overwrite', overwrite, ...
                    'fancy', fancy, ...
                    'progress_step_pct', 20, ...
                    'verbose', false) ;
            end
            clear out_y1
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Prepare for next iteration %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clear out_y1_past_YX
        if y < length(years)
            in_y0 = in_y1 ;
            in_y0_2deg = in_y1_2deg ;
            in_y0_agri_YXv = in_y1_agri_YXv ;
            in_y0_2deg_agri_YXv = in_y1_2deg_agri_YXv ;
            [~,IA,IB] = intersect(out_y0.varNames,LUnames_agri,'stable') ;
            out_y0.maps_YXv(:,:,IA) = out_y1_agri_YXv(:,:,IB) ;
            clear IA IB
            out_y0.maps_YXv(:,:,strcmp(LUnames,'NATURAL')) = out_y1_ntrl_YX ;
            out_y0.maps_YXv(:,:,strcmp(LUnames,'BARREN')) = out_y1_bare_YX ;
            out_y0_2deg = out_y1_2deg ;
            out_y0_2deg_agri_YXv = out_y1_2deg_agri_YXv ;
            out_y0_agri_YXv = out_y1_agri_YXv ;
            bareFrac_y0_YX = in_y1_bareFrac_YX ;
            out_y0_nfert_2deg.maps_YXv = out_y1_2deg_nfert_YXv ;
            out_y0_irrig_2deg.maps_YXv = out_y1_2deg_irrig_YXv ;
            in_y0_2deg_vegd_YX = in_y1_2deg_vegd_YX ;
            in_y0_nfert_2deg = in_y1_nfert_2deg ;
            in_y0_irrig_2deg = in_y1_irrig_2deg ;
            max_orig_nfert_y0 = max_orig_nfert_y1 ;
            max_nfert_y0 = max_nfert_y1 ;
            out_y0_nfert_YXv = out_y1_nfert_YXv ;
            out_y0_irrig_YXv = out_y1_irrig_YXv ;
            if do_debug
                in_y0orig_2deg = in_y1orig_2deg ;
            end
            clear *y1*
        end
        
        % Save full-precision outputs for use in restarting
        save([PLUM_out_top num2str(thisYear) 'post.base' num2str(base_year) '.mat'], ...
            '*y0*','latestPLUMin_*','-v7.3') ;

        disp(['  Done (' toc_hms(toc) ').'])


    end % years loop

    disp('Done')

end






