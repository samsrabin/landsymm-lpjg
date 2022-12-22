%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process PLUM outputs for LPJ-GUESS, as PLUM types %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisVer = '20210217' ;
% Like 20181031, but removed all "pasture-only" stuff, which I think had
% been intended to specify what kind of grass should grow on pastures.

exclude_nonlpjg_cells = true ;

% List of directories to process.
% Specify paths relative to parent directory:
%   - /Users/Shared/PLUM/PLUM_outputs_for_LPJG
%   - /Volumes/*/PLUM/PLUM_outputs_for_LPJG
%   - /Volumes/Reacher/LandSyMM/inputs/LU/
% If none of those are the parent directory, modify find_PLUM2LPJG_inputs()
dirList = {...
%     'SSP1.v11.s1';
%     'SSP3.v11.s1';
%     'SSP4.v11.s1';
%     'SSP5.v11.s1';
%     'HE_ensemble/baseline/s1';
%     'halfearth/HEoct/baseline/s1';
    'halfearth/HEoct/halfearth/s1';
    };

% Can turn some of these off if not needed, to spare memory.
save_gl = true ; % gridlist
save_lu = true ;
save_cropfracs = true ;
save_nfert = true ;
save_irr = true ;
save_yield = true ;

% Trying to avoid new crop spinup time
y1_pre = Inf ;    % Will repeat first PLUMout year for y1_pre:(y1-1)
someofall = false ; % Make it so that each gridcell always has at least some tiny amount of every crop


%% Setup

disp('Setting up...')

% Unit conversion factors
cf_ktNha_kgNm2 = 1e-4 ;
cf_tPha_kgPm2 = 0.1 ;

% Settings for outputs of this script
outPrec_LC = 6 ;
outPrec_mgmtInputs = outPrec_LC ;
force_overwrite = true ;
fclose_every = 1000 ;
pct_progress = 25 ;

% Path to this script (so we can copy it to the output directory)
scripts_dir = '/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/MATLAB_work' ;
addpath(scripts_dir)
thisScript = sprintf('%s/process_PLUMout_forLPJGin_%s.m', ...
    scripts_dir, thisVer) ;
if ~exist(thisScript, 'file')
    error('thisScript not found: %s', thisScript)
end

PLUMcrops = {'wheat';
             'maize';
             'rice';
             'oilcrops';
             'pulses';
             'starchyRoots';
             'energycrops';
             'setaside';
             'fruitveg';
             'sugar'} ;
LPJGcrops = {'CerealsC3';
             'CerealsC4';
             'Rice';
             'Oilcrops';
             'Pulses';
             'StarchyRoots';
             'Miscanthus';
             'ExtraCrop';
             'FruitAndVeg';
             'Sugar'} ;
LPJGcrops_mid = [strcat(LPJGcrops,'i') ; LPJGcrops] ;
Ncfts_plum = length(PLUMcrops) ;
Ncfts_lpjg = length(LPJGcrops) ;
PLUMtoLPJG = {} ;
PLUMtoLPJG{strcmp(LPJGcrops,'CerealsC3')} = {'wheat'} ;
PLUMtoLPJG{strcmp(LPJGcrops,'CerealsC4')} = {'maize'} ;
PLUMtoLPJG{strcmp(LPJGcrops,'Rice')} = {'rice'} ;
PLUMtoLPJG{strcmp(LPJGcrops,'Oilcrops')} = {'oilcrops'} ;
PLUMtoLPJG{strcmp(LPJGcrops,'Pulses')} = {'pulses'} ;
PLUMtoLPJG{strcmp(LPJGcrops,'StarchyRoots')} = {'starchyRoots'} ;
PLUMtoLPJG{strcmp(LPJGcrops,'Miscanthus')} = {'energycrops'} ;
PLUMtoLPJG{strcmp(LPJGcrops,'ExtraCrop')} = {'setaside'} ;
PLUMtoLPJG{strcmp(LPJGcrops,'FruitAndVeg')} = {'fruitveg'} ;
PLUMtoLPJG{strcmp(LPJGcrops,'Sugar')} = {'sugar'} ;

% Set up anonymous functions
getPi = @(x) find(strcmp(PLUMcrops,x)) ;
getLi = @(x) find(strcmp(LPJGcrops,x)) ;

% Get cells present in previous LPJ-GUESS output
lpjg_in = lpjgu_matlab_readTable('/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851/yield.out', ...
    'dont_save_MAT', true, 'do_save_MAT', false, ...
    'verboseIfNoMat',false) ;
lons_lpjg = lpjg_in.Lon ;
lats_lpjg = lpjg_in.Lat ;
clear lpjg_in

% Get land uses from previous LPJ-GUESS run
LUfile = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
lpj_lumap = lpjgu_matlab_readTable_then2map(LUfile) ;
barren_lpj_map = lpj_lumap.maps_YXvy(:,:,strcmp(lpj_lumap.varNames,'BARREN'),end) ;
vegd_lpj_map = sum(lpj_lumap.maps_YXvy(:,:,strcmp(lpj_lumap.varNames,'NATURAL')|strcmp(lpj_lumap.varNames,'CROPLAND')|strcmp(lpj_lumap.varNames,'PASTURE'),end),3) ;
lons4map = -0.25 + -179.75:0.5:179.75 ;
lats4map = -0.25 + -89.75:0.5:89.75 ;
lons_map = repmat(lons4map,[length(lats4map) 1]) ;
lats_map = repmat(lats4map',[1 length(lons4map)]) ;
clear lpj_lumap

disp('Done setting up.')


%% Check mapping

% Check that every PLUM crop is mapped
wheatfound2x = false ;
setasidefound2x = false ;
for cP = 1:Ncfts_plum
    thisCrop_plum = PLUMcrops{cP} ;
    found1x = false ;
    for cL = 1:Ncfts_lpjg
        Nfound = length(find(strcmp(PLUMtoLPJG{cL},thisCrop_plum))) ;
        if Nfound==1
            if ~found1x
                found1x = true ;
            elseif strcmp(thisCrop_plum,'wheat')
                if wheatfound2x
                    error('wheat found in more than two elements of PLUMtoLPJG!') ;
                else
                    wheatfound2x = true ;
                end
            elseif strcmp(thisCrop_plum,'setaside')
                if setasidefound2x
                    error('setaside found in more than two elements of PLUMtoLPJG!') ;
                else
                    setasidefound2x = true ;
                end
            else
                error([thisCrop_plum ' found in more than one element of PLUMtoLPJG!']) ;
            end
        elseif Nfound > 1
            error([thisCrop_plum ' found more than once in PLUMtoLPJG{' num2str(cL) '}!'])
        end
    end ; clear cL
    if ~found1x
        error([thisCrop_plum ' not found in PLUMtoLPJG!']) ;
    end
    clear thisCrop_plum found
end ; clear cP
if wheatfound2x && setasidefound2x
    disp('Every PLUM crop is mapped once (except for wheat and setaside each being mapped twice, correctly).')
elseif wheatfound2x
    disp('Every PLUM crop is mapped once (except for wheat being mapped twice, correctly).')
elseif setasidefound2x
    disp('Every PLUM crop is mapped once (except for setaside being mapped twice, correctly).')
else
    disp('Every PLUM crop is mapped once.')
end

% Check that each LPJGcrop has a corresponding element of PLUMtoLPJG
if length(PLUMtoLPJG) ~= Ncfts_lpjg
    error('length(PLUMtoLPJG) ~= Ncfts_lpjg') ;
end
disp('Each LPJG crop has a corresponding member of PLUMtoLPJG.')


%% Do it

for d = 1:length(dirList)
    inDir = find_PLUM2LPJG_inputs(dirList{d}) ;
    disp(inDir)
    
    %%%%%%%%%%%%%
    %%% Setup %%%
    %%%%%%%%%%%%%
    
    % Get directories
    outDir = addslashifneeded([removeslashifneeded(inDir) '.forLPJG.MATLAB.' thisVer]) ;
    outDir = strrep(outDir,' ','\ ') ;
    unix(['mkdir -p ' outDir]) ;
    unix(['cp ' thisScript ' ' outDir]) ;
        
    inDir_contents = dir(inDir) ;
    inDir_contents = {inDir_contents([inDir_contents.isdir]).name} ;
    yearList = inDir_contents(~cellfun(@isempty, regexp(inDir_contents, '\d\d\d\d'))) ;
    yearList = cellfun(@str2num, yearList) ;
    Nyears = length(yearList) ;
    y1 = yearList(1) ;
    if y1_pre<y1
        yearList_xtra = y1_pre:(y1-1) ;
    else
        yearList_xtra = [] ;
    end
    Nyears_xtra = length(yearList_xtra) ;
    
    %%%%%%%%%%%%%%
    %%% Import %%%
    %%%%%%%%%%%%%%
    
    for y = 1:Nyears
        thisYear = yearList(y) ;
        disp(['Reading ' num2str(thisYear) '...'])
        
        landcover_in_tmp = lpjgu_matlab_readTable([inDir num2str(thisYear) ...
            '/LandCoverFract.txt'], ...
            'dont_save_MAT', true, 'do_save_MAT', false, ...
            'verboseIfNoMat',false) ;
        details_in_tmp = lpjgu_matlab_readTable([inDir num2str(thisYear) ...
            '/LandUse.txt'], ...
            'dont_save_MAT', true, 'do_save_MAT', false, ...
            'verboseIfNoMat',false) ;
        cropland_tmp = landcover_in_tmp.CROPLAND ;
        cropfracs2_in_tmp = lpjgu_matlab_readTable([inDir num2str(thisYear) ...
            '/CropFract.txt'], ...
            'dont_save_MAT', true, 'do_save_MAT', false, ...
            'verboseIfNoMat',false) ;
        
        if y==1
            % Get variable names, part 1
            landcover_cols = landcover_in_tmp.Properties.VariableNames ;
            details_cols = details_in_tmp.Properties.VariableNames ;
            
            % Deal with cropfracs2 variables
            cropfracs2_cols = cropfracs2_in_tmp.Properties.VariableNames ;
            cftList_cropfracs2 = strrep(cropfracs2_cols(4:end),'irr','') ;
            Ncfts_cropfracs2 = length(cftList_cropfracs2) ;
            
            % Get columns indices from detailed table
            inds_PLUMcrops_cropfracs = find(...
                contains(details_cols,'_A') ...
                & ~contains(details_cols,'ruminants') ...
                & ~contains(details_cols,'monogastrics') ...
                & ~contains(details_cols,'pasture') ...
                ) ;
            inds_PLUMcrops_nfert = find(...
                contains(details_cols,'_FQ') ...
                & ~contains(details_cols,'ruminants') ...
                & ~contains(details_cols,'monogastrics') ...
                & ~contains(details_cols,'pasture') ...
                ) ;
            inds_PLUMcrops_irrig = find(...
                contains(details_cols,'_II') ...
                & ~contains(details_cols,'ruminants') ...
                & ~contains(details_cols,'monogastrics') ...
                & ~contains(details_cols,'pasture') ...
                ) ;
            inds_PLUMcrops_yield = find(...
                contains(details_cols,'_Y') ...
                & ~contains(details_cols,'ruminants') ...
                & ~contains(details_cols,'monogastrics') ...
                & ~contains(details_cols,'pasture') ...
                ) ;
            
            % Get temporary tables from detailed table
            cropfracs_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_cropfracs]) ;
            if save_nfert
                nfert_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_nfert]) ;
            end
            if save_irr
                irrig_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_irrig]) ;
            end
            if save_yield
                yield_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_yield]) ;
            end
            
            % Get variable names, part 2
            cropfracs_cols = cropfracs_in_tmp.Properties.VariableNames ;
            cropfracs_cols2 = cropfracs_cols(3:end) ;
            if save_nfert
                nfert_cols = nfert_in_tmp.Properties.VariableNames ;
                nfert_cols2 = nfert_cols(3:end) ;
            end
            if save_irr
                irrig_cols = irrig_in_tmp.Properties.VariableNames ;
                irrig_cols2 = irrig_cols(3:end) ;
            end
            if save_yield
                yield_cols = yield_in_tmp.Properties.VariableNames ;
                yield_cols2 = yield_cols(3:end) ;
            end
            
            % Convert crop names
            for c = 1:Ncfts_plum
                thisCrop_plum = PLUMcrops{c} ;
                cropfracs_cols(strcmp(cropfracs_cols,[thisCrop_plum '_A'])) = {[thisCrop_plum 'i']} ;
                cropfracs_cols2(strcmp(cropfracs_cols2,[thisCrop_plum '_A'])) = {[thisCrop_plum 'i']} ;
                if save_nfert
                    nfert_cols(strcmp(nfert_cols,[thisCrop_plum '_FQ'])) = {[thisCrop_plum 'i']} ;
                    nfert_cols2(strcmp(nfert_cols2,[thisCrop_plum '_FQ'])) = {[thisCrop_plum 'i']} ;
                end
                if save_irr
                    irrig_cols(strcmp(irrig_cols,[thisCrop_plum '_II'])) = {[thisCrop_plum 'i']} ;
                    irrig_cols2(strcmp(irrig_cols2,[thisCrop_plum '_II'])) = {[thisCrop_plum 'i']} ;
                end
                if save_yield
                    yield_cols(strcmp(yield_cols,[thisCrop_plum '_Y'])) = {[thisCrop_plum 'i']} ;
                    yield_cols2(strcmp(yield_cols2,[thisCrop_plum '_Y'])) = {[thisCrop_plum 'i']} ;
                end
            end
            
            % Remove extraneous crop types
            PLUMcrops_toRemove = [] ;
            for c = 1:length(PLUMcrops)
                thisCrop = PLUMcrops{c} ;
                if ~any(strcmp(cropfracs_cols2, thisCrop)) ...
                && ~any(strcmp(cropfracs_cols2, [thisCrop 'i']))
                    warning('%s not found in PLUMcrops; removing.', thisCrop)
                    PLUMcrops_toRemove = [PLUMcrops_toRemove c] ;
                    found = false ;
                    for cc = 1:length(PLUMtoLPJG)
                        thisCell = PLUMtoLPJG{cc} ;
                        for ccc = 1:length(thisCell)
                            thisCellCrop = thisCell{ccc} ;
                            fprintf('%s = %s?\n', thisCellCrop, thisCrop)
                            if strcmp(thisCellCrop, thisCrop)
                                found = true ;
                                thisCell(ccc) = [] ;
                                if isempty(thisCell)
                                    PLUMtoLPJG(cc) = [] ;
                                    LPJGcrops(cc) = [] ;
                                else
                                    PLUMtoLPJG{cc} = thisCell ;
                                end
                                break
                            end
                        end
                        if found; break; end
                    end
                    if ~found; error('???'); end
                    LPJGcrops_mid = [strcat(LPJGcrops,'i') ; LPJGcrops] ;
                    Ncfts_lpjg = length(LPJGcrops) ;
                end
            end
            if ~isempty(PLUMcrops_toRemove)
                PLUMcrops(PLUMcrops_toRemove) = [] ;
                Ncfts_plum = length(PLUMcrops) ;
            end
            
            % Get # variables
            Nvar_landcover = length(landcover_cols) - 3 ;
            Nvar_cropfracs = length(cropfracs_cols) - 2 ;
            if save_nfert
                Nvar_nfert = length(nfert_cols) - 2 ;
            end
            if save_irr
                Nvar_irrig = length(irrig_cols) - 2 ;
            end
            if save_yield
                Nvar_yield = length(yield_cols) - 2 ;
            end
            
            % Get lat/lons
            lons = landcover_in_tmp.Lon ;
            lats = landcover_in_tmp.Lat ;
            Ncells = length(lons) ;
            lonlats_randomized = false ;
            
        else
            % Get temporary tables from detailed table
            cropfracs_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_cropfracs]) ;
            if save_nfert
                nfert_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_nfert]) ;
            end
            if save_irr
                irrig_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_irrig]) ;
            end
            if save_yield
                yield_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_yield]) ;
            end
        end
        
        % Make sure that lon/lat columns of files are identical
        if ~(isequal(cropfracs_in_tmp.Lon,cropfracs2_in_tmp.Lon) && isequal(cropfracs_in_tmp.Lat,cropfracs2_in_tmp.Lat))
            error('Lon/Lat columns of cropfracs_in_tmp and cropfracs2_in_tmp not identical!')
        end
        
        % (Find and) remove rows that sum to zero
        if y==1
            i_plum = 1:Ncells ;
            remove_these_cells = find(sum(table2array(landcover_in_tmp(:,4:end)),2)==0) ;
            lons(remove_these_cells) = [] ;
            lats(remove_these_cells) = [] ;
            Ncells = length(lons) ;
            i_plum(remove_these_cells) = [] ;
        end
        landcover_in_tmp = landcover_in_tmp(i_plum,:) ;
        details_in_tmp = details_in_tmp(i_plum,:) ;
        cropland_tmp = cropland_tmp(i_plum) ;
        cropfracs2_in_tmp = cropfracs2_in_tmp(i_plum,:) ;
        cropfracs_in_tmp = cropfracs_in_tmp(i_plum,:) ;
        if save_nfert
            nfert_in_tmp = nfert_in_tmp(i_plum,:) ;
        end
        if save_irr
            irrig_in_tmp = irrig_in_tmp(i_plum,:) ;
        end
        if save_yield
            yield_in_tmp = yield_in_tmp(i_plum,:) ;
        end

        if any(sum(table2array(landcover_in_tmp(:,4:end)),2)==0)
            error('At least one row in landcover_in_tmp sums to zero!')
        end
        
        % If cropfracs sums to zero, make sure crop area is set to zero.
        % No need to do this with cropfracs2, as it does not have this issue.
        cropfracs_in_sum = sum(table2array(cropfracs_in_tmp(:,3:end)),2) ;
        if any(cropfracs_in_sum==0 & cropland_tmp>0)
            arebad = find(cropfracs_in_sum==0 & cropland_tmp>0) ;
            nbad = length(arebad) ;
            max_cropFrac_inBad = max(cropland_tmp(arebad)) ;
            warning(['cropfracs_in_sum (i.e., sum of all crops_A in detailed table) is zero in ' num2str(nbad) ' row(s) that contain cropland according to LandCoverFract.txt (less-detailed file). Setting to zero. (Max crop frac = ' num2str(max_cropFrac_inBad) ')'])
            if max_cropFrac_inBad>1e-4
                error('max_cropFrac_inBad is too large for me to be comfortable with this!')
            end
            landcover_in_tmp.CROPLAND(arebad) = 0 ;
            cropland_tmp = landcover_in_tmp.CROPLAND ;
            clear arebad nbad max_cropFrac_inBad
        end
        
        % If cropfracs2 sums to zero, make sure crop area is set to zero.
        cropfracs2_in_sum = sum(table2array(cropfracs2_in_tmp(:,3:end)),2) ;
        if any(cropfracs2_in_sum==0 & cropland_tmp>0)
            arebad = find(cropfracs2_in_sum==0 & cropland_tmp>0) ;
            nbad = length(arebad) ;
            max_cropFrac_inBad = max(cropland_tmp(arebad)) ;
            warning(['cropfracs2_in_sum (i.e., sum of all crops in CropFract.txt) is zero in ' num2str(nbad) ' row(s) that contain cropland according to LandCoverFract.txt (less-detailed file). Setting to zero. (Max crop frac = ' num2str(max_cropFrac_inBad) ')'])
            if max_cropFrac_inBad>1e-4
                error('max_cropFrac_inBad is too large for me to be comfortable with this!')
            end
            landcover_in_tmp.CROPLAND(arebad) = 0 ;
            cropland_tmp = landcover_in_tmp.CROPLAND ;
            clear arebad nbad max_cropFrac_inBad
        end
        
        
        % Make sure lat/lon cols are equal
        if ~isequal(lons,landcover_in_tmp.Lon) ...
                || ~isequal(lons,cropfracs_in_tmp.Lon) ...
                || ~isequal(lons,cropfracs2_in_tmp.Lon) ...
                || (save_nfert && ~isequal(lons,nfert_in_tmp.Lon)) ...
                || (save_irr && ~isequal(lons,irrig_in_tmp.Lon)) ...
                || (save_yield && ~isequal(lons,yield_in_tmp.Lon)) ...
                || (save_yield && ~isequal(lons,yield_in_tmp.Lon))
            error('Lons do not match!')
        end
        if ~isequal(lats,landcover_in_tmp.Lat) ...
                || ~isequal(lats,cropfracs_in_tmp.Lat) ...
                || ~isequal(lats,cropfracs2_in_tmp.Lat) ...
                || (save_nfert && ~isequal(lats,nfert_in_tmp.Lat)) ...
                || (save_irr && ~isequal(lats,irrig_in_tmp.Lat)) ...
                || (save_yield && ~isequal(lats,yield_in_tmp.Lat))
            error('Lats do not match!')
        end
        
        if y==1
            % Set up empty Xvy arrays
            landcover_in_Xvy = nan(Ncells,Nvar_landcover,Nyears) ;
            cropfracs_in_Xvy = nan(Ncells,Ncfts_lpjg,Nyears) ;
            if save_nfert
                nfert_in_Xvy = nan(Ncells,Ncfts_lpjg,Nyears) ;
            end
            if save_irr
                irrig_in_Xvy = nan(Ncells,Ncfts_lpjg,Nyears) ;
            end
            if save_yield
                yield_in_Xvy = nan(Ncells,Ncfts_lpjg,Nyears) ;
            end
        end
        
        % If doing so, make sure that every grid cell has at least some cropland
        if someofall
            mincropfrac = 10^(-outPrec_LC) ;
            warning('Should no_cropland be calculated as landcover_in_tmp.CROPLAND<mincropfrac instead??')
%             no_cropland = landcover_in_tmp.CROPLAND==0 ;
            no_cropland = landcover_in_tmp.CROPLAND<mincropfrac ;
            involved = landcover_in_tmp.BARREN>=mincropfrac & no_cropland ;
            landcover_in_tmp.BARREN(involved) = landcover_in_tmp.BARREN(involved) - mincropfrac ;
            landcover_in_tmp.CROPLAND(involved) = landcover_in_tmp.CROPLAND(involved) + mincropfrac ;
            no_cropland = landcover_in_tmp.CROPLAND==0 ;
            if any(no_cropland)
                involved = landcover_in_tmp.URBAN>=mincropfrac & no_cropland ;
                landcover_in_tmp.URBAN(involved) = landcover_in_tmp.URBAN(involved) - mincropfrac ;
                landcover_in_tmp.CROPLAND(involved) = landcover_in_tmp.CROPLAND(involved) + mincropfrac ;
                no_cropland = landcover_in_tmp.CROPLAND==0 ;
                if any(no_cropland)
                    involved = landcover_in_tmp.PASTURE>landcover_in_tmp.NATURAL & no_cropland ;
                    landcover_in_tmp.PASTURE(involved) = landcover_in_tmp.PASTURE(involved) - mincropfrac ;
                    landcover_in_tmp.CROPLAND(involved) = landcover_in_tmp.CROPLAND(involved) + mincropfrac ;
                    no_cropland = landcover_in_tmp.CROPLAND==0 ;
                    if any(no_cropland)
                        involved = landcover_in_tmp.PASTURE<landcover_in_tmp.NATURAL & no_cropland ;
                        landcover_in_tmp.NATURAL(involved) = landcover_in_tmp.NATURAL(involved) - mincropfrac ;
                        landcover_in_tmp.CROPLAND(involved) = landcover_in_tmp.CROPLAND(involved) + mincropfrac ;
                        no_cropland = landcover_in_tmp.CROPLAND==0 ;
                        if any(no_cropland)
                            error('GET CROPLAND FROM SOMEWHERE')
                        end
                    end
                end
            end
        end
        
        % Convert to arrays and set up empty Xv arrays for LPJG crops
        landcover_in_tmp = table2array(landcover_in_tmp(:,4:end)) ;
        cropfracs_in_tmp = table2array(cropfracs_in_tmp(:,3:end)) ;
        cropfracs_in_Xv = nan(Ncells,Ncfts_lpjg) ;
        if save_nfert
            nfert_in_tmp = table2array(nfert_in_tmp(:,3:end)) ;
            nfert_in_Xv = nan(Ncells,Ncfts_lpjg) ;
        end
        if save_irr
            irrig_in_tmp = table2array(irrig_in_tmp(:,3:end)) ;
            irrig_in_Xv = nan(Ncells,Ncfts_lpjg) ;
        end
        if save_yield
            yield_in_tmp = table2array(yield_in_tmp(:,3:end)) ;
            yield_in_Xv = nan(Ncells,Ncfts_lpjg) ;
        end
        
        % Convert PLUM cropfrac to LPJG cropfrac
        for c = 1:Ncfts_lpjg
            thisLPJGcrop = LPJGcrops{c} ;
            thisMap = PLUMtoLPJG{c} ;
            [~,IA] = intersect(cropfracs_cols2,strcat(thisMap,'i')) ;
            thisSum = sum(cropfracs_in_tmp(:,IA),2) ;
            if strcmp(thisLPJGcrop,'CerealsC3w')
                thisSum(thisSum>0 & cropfracs2_in_tmp.TeWWirr==0) = 0 ;
            elseif strcmp(thisLPJGcrop,'CerealsC3s')
                thisSum(thisSum>0 & cropfracs2_in_tmp.TeWWirr>0) = 0 ;
            end
            cropfracs_in_Xv(:,c) = thisSum ;
            clear thisLPJGcrop thisMap thisSum
        end
        
        % If doing so, make sure that every grid cell has at least some of
        % each crop.
        % Make sure landcover and cropfracs sum to 1
        cropfracs_in_sum_whereCrop = sum(cropfracs_in_Xv(cropland_tmp>0,:),2) ;
        if any(abs(cropfracs_in_sum_whereCrop - 1) > 0.01)
            error('Cropfrac difference too big to be comfortable!')
        end
        while any(abs(cropfracs_in_sum_whereCrop - 1) > 10^-(outPrec_LC+2))
            cropfracs_in_Xv(cropland_tmp>0,:) = cropfracs_in_Xv(cropland_tmp>0,:) ./ repmat(cropfracs_in_sum_whereCrop,[1 Ncfts_lpjg]) ;
            if any(isnan(cropfracs_in_Xv))
                error('While ensuring sum to 1 of cropfracs_Xv, NaN created!')
            end
            cropfracs_in_sum_whereCrop = sum(cropfracs_in_Xv(cropland_tmp>0,:),2) ;
            if any(abs(cropfracs_in_sum_whereCrop - 1) > 0.01)
                error('Cropfrac difference too big to be comfortable!')
            end
        end
        if someofall
            mincropfrac = 10^(-outPrec_LC) ;
            % First fill in all-zero rows
            nocropfracs = sum(cropfracs_in_Xv,2) ;
            cropfracs_in_Xv(repmat(nocropfracs,[1 Ncfts_lpjg])==0) = 1/Ncfts_lpjg ;
            % Now fill in remaining
            for c = 1:Ncfts_lpjg
                % Find rows with 0 for this crop
                thiscrop = cropfracs_in_Xv(:,c) ;
                iszerothiscrop = thiscrop==0 ;
                % Find the crop that currently has the greatest area
                maxcropfrac_Xv = repmat(max(cropfracs_in_Xv,[],2),[1 Ncfts_lpjg]) ;
                ismaxcropfrac_Xv = cropfracs_in_Xv==maxcropfrac_Xv ;
                % Make sure there's only one ismaxcropfrac in each row
                for i = fliplr(2:Ncfts_lpjg)
                    tmp = ismaxcropfrac_Xv(:,i) ;
                    sumtoleft = sum(ismaxcropfrac_Xv(:,1:(i-1)),2) ;
                    tmp(sumtoleft>0) = false ;
                    ismaxcropfrac_Xv(:,i) = tmp ;
                end
                cropfracs_in_Xv(ismaxcropfrac_Xv & iszerothiscrop) = cropfracs_in_Xv(ismaxcropfrac_Xv & iszerothiscrop) - mincropfrac ;
                cropfracs_in_Xv(iszerothiscrop,c) = cropfracs_in_Xv(iszerothiscrop,c) + mincropfrac ;
            end
        end
        
        % Convert PLUM irrigation, fertilizer, and yield to LPJG crop
        % types (unit conversion happens later)
        for cL = 1:Ncfts_lpjg
            thisMap = PLUMtoLPJG{cL} ;
            if save_nfert
                nfert_in_Xv(:,cL) = weighted_average_for_PLUM2LPJG(thisMap,nfert_cols2,nfert_in_tmp,cropfracs_cols2,cropfracs_in_tmp) ;            
            end
            if save_irr
                irrig_in_Xv(:,cL) = weighted_average_for_PLUM2LPJG(thisMap,irrig_cols2,irrig_in_tmp,cropfracs_cols2,cropfracs_in_tmp) ;
            end
            if save_yield
                thisYield = weighted_average_for_PLUM2LPJG(thisMap,yield_cols2,yield_in_tmp,cropfracs_cols2,cropfracs_in_tmp) ;
                % Setaside should be handled gracefully
                % Special operation if specifying winter/spring wheat
                thisLPJGcrop = LPJGcrops{cL} ;
                if strcmp(thisLPJGcrop,'CerealsC3w')
                    thisYield(thisYield>0 & cropfracs2_in_tmp.TeWWirr==0) = 0 ;
                elseif strcmp(thisLPJGcrop,'CerealsC3s')
                    thisYield(thisYield>0 & cropfracs2_in_tmp.TeWWirr>0) = 0 ;
                end
                yield_in_Xv(:,cL) = thisYield ;
                clear thisLPJGcrop thisYield
            end
            clear thisMap
        end
        clear cropfracs_in_tmp 
        clear nfert_in_tmp irrig_in_tmp
        
        % Check for NaNs
        if any(isnan(landcover_in_tmp(:)))
            error('Some value of landcover_in_tmp is NaN (before ensuring sum to 1).')
        end
        if any(isnan(cropfracs_in_Xv(:)))
            error('Some value of cropfracs_in_Xv is NaN (before ensuring sum to 1).')
        end
        if save_nfert && any(isnan(nfert_in_Xv(:)))
            error('Some value of nfert_in_Xv is NaN.')
        end
        if save_irr && any(isnan(irrig_in_Xv(:)))
            error('Some value of irrig_in_Xv is NaN.')
        end
        if save_yield && any(isnan(yield_in_tmp(:)))
            error('Some value of yield_in_tmp is NaN.')
        end
        
        % Make sure there is no grid cell with both CerealsC3w and CerealsC3s
        if ~isempty(strcmp(LPJGcrops,'CerealsC3w')) && ~isempty(strcmp(LPJGcrops,'CerealsC3s'))
            if any(cropfracs_in_Xv(:,strcmp(LPJGcrops,'CerealsC3w'))>0 & cropfracs_in_Xv(:,strcmp(LPJGcrops,'CerealsC3s'))>0)
                error('At least one cell has both CerealsC3w and CerealsC3s!')
            end
        end
        
        % Make sure landcover and cropfracs sum to 1
        % (Previously checked to avoid rows that sum to 0)
        while any(abs(sum(landcover_in_tmp,2) - 1) > 10^-(outPrec_LC+2))
            if y==1
                landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover]) ;
            else
                % Because got rid of URBAN after doing this on y==1
                landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover+1]) ;
            end
            if any(isnan(landcover_in_tmp))
                error('While ensuring sum to 1 of landcover_in_tmp, NaN created!')
            end
        end
        cropfracs_in_sum_whereCrop = sum(cropfracs_in_Xv(cropland_tmp>0,:),2) ;
        while any(abs(cropfracs_in_sum_whereCrop - 1) > 10^-(outPrec_LC+2))
            cropfracs_in_Xv(cropland_tmp>0,:) = cropfracs_in_Xv(cropland_tmp>0,:) ./ repmat(cropfracs_in_sum_whereCrop,[1 Ncfts_lpjg]) ;
            if any(abs(cropfracs_in_sum_whereCrop - 1) > 0.01)
                error('Cropfrac difference too big to be comfortable!')
            end
            if any(isnan(cropfracs_in_Xv))
                error('While ensuring sum to 1 of cropfracs_Xv, NaN created!')
            end
            cropfracs_in_sum_whereCrop = sum(cropfracs_in_Xv(cropland_tmp>0,:),2) ;
        end
        
        % Deal with NATURAL/URBAN/BARREN
        if y==1
            plum_list2map = nan(Ncells,1) ;
            for c = 1:Ncells
                thisLon = lons(c) ;
                thisLat = lats(c) ;
                plum_list2map(c) = find(lons_map==thisLon & lats_map==thisLat) ;
            end
            barren_lpj = barren_lpj_map(plum_list2map) ;
            vegd_lpj = vegd_lpj_map(plum_list2map) ;
        end
        barren_old = landcover_in_tmp(:,find(strcmp(landcover_cols,'BARREN'))-3) ;
        natural_old = landcover_in_tmp(:,find(strcmp(landcover_cols,'NATURAL'))-3) ;
        urban_old = landcover_in_tmp(:,find(strcmp(landcover_cols,'URBAN'))-3) ;
        vegd_old = sum(landcover_in_tmp(:,find(strcmp(landcover_cols,'NATURAL')|strcmp(landcover_cols,'CROPLAND')|strcmp(landcover_cols,'PASTURE'))-3),2) ;
        mgmt_old = vegd_old - natural_old ;
        [barren_new, natural_new, ~] = ...
            BarNatUrb(barren_old, mgmt_old, natural_old, urban_old, ...
                      barren_lpj, vegd_lpj, outPrec_LC) ;
        
        % Save
        landcover_in_tmp(:,find(strcmp(landcover_cols,'BARREN'))-3) = barren_new ;
        landcover_in_tmp(:,find(strcmp(landcover_cols,'NATURAL'))-3) = natural_new ;
        landcover_in_tmp(:,find(strcmp(landcover_cols,'URBAN'))-3) = [] ;
        if y==1
            % Because we got rid of URBAN
            landcover_in_Xvy(:,end,:) = [] ;
            Nvar_landcover = Nvar_landcover - 1 ;
        end
        
        % Again make sure landcover and cropfracs sum to 1
        % (Previously checked to avoid rows that sum to 0)
        while any(abs(sum(landcover_in_tmp,2) - 1) > 10^-(outPrec_LC+2))
            landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover]) ;
            if any(isnan(landcover_in_tmp))
                error('While ensuring sum to 1 of landcover_in_tmp, NaN created!')
            end
        end
                
        % Save to array
        landcover_in_Xvy(:,:,y) = landcover_in_tmp ;
        cropfracs_in_Xvy(:,:,y) = cropfracs_in_Xv ;
        if save_nfert
            nfert_in_Xvy(:,:,y) = nfert_in_Xv ;
        end
        if save_irr
            irrig_in_Xvy(:,:,y) = irrig_in_Xv ;
        end
        if save_yield
            yield_in_Xvy(:,:,y) = yield_in_Xv ;
        end
        
        clear *tmp *_Xv
    end
    
    % Check for NaNs
    if any(isnan(landcover_in_Xvy(:)))
        error('Some value of landcover_in_Xvy is NaN.')
    elseif any(isnan(cropfracs_in_Xvy(:)))
        error('Some value of cropfracs_in_Xvy is NaN.')
    elseif save_nfert && any(isnan(nfert_in_Xvy(:)))
        error('Some value of nfert_in_Xvy is NaN.')
    elseif save_irr && any(isnan(irrig_in_Xvy(:)))
        error('Some value of irrig_in_Xvy is NaN.')
    elseif save_yield && any(isnan(yield_in_Xvy(:)))
        error('Some value of yield_in_Xvy is NaN.')
    end
    
    % Check for nonsensical values
    if save_nfert && min(nfert_in_Xvy(:))<0
        warning('min(nfert_in_Xvy)<0')
    elseif save_irr && min(irrig_in_Xvy(:))<0
        warning('min(irrig_in_Xvy)<0')
    elseif save_yield && min(yield_in_Xvy(:))<0
        warning('min(yield_in_Xvy)<0')
    elseif save_irr && max(irrig_in_Xvy(:))>1
        warning('max(irrig_in_Xvy)>1')
    end
    
    % Change lon/lat from lower-left to center
    lons = lons + 0.25 ;
    lats = lats + 0.25 ;
    
    % Convert ktN/ha to kgN/m2
    if save_nfert
        nfert_in_Xvy = nfert_in_Xvy*cf_ktNha_kgNm2 ;
    end
    if save_yield
        yield_in_Xvy = yield_in_Xvy*cf_tPha_kgPm2 ;
    end
    
    % Add rainfed crops (zeros)
    cropfracs_in_Xvy = cat(2,cropfracs_in_Xvy,zeros(Ncells,Ncfts_lpjg,Nyears)) ;
    if save_nfert
        nfert_in_Xvy = cat(2,nfert_in_Xvy,zeros(Ncells,Ncfts_lpjg,Nyears)) ;
    end
    if save_irr
        irrig_in_Xvy = cat(2,irrig_in_Xvy,zeros(Ncells,Ncfts_lpjg,Nyears)) ;
    end
    
    % Swap ExtraCropi to rainfed and remove, because this has no 
    % irrigated area. If in the future PLUM (or something else) specifies
    % irrigation for this, remove this part.
    cropfracs_in_Xvy(:,strcmp(LPJGcrops_mid,'ExtraCrop'),:) = ...
        cropfracs_in_Xvy(:,strcmp(LPJGcrops_mid,'ExtraCropi'),:) ;
    cropfracs_in_Xvy(:,contains(LPJGcrops_mid,{'ExtraCropi'}),:) = [] ;
    if save_nfert
        nfert_in_Xvy(:,strcmp(LPJGcrops_mid,'ExtraCrop'),:) = ...
            nfert_in_Xvy(:,strcmp(LPJGcrops_mid,'ExtraCropi'),:) ;
        nfert_in_Xvy(:,contains(LPJGcrops_mid,{'ExtraCropi'}),:) = [] ;
    end
    if save_irr
        irrig_in_Xvy(:,strcmp(LPJGcrops_mid,'ExtraCrop'),:) = ...
            irrig_in_Xvy(:,strcmp(LPJGcrops_mid,'ExtraCropi'),:) ;
        irrig_in_Xvy(:,contains(LPJGcrops_mid,{'ExtraCropi'}),:) = [] ;
    end
    LPJGcrops_out = LPJGcrops_mid ;
    LPJGcrops_out(contains(LPJGcrops_out,{'ExtraCropi'})) = [] ;
    Ncfts_lpjg_out = length(LPJGcrops_out) ;
    
    % Add extra years, if doing so
    if Nyears_xtra>0
        landcover_in_Xvy = cat(3,repmat(landcover_in_Xvy(:,:,1),[1 1 Nyears_xtra]),landcover_in_Xvy) ;
        cropfracs_in_Xvy = cat(3,repmat(cropfracs_in_Xvy(:,:,1),[1 1 Nyears_xtra]),cropfracs_in_Xvy) ;
        if save_nfert
            nfert_in_Xvy = cat(3,repmat(nfert_in_Xvy(:,:,1),[1 1 Nyears_xtra]),nfert_in_Xvy) ;
        end
        if save_irr
            irrig_in_Xvy = cat(3,repmat(irrig_in_Xvy(:,:,1),[1 1 Nyears_xtra]),irrig_in_Xvy) ;
        end
        if save_yield
            yield_in_Xvy = cat(3,repmat(yield_in_Xvy(:,:,1),[1 1 Nyears_xtra]),yield_in_Xvy) ;
        end
    end
    
    disp('Done.')
    
    % Exclude cells that aren't in LPJ-GUESS output, if doing so
    [~,not_in_lpjg] = setdiff([lons lats],[lons_lpjg lats_lpjg],'rows') ;
    if exclude_nonlpjg_cells
        lons(not_in_lpjg) = [] ;
        lats(not_in_lpjg) = [] ;
        Ncells = length(lons) ;
    end
             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write out: Gridlist %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if save_gl
        % Find cells not in LPJ-GUESS output
        % (including these in other outputs because they might be needed at some
        % point?)
        lons_4gl = lons ;
        lats_4gl = lats ;
        if ~exclude_nonlpjg_cells
            lons_4gl(not_in_lpjg) = [] ;
            lats_4gl(not_in_lpjg) = [] ;
            Ncells_4gl = length(lons_4gl) ;
        else
            Ncells_4gl = Ncells ;
        end
        
        % Get random order for output
        rng(20171116) ;
        rdmsam = randsample(Ncells_4gl,Ncells_4gl) ;
        
        % Save gridlist
        outFile_gridlist = [outDir 'PLUMout_gridlist.txt'] ;
        if exist(strrep(outFile_gridlist,'\ ',' '),'file') && ~force_overwrite
            error(['Outfile already exists! (' strrep(outFile_gridlist,'\ ',' ') ')']) ;
        end
        out_formatSpec_gridlist = '%4.2f %4.2f\n' ;
        fid1_gridlist = fopen(strrep(outFile_gridlist,'\ ',' '),'w') ;
        fprintf(fid1_gridlist,out_formatSpec_gridlist,[lons_4gl(rdmsam) lats_4gl(rdmsam)]') ;
        fclose(fid1_gridlist) ;
        fid1_gridlist = fopen(strrep(outFile_gridlist,'\ ',' '),'a+') ;
        fprintf(fid1_gridlist,'%s','') ;
        fclose(fid1_gridlist) ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write out: Landcover %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if save_lu
        % Get header and formatSpec
        out_header_landcover = 'Lon Lat Year' ;
        for c = 1:Nvar_landcover
            out_header_landcover = [out_header_landcover ' ' landcover_cols{c+3}] ;
        end
        out_formatSpec_landcover = ['%4.2f %4.2f %4.0f' repmat([' %5.' num2str(outPrec_LC) 'f'],[1 Nvar_landcover])] ;
        out_formatSpec_landcover = [out_formatSpec_landcover '\n'] ;

        % Set up output files
        outFile_landcover = [outDir 'landcover.txt'] ;
        if exist(outFile_landcover,'file') && ~force_overwrite
            error(['Outfile already exists! (' outFile_landcover ')']) ;
        end
        fid1_landcover = fopen(strrep(outFile_landcover,'\ ',' '), 'w') ;
        fprintf(fid1_landcover,'%s \n',out_header_landcover);

        disp('Saving landcover...')
        for c = 1:Ncells
            thisLon = lons(c) ;
            thisLat = lats(c) ;

            landcover_out = [repmat(thisLon,[Nyears+Nyears_xtra 1]) repmat(thisLat,[Nyears+Nyears_xtra 1]) [yearList_xtra yearList]' squeeze(landcover_in_Xvy(c,:,:))'] ;

            if rem(c-1,fclose_every)==0
                fid1_landcover = fopen(strrep(outFile_landcover,'\ ',' '), 'a+') ;
            end

            fprintf(fid1_landcover,out_formatSpec_landcover,landcover_out') ;

            if rem(c,fclose_every)==0 || c==Ncells
                fclose(fid1_landcover) ;
            end

            if rem(c,pct_progress*floor(Ncells/100))==0
                compltn = 100*(c)/Ncells ;
                if compltn>0
                    disp(['   ' num2str(round(compltn)) '% complete...'])
                end
            end
        end
        disp('gzipping...')
        if exist([strrep(outFile_landcover,'\ ',' ') '.gz'],'file')
            unix(['rm ' outFile_landcover '.gz']) ;
        end
        unix(['gzip ' outFile_landcover]) ;
        disp('Done.')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write out: Cropfracs %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if save_cropfracs
        % Get header and formatSpec
        out_header_cropfracs = 'Lon Lat Year' ;
        for c = 1:Ncfts_lpjg_out
            out_header_cropfracs = [out_header_cropfracs ' ' LPJGcrops_out{c}] ;
        end
        out_formatSpec_cropfracs = ['%4.2f %4.2f %4.0f' repmat([' %5.' num2str(outPrec_LC) 'f'],[1 Ncfts_lpjg_out])] ;
        out_formatSpec_cropfracs = [out_formatSpec_cropfracs '\n'] ;

        % Set up output files
        outFile_cropfracs = [outDir 'cropfractions.txt'] ;
        if exist(outFile_cropfracs,'file') && ~force_overwrite
            error(['Outfile already exists! (' outFile_cropfracs ')']) ;
        end
        fid1_cropfracs = fopen(strrep(outFile_cropfracs,'\ ',' '), 'w') ;
        fprintf(fid1_cropfracs,'%s \n',out_header_cropfracs);

        disp('Saving cropfracs...')
        for c = 1:Ncells
            thisLon = lons(c) ;
            thisLat = lats(c) ;

            cropfracs_out = [repmat(thisLon,[Nyears+Nyears_xtra 1]) repmat(thisLat,[Nyears+Nyears_xtra 1]) [yearList_xtra yearList]' squeeze(cropfracs_in_Xvy(c,:,:))'] ;

            if rem(c-1,fclose_every)==0
                fid1_cropfracs = fopen(strrep(outFile_cropfracs,'\ ',' '), 'a+') ;
            end

            fprintf(fid1_cropfracs,out_formatSpec_cropfracs,cropfracs_out') ;

            if rem(c,fclose_every)==0 || c==Ncells
                fclose(fid1_cropfracs) ;
            end

            if rem(c,pct_progress*floor(Ncells/100))==0
                compltn = 100*(c)/Ncells ;
                if compltn>0
                    disp(['   ' num2str(round(compltn)) '% complete...'])
                end
            end
        end
        disp('gzipping...')
        if exist([strrep(outFile_cropfracs,'\ ',' ') '.gz'],'file')
            unix(['rm ' outFile_cropfracs '.gz']) ;
        end
        unix(['gzip ' outFile_cropfracs]) ;
        disp('Done.')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write out: Nfert %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    if save_nfert
        % Get header and formatSpec
        out_header_nfert = 'Lon Lat Year' ;
        for c = 1:Ncfts_lpjg_out
            out_header_nfert = [out_header_nfert ' ' LPJGcrops_out{c}] ;
        end
        out_formatSpec_nfert = ['%4.2f %4.2f %4.0f' repmat([' %5.' num2str(outPrec_mgmtInputs) 'f'],[1 Ncfts_lpjg_out])] ;
        out_formatSpec_nfert = [out_formatSpec_nfert '\n'] ;

        % Set up output files
        outFile_nfert = [outDir 'nfert.txt'] ;
        if exist(strrep(outFile_nfert,'\ ',' '),'file') && ~force_overwrite
            error(['Outfile already exists! (' outFile_nfert ')']) ;
        end
        fid1_nfert = fopen(strrep(outFile_nfert,'\ ',' '), 'w') ;
        fprintf(fid1_nfert,'%s \n',out_header_nfert);

        disp('Saving Nfert...')
        for c = 1:Ncells
            thisLon = lons(c) ;
            thisLat = lats(c) ;

            nfert_out = [repmat(thisLon,[Nyears+Nyears_xtra 1]) repmat(thisLat,[Nyears+Nyears_xtra 1]) [yearList_xtra yearList]' squeeze(nfert_in_Xvy(c,:,:))'] ;

            if rem(c-1,fclose_every)==0
                fid1_nfert = fopen(strrep(outFile_nfert,'\ ',' '), 'a+') ;
            end

            fprintf(fid1_nfert,out_formatSpec_nfert,nfert_out') ;

            if rem(c,fclose_every)==0 || c==Ncells
                fclose(fid1_nfert) ;
            end

            if rem(c,pct_progress*floor(Ncells/100))==0
                compltn = 100*(c)/Ncells ;
                if compltn>0
                    disp(['   ' num2str(round(compltn)) '% complete...'])
                end
            end
        end
        disp('gzipping...')
        if exist([strrep(outFile_nfert,'\ ',' ') '.gz'],'file')
            unix(['rm ' outFile_nfert '.gz']) ;
        end
        unix(['gzip ' outFile_nfert]) ;
        disp('Done.')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write out: Irrigation %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if save_irr
        % Get header and formatSpec
        out_header_irrig = 'Lon Lat Year' ;
        for c = 1:Ncfts_lpjg_out
            out_header_irrig = [out_header_irrig ' ' LPJGcrops_out{c}] ;
        end
        out_formatSpec_irrig = ['%4.2f %4.2f %4.0f' repmat([' %1.' num2str(outPrec_mgmtInputs) 'f'],[1 Ncfts_lpjg_out])] ;
        out_formatSpec_irrig = [out_formatSpec_irrig '\n'] ;
        
        % Set up output files
        outFile_irrig = [outDir 'irrig.txt'] ;
        if exist(strrep(outFile_irrig,'\ ',' '),'file') && ~force_overwrite
            error(['Outfile already exists! (' outFile_irrig ')']) ;
        end
        fid1_irrig = fopen(strrep(outFile_irrig,'\ ',' '), 'w') ;
        fprintf(fid1_irrig,'%s \n',out_header_irrig);
        
        disp('Saving irrig...')
        for c = 1:Ncells
            thisLon = lons(c) ;
            thisLat = lats(c) ;
            
            irrig_out = [repmat(thisLon,[Nyears+Nyears_xtra 1]) repmat(thisLat,[Nyears+Nyears_xtra 1]) [yearList_xtra yearList]' squeeze(irrig_in_Xvy(c,:,:))'] ;
            
            if rem(c-1,fclose_every)==0
                fid1_irrig = fopen(strrep(outFile_irrig,'\ ',' '), 'a+') ;
            end
            
            fprintf(fid1_irrig,out_formatSpec_irrig,irrig_out') ;
            
            if rem(c,fclose_every)==0 || c==Ncells
                fclose(fid1_irrig) ;
            end
            
            if rem(c,pct_progress*floor(Ncells/100))==0
                compltn = 100*(c)/Ncells ;
                if compltn>0
                    disp(['   ' num2str(round(compltn)) '% complete...'])
                end
            end
        end
        disp('gzipping...')
        if exist([strrep(outFile_irrig,'\ ',' ') '.gz'],'file')
            unix(['rm ' outFile_irrig '.gz']) ;
        end
        unix(['gzip ' outFile_irrig]) ;
        disp('Done.')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write out: Expected yield %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if save_yield
        % Get header and formatSpec
        out_header_yield = 'Lon Lat Year' ;
        % Only include "rainfed" crops
        for c = 1:Ncfts_lpjg
            out_header_yield = [out_header_yield ' ' LPJGcrops{c}] ;
        end
        out_formatSpec_yield = ['%4.2f %4.2f %4.0f' repmat([' %1.' num2str(outPrec_mgmtInputs) 'f'],[1 Ncfts_lpjg])] ;
        out_formatSpec_yield = [out_formatSpec_yield '\n'] ;
        
        % Set up output files
        outFile_yield = [outDir 'yield.txt'] ;
        if exist(strrep(outFile_yield,'\ ',' '),'file') && ~force_overwrite
            error(['Outfile already exists! (' outFile_yield ')']) ;
        end
        fid1_yield = fopen(strrep(outFile_yield,'\ ',' '), 'w') ;
        fprintf(fid1_yield,'%s \n',out_header_yield);
        
        disp('Saving yield...')
        for c = 1:Ncells
            thisLon = lons(c) ;
            thisLat = lats(c) ;
            
            yield_out = [repmat(thisLon,[Nyears+Nyears_xtra 1]) repmat(thisLat,[Nyears+Nyears_xtra 1]) [yearList_xtra yearList]' squeeze(yield_in_Xvy(c,:,:))'] ;
            
            if rem(c-1,fclose_every)==0
                fid1_yield = fopen(strrep(outFile_yield,'\ ',' '), 'a+') ;
            end
            
            fprintf(fid1_yield,out_formatSpec_yield,yield_out') ;
            
            if rem(c,fclose_every)==0 || c==Ncells
                fclose(fid1_yield) ;
            end
            
            if rem(c,pct_progress*floor(Ncells/100))==0
                compltn = 100*(c)/Ncells ;
                if compltn>0
                    disp(['   ' num2str(round(compltn)) '% complete...'])
                end
            end
        end
        disp('gzipping...')
        if exist([strrep(outFile_yield,'\ ',' ') '.gz'],'file')
            unix(['rm ' outFile_yield '.gz']) ;
        end
        unix(['gzip ' outFile_yield]) ;
        disp('Done.')
    end
        
end



