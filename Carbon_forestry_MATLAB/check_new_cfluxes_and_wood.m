%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing outputs after giving PLUM all wood %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/landsymm-dev-forestry/all_wood_to_plum' ;
testDir1 = sprintf('%s/out.before.latest', topDir) ;
testDir2 = sprintf('%s/out.after.latest', topDir) ;

% Set to true if you ever do wood harvest (real or potential) as part of
% forestry rather than just during land-use change.
woodharvest_outside_LUC = false ;

files_identical = {'aaet.out', 'agestruct_dens_forest.out', 'agestruct_dens_natural.out', 'agpp.out', ...
    'annual_burned_area.out', 'anpp.out', 'anpp_forest.out', 'anpp_natural.out', 'anpp_sts.out', ...
    'clitter.out', 'clitter_sts.out', 'cmass.out', 'cmass_forest.out', 'cmass_killed_harv_sts.out', ...
    'cmass_natural.out', 'cmass_sts.out', 'cmass_tree_mort_sts.out', 'cmass_tree_sts.out', ...
    'cmass_wood_harv_sts.out', 'cmass_wood_potharv_sts.out', 'cmass_wood_sts.out', ...
    'csink_sts.out', 'csoil_sts.out', 'cutinterval_sts.out', 'dens.out', ...
    'dens_forest.out', 'dens_natural.out', 'dens_sts.out', 'diam.out', 'diam_forest.out', 'diam_g_sts.out', ...
    'diam_natural.out', 'diamstruct_cmass_forest.out', 'diamstruct_cmass_natural.out', ...
    'diamstruct_dens_forest.out', 'diamstruct_dens_natural.out', 'fpc.out', 'fpc_forest.out', ...
    'fpc_natural.out', 'height.out', 'height_forest.out', 'height_natural.out', 'lai.out', 'lai_forest.out', ...
    'lai_natural.out', 'landsymm_ctree_ag_sts.out', 'landsymm_ctree_sts.out', ...
    'landsymm_plutC_from_crop.out', 'landsymm_plutC_from_past.out', ...
    'landsymm_plutW_from_crop.out', 'landsymm_plutW_from_past.out', ...
    'simfireanalysis.out', ...
    'cmass_wood_harv_leaf_sts.out', 'cmass_wood_harv_twig_sts.out', 'cmass_wood_harv_root_sts.out'} ;

files_cflux_nonsts = {'cflux.out', 'cflux_cropland.out', 'cflux_forest.out', 'cflux_natural.out', 'cflux_pasture.out'} ;
files_cpool_nonsts = {'cpool.out', 'cpool_cropland.out', 'cpool_forest.out', 'cpool_natural.out', 'cpool_pasture.out'} ;

files_zero_in_before = {'cmass_wood_harv_toplum_sts.out', 'cmass_wood_potharv_toplum_sts.out'} ;
files_zero_in_after = {'cflux_slowh_sts.out', 'cmass_wood_harv_toprod_sts.out', 'cmass_wood_potharv_toprod_sts.out', ...
    'cslowprod_sts.out'} ;

files_pool_diff_bc_slowpool_vs_2plum = {'ctotal_sts.out'} ;

files_landsymm_plutC = {'landsymm_plutC_from_forC.out', 'landsymm_plutC_from_ntrl.out'} ;

% files_harvPoolDiff = {'cpool.out', 'ctotal_sts.out'} ;
files_harvPoolDiff = {} ; % Because I'm not actually handling these yet

cropList = {'CerealsC3', 'CerealsC3i', 'CerealsC4', 'CerealsC4i', 'Rice', 'Ricei', 'OilNfix', 'OilNfixi', ...
    'OilOther', 'OilOtheri', 'Pulses', 'Pulsesi', 'StarchyRoots', 'StarchyRootsi', 'FruitAndVeg', ...
    'FruitAndVegi', 'Sugarbeet', 'Sugarbeeti', 'Sugarcane', 'Sugarcanei', 'ExtraCrop'} ;

stem_harvest_slow_frac = 0.33 ;


%% Make sure that both dirs have same output files

[outFiles1, outFilenames1] = get_outFiles(testDir1) ;
[outFiles2, outFilenames2] = get_outFiles(testDir2) ;

diffFilenames_in1not2 = setdiff(outFilenames1, outFilenames2) ;
diffFilenames_in2not1 = setdiff(outFilenames2, outFilenames1) ;

files_ok_in1not2 = {} ;
files_ok_in2not1 = {} ;

ok1 = compare_file_lists(diffFilenames_in1not2, files_ok_in1not2, 1) ;
ok2 = compare_file_lists(diffFilenames_in2not1, files_ok_in2not1, 2) ;
if ~ok1 || ~ok2
    error('File list mismatch')
else
    fprintf('%d files found in each directory; all match\n', length(outFilenames1))
end


%% Finish categorizing files

files_landsymm_plutW = strrep(files_landsymm_plutC, 'plutC', 'plutW') ;
files_landsymm_plut = [files_landsymm_plutC files_landsymm_plutW] ;

% Will be identical if wood harvest only occurs as part of LUC.
files_identical_if_woodharvest_only_in_luc = { ...
    'cflux_landsymm_sts.out', 'cflux_noslowh_sts.out'
    } ;
if ~woodharvest_outside_LUC
    files_identical = [files_identical files_identical_if_woodharvest_only_in_luc] ; 
end

unclassifiedFiles = setdiff(outFilenames1, ...
    [files_identical files_identical_if_woodharvest_only_in_luc files_harvPoolDiff files_cflux_nonsts ...
    files_cpool_nonsts files_zero_in_before files_zero_in_after files_landsymm_plut]) ;
if ~isempty(unclassifiedFiles)
    fprintf('%d unclassified files:\n', length(unclassifiedFiles))
    disp(unclassifiedFiles')
end

files_NOTidentical = setdiff(outFilenames1, files_identical) ;


%% Check that wood harvest components all add up

% Amount transferred from veg to slow-product or to-PLUM pools
%%% Columns = DESTINATION land use
[wood_harv_1, wood_harv_2] = ...
    read_two_tables(testDir1, testDir2, 'cmass_wood_harv_sts.out', cropList) ;
[prod_harv_1, prod_harv_2] = ...
    read_two_tables(testDir1, testDir2, 'cmass_wood_harv_toprod_sts.out', cropList) ;
[plum_harv_1, plum_harv_2] = ...
    read_two_tables(testDir1, testDir2, 'cmass_wood_harv_toplum_sts.out', cropList) ;

% C "harvested" from stems, "leaves attached to twigs," twigs, and coarse roots
%%% Columns = DESTINATION land use
[hwstem1, hwstem2] = ...
    read_two_tables(testDir1, testDir2, 'cmass_wood_harv_sts.out', cropList) ;
[hwleaf1, hwleaf2] = ...
    read_two_tables(testDir1, testDir2, 'cmass_wood_harv_leaf_sts.out', cropList) ;
[hwtwig1, hwtwig2] = ...
    read_two_tables(testDir1, testDir2, 'cmass_wood_harv_twig_sts.out', cropList) ;
[hwroot1, hwroot2] = ...
    read_two_tables(testDir1, testDir2, 'cmass_wood_harv_root_sts.out', cropList) ;

% Test that sums work out
t2a = @(x) table2array(x(:,4:end)) ;
prec = 1e-4 ;
if max(max(abs(t2a(prod_harv_1) - t2a(hwstem1)*stem_harvest_slow_frac))) > prec
    warning('Sums in "before" don''t match expectation')
end
if max(max(abs(t2a(plum_harv_2) - (t2a(hwstem2) + t2a(hwtwig2) + t2a(hwroot2))))) > prec
    warning('Sums in "after" don''t match expectation')
end

t2a(plum_harv_2) ./ t2a(prod_harv_1)

disp('Done checking: Wood harvest components add up.')


%%

[wood_potharv_1, wood_potharv_2] = ...
    read_two_tables(testDir1, testDir2, 'cmass_wood_potharv_sts.out', cropList) ;
[prod_potharv_1, prod_potharv_2] = ...
    read_two_tables(testDir1, testDir2, 'cmass_wood_potharv_toprod_sts.out', cropList) ;
[plum_potharv_1, plum_potharv_2] = ...
    read_two_tables(testDir1, testDir2, 'cmass_wood_potharv_toplum_sts.out', cropList) ;

head(plum_potharv_2)
t2a(plum_potharv_2) ./ t2a(prod_potharv_1)

(1) Why are these not constant?
(2) Why are these not the same for forest and natural?
(3) You need to also provide Bart a potential CUT file! This can probably come from a potharv file.


%% Check that landsymm_plut files match expectations

t2a = @(x) table2array(x(:,4:end)) ;
prec = 3 ; % In code: landsymm_outprec

for f = 1:length(files_landsymm_plutC)
    thisFile_C = files_landsymm_plutC{f} ;
    [plutC_1, plutC_2] = ...
        read_two_tables(testDir1, testDir2, thisFile_C, cropList) ;
    thisFile_W = files_landsymm_plutW{f} ;
    [plutW_1, plutW_2] = ...
        read_two_tables(testDir1, testDir2, thisFile_W, cropList) ;

    thisLU = regexp(thisFile_C,'from_[a-z]+','match') ;
    if length(thisLU) ~= 1
        error('Expected to find one match in file name; found %d', length(thisLU))
    end
    thisLU = strrep(thisLU{1}, 'from_', '') ;

    % Total removed in LU transition should not change
    if max(max(abs(t2a(plutC_1) + t2a(plutW_1) - (t2a(plutC_2) + t2a(plutW_2))))) > prec
        warning('LU transition C+W differ from "before" to "after" for %s', thisLU)
    end

    % ntrl-to-ntrl transition should be all zeros (because no change
    % actually happens)
    if strcmp(thisLU, 'ntrl')
        if any(plutC_1.('to_ntrl') > 0)
            warning('ntrl-to-ntrl plutC > 0 in "before"')
        end
        if any(plutW_1.('to_ntrl') > 0)
            warning('ntrl-to-ntrl plutW > 0 in "before"')
        end
        if any(plutC_2.('to_ntrl') > 0)
            warning('ntrl-to-ntrl plutC > 0 in "after"')
        end
        if any(plutW_2.('to_ntrl') > 0)
            warning('ntrl-to-ntrl plutW > 0 in "after"')
        end
    end

end

disp('Done checking that landsymm_plut files match expectations.')

%% Check that various files are identical

% Show up to $nrows of differences, optionally stopping after first different file
nrows = 3 ;
stop_after_first_file = false ;

check_tables_identical_loop(true, files_identical, testDir1, testDir2, stop_after_first_file, nrows, cropList);

fprintf('\nDone checking: Files that should be identical are.\n')


%% Check that various files are NOT identical

% Optionally stop after first different file
stop_after_first_file = false ;

check_tables_identical_loop(false, files_NOTidentical, testDir1, testDir2, stop_after_first_file, 0, cropList);

fprintf('\nDone checking: Files that should NOT be identical aren''t.\n')


%% Check cflux files (columns of flux type, not stand)

cflux_type_cols_equal_all = {'Veg', 'Repr', 'Soil', 'Fire', 'Est', 'Seed', 'Harvest', 'twigleaf', 'krv', 'bnl'} ;
% "Harvest" should never include wood, because wood harvest is all lumped
% into LU_ch in my setup.
cflux_type_cols_unequal_all = {} ;

cflux_type_cols_equal_lc = ...
    {'cropland', {'LU_ch', 'LUch_wha', 'LUch_clr', 'LUch_luc', 'hrv'} ; % Because no wood on cropland
     'pasture', {'LU_ch', 'LUch_wha', 'LUch_clr', 'LUch_luc', 'hrv'} ; % Because no wood on pasture
     'forest', {'LU_ch', 'LUch_wha', 'LUch_clr', 'LUch_luc', 'hrv'} ;   % Because mgd. forest is never an area donor in this setup
     'natural', {'LUch_wha', 'LUch_luc'} ; % Because wood harvest is "clearing" flux
     'total', {'LUch_wha', 'LUch_luc'} ; % Because wood harvest is "clearing" flux
    } ;

cflux_type_cols_unequal_lc = ...
    {'natural', {'LUch_clr', 'hrv'} ; % Because wood harvest is "clearing" flux
     'total', {'LUch_clr', 'hrv'} ; % Because wood harvest is "clearing" flux
    } ;

for f = 1:length(files_cflux_nonsts)
    thisFile = files_cflux_nonsts{f} ;
    disp(thisFile)
    if strcmp(thisFile, 'cflux.out')
        thisLU = 'total' ;
    else
        thisLU = strrep(strrep(thisFile, 'cflux_', ''), '.out', '') ;
    end
    [t1, t2] = read_two_tables(testDir1, testDir2, thisFile, cropList) ;

    varList = setdiff(t1.Properties.VariableNames, {'Lon', 'Lat', 'Year'}, 'stable') ;
    for v = 1:length(varList)
        thisVar = varList{v} ;
        v1 = t1.(thisVar) ;
        v2 = t2.(thisVar) ;

        is_special_equal_lu = any(strcmp(cflux_type_cols_equal_lc(:,1), thisLU)) ;
        is_special_equal_var = is_special_equal_lu && any(strcmp(thisVar, cflux_type_cols_equal_lc{strcmp(cflux_type_cols_equal_lc(:,1), thisLU),2})) ;
        is_special_unequal_lu = any(strcmp(cflux_type_cols_unequal_lc(:,1), thisLU)) ;
        is_special_unequal_var = is_special_unequal_lu && any(strcmp(thisVar, cflux_type_cols_unequal_lc{strcmp(cflux_type_cols_unequal_lc(:,1), thisLU),2})) ;
        if is_special_equal_var && is_special_unequal_var
            error('is_special_equal_var && is_special_unequal_var')
        elseif any(strcmp(thisVar, cflux_type_cols_equal_all)) || is_special_equal_var
            if ~isequal(v1, v2)
                warning('%s unexpectedly unequal', thisVar)
            end
        elseif any(strcmp(thisVar, cflux_type_cols_unequal_all)) || is_special_unequal_var
            if isequal(v1, v2)
                warning('%s unexpectedly equal', thisVar)
            end
        elseif strcmp(thisVar, 'Slow_h')
            if any(v2 > 0)
                warning('Slow_h unexpectedly nonzero in "after"')
            end
        elseif strcmp(thisVar, 'NEE')
            if any(v2 > v1)
                warning('NEE unexpectedly higher in "after"')
            end
        elseif strcmp(thisVar, 'LU_ch')
%             check_twigleaf_krv_bnl(v1, t1, 'before')
            check_twigleaf_krv_bnl(v2, t2, 'after')
        else
            fprintf('Skipping column %s (equal? %d)\n', thisVar, isequal(v1, v2))
        end
    end

    if any(t2.PLUMwood)
        warning('"PLUMwood" in "after" should always be zero (b/c no flux to atmosphere); it''s not')
    end

    disp(' ')
end


%% Check cpool files (columns of pool type, not stand)

cpool_type_cols_equal_all = {'VegC', 'LitterC', 'SoilC'} ;
cpool_type_cols_unequal_all = {'HarvSlowC/PLUMwoodC', 'Total'} ;
% These are unequal for all because they go to the destination LU, so the
% initial transition from Natural will result in them always being unequal
% starting 1850 (or whenever first LUC happens).

% cpool_type_cols_equal_lc = ...
%     {'cropland', {'LU_ch', 'LUch_wha', 'LUch_clr', 'LUch_luc', 'hrv'} ; % Because no wood on cropland
%      'pasture', {'LU_ch', 'LUch_wha', 'LUch_clr', 'LUch_luc', 'hrv'} ; % Because no wood on pasture
%      'forest', {'LU_ch', 'LUch_wha', 'LUch_clr', 'LUch_luc', 'hrv'} ;   % Because mgd. forest is never an area donor in this setup
%      'natural', {'LUch_wha', 'LUch_luc'} ; % Because wood harvest is "clearing" pool
%      'total', {'LUch_wha', 'LUch_luc'} ; % Because wood harvest is "clearing" pool
%     } ;
% 
% cpool_type_cols_unequal_lc = ...
%     {'natural', {'LUch_clr', 'hrv'} ; % Because wood harvest is "clearing" pool
%      'total', {'LUch_clr', 'hrv'} ; % Because wood harvest is "clearing" pool
%     } ;

cpool_type_cols_equal_lc = {} ;
cpool_type_cols_unequal_lc = {} ;

for f = 1:length(files_cpool_nonsts)
    thisFile = files_cpool_nonsts{f} ;
    disp(thisFile)
    if strcmp(thisFile, 'cpool.out')
        thisLU = 'total' ;
    else
        thisLU = strrep(strrep(thisFile, 'cpool_', ''), '.out', '') ;
    end
    [t1, t2] = read_two_tables(testDir1, testDir2, thisFile, cropList) ;
    
    varList = setdiff(t1.Properties.VariableNames, {'Lon', 'Lat', 'Year'}, 'stable') ;
    for v = 1:length(varList)
        thisVar = varList{v} ;

        v1 = t1.(thisVar) ;
        if strcmp(thisVar, 'HarvSlowC')
            thisVar = 'HarvSlowC/PLUMwoodC' ;
            v2 = t2.('PLUMwoodC') ;
        else
            v2 = t2.(thisVar) ;
        end

        is_special_equal_lu = ~isempty(cpool_type_cols_equal_lc) && any(strcmp(cpool_type_cols_equal_lc(:,1), thisLU)) ;
        is_special_equal_var = is_special_equal_lu && any(strcmp(thisVar, cpool_type_cols_equal_lc{strcmp(cpool_type_cols_equal_lc(:,1), thisLU),2})) ;
        is_special_unequal_lu = ~isempty(cpool_type_cols_unequal_lc) && any(strcmp(cpool_type_cols_unequal_lc(:,1), thisLU)) ;
        is_special_unequal_var = is_special_unequal_lu && any(strcmp(thisVar, cpool_type_cols_unequal_lc{strcmp(cpool_type_cols_unequal_lc(:,1), thisLU),2})) ;
        if is_special_equal_var && is_special_unequal_var
            error('is_special_equal_var && is_special_unequal_var')
        elseif any(strcmp(thisVar, cpool_type_cols_equal_all)) || is_special_equal_var
            if ~isequal(v1, v2)
                warning('%s unexpectedly unequal', thisVar)
            end
        elseif any(strcmp(thisVar, cpool_type_cols_unequal_all)) || is_special_unequal_var
            if isequal(v1, v2)
                warning('%s unexpectedly equal', thisVar)
            end
        else
            fprintf('Skipping column %s (equal? %d)\n', thisVar, isequal(v1, v2))
        end
    end

    if any(t1.HarvSlowC > t2.PLUMwoodC)
        error('any(t1.HarvSlowC > t2.PLUMwoodC)')
    end

    disp(' ')
end


%% For files that should be zero in "before" XOR "after," make sure they are

% Zero in "before"
for f = 1:length(files_zero_in_before)
    thisFile = files_zero_in_before{f} ;
    [t1, t2] = read_two_tables(testDir1, testDir2, thisFile, cropList) ;
    get_sum = @(x) sum(sum(table2array(x(:,4:end)))) ;
    if get_sum(t1) >= 0
        warning('%s "before" unexpectedly nonzero')
    end
    if get_sum(t2) == 0
        warning('%s "after" unexpectedly zero')
    end
end
disp('Done checking files that should be zero in "after"')

% Zero in "after"
for f = 1:length(files_zero_in_after)
    thisFile = files_zero_in_after{f} ;
    [t1, t2] = read_two_tables(testDir1, testDir2, thisFile, cropList) ;
    get_sum = @(x) sum(sum(table2array(x(:,4:end)))) ;
    if get_sum(t1) == 0
        warning('%s "before" unexpectedly zero')
    end
    if get_sum(t2) > 0
        warning('%s "after" unexpectedly nonzero')
    end
end
disp('Done checking files that should be zero in "after"')


%% after > before by ∆HarvSlowC+HarvSlowC_flux2atm

% cpool LitterC gets    patchpft.litter_leaf + patchpft.litter_root + patchpft.litter_sap + patchpft.litter_heart + patchpft.litter_repr + patch.soil.sompool[r].cmass
% ctotal_sts    gets    patchpft.litter_leaf + patchpft.litter_root + patchpft.litter_sap + patchpft.litter_heart + patchpft.litter_repr


%% FUNCTIONS

function check_twigleaf_krv_bnl(v, t, which_run) %#ok<INUSL> 

the_expr = 't.twigleaf + t.krv + t.bnl + t.hrv ;' ;
the_sum = eval(the_expr) ; %#ok<EVLCS> 
prec = 1e-6 ;

if any(abs(v - the_sum) - prec > 1e-18)
    if all(v >= the_sum)
        warning('LU_ch "%s" >= %s. Why?', which_run, the_expr)
    elseif all(v <= the_sum)
        warning('LU_ch "%s" <= %s. Why?', which_run, the_expr)
    elseif any(v ~= the_sum)
        warning('LU_ch "%s" ~= %s. Why?', which_run, the_expr)
    end
    keyboard
    % array2table([t.Year, v, the_sum, t.twigleaf, t.krv, t.bnl, t.hrv, v - the_sum])
%     t(t.Year>=1898 & t.Year<=1902, :)
end

end


function [outFiles, outFilenames] = get_outFiles(whichDir)
pattern = sprintf('%s/*out', whichDir) ;
outFiles = dir(pattern) ;
if isempty(outFiles)
    error('No files found matching pattern: %s', pattern)
end
outFilenames = {outFiles.name} ;
end


function ok = compare_file_lists(diffFilenames_inAnotB, files_ok_inAnotB, which_file_A)
ok = true ;
if ~isempty(diffFilenames_inAnotB)
    for f = 1:length(diffFilenames_inAnotB)
        thisFile = diffFilenames_inAnotB{f} ;
        if ~any(strcmp(files_ok_inAnotB, thisFile))
            if which_file_A == 1
                fprintf('Found in before but not after: %s\n', thisFile)
            else
                fprintf('Found in after but not before: %s\n', thisFile)
            end
            ok = false ;
        end
    end
end
end


function [errMsg, t1, t2] = check_tables_identical(t1, t2, thisFile, testDir1, testDir2, noteWhenDiff)

errMsg = '' ;
print_errMsg = true ;

if noteWhenDiff
    if ~isequaln(t1, t2)
        cols1 = t1.Properties.VariableNames ;
        cols2 = t2.Properties.VariableNames ;
        if length(cols1) ~= length(cols2)
            errMsg = sprintf('Different # columns (%d vs %d)', length(cols1), length(cols2)) ;
        elseif ~isequaln(cols1, cols2)
            if any(strcmp(cols1, 'HarvSlowC')) && any(strcmp(cols2, 'PLUMwoodC'))
                t1.Properties.VariableNames(strcmp(cols1, 'HarvSlowC')) = {'PLUMwoodC'} ;
                errMsg = check_tables_identical(t1, t2, testDir1, testDir2, thisFile) ;
                print_errMsg = false ; % To avoid duplication when recursing
            else
                errMsg = 'Same # columns but different names' ;
            end
        else
            errMsg = 'Other' ;
        end
        if ~strcmp(errMsg, '') && print_errMsg
            fprintf('%s do not match: %s\n', thisFile, errMsg)
        end
    end
elseif isequaln(t1, t2)
    errMsg = sprintf('%s are unexpectedly equal', thisFile) ;

    % Could be identical if, in "before", `harvest_slow_frac` is 0
%     if strcmp(thisFile, 'cflux_landsymm_sts.out')
    if strcmp(thisFile, 'cflux_noslowh_sts.out')
        x1 = lpjgu_matlab_readTable(sprintf('%s/cmass_wood_harv_toprod_sts.out', testDir1), ...
            'dont_save_MAT', true, 'do_save_MAT', false, 'verboseIfNoMat', false) ;
        x2 = lpjgu_matlab_readTable(sprintf('%s/cmass_wood_harv_toprod_sts.out', testDir2), ...
            'dont_save_MAT', true, 'do_save_MAT', false, 'verboseIfNoMat', false) ;
        t1(t1.Year>=1848 & t1.Year<=1852,:)
        x1(t1.Year>=1848 & t1.Year<=1852,:)
        x2(t1.Year>=1848 & t1.Year<=1852,:)
        stop
        [~, lu_inds] = setdiff(x1.Properties.VariableNames, {'Lon', 'Lat', 'Year'}) ;
        if ~any(table2array(x1(:,lu_inds)) > 0, 1)
            errMsg = '' ;
        else
            errMsg = [errMsg ' (cmass_wood_harv_toprod_sts.out has nonzeros)'] ;
        end
    end
    if ~strcmp(errMsg, '')
        disp(errMsg)
    end
        
end

end


function check_tables_identical_loop(noteWhenDiff, fileList, testDir1, testDir2, stop_after_first_file, nrows, cropList)
for f = 1:length(fileList)
    thisFile = fileList{f} ;
    [t1, t2] = read_two_tables(testDir1, testDir2, thisFile, cropList) ;
    errMsg = check_tables_identical(t1, t2, thisFile, testDir1, testDir2, noteWhenDiff) ;
    if stop_after_first_file && ~strcmp(errMsg, '')
        break
    end
    if noteWhenDiff
        if strcmp(errMsg, 'Other')
            differences = table2array(t1) ~= table2array(t2) ;
            
            % Find mismatched columns
            badCols = find(any(differences, 1)) ;
            if isempty(badCols)
                error('???')
            end
            
            % Find mismatched rows
            badRows = find(all(differences(:,badCols), 2)) ;
            if isempty(badRows)
                disp('No rows found where all badCols are bad; falling back to rows where only some are')
                badRows = any(differences, 2) ;
            end
            
            % Display
            bad1 = t1(badRows,:) ;
            bad2 = t2(badRows,:) ;
            fprintf('Showing only %d bad rows:\n', nrows')
            disp('  Before:')
            disp(head(bad1, nrows))
            disp('  After:')
            disp(head(bad2, nrows))
    
        elseif ~strcmp(errMsg, '')
            stop
        end
    end
end
end


function [t1, t2] = read_two_tables(testDir1, testDir2, thisFile, cropList)

t1 = lpjgu_matlab_readTable(sprintf('%s/%s', testDir1, thisFile), ...
    'dont_save_MAT', true, 'do_save_MAT', false, 'verboseIfNoMat', false) ;
t2 = lpjgu_matlab_readTable(sprintf('%s/%s', testDir2, thisFile), ...
    'dont_save_MAT', true, 'do_save_MAT', false, 'verboseIfNoMat', false) ;

% Get one representative crop
varNames = t1.Properties.VariableNames ;
[crops_found, IA] = intersect(varNames, cropList, 'stable') ;
other_cols = setdiff(1:length(varNames), IA) ;
if ~isempty(crops_found)
    a12_sum = table2array(t1(:,IA)) + table2array(t2(:,IA)) ;
    a12_colsums = sum(a12_sum, 1) ;
    if any(a12_colsums)
        crop_inds_incl = find(a12_colsums) ;
        crop_ind_keep = IA(crop_inds_incl(1)) ;
    else
        crop_ind_keep = IA(1) ;
    end
    cols_keep = [other_cols crop_ind_keep] ;
    t1 = t1(:,cols_keep) ;
    t2 = t2(:,cols_keep) ;
    t1.Properties.VariableNames{end} = 'crop' ;
    t2.Properties.VariableNames{end} = 'crop' ;
end

end






