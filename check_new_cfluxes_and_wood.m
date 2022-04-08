%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing outputs after giving PLUM all wood %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/landsymm-dev-forestry/all_wood_to_plum' ;
testDir1 = sprintf('%s/out.before.latest', topDir) ;
testDir2 = sprintf('%s/out.after.latest', topDir) ;

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
    'lai_natural.out', 'landsymm_ctree_ag_sts.out', 'landsymm_ctree_sts.out', 'simfireanalysis.out'} ;

files_harvPoolDiff = {'cpool.out', 'ctotal_sts.out'} ;

% cropList = {'CerealsC3', 'CerealsC3i', 'CerealsC4', 'CerealsC4i', 'Rice', 'Ricei', 'OilNfix', 'OilNfixi', ...
%     'OilOther', 'OilOtheri', 'Pulses', 'Pulsesi', 'StarchyRoots', 'StarchyRootsi', 'FruitAndVeg', ...
%     'FruitAndVegi', 'Sugarbeet', 'Sugarbeeti', 'Sugarcane', 'Sugarcanei', 'ExtraCrop'} ;


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

unclassifiedFiles = setdiff(outFilenames1, [files_identical]) ;
if ~isempty(unclassifiedFiles)
    disp('Unclassified files:')
    disp(unclassifiedFiles')
end

files_NOTidentical = setdiff(outFilenames1, files_identical) ;


%% Check that various files are identical

% Show up to $nrows of differences, optionally stopping after first different file
nrows = 3 ;
stop_after_first_file = false ;

check_tables_identical_loop(true, files_identical, testDir1, testDir2, stop_after_first_file, nrows);

fprintf('\nDone checking: Files that should be identical are.\n')


%% Check that various files are NOT identical

% Optionally stop after first different file
stop_after_first_file = false ;

check_tables_identical_loop(false, files_NOTidentical, testDir1, testDir2, stop_after_first_file, 0);

fprintf('\nDone checking: Files that should NOT be identical aren''t.\n')


%% FUNCTIONS

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


function [errMsg, t1, t2] = check_tables_identical(t1, t2, thisFile, noteWhenDiff)

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
                errMsg = check_tables_identical(t1, t2, thisFile) ;
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
    errMsg = fprintf('%s are unexpectedly equal\n', thisFile) ;
end

end


function check_tables_identical_loop(noteWhenDiff, fileList, testDir1, testDir2, stop_after_first_file, nrows)
for f = 1:length(fileList)
    thisFile = fileList{f} ;
    t1 = lpjgu_matlab_readTable(sprintf('%s/%s', testDir1, thisFile), ...
        'dont_save_MAT', true, 'do_save_MAT', false, 'verboseIfNoMat', false) ;
    t2 = lpjgu_matlab_readTable(sprintf('%s/%s', testDir2, thisFile), ...
        'dont_save_MAT', true, 'do_save_MAT', false, 'verboseIfNoMat', false) ;
    errMsg = check_tables_identical(t1, t2, thisFile, noteWhenDiff) ;
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