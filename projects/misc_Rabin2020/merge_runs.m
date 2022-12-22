%%%%%%%%%%%%%%%%%%%%
%%% Merge 2 runs %%%
%%%%%%%%%%%%%%%%%%%%

% inDir_1 = 'PLUM2LPJGblIrr_26_s1_1850-2005_moreProcs/output-2017-12-23-115308' ;
% inDir_2 = 'PLUM2LPJGblIrr_SSP1v2_RCP45/output-2017-12-25-140114' ;
% outDir = '/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/testOut' ;

% inDir_1 = 'PLUM2LPJGblIrrNoFG_26_s1_1850-2005_moreProcs/output-2018-01-06-071106' ;
% inDir_2 = 'PLUM2LPJGblIrrNoFG_26_s1_2006-2010_moreProcs/output-2018-01-08-111859' ;
% outDir = '/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/PLUM2LPJGblIrrNoFG_26_s1_1850-2010_moreProcs' ;

% inDir_1 = 'PLUM2LPJG_PLUM7_1850-2005/output-2018-01-23-195509' ;
% inDir_2 = 'PLUM2LPJG_PLUM7_2006-2010/output-2018-01-23-204559' ;
% % outDir = '/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/PLUM2LPJG_PLUM7_1850-2010' ;
% outDir = '/Users/Shared/PLUM/trunk_runs/PLUM2LPJG_PLUM7_1850-2010' ;

% inDir_1 = 'LPJGPLUM_1850-2005_PLUM6xtra/output-2018-04-24-042556' ;
% inDir_2 = 'LPJGPLUM_2006-2010_PLUM6xtra/output-2018-04-21-061836' ;
% outDir = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_PLUM6xtra' ;

% inDir_1 = 'LPJGPLUM_1850-2005_PLUM6xtraMisc/output-2018-04-29-103725' ;
% inDir_2 = 'LPJGPLUM_2006-2010_PLUM6xtraMisc/output-2018-05-02-104249' ;
% outDir = '/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_PLUM6xtraMisc' ;
% outDir_full = addslashifneeded('/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180502104840') ;

% inDir_1 = 'LPJGPLUM_1850-2005_PLUM6xtraMisc/output-2018-06-05-022332' ;
% inDir_2 = 'LPJGPLUM_2006-2010_PLUM6xtraMisc/output-2018-06-05-133701' ;
% outDir = '/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_PLUM6xtraMisc' ;

inDir_1 = 'LPJGPLUM_1850-2005_PLUM6xtraMisc_cg/output-2018-07-31-055344' ;
inDir_2 = 'LPJGPLUM_2006-2010_PLUM6xtraMisc_cg/output-2018-08-01-082312' ;
outDir = '/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_PLUM6xtraMisc_cg' ;



%% Setup

addpath(genpath(landsymm_lpjg_path()))

% Get directories
inDir_1_full = find_PLUM2LPJG_run(inDir_1) ;
inDir_2_full = find_PLUM2LPJG_run(inDir_2) ;
if ~exist('outDir_full','var')
    outDir_full = addslashifneeded([addslashifneeded(outDir) 'matlab_merge_' datestr(now,'yyyymmddHHMMSS')]) ;
else
    disp(['Using outDir_full as specified: ' outDir_full])
end
if ~exist(outDir_full,'dir')
    mkdir(outDir_full)
end

% Get list of files
files_inDir_1 = dir([inDir_1_full '*.out.gz']) ;
files_inDir_2 = dir([inDir_2_full '*.out.gz']) ;
files_inDir_1 = {files_inDir_1.name} ;
files_inDir_2 = {files_inDir_2.name} ;
files_in = intersect(files_inDir_1,files_inDir_2) ;
Nfiles = length(files_in) ;
in1not2 = setdiff(files_inDir_1,files_inDir_2) ;
in2not1 = setdiff(files_inDir_2,files_inDir_1) ;
disp('Files to be processed:')
for f = 1:Nfiles
    disp(['   ' files_in{f}])
end
if ~isempty(in1not2) || ~isempty(in2not1)
    disp('Files to be ignored:')
    for f = 1:length(in1not2)
        disp(['   ' inDir_1 ': ' in1not2{f}])
    end
    for f = 1:length(in2not1)
        disp(['   ' inDir_2 ': ' in2not1{f}])
    end
else
    disp('No files being ignored.')
end
disp(' ')


%% Process

unix(['echo ' inDir_1 ' > ' outDir_full 'in_dirs.txt']) ;
unix(['echo ' inDir_2 ' >> ' outDir_full 'in_dirs.txt']) ;

for f = 1:Nfiles
    thisFile = files_in{f} ;
    disp([thisFile ':'])
    
    % Check for already-existing out_file
    out_file = [outDir_full strrep(thisFile,'.gz','.mat')] ;
    if exist(out_file,'file')
        warning('   out_file already exists; skipping.')
        continue
    end
    
    % Import
    disp('   Importing table 1...')
    inTable_1 = lpjgu_matlab_readTable([inDir_1_full thisFile],'verboseIfNoMat',true,'do_save_MAT',true,'dispPrefix','   ') ;
    disp('   Importing table 2...')
    inTable_2 = lpjgu_matlab_readTable([inDir_2_full thisFile],'verboseIfNoMat',true,'do_save_MAT',true,'dispPrefix','   ') ;
    
    % Sanity checks
    disp('   Checking...')
    %%% Equal gridlist
    if ~isequal(unique([inTable_1.Lon inTable_1.Lat],'rows'),unique([inTable_2.Lon inTable_2.Lat],'rows'))
        error('Gridlists do not match!') ;
    end
    %%% Contiguous yearLists
    if max(inTable_1.Year)+1 ~= min(inTable_2.Year)
        error('Year lists are not contiguous!')
    end
    
    % Create out_table
    disp('   Making out_table...')
    out_table = [inTable_1 ; inTable_2] ;
    clear inTable_*
    out_table = sortrows(out_table) ;
    
    % Save out_table
    disp('   Saving...')
    
    save(out_file,'out_table','-v7.3') ;
    
    clear out_table
end
disp('Copying ins-files from first run...')
['cp ' inDir_1_full '*.ins ' outDir_full]
disp('Done.')
disp(' ')

