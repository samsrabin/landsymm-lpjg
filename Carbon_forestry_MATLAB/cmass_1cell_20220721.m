%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Look at cmass for individual stands %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import

% thisDir = '/Users/Shared/PLUM/forestry_tests_sam/outForPLUM-2022-06-30-040356/1850pot_1850-2014' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot-1949_20220628/out.distoff.nofire.10xlongev' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot-1949_20220628/out.20220721.nosuppress' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot-1949_20220628/out.20220721.suppress' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot/out.20220721' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot/out.20220928' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot.suppress/out' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot_manageforest_olddefault/out' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot_stripped/out' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot_tdf_simple/out.af593ab' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot_tdf_simple/out.36f38c0' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot_tdf_simple/out.462abe3' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot_tdf_simple/out.6273f4c' ;
thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot_tdf_simple/out.4be89cc' ;
% thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/1850pot_ts_20220721/pot_tdf_simple/out' ;

% thisFile = 'cmass_sts' ; units = 'kgC/m2' ; theseCols = {'ntrl'} ;
% thisFile = 'landsymm_pcutW_sts' ; units = 'kgC/m2' ; thisCol = {'ntrl'} ;
% thisFile = 'height_natural' ; units = 'm' ; theseCols = {'BNE', 'BINE', 'BNS', 'TeNE', 'TeBS', 'IBS', 'TeBE'} ;
thisFile = 'cmass_natural' ; units = 'kgC/m2' ; theseCols = {'BNE', 'BINE', 'BNS', 'TeNE', 'TeBS', 'IBS', 'TeBE'} ;

inFile = sprintf('%s/%s.out', thisDir, thisFile) ;
S = lpjgu_matlab_read2geoArray(inFile, 'force_mat_save', false, 'force_mat_nosave', true, ...
    'xres', 0.5, 'yres', 0.5) ;
Ncells = length(S.list2map) ;

thisLegend = {} ;
for c = 1:Ncells
    thisLegend{c} = sprintf('%0.2f, %0.2f', S.lonlats(c,1), S.lonlats(c,2)) ; %#ok<SAGROW> 
end


%% Plot

Ncols = length(theseCols) ;
if Ncols == 1
    ny = 1 ;
    nx = 1 ;
    thisPos = [421   652   890   464] ;
    spacing = [0.1 0.05] ; % y x
elseif Ncols >= 7 && Ncols <= 8
    ny = 2 ;
    nx = 4 ;
    thisPos = figurePos ;
    spacing = [0.1 0.05] ; % y x
else
    error('Specify figure properties for Ncols %d', Ncols)
end

figure('Color', 'w', 'Position', thisPos)

for c = 1:length(theseCols)

    subplot_tight(ny, nx, c, spacing)
    
    thisCol = theseCols{c} ;
    isThisCol = find(strcmp(S.varNames, thisCol)) ;
    if length(isThisCol) ~= 1
        error('Expected 1 matching stand; found %d', length(isThisCol))
    end
    plot(S.yearList, squeeze(S.garr_xvy(:,isThisCol,:)))
    
    ylabel(sprintf('%s', units))
    title(strrep(sprintf('%s: %s', thisFile, thisCol), '_', ', '))
    set(gca, 'FontSize', 14)
    
    hl = legend(gca, thisLegend, 'Location', 'Best') ;
    title(hl, 'Gridcells (lon, lat)')
end

export_fig(sprintf('%s/%s.png', thisDir, thisFile), '-r75')