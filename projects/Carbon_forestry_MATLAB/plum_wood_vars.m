%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Is what I need already in outputs? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisDir = '/Users/Shared/lpj-guess_git-svn_20190828_20211202/ssr/landsymm-dev-forestry/all_wood_to_plum/out' ;


%% Import

% From original outputs
cmass_wood_potharv_sts = lpjgu_matlab_readTable(sprintf('%s/cmass_wood_potharv_sts.out', thisDir), ...
    'do_save_MAT', false, 'dont_save_MAT', true, 'verboseIfNoMat', false) ;
cmass_wood_potharv_toprod_sts = lpjgu_matlab_readTable(sprintf('%s/cmass_wood_potharv_toprod_sts.out', thisDir), ...
    'do_save_MAT', false, 'dont_save_MAT', true, 'verboseIfNoMat', false) ;
forest_vegc = lpjgu_matlab_readTable(sprintf('%s/forest_vegc.out', thisDir), ...
    'do_save_MAT', false, 'dont_save_MAT', true, 'verboseIfNoMat', false) ;
orig2 = cmass_wood_potharv_sts.forC + forest_vegc.for_fuel ;
orig3 = orig2 + cmass_wood_potharv_toprod_sts.forC ;

% From my new output
cmass_wood_potharv_toplum_sts = lpjgu_matlab_readTable(sprintf('%s/cmass_wood_potharv_toplum_sts.out', thisDir), ...
    'do_save_MAT', false, 'dont_save_MAT', true, 'verboseIfNoMat', false) ;
mine = cmass_wood_potharv_toplum_sts.forC ;

yearList = cmass_wood_potharv_toplum_sts.Year ;


%% Plot

figure('Color', 'w', 'Position', [433   617   757   351]) ;
plot(yearList, [mine orig2 orig3])
ylabel('kgC m^{-2}')
set(gca, 'FontSize', 14)

thisLegend = {'mine' ;
              'orig2: cmass_wood_potharv + for_fuel' ;
              'orig3: orig2 + cmass_wood_potharv_toprod' ;
              } ;
thisLegend = strrep(thisLegend, '_', '\_') ;
legend(thisLegend, 'Location', 'Best')


%% Plot leaf mass as % of wood mass

leafmass = orig2 - mine ;
woodmass = orig2 - leafmass ;

figure('Color', 'w', 'Position', [433   617   757   351]) ;
plot(yearList, 100 * leafmass ./ woodmass)
ylabel('%')
set(gca, 'FontSize', 14)
title('Leaf mass as % of wood mass')

