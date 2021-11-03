%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check for missing irrigation data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% topDir = '/Users/sam/Documents/Dropbox/2016_KIT/GGCMI/GGCMI2PLUM_DB/CMIP_emulated_forPLUM_tars/CMIP_emulated_forPLUM_20211005/A1_v2.5_UKESM1-0-LL_ssp126_20210923_ph2bl_intpinfs_rmolend_fake1k/sim_LPJ-GUESS' ;
topDir = '/Users/sam/Documents/Dropbox/2016_KIT/GGCMI/GGCMI2PLUM_DB/CMIP_emulated_forPLUM_tars/sim_LPJ-GUESS' ;


%%

ny = 3 ;
nx = 3 ;
yLim = 65:360 ;
spacing = [0.025 0.025] ; % v h
fontSize = 14 ;

dirList = dir(sprintf('%s/*/tot_runoff.out*', topDir)) ;
figure('Color', 'w', 'Position', figurePos) ;

for d = 1:length(dirList)
    thisFile = sprintf('%s/%s', dirList(d).folder, dirList(d).name) ;
    runoff = lpjgu_matlab_read2geoArray(thisFile, ...
        'force_mat_save', false, ...
        'force_mat_nosave', true, ...
        'verboseIfNoMat', false) ;
    surf_YX = lpjgu_vector2map( ...
        runoff.garr_xv(:,strcmp(runoff.varNames, 'Surf')), ...
        [360 720], runoff.list2map) ;
    
    
    subplot_tight(ny, nx, d, spacing) ;
    pcolor(log10(surf_YX(yLim,:))) ; shading flat ; axis equal tight off
    
    hcb = colorbar ;
    ylabel(hcb, 'log10(mm)')
    
    set(gca, 'FontSize', 14)
    
    tmp = strsplit(dirList(d).folder, '/') ;
    title(tmp{end})
end