n = 60 ;

spacing = [0.05 0.05] ; % y x
fontSize = 18 ;
clims = [0 4] ;


figure('Color','w','Position',figurePos) ;

for adapt = 0:1
    
    file_ph2 = sprintf(['/Volumes/Reacher/GGCMI/AgMIP.output/EPIC-TAMU/phase2/spring_wheat/A%d/yield/' ...
        'epic-tamu_agmerra_fullharm_yield_swh_global_annual_1980_2010_C360_T0_Winf_N%d_A%d.nc4'], ...
        adapt, n, adapt) ;
    file_emu = sprintf(['/Volumes/Reacher/GGCMI/AgMIP.output/CMIP_emulated/'...
        'yields/CMIP6/A%d_N%d/ssp126/EPIC-TAMU/'...
        'cmip6_ssp126_UKESM1-0-LL_r1i1p1_spring_wheat_EPIC-TAMU_A%d_N%d_emulated_yield_baseline_1980_2010_average_v2.5.nc4'], ...
        adapt, n, adapt, n) ;
    
    ph2_YXy = flip(permute(ncread(file_ph2, 'yield_swh'), [2 1 3]), 1) ;
    ph2_YX = nanmean(ph2_YXy(:,:,1:end-1), 3) ;
    
    emu_YX = flip(permute(ncread(file_emu, 'yield_ir_swh'), [2 1 3]), 1) ;

    ii = (adapt*2)+1 ;
    subplot_tight(2, 2, ii, spacing)
    pcolor(ph2_YX); shading flat
    axis equal tight off
    title(sprintf('Phase 2 1980-2009 mean yield (A%d N%d)', ...
        adapt, n))
    hcb = colorbar ;
    caxis(clims); hcb.TickLabels{end} = sprintf('%s%s','≥',hcb.TickLabels{end}) ;
    set(gca, 'FontSize', fontSize)
    
    ii = (adapt*2)+2 ;
    subplot_tight(2, 2, ii, spacing)
    pcolor(emu_YX); shading flat
    axis equal tight off
    title(sprintf('Emulated baseline yield (A%d N%d)', ...
        adapt, n))
    hcb = colorbar ;
    caxis(clims); hcb.TickLabels{end} = sprintf('%s%s','≥',hcb.TickLabels{end})
    set(gca, 'FontSize', fontSize)
    
end