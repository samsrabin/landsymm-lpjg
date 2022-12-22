%%

% bl = repmat(yield_bl_lpj.garr_xv, [1 1 8]) ;
% fu = yield_fu_lpj.garr_xvt ;
bl = yield_bl_lpj.garr_xv ;
fu = yield_fu_lpj.garr_xvt(:,:,end) ;

bl = bl(:) ;
fu = fu(:) ;
delta =  fu ./ bl ;
delta(bl==0) = [] ;
bl(bl==0) = [] ;

figure; hist(delta, 100)
% figure; hist(log10(delta_xvt), 100)

thresh_tph = 0.1 ;
figure; hist(delta(bl>thresh_tph*0.1), 100)




%%

orig1 = lpjgu_matlab_readTable_then2map( ...
    '/Users/Shared/GGCMI2PLUM_sh/lpj-guess_runs/GGCMIPLUM_2001-2100_remap6p7_forPotYields_rcp45_forED_20191104133630/2011-2015/yield.out.gz', ...
    'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat', false) ;
orig2 = lpjgu_matlab_readTable_then2map( ...
    '/Users/Shared/GGCMI2PLUM_sh/lpj-guess_runs/GGCMIPLUM_2001-2100_remap6p7_forPotYields_rcp45_forED_20191104133630/2016-2020/yield.out.gz', ...
    'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat', false) ;
orig = orig1 ;
orig.maps_YXv = mean(cat(4, orig1.maps_YXv, orig2.maps_YXv), 4) ;

new = lpjgu_matlab_readTable_then2map( ...
    '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_GGCMIcrops_20191008/IPSL-CM5A-MR/LPJmL/rcp45/2011-2020/yield_deltad.out', ...
    'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat', false) ;


emubl = lpjgu_matlab_readTable_then2map( ...
    '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_GGCMIcropsBaseline_20191008/LPJmL/yield.out.gz', ...
    'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat', false) ;


%%

thisVar = 'Pulsesi060' ;
    
shademap(orig.maps_YXv(:,:,strcmp(orig.varNames, thisVar))) ;
set(gcf,'Position',[2721         170        1280         635])
this_caxis = caxis ;
title(sprintf('Original: %s', thisVar))

shademap(new.maps_YXv(:,:,strcmp(new.varNames, thisVar))) ;
set(gcf,'Position',[2721        -539        1280         635])
% caxis(this_caxis)
title(sprintf('Emulated: %s', thisVar))

%%

thisVar = 'soyi060' ;
tmp = emubl.maps_YXv(:,:,strcmp(emubl.varNames,thisVar)) < 0.0001 ;
tmp = double(tmp) ;
tmp(isnan(emubl.maps_YXv(:,:,strcmp(emubl.varNames,thisVar)))) = NaN ;
shademap(tmp)



%%

for v = 1:length(new.varNames)
    thisVar = new.varNames{v} ;
    
    shademap(orig.maps_YXv(:,:,strcmp(orig.varNames, thisVar))) ;
    set(gcf,'Position',[2721         170        1280         635])
    this_caxis = caxis ;
    title(sprintf('Original: %s', thisVar))
    
    shademap(new.maps_YXv(:,:,strcmp(new.varNames, thisVar))) ;
    set(gcf,'Position',[2721        -539        1280         635])
%     caxis(this_caxis)
    title(sprintf('Emulated: %s', thisVar))
    
    pause(5)
    close all
end



%% Check order of N

new = lpjgu_matlab_read2geoArray( ...
    '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_GGCMIcrops_20191008/IPSL-CM5A-MR/LPJmL/rcp45/2011-2020/yield_deltad.out', ...
    'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat', false) ;

Nvars = length(new.varNames) ;
for c = 1:Nvars/3
    cc = (c-1)*3 + 1 ;
    tmp = new.garr_xv(:,cc:cc+2) ;
    is_123 = tmp(:,1) <= tmp(:,2) & tmp(:,2) <= tmp(:,3) ;
    is_xx1 = tmp(:,1) >= tmp(:,2) & tmp(:,1) >= tmp(:,3) ;
    is_x2x = tmp(:,1) <= tmp(:,2) & tmp(:,2) >= tmp(:,3) ;
    is_321 = tmp(:,1) >= tmp(:,2) & tmp(:,2) >= tmp(:,3) ;
    
    thisCrop = new.varNames{cc} ;
    thisCrop = thisCrop(1:end-3) ;
    fprintf('%s: %d cells not 1-2-3 or 3-2-1\n', thisCrop, length(find(~is_123 & ~is_xx1))) ;
end


