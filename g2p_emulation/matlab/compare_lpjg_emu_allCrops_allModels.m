%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare LPJ-GUESS and emulator outputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

emu_list = {'LPJ-GUESS', 'LPJmL', 'pDSSAT', 'EPIC-TAMU'} ;

incl_cf = true ;

top_dir = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/potYields_20200310/IPSL-CM5A-MR_r1i1p1' ;
fig_dir = '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/outputs_figs/potYields_20200310_IPSL-CM5A-MR' ;


%% Setup

crop_list = {'CerealsC3','CerealsC4','Rice','Oilcrops','StarchyRoots','Pulses'} ;
Ncrops = length(crop_list) ;
irr_list = {'','i'} ;
N_list_lpj = {'0','0200'} ;
N_list_emu = {'010','200'} ;
Nn = length(N_list_lpj) ;
Nemu = length(emu_list) ;

if incl_cf
    cf_dir = '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/outputs_figs/calibration/calibration_tables' ;
    fig_dir = strrep(fig_dir, 'potYields', 'potYields_cf') ;
end

if ~exist(fig_dir,'dir')
    mkdir(fig_dir) ;
end


%% Read files: LPJ-GUESS (actual run)

yield_lpj1 = lpjgu_matlab_read2geoArray('/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/plum/yield.2011-2015.out.gz') ;
yield_lpj2 = lpjgu_matlab_read2geoArray('/Users/Shared/GGCMI2PLUM_sh/emulation/inputs/plum/yield.2016-2020.out.gz') ;
yield_lpj = yield_lpj1 ;
yield_lpj.garr_xv = nanmean(cat(3,yield_lpj1.garr_xv,yield_lpj2.garr_xv),3) ;
clear yield_lpj1 yield_lpj2
[~, IA] = intersect(yield_lpj.varNames, crop_list) ;
yield_lpj.garr_xv(:,IA) = [] ;
yield_lpj.varNames(IA) = [] ;
[~, IA] = intersect(yield_lpj.varNames, strcat(crop_list,'i')) ;
yield_lpj.garr_xv(:,IA) = [] ;
yield_lpj.varNames(IA) = [] ;

% Convert from kg/m2 to tons/ha
yield_lpj.garr_xv = 10 * yield_lpj.garr_xv ;

% Apply calibration factors, if needed
if incl_cf
    % Used in ES paper
    cf_lpj = nan(Ncrops,1) ;
    cf_lpj(strcmp(crop_list, 'CerealsC3')) = 1.056 ;
    cf_lpj(strcmp(crop_list, 'CerealsC4')) = 0.738 ;
    cf_lpj(strcmp(crop_list, 'Rice')) = 1.052 ;
    cf_lpj(strcmp(crop_list, 'Oilcrops')) = 0.687 ;
    cf_lpj(strcmp(crop_list, 'Pulses')) = 0.865 ;
    cf_lpj(strcmp(crop_list, 'StarchyRoots')) = 5.443 ;
    for c = 1:length(crop_list)
        thisCrop = crop_list{c} ;
        isThisCrop = contains(yield_lpj.varNames, thisCrop) ;
        if length(find(isThisCrop)) ~= 6
            error('length(find(isThisCrop)) ~= 6')
        end
        fprintf('%s lpj: %0.3f\n', thisCrop, cf_lpj(c)) ;
        yield_lpj.garr_xv(:,isThisCrop) = ...
            cf_lpj(c) * yield_lpj.garr_xv(:,isThisCrop) ;
    end
end


%% Read files: Emulators

for m = 1:Nemu
    thisEmu = emu_list{m} ;
    emu_dir = sprintf('%s/%s/rcp45/2011-2020', top_dir, thisEmu) ;
    yield_emu_tmp = lpjgu_matlab_read2geoArray(...
        sprintf('%s/yield.out', emu_dir)) ;
    
    % Set up big array
    if m==1
        yield_emu = yield_emu_tmp ;
        yield_emu = rmfield(yield_emu, 'garr_xv') ;
        yield_emu.garr_xvm = nan([size(yield_emu_tmp.garr_xv) Nemu]) ;
    end
    
    % Read calibration factors, if needed
    if incl_cf
        cf_file_emu = dir(sprintf('%s/*%s*', cf_dir, thisEmu)) ;
        if length(cf_file_emu) ~= 1
            error('%d possible cf_file_emu found', length(cf_file_emu))
        end
        cf_file_emu = sprintf('%s/%s', cf_dir, cf_file_emu.name) ;
        cf_table_emu = readtable(cf_file_emu) ;
        [C, IA, IB] = intersect(crop_list, cf_table_emu.PLUM_crop, 'stable') ;
        if length(IB) ~= length(crop_list)
            error('Mismatch between crops in crop_list and cf_table_emu')
        end
        cf_emu = cf_table_emu.calib_factor(IB) ;
        for c = 1:length(crop_list)
            thisCrop = crop_list{c} ;
            isThisCrop = contains(yield_emu_tmp.varNames, thisCrop) ;
            if length(find(isThisCrop)) ~= 6
                error('length(find(isThisCrop)) ~= 6')
            end
            if m==1
                fprintf('%s emu: %0.3f\n', thisCrop, cf_emu(c)) ;
            end
            yield_emu_tmp.garr_xv(:,isThisCrop) = ...
                cf_emu(c) * yield_emu_tmp.garr_xv(:,isThisCrop) ;
        end
    end
    
    % Save into big array, converting from kg/m2 to tons/ha
    yield_emu.garr_xvm(:,:,m) = 10 * yield_emu_tmp.garr_xv ;
    clear yield_emu_tmp
    
end


%% Compare for all crops: Maps

%%%% Options
spacing = [0.04 0.015] ; % v h
yrange = 70:360 ;
% thisPos = [1    34   720   771] ;
thisPos = figurePos ;
% thisPos = [1 375 1440 430] ;
%%%%

for c = 1:length(crop_list)
    thisCrop = crop_list{c} ;
    
    % Get data
    tmp_lpj_YXi = nan(360,720,2) ;
    tmp_emu_YXim = nan(360,720,2,Nemu) ;
    for ii = 1:2
        thisIrr = irr_list{ii} ;
        thisVar_lpj_Nmin = sprintf('%s%s%s',thisCrop,thisIrr,N_list_lpj{1}) ;
        thisVar_lpj_Nmax = sprintf('%s%s%s',thisCrop,thisIrr,N_list_lpj{end}) ;
        tmp_min = yield_lpj.garr_xv(:,strcmp(yield_lpj.varNames,thisVar_lpj_Nmin)) ;
        tmp_max = yield_lpj.garr_xv(:,strcmp(yield_lpj.varNames,thisVar_lpj_Nmax)) ;
        tmp = 100*(tmp_max-tmp_min)./tmp_min ;
        tmp(tmp_min==0 & tmp_max==0) = 0 ;
        tmp_YX = nan(360,720) ;
        tmp_YX(yield_lpj.list2map) = tmp ;
        tmp_lpj_YXi(:,:,ii) = tmp_YX ;
        
        for m = 1:Nemu
            thisVar_emu_Nmin = sprintf('%s%s%s',thisCrop,thisIrr,N_list_emu{1}) ;
            thisVar_emu_Nmax = sprintf('%s%s%s',thisCrop,thisIrr,N_list_emu{end}) ;
            tmp_min = yield_emu.garr_xvm(:,strcmp(yield_emu.varNames,thisVar_emu_Nmin),m) ;
            tmp_max = yield_emu.garr_xvm(:,strcmp(yield_emu.varNames,thisVar_emu_Nmax),m) ;
            tmp = 100*(tmp_max-tmp_min)./tmp_min ;
            tmp(tmp_min==0 & tmp_max==0) = 0 ;
            tmp_YX = nan(360,720) ;
            tmp_YX(yield_emu.list2map) = tmp ;
            tmp_emu_YXim(:,:,ii,m) = tmp_YX ;
        end
    end
    
    % Make figure
    figure('Color','w','Position',thisPos)
    excl_pctile = 2.5 ; new_caxis = [ ...
        prctile(cat(1,tmp_lpj_YXi(~isnan(tmp_lpj_YXi) & ~isinf(tmp_lpj_YXi)),tmp_emu_YXim(~isnan(tmp_emu_YXim) & ~isinf(tmp_emu_YXim))),excl_pctile) ...
        prctile(cat(1,tmp_lpj_YXi(~isnan(tmp_lpj_YXi) & ~isinf(tmp_lpj_YXi)),tmp_emu_YXim(~isnan(tmp_emu_YXim) & ~isinf(tmp_emu_YXim))),100-excl_pctile) ...
        ] ;
    units = sprintf('%% (excluding top/bottom %0.1f percentiles)', excl_pctile) ;
    outfile = sprintf('%s/LPJGcomp_map_ALL_%s_pctDiff.png', ...
        fig_dir, thisCrop) ;
    
    ny = Nemu/2 + 1 ;
    nx = 4 ;
    pltind = nx ;
    for ii = 1:2
        thisIrr = irr_list{ii} ;
        
        subplot_tight(ny,nx,ii+1,spacing)
        pcolor(tmp_lpj_YXi(yrange,:,ii)); shading flat; axis equal tight off
        caxis(gca,new_caxis) ; colormap(gca,'jet');
        title(sprintf('LPJ-GUESS: %s%s (N%s vs N%s)', thisCrop, thisIrr, N_list_lpj{end}, N_list_lpj{1}))
        hcb = colorbar('Location','SouthOutside') ; xlabel(hcb,units)
        
        % Line separating subplots
        if ii == 1
            annotation('line', [0 1],(ny-1)/ny*[1 1])
        end
        
        for m = 1:Nemu
            thisEmu = emu_list{m} ;
            pltind = pltind+1 ;
            subplot_tight(ny,nx,pltind,spacing)
            pcolor(tmp_emu_YXim(yrange,:,ii,m)); shading flat; axis equal tight off
            caxis(gca,new_caxis) ; colormap(gca,'jet');
            title(sprintf('%s emulator: %s%s (N%s vs N%s)', thisEmu, thisCrop, thisIrr, N_list_emu{end}, N_list_emu{1}))
            hcb = colorbar('Location','SouthOutside') ; xlabel(hcb,units)
            
            % Lines separating subplots
            if ii == 1 && m < Nemu
                annotation('line', m*[0.25 0.25],[0 (ny-1)/ny])
            end
        end

    end
    
    export_fig(outfile, '-r150') ;
    close
end

disp('Done')


%% FUNCTIONS

function [hcb, hyl] = add_big_colorbar(h, fontSize, units_map, bins_lowBnds)

thisPos = get(h,'Position') ;
hcb = colorbar(h,'Location','SouthOutside') ;
set(h,'Position',thisPos)
hcb.Position(2) = 0.075 ;
hyl = ylabel(hcb, sprintf('(%s)',units_map)) ;
hcb.FontSize = fontSize ;

% Mess with ticks
hcb.TickDirection = 'out' ;
if ~isempty(bins_lowBnds)
    Nbins = length(bins_lowBnds) ;
    hcb.Ticks = 1:(Nbins+1) ;
    bins_lowBnds_str = strrep(cellstr(num2str(bins_lowBnds')), ' ', '') ;
    bins_lowBnds_str = [bins_lowBnds_str ; {'100'}] ;
    hcb.TickLabels = bins_lowBnds_str ;
    for b = 1:length(bins_lowBnds_str)
        thisTick = str2num(bins_lowBnds_str{b}) ;
        if thisTick > 0
            hcb.TickLabels{b} = ['+' hcb.TickLabels{b}] ;
        end
    end
end

% Move units label
hyl.Units = 'normalized' ;
hyl.Position(1) = 1.09 ;
hyl.Position(2) = -0.5 ;


end
