%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare LPJ-GUESS and emulator outputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisEmu = 'LPJ-GUESS' ;
% thisEmu = 'LPJmL' ;
% thisEmu = 'pDSSAT' ;
% thisEmu = 'EPIC-TAMU' ;

incl_cf = true ;

top_dir = '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/potYields_20200310/IPSL-CM5A-MR_r1i1p1' ;
fig_dir = '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/outputs_figs/potYields_20200310_IPSL-CM5A-MR' ;


%% Setup

crop_list = {'CerealsC3','CerealsC4','Rice','Oilcrops','StarchyRoots','Pulses'} ;
irr_list = {'','i'} ;
N_list_lpj = {'0','0200'} ;
N_list_emu = {'010','200'} ;

if incl_cf
    cf_dir = '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/outputs_figs/calibration/calibration_tables' ;
    cf_file_emu = dir(sprintf('%s/*%s*', cf_dir, thisEmu)) ;
    if length(cf_file_emu) ~= 1
        error('%d possible cf_file_emu found', length(cf_file_emu))
    end
    cf_file_emu = sprintf('%s/%s', cf_dir, cf_file_emu.name) ;
    fig_dir = strrep(fig_dir, 'potYields', 'potYields_cf') ;
end

if ~exist(fig_dir,'dir')
    mkdir(fig_dir) ;
end


%% Read files

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

emu_dir = sprintf(...
    '%s/%s/rcp45/2011-2020', ...
    top_dir, thisEmu) ;
yield_emu = lpjgu_matlab_read2geoArray(...
    sprintf('%s/yield.out', emu_dir)) ;

% Convert from kg/m2 to tons/ha
yield_lpj.garr_xv = 10 * yield_lpj.garr_xv ;
yield_emu.garr_xv = 10 * yield_emu.garr_xv ;

% Read calibration factors, if needed
if incl_cf
    cf_table_emu = readtable(cf_file_emu) ;
    [C, IA, IB] = intersect(crop_list, cf_table_emu.PLUM_crop, 'stable') ;
    if length(IB) ~= length(crop_list)
        error('Mismatch between crops in crop_list and cf_table_emu')
    end
    cf_emu = cf_table_emu.calib_factor(IB) ;
    % Used in ES paper
    cf_lpj = nan(size(cf_emu)) ;
    cf_lpj(strcmp(crop_list, 'CerealsC3')) = 1.056 ;
    cf_lpj(strcmp(crop_list, 'CerealsC4')) = 0.738 ;
    cf_lpj(strcmp(crop_list, 'Rice')) = 1.052 ;
    cf_lpj(strcmp(crop_list, 'Oilcrops')) = 0.687 ;
    cf_lpj(strcmp(crop_list, 'Pulses')) = 0.865 ;
    cf_lpj(strcmp(crop_list, 'StarchyRoots')) = 5.443 ;
    for c = 1:length(crop_list)
        thisCrop = crop_list{c} ;
        isThisCrop = contains(yield_emu.varNames, thisCrop) ;
        if length(find(isThisCrop)) ~= 6
            error('length(find(isThisCrop)) ~= 6')
        end
        fprintf('%s emu: %0.3f\n', thisCrop, cf_emu(c)) ;
        yield_emu.garr_xv(:,isThisCrop) = ...
            cf_emu(c) * yield_emu.garr_xv(:,isThisCrop) ;
        isThisCrop = contains(yield_lpj.varNames, thisCrop) ;
        if length(find(isThisCrop)) ~= 6
            error('length(find(isThisCrop)) ~= 6')
        end
        fprintf('%s lpj: %0.3f\n', thisCrop, cf_lpj(c)) ;
        yield_lpj.garr_xv(:,isThisCrop) = ...
            cf_lpj(c) * yield_lpj.garr_xv(:,isThisCrop) ;
    end
end


%% Compare for all crops: Maps

%%%% Options
spacing = [0.04 0.015] ; % v h
yrange = 70:360 ;
% thisPos = [1    34   720   771] ;
% thisPos = figurePos ;
thisPos = [1 375 1440 430] ;
%%%%

Nn = length(N_list_lpj) ;

for c = 1:length(crop_list)
    thisCrop = crop_list{c} ;
    
    tmp_lpj_YXin = nan(360,720,2,Nn) ;
    tmp_emu_YXin = nan(360,720,2,Nn) ;
    for ii = 1:2
        thisIrr = irr_list{ii} ;
        for n = 1:Nn
            thisVar_lpj = sprintf('%s%s%s',thisCrop,thisIrr,N_list_lpj{n}) ;
            thisVar_emu = sprintf('%s%s%s',thisCrop,thisIrr,N_list_emu{n}) ;
            
            tmp_lpj = yield_lpj.garr_xv(:,strcmp(yield_lpj.varNames,thisVar_lpj)) ;
            tmp_emu = yield_emu.garr_xv(:,strcmp(yield_emu.varNames,thisVar_emu)) ;
            
            tmp_lpj_YX = nan(360,720) ;
            tmp_lpj_YX(yield_lpj.list2map) = tmp_lpj ;
            tmp_lpj_YXin(:,:,ii,n) = tmp_lpj_YX ;
            tmp_emu_YX = nan(360,720) ;
            tmp_emu_YX(yield_emu.list2map) = tmp_emu ;
            tmp_emu_YXin(:,:,ii,n) = tmp_emu_YX ;
        end
        
    end
    
    figure('Color','w','Position',thisPos)
    new_caxis = [0 prctile(cat(1,tmp_lpj_YXin(~isnan(tmp_lpj_YXin)),tmp_emu_YXin(~isnan(tmp_emu_YXin))),99.9)] ;
    outfile = sprintf('%s/LPJGcomp_map_%s_%s.png', ...
        fig_dir, thisEmu, thisCrop) ;
    
    nx = 2*Nn ;
    ny = 2 ;
    pltind1 = 0;
    
    for ii = 1:2
        thisIrr = irr_list{ii} ;
        for n = 1:Nn
            thisVar_lpj = sprintf('%s%s%s',thisCrop,thisIrr,N_list_lpj{n}) ;
            thisVar_emu = sprintf('%s%s%s',thisCrop,thisIrr,N_list_emu{n}) ;
            
            pltind1 = pltind1 + 1 ;
            pltind2 = pltind1 + 2*Nn ;
            
            subplot_tight(ny,nx,pltind1,spacing)
            pcolor(tmp_lpj_YXin(yrange,:,ii,n)); shading flat; axis equal tight off
            caxis(new_caxis) ; colormap(gca,'jet'); hcb = colorbar('Location','SouthOutside') ;
            title(sprintf('LPJ-GUESS: %s (N%s)', thisVar_lpj, N_list_lpj{n}))
            xlabel(hcb,'tons/ha')
            
            subplot_tight(ny,nx,pltind2,spacing)
            pcolor(tmp_emu_YXin(yrange,:,ii,n)); shading flat; axis equal tight off
            caxis(new_caxis) ; colormap(gca,'jet'); hcb = colorbar('Location','SouthOutside') ;
            title(sprintf('%s emulator: %s (N%s)', thisEmu, thisVar_emu, N_list_emu{n}))
            xlabel(hcb,'tons/ha')
        end
    end
    
    export_fig(outfile, '-r150') ;
    close
end

disp('Done')


%% Compare for all crops: Scatter

%%%% Options
spacing = [0.1 0.025] ; % v h
thisPos = [254    33   792   772] ;
alfa = 0.1 ;
fontSize = 14 ;
%%%%

for c = 1:length(crop_list)
    thisCrop = crop_list{c} ;
    figure('Color','w','Position',thisPos) ;
    ht = sgtitle(sprintf('%s (tons/ha)',thisCrop)) ;
    set(ht, 'FontSize', fontSize*1.5, 'FontWeight', 'bold')
    for ii = 1:2
        thisIrr = irr_list{ii} ;
        for n = 1:length(N_list_lpj)
            thisVar_lpj = sprintf('%s%s%s',thisCrop,thisIrr,N_list_lpj{n}) ;
            thisVar_emu = sprintf('%s%s%s',thisCrop,thisIrr,N_list_emu{n}) ;
            
            tmp_lpj = yield_lpj.garr_xv(:,strcmp(yield_lpj.varNames,thisVar_lpj)) ;
            tmp_emu = yield_emu.garr_xv(:,strcmp(yield_emu.varNames,thisVar_emu)) ;
            
            tmp_lpj_YX = nan(360,720) ;
            tmp_lpj_YX(yield_lpj.list2map) = tmp_lpj ;
            tmp_emu_YX = nan(360,720) ;
            tmp_emu_YX(yield_emu.list2map) = tmp_emu ; 
            
            ok_YX = ~isnan(tmp_lpj_YX) & ~isnan(tmp_emu_YX) ;
            tmp2_lpj = tmp_lpj_YX(ok_YX) ;
            tmp2_emu = tmp_emu_YX(ok_YX) ;
            
            % Make plot
            thisInd = (ii-1)*2 + n ;
            ha = subplot_tight(2,2,thisInd,spacing) ;
            hp = scatter(tmp2_lpj,tmp2_emu,'.b') ;
            set(hp,'MarkerFaceAlpha',alfa,'MarkerEdgeAlpha',alfa)
            
            % Equalize axes
            new_lims = [0 max(ha.XLim(2),ha.YLim(2))] ;
            set(ha,'XLim',new_lims,'YLim',new_lims) ;
            axis equal
            set(ha,'XLim',new_lims,'YLim',new_lims) ;
            
            % Add 1:1 line
            hold on
            plot(new_lims,new_lims,'--k')
            hold off
            
            % Add text
            xlabel(sprintf('LPJ-GUESS: %s', thisVar_lpj))
            ylabel(sprintf('%s emulator: %s', thisEmu, thisVar_emu))
            set(ha,'FontSize',fontSize)
            xticks = ha.XTick ;
            yticks = ha.YTick ;
            if ~isequal(xticks,yticks)
                if length(xticks) < length(yticks)
                    new_ticks = xticks ;
                    new_tickLabels = ha.XTickLabel ;
                else
                    new_ticks = yticks ;
                    new_tickLabels = ha.YTickLabel ;
                end
                set(ha,'XTick',new_ticks,'XTickLabel',new_tickLabels) ;
                set(ha,'YTick',new_ticks,'YTickLabel',new_tickLabels) ;
            end
                        
            

        end
    end
    outfile = sprintf('%s/LPJGcomp_scatter_%s_%s.png', ...
                fig_dir, thisEmu, thisCrop) ;
    export_fig(outfile, '-r150') ;
    close
end

disp('Done')

