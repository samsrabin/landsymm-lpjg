%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Does the GGCMI emulator reproduce the Phase II results? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% thisModel = 'LPJmL' ;
% thisModel = 'LPJ-GUESS' ;
thisModel = 'pDSSAT' ;

emulator_file = sprintf('/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_recreate_phase2/%s.out.gz', thisModel) ;
phase2_topdir = sprintf('/Volumes/WDMPP_Storage/GGCMI/AgMIP.output/%s/phase2', thisModel) ;
outdir = '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/outputs_figs/recreate_phase2_test' ;


%% Setup

cropList_long = {'winter_wheat','spring_wheat','maize'} ;
cropList = {'wwh','swh','mai'} ;
if ~strcmp(thisModel, 'LPJ-GUESS')
    cropList_long = [cropList_long {'soy','rice'}] ;
    cropList = [cropList {'soy','ric'}] ;
end
irrList_emu = {'','i'} ;
irrList_ph2 = {'0','inf'} ;
N_list = {} ;

lats_YX = repmat(0.25 + (-90:0.5:89.5)',[1 720]) ;


%% Import

% Emulator outputs
emu = lpjgu_matlab_read2geoArray(emulator_file,'force_mat_save',false,'force_mat_nosave',true) ;

% Phase II results
ph2_xv = nan(size(emu.garr_xv)) ;
for v = 1:length(emu.varNames)
    thisVar = emu.varNames{v} ;
    disp(thisVar)
    thisCropi = thisVar(1:end-3) ;
    thisCrop = thisCropi(1:3) ;
    thisCrop_long = cropList_long{find(strcmp(cropList, thisCrop))} ;
    is_irr = length(thisCropi)==4 ;
    thisN = thisVar(end-2:end) ;
    if isempty(find(strcmp(N_list, thisN)))
        N_list{end+1} = thisN ;
    end
    filename = sprintf(['%s/%s/A0/yield/', ...
        '%s_agmerra_fullharm_yield_%s_global_annual_1980_2010_', ...
        'C360_T0_W%s_N%s_A0.nc4'], ...
        phase2_topdir, thisCrop_long, lower(thisModel), thisCrop, irrList_ph2{is_irr+1}, ...
        num2str(str2double(thisN))) ;
    if ~exist(filename,'file')
        error('File not found: %s', filename)
    end
    tmp_YX = flipud(nanmean(ncread(filename, sprintf('yield_%s', thisCrop)),3)') ;
    ph2_xv(:,v) = tmp_YX(emu.list2map) ;
    clear tmp_YX
end
disp('Done')


%% Compare: Scatter

%%%% Options
spacing = [0.1 0.075] ; % v h
thisPos = [254    33   882   772] ;
alfa = 0.1 ;
fontSize = 14 ;
%%%%

for c = 1:length(cropList)
    thisCrop = cropList{c} ;
    figure('Color','w','Position',thisPos) ;
    ht = sgtitle(sprintf('%s: %s',thisModel,strrep(cropList_long{c},'_','\_'))) ;
    set(ht, 'FontSize', fontSize*1.5, 'FontWeight', 'bold')
    for ii = 1:2
        thisIrr = irrList_emu{ii} ;
        
        for N = 1:length(N_list)
            thisN = N_list{N} ;
            thisVar = [thisCrop thisIrr thisN] ;
            thisIndex = find(strcmp(emu.varNames,thisVar)) ;
            
            emu_YX = nan(360,720) ;
            ph2_YX = nan(360,720) ;
            emu_YX(emu.list2map) = emu.garr_xv(:,thisIndex) ;
            ph2_YX(emu.list2map) = ph2_xv(:,thisIndex) ;
            
            ok_YX = ~isnan(emu_YX) & ~isnan(ph2_YX) ;
            tmp2_emu = emu_YX(ok_YX) ;
            tmp2_ph2 = ph2_YX(ok_YX) ;
            
            % Make plot
            thisInd = (ii-1)*3 + N ;
            ha = subplot_tight(2,3,thisInd,spacing) ;
            hp = scatter(tmp2_ph2,tmp2_emu,'.b') ;
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
            title(thisVar)
            title(sprintf('%s%s',thisIrr, thisN))
            xlabel('Phase II')
            ylabel('Emulator')
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
    filename = sprintf('%s/%s_%s_scatter.png', outdir, thisModel, thisCrop) ;
    export_fig(filename, '-r150')
    close
end


%% Compare: Maps

spacing = [0.025 0.025] ;
% thisPos = [1    34   628   771] ;
thisPos = figurePos ;
hs = [] ;
min_yield = 1 ; % t/ha
use_pctDiff = true ;

emu_YX_old = [] ;
ph2_YX_old = [] ;
for c = 1:length(cropList)
    thisCrop = cropList{c} ;
    for ii = 1:2
        thisIrr = irrList_emu{ii} ;
        fprintf('%s%s (%d of %d)...\n', thisCrop, thisIrr, (c-1)*2 + ii, length(cropList)*2) ;
        figure('Color','w','Position',thisPos) ;
        for N = 1:length(N_list)
            thisN = N_list{N} ;
            thisVar = [thisCrop thisIrr thisN] ;
            thisIndex = find(strcmp(emu.varNames,thisVar)) ;
            
            emu_YX = nan(360,720) ;
            ph2_YX = nan(360,720) ;
            emu_YX(emu.list2map) = emu.garr_xv(:,thisIndex) ;
            ph2_YX(emu.list2map) = ph2_xv(:,thisIndex) ;
            
            if isequaln(emu_YX,emu_YX_old)
                stop
            elseif isequaln(ph2_YX,ph2_YX_old)
                stop
            end
            
            new_caxis = [0 max(max(ph2_YX(:)), max(emu_YX(:)))] ;

            hs(end+1) = subplot_tight(3,3,N*3-2,spacing) ;
            pcolor(ph2_YX(30:end,:)) ;
            shading flat; axis equal tight off
            caxis(new_caxis)
            colormap(gca,'jet')
            colorbar
            title(sprintf('%s: Phase II', thisVar))
            set(gca,'FontSize',14)
            
            hs(end+1) = subplot_tight(3,3,N*3-1,spacing) ;
            pcolor(emu_YX(30:end,:)) ;
            shading flat; axis equal tight off
            caxis(new_caxis)
            colormap(gca,'jet')
            colorbar
            title(sprintf('%s: Emulated', thisVar))
            set(gca,'FontSize',14)
            
            subplot_tight(3,3,N*3,spacing) ;
            if use_pctDiff
                diff_YX = 100*(emu_YX - ph2_YX) ./ ph2_YX ;
                diff_YX(ph2_YX<min_yield & emu_YX<min_yield) = NaN ;
                pcolor(diff_YX(30:end,:)) ;
            else
                diff_YX = emu_YX - ph2_YX ;
%                 diff_YX(ph2_YX<3) = NaN ;
            end
            pcolor(diff_YX(30:end,:)) ;
            shading flat; axis equal tight off
            colormap(gca,'jet')
            if use_pctDiff
                caxis(100*[-1 1])
            else
%                 caxis(max(abs(caxis))*[-1 1])
                caxis([-1 1])
                colormap(gca,brewermap(64,'RdBu_ssr'))
            end
            hcb = colorbar ;
            if any(any(diff_YX > max(caxis)))
                hcb.TickLabels{end} = ['\geq ' hcb.TickLabels{end}] ;
            end
            if any(any(diff_YX < min(caxis)))
                hcb.TickLabels{1} = ['\leq ' hcb.TickLabels{1}] ;
            end
            title(sprintf('%s: emu. - Ph.2 (t/ha)', thisVar))
            set(gca,'FontSize',14)
            
            emu_YX_old = emu_YX ;
            ph2_YX_old = ph2_YX ;
        end
        pause(0.1)
        if use_pctDiff
            filename = sprintf('%s/%s_%s%s_pctDiff.png', outdir, thisModel, thisCrop, thisIrr) ;
        else
            filename = sprintf('%s/%s_%s%s_diff.png', outdir, thisModel, thisCrop, thisIrr) ; 
        end
        export_fig(filename, '-r150')
        close
    end
end

disp('Done')


%% Maps of Phase II results

spacing = [0.025 0.025] ;
% thisPos = [1    34   628   771] ;
thisPos = [1         373        1440         432] ;
hs = [] ;
min_yield = 1 ; % t/ha

emu_YX_old = [] ;
ph2_YX_old = [] ;
for c = 1:length(cropList)
    thisCrop = cropList{c} ;
    figure('Color','w','Position',thisPos) ;
    for ii = 1:2
        thisIrr = irrList_emu{ii} ;
        for N = 1:length(N_list)
            thisN = N_list{N} ;
            thisVar = [thisCrop thisIrr thisN] ;
            thisIndex = find(strcmp(emu.varNames,thisVar)) ;
            
            ph2_YX = nan(360,720) ;
            ph2_YX(emu.list2map) = ph2_xv(:,thisIndex) ;
            
            hs(end+1) = subplot_tight(2,3,(ii-1)*3+N,spacing) ;
            pcolor(ph2_YX(30:end,:)) ;
            shading flat; axis equal tight off
            colormap(gca,'jet')
            colorbar
            title(sprintf('%s: Phase II', thisVar))
            set(gca,'FontSize',14)
                        
        end
    end
    pause(0.1)
    filename = sprintf('%s/%s_%s_phase2.png', outdir, thisModel, thisCrop) ;
    export_fig(filename, '-r150')
    close
end

disp('Done.')



%% Spring wheat vs. rice?

spacing = [0.025 0.025] ;
var_swh = 'swh200' ;
var_ric = strrep(var_swh,'swh','ric') ;
% useRice_cond = 'abs(lats_YX)<23.5' ;
useRice_cond = 'ric_YX>swh_YX' ;

swh = emu.garr_xv(:,strcmp(emu.varNames,var_swh)) ;
ric = emu.garr_xv(:,strcmp(emu.varNames,var_ric)) ;

% Normalize
swh = swh / mean(swh) ;
ric = ric / mean(ric) ;

new_caxis = [0 prctile(cat(1,swh,ric),99.9)] ;

swh_YX = nan(360,720) ;
ric_YX = nan(360,720) ;
swh_YX(emu.list2map) = swh ;
ric_YX(emu.list2map) = ric ;

eval(sprintf('useRice_YX = %s ;', useRice_cond)) ;

tmp_YX = swh_YX ;
tmp_YX(useRice_YX) = ric_YX(useRice_YX) ;

useRice_YX = double(useRice_YX) ;
useRice_YX(isnan(tmp_YX)) = NaN ;

figure('Color','w','Position',thisPos) ;

subplot_tight(2,2,1,spacing)
pcolor(swh_YX); shading flat; axis equal tight off
colormap(gca,'jet'); caxis(new_caxis); colorbar
title(var_swh)

subplot_tight(2,2,2,spacing)
pcolor(ric_YX); shading flat; axis equal tight off
colormap(gca,'jet'); caxis(new_caxis); colorbar
title(var_ric)

subplot_tight(2,2,3,spacing)
pcolor(tmp_YX); shading flat; axis equal tight off
colormap(gca,'jet'); caxis(new_caxis); colorbar
title('Combined')

subplot_tight(2,2,4,spacing)
pcolor(useRice_YX); shading flat; axis equal tight off
colorbar
title(sprintf('Use rice if %s',strrep(useRice_cond,'_','\_')))







