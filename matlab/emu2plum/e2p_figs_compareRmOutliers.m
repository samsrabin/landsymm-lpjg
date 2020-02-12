%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare yields with and without outliers removed %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

which_file = 'gsirrigation' ;

in_dir_top = '/Users/Shared/GGCMI2PLUM_sh/send_to_plum/emulator_test_20191106' ;
ggcm = 'LPJmL' ;


%% Import and map

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spacing = [0.025 0.025] ; % v h
% thisPos = figurePos ;
thisPos = [1 137 1440 668] ;
fontSize = 14 ;
ylims = 60:360 ;
textpos = [120 0] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

getbasenamei = @(x) regexprep(x,'\d\d\d$','') ;

in_dir_lpj = sprintf('%s/LPJ-GUESS', in_dir_top) ;
in_dir_emu = sprintf('%s/emul_%s', in_dir_top, ggcm) ;

dirList_lpj = dir(sprintf('%s/*-*', in_dir_lpj)) ;
tpers = {dirList_lpj.name} ;
Ntpers = length(tpers) ;

for t = 1:Ntpers
    
    thisTper = tpers{t} ;
    
    lpj_with = lpjgu_matlab_read2geoArray( ...
        sprintf('%s/%s/%s_intpinfs.out', in_dir_lpj, thisTper, which_file), ...
        'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat',false) ;
    
    lpj_rmed = lpjgu_matlab_read2geoArray( ...
        sprintf('%s/%s/%s_intpinfs_rmol.out', in_dir_lpj, thisTper, which_file), ...
        'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat',false) ;
    
    emu_with = lpjgu_matlab_read2geoArray( ...
        sprintf('%s/%s/%s_intpinfs.out', in_dir_emu, thisTper, which_file), ...
        'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat',false) ;
    
    emu_rmed = lpjgu_matlab_read2geoArray( ...
        sprintf('%s/%s/%s_intpinfs_rmol.out', in_dir_emu, thisTper, which_file), ...
        'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat',false) ;
    
    if ~isequal(lpj_with.varNames, lpj_rmed.varNames) ...
            || ~isequal(lpj_with.varNames, emu_with.varNames) ...
            || ~isequal(lpj_with.varNames, emu_rmed.varNames)
        error('Mismatch of variable names')
    end
    if strcmp(which_file,'gsirrigation')
        Nvars = length(lpj_with.varNames) / 2 ;
    else
        Nvars = length(lpj_with.varNames) ;
    end
    
    n = 0 ;
    for v = 1:Nvars
        
        ii = (t-1)*Nvars + v ;
        thisVar = lpj_with.varNames{v} ;
        thisCrop = getbasenamei(thisVar) ;
        if strcmp(which_file,'gsirrigation') && ~strcmp(thisCrop(end),'i')
            continue
        end
        n = n + 1 ;
        fprintf('%s (%d/%d)...\n', thisVar, n, Ntpers*Nvars)
        
        lpj_rmed_YX = make_YX(lpj_rmed.garr_xv(:,v), lpj_rmed.list2map) ;
        lpj_with_YX = make_YX(lpj_with.garr_xv(:,v), lpj_rmed.list2map) ;
        emu_rmed_YX = make_YX(emu_rmed.garr_xv(:,v), emu_rmed.list2map) ;
        emu_with_YX = make_YX(emu_with.garr_xv(:,v), emu_with.list2map) ;
        
        figure('Color','w','Position',thisPos) ;
        
        new_caxis = [0 max(max([lpj_with_YX emu_with_YX]))] ;
        
        subplot_tight(2,2,1,spacing) ;
        pcolor(lpj_with_YX(ylims,:)); shading flat; axis equal tight off
        caxis(new_caxis)
        colorbar
        title('With outliers', 'FontSize', fontSize*1.25)
        text(textpos(2), textpos(1), 'LPJ-GUESS', ...
            'FontSize', fontSize, 'FontWeight', 'bold') ;
        
        
        subplot_tight(2,2,3,spacing) ;
        pcolor(emu_with_YX(ylims,:)); shading flat; axis equal tight off
        caxis(new_caxis)
        colorbar
        text(textpos(2), textpos(1), sprintf('%s emulator', ggcm), ...
            'FontSize', fontSize, 'FontWeight', 'bold') ;
        
        new_caxis = [0 max(max([lpj_rmed_YX emu_rmed_YX]))] ;
        
        subplot_tight(2,2,2,spacing) ;
        pcolor(lpj_rmed_YX(ylims,:)); shading flat; axis equal tight off
        caxis(new_caxis)
        colorbar
        title('Outliers removed', 'FontSize', fontSize*1.25)
        text(textpos(2), textpos(1), 'LPJ-GUESS', ...
            'FontSize', fontSize, 'FontWeight', 'bold') ;
        
        subplot_tight(2,2,4,spacing) ;
        pcolor(emu_rmed_YX(ylims,:)); shading flat; axis equal tight off
        caxis(new_caxis)
        colorbar
        text(textpos(2), textpos(1), sprintf('%s emulator', ggcm), ...
            'FontSize', fontSize, 'FontWeight', 'bold') ;
                
        % Add overall title
        a = axes;
        a.Position = [0 0 1 0.95] ;
        plot(0.5*[1 1], [0.05 0.98], '-k')
        a.XLim = [0 1];
        axis off
        t1 = title(sprintf('%s (%s)', thisVar, thisTper));
        t1.Visible = 'on'; % set(t1,'Visible','on');
        a.FontSize = fontSize*1.5 ;
                
        % Save
        outDir = sprintf('%s/compareRmOutliers/%s', in_dir_top, which_file) ;
        if ~exist(outDir, 'dir')
            s = unix(sprintf('mkdir -p %s', outDir)) ;
            if s~=0
                error('Error in mkdir -p')
            end
            clear s
        end
        out_file = sprintf('%s/%s_%s_%s.png', outDir, ggcm, thisVar, thisTper) ;
        export_fig(out_file, '-r75') ;
        
        close
        
    end
    
    stop
    
end



%% FUNCTIONS

function out_YX = make_YX(x, list2map)

out_YX = nan(360,720) ;
out_YX(list2map) = x ;


end



