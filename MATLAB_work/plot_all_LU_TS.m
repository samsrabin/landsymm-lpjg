%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot range of LU trajectories for each SSP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_save = true ;

%% Import

topDir = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/ssp12' ;
SSPlist = [1 3 4 5] ;
Nssp = length(SSPlist) ;

for s = 1:Nssp
    thisSSP = SSPlist(s) ;
    thisFile = sprintf('%s/SSP%d/LU_timeseries.mat', topDir, thisSSP) ;
    tmp = load(thisFile) ;
    if s==1
        LUts_yvds = nan([size(tmp.LUts_yvd) Nssp]) ;
        LUnames = tmp.LUnames ;
        yearList = tmp.yearList ;
    else
        if ~strcmp(LUnames, tmp.LUnames)
            error('LUnames mismatch')
        elseif ~isequal(yearList, tmp.yearList)
            error('yearList mismatch')
        end
    end
    
    [badYr, badDir] = find(squeeze(sum(tmp.LUts_yvd,2))<1.2616e14) ;
    if ~isempty(badYr)
%         warnMsg = '' ;
        warnMsg = sprintf('\nSSP%d', thisSSP) ;
        for i = 1:length(badYr)
            tmpMsg = sprintf('\n s%d %d\n', badDir(i), yearList(badYr(i))) ;
            warnMsg = [warnMsg tmpMsg] ;
            tmp.LUts_yvd(badYr(i),:,badDir(i)) = NaN ;
        end
        warnMsg = ['Something went wrong in the following (setting to NaN):' warnMsg] ;
        warning(warnMsg)
    end

    
    LUts_yvds(:,:,:,s) = tmp.LUts_yvd ;
    clear tmp
end


%% Spaghetti plot

thisLU = 'NATURAL' ;

cmap_lines = colormap('lines') ;
figure('Color','w') ;

for s = 1:Nssp
    hold on
    thisColor = cmap_lines(s,:) ;
    h = plot(yearList, squeeze(LUts_yvds(:,strcmp(LUnames,thisLU),:,s)), ...
        'Color', thisColor, ...
        'LineWidth', 3) ;
    for i = 1:length(h)
        h(i).Color = [h(i).Color 0.1] ;
    end
    hold off
end


%% SEM shading

% Options %%%%%%%%%%%%
thisLU = 'NATURAL' ;
fontSize = 20 ;
convFact = 1e-12 ;   % m2 to Mkm2
thisPos = [294   187   867   422] ;
%%%%%%%%%%%%%%%%%%%%%%

cmap_lines = colormap('lines') ; close
figure('Color', 'w', 'Position', thisPos) ;

hls = {} ;
for s = 1:Nssp
    hold on
    thisColor = cmap_lines(s,:) ;
    data_yd = convFact*squeeze(LUts_yvds(:,strcmp(LUnames,thisLU),:,s)) ;
    data_mean = nanmean(data_yd,2) ;
    data_unc = std(data_yd,0,2,'omitnan') ;
%     data_unc = 1.96 * std(data_yd,0,2,'omitnan') / sqrt(size(data_yd,2)) ;
    [hl, hp] = boundedline(yearList, data_mean, data_unc, ...
        'cmap', thisColor) ;
    hl.LineWidth = 3 ;
    hp.FaceAlpha = 0.3 ;
    hout = outlinebounds(hl, hp) ;
    hold off
    hls = [hls hl] ;
end

set(gca, ...
    'FontSize', fontSize, ...
    'XLim', yearList([1 end]))
legend(hls, strcat('SSP', strsplit(num2str(SSPlist))),'Location','Southwest')
title('Natural land area')
ylabel('Million km^2')

if do_save
    export_fig(sprintf('%s/LU_ts.png', topDir),'-r300')
end


%% Just means

% Options %%%%%%%%%%%%
thisLU = 'NATURAL' ;
fontSize = 20 ;
convFact = 1e-12 ;   % m2 to Mkm2
thisPos = [294   187   867   422] ;
%%%%%%%%%%%%%%%%%%%%%%

cmap_lines = colormap('lines') ; close

figure('Color', 'w', 'Position', thisPos) ;
data_yr = convFact*squeeze(nanmean(LUts_yvds(:,strcmp(LUnames,thisLU),:,:),3)) ;
plot(yearList, data_yr, 'LineWidth', 3)

set(gca, ...
    'FontSize', fontSize, ...
    'XLim', yearList([1 end]), ...
    'YLim', [62 76])
legend(strcat('SSP', strsplit(num2str(SSPlist))),'Location','Southwest')
title('Natural land area')
ylabel('Million km^2')

if do_save
    export_fig(sprintf('%s/LU_tsjustMeans.png', topDir),'-r300')
end
