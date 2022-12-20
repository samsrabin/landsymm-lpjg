thisCrop = 'spring_wheat' ;
adaptation = 1 ;
N = 200 ;
ssp = 126 ;
emuVer = 'v2.5' ;
ggcm = 'LPJmL' ;
period = 'baseline' ;


%% Setup

irrList = {'rf', 'ir'} ;

thisDir = sprintf('/Volumes/Reacher/GGCMI/CMIP_emulated/yields/CMIP6/A%d_N%d_%s/ssp%d/%s', ...
    adaptation, N, emuVer, ssp, ggcm) ;
cd(thisDir)

outDir = [thisDir '/figs'] ;
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

fileList = dir(sprintf('*_%s_*_%s_*', thisCrop, period)) ;
fileList = {fileList.name}' ;
Ngcm = length(fileList) ;
gcmList = cell(Ngcm,1) ;
for g = 1:Ngcm
    thisFile = fileList{g} ;
    parts = strsplit(thisFile, '_') ;
    gcmList{g} = parts{3} ;
end


%% Import and plot

spacing = [0.05 0.05] ; % v h
fontSize = 14 ;
thisPos = [1         139        1440         420] ;

switch thisCrop
    case 'spring_wheat'
        thisCrop_short = 'swh' ;
    otherwise
        error('thisCrop %s not recognized', thisCrop)
end

for g = 1:Ngcm
    thisFile = sprintf('%s/%s', thisDir, fileList{g}) ;
    thisGCM = gcmList{g} ;
    fprintf('%s (%d/%d)...\n', thisGCM, g, Ngcm) ;
    figure('Color', 'w', 'Position', thisPos);
    for ii = 1:2
        % Import
        thisIrr = irrList{ii} ;
        thisVar = sprintf('yield_%s_%s', thisIrr, thisCrop_short) ;
        tmp_YX = 0<flipud(permute(ncread(thisFile, thisVar), [2 1 3])) ;
        
        % Plot
        subplot_tight(1, 2, ii, spacing)
        pcolor(tmp_YX); shading flat; axis equal tight off
        title(thisIrr)
        set(gca, 'FontSize', fontSize) ;
    end
    
    % Add main title
    hold on
    hat = axes ;
    hat.Position = [0 0 1 1] ;
    ht = text(hat, 100, 100, ...
        strrep(fileList{g}, '_', '\_'), ...
        'FontSize', fontSize+4, 'FontWeight', 'bold') ;
    ht.Units = 'normalized' ;
    ht.HorizontalAlignment = 'center' ;
    ht.Position(1) = 0.5 ;
    ht.Position(2) = 0.95 ;
    hat.Visible = 'off' ;
    hold off
    
    % Save
    outFile = sprintf('%s/%s_%s_gtZero.png', outDir, thisCrop_short, period) ;
    export_fig(outFile, '-r100') ;
    close
end
disp('Done')





