%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Playing around with GGCMI emulator outputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inDir = 'outputs_20190802' ;
% inDir = 'outputs_20190930' ;
% inDir = 'outputs_GGCMIcrops_20191008' ;
% inDir = 'outputs_GGCMIcrops_20200204' ;
inDir = 'outputs_GGCMIcrops_20200310' ;
thisVar = 'yield' ;
% thisVar = 'gsirrigation' ;


%% Setup

cd '/Users/Shared/GGCMI2PLUM_sh/emulation/outputs' ;

% Get GCMs, RCPs, GGCMs, time periods
list_gcms = dir_noInvis(inDir) ;
list_gcms = {list_gcms.name} ;
list_gcms(strcmp(list_gcms,'figures')) = [] ;
Ngcms = length(list_gcms) ;
list_ggcms = dir_noInvis(sprintf('%s/%s',inDir, list_gcms{1})) ; 
list_ggcms = {list_ggcms.name} ; Nggcms = length(list_ggcms) ;
list_rcps = dir_noInvis(sprintf('%s/%s/%s',inDir, list_gcms{1}, list_ggcms{1})) ; 
list_rcps = {list_rcps.name} ; Nrcps = length(list_rcps) ;
list_tpers = dir_noInvis(sprintf('%s/%s/%s/%s',inDir, list_gcms{1}, list_ggcms{1}, list_rcps{1})) ; 
list_tpers = {list_tpers.name} ; Ntpers = length(list_tpers) ;

% Get crops and N treatments
list_cropLists = cell(Nggcms,1) ;
list_varNames = cell(Nggcms,1) ;
list_allCrops = {} ;
for j = 1:Nggcms
    thisFile = sprintf('%s/%s/%s/%s/%s/%s.out', ...
        inDir, list_gcms{1}, list_ggcms{j}, list_rcps{1}, list_tpers{1}, thisVar) ;
    S = lpjgu_matlab_read2geoArray(thisFile, 'verboseIfNoMat', false) ;
    list_cropLists{j} = unique(regexprep(S.varNames','\d\d\d$','')) ;
    if j==1
        getlast3 = @(x) x(end-2:end) ;
        list_n = unique(cellfun(getlast3, S.varNames', 'UniformOutput', false)) ;
        Nn = length(list_n) ;
    end
    list_varNames{j} = S.varNames ;
    list_allCrops = unique([list_allCrops;list_cropLists{j}]) ;
    clear S
end

outDir = sprintf('/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/outputs_figs/%s', inDir) ;
if ~exist(outDir,'dir')
    mkdir(outDir)
end

firstCap = @(x) [upper(x(1)) x(2:end)] ;
list_ops = {'max','min','mean'} ;

list_gcms
list_ggcms
list_rcps
list_tpers
list_allCrops
list_n

if strcmp(thisVar,'yield')
    units = 't/ha' ;
elseif strcmp(thisVar,'gsirrigation')
    units = 'mm' ;
else
    error('thisVar (%s) not recognized for units', thisVar) ;
end


%% See max/mean/etc. projected yield for a given GCM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
just_summaries = true ;
spacing = [0.035 0.010] ; % v h
ny = 3 ;
nx = 4 ;
ylims = 65:360 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gcm = 1:Ngcms
    thisGcm = list_gcms{gcm} ;
    for ggcm = 1:Nggcms
        thisGgcm = list_ggcms{ggcm} ;
        list_crops = list_cropLists{ggcm} ;
        varNames = list_varNames{ggcm} ;
        Ncrops = length(list_crops) ;
        if ny*nx < Ncrops
            error('ny*nx (%d) < Ncrops (%d)', ny*nx, Ncrops)
        end
        for rcp = 1:Nrcps
            thisRcp = list_rcps{rcp} ;
            fprintf('%s: %s %s\n',thisGgcm, thisGcm, thisRcp)
            
            % Read data
            for t = 1:Ntpers
                thisTper = list_tpers{t} ;
                fprintf('    %s... ', thisTper)
                thisFile = sprintf('%s/%s/%s/%s/%s/%s.out', ...
                    inDir, thisGcm, thisGgcm, thisRcp, thisTper, thisVar) ;
                if rcp==1 && t==1
                    S = lpjgu_matlab_read2geoArray(thisFile, ...
                        'verboseIfNoMat', false) ;
                    target = {S.lonlats, S.list2map} ;
                else
                    S = lpjgu_matlab_read2geoArray(thisFile, ...
                        'target', target, 'verboseIfNoMat', false) ;
                end
                
                if strcmp(thisVar,'yield')
                    % Convert from kg/m2 to tons/ha
                    S.garr_xv = 10 * S.garr_xv ;
                end
                
                % Make figure for this time period
                if just_summaries
                    fprintf('\n')
                else
                    for n = 1:Nn
                        thisN = list_n{n} ;
                        make_fig(S.garr_xv, target, list_crops, list_allCrops, varNames, ...
                            figurePos, outDir, ...
                            thisGgcm, thisGcm, thisRcp, thisVar, thisTper, thisN, ...
                            nx, ny, spacing, ylims, units)
                    end
                end
                
                if t==1
                    tmp_garr_gvt = nan([size(S.garr_xv) Ntpers]) ;
                end
                tmp_garr_gvt(:,:,t) = S.garr_xv ;
                clear S
            end
                        
            % Make summary-statistic figures
            for o = 1:length(list_ops)
                thisOperation = list_ops{o} ;
                if strcmp(thisOperation,'max')
                    tmp_garr_gv = max(tmp_garr_gvt,[],3) ;
                elseif strcmp(thisOperation,'min')
                    tmp_garr_gv = min(tmp_garr_gvt,[],3) ;
                elseif strcmp(thisOperation,'mean')
                    tmp_garr_gv = mean(tmp_garr_gvt,3) ;
                else
                    error('\nthisOperation (%s) not recognized', thisOperation)
                end
                for n = 1:Nn
                    thisN = list_n{n} ;
                    fprintf('    %s %s... ', thisOperation, thisN)
                    make_fig(tmp_garr_gv, target, list_crops, list_allCrops, varNames, ...
                        figurePos, outDir, ...
                        thisGgcm, thisGcm, thisRcp, thisVar, thisOperation, thisN, ...
                        nx, ny, spacing, ylims, units)
                end
            end
                        
        end
    end
end
disp('Done!')


%% FUNCTIONS

function S = dir_noInvis(char_in)

S = dir(char_in) ;
names = {S.name} ;
Nfound = length(names) ;
bad = false(Nfound,1) ;
for f = 1:Nfound
    if strcmp(S(f).name(1),'.')
        bad(f) = true ;
    end
end
S(bad) = [] ;


end

function make_fig(tmp_garr_gv, target, list_crops, list_allCrops, varNames, ...
    figurePos, outDir, ...
    thisGgcm, thisGcm, thisRcp, thisVar, thisTp, thisN, ...
    nx, ny, spacing, ylims, units)

out_file = sprintf('%s/%s_%s_%s_N%s_%s_%s.png', outDir, thisGgcm, thisGcm, thisRcp, thisN, thisVar, thisTp) ;
if exist(out_file,'file')
    fprintf('Skipping (output file exists).\n')
    return
end

this_sgtitle = sprintf('%s N%s: %s %s (%s)', thisGgcm, thisN, thisGcm, thisRcp, thisTp) ;

figure('Position',figurePos,'Color','w')
clear hsp_prev thisMap_YX_prev
for c = 1:length(list_crops)
    
    thisCrop = list_crops{c} ;
    thisFigIndex = find(strcmp(list_allCrops,thisCrop)) ;
    if length(thisFigIndex) ~= 1
        error('\nError finding thisFigIndex')
    end
    thisDataIndex = find(strcmp(varNames, ...
        sprintf('%s%s', thisCrop, thisN))) ;
    if length(thisDataIndex) ~= 1
        error('\nError finding thisDataIndex')
    end
    
    hsp = subplot_tight(ny, nx, thisFigIndex, spacing) ;
    thisMap_YX = nan(360,720) ;
    thisMap_YX(target{2}) = tmp_garr_gv(:,thisDataIndex) ; %#ok<FNDSB>
    thisPrctile = prctile(thisMap_YX(~isnan(thisMap_YX)),99) ;
    thisMap_YX(thisMap_YX>thisPrctile) = thisPrctile ;
    pcolor(thisMap_YX(ylims,:)) ; shading flat; axis equal tight off
    title(sprintf('%s N%s', strrep(thisCrop,'_','\_'), thisN))
    hcb = colorbar('Location','SouthOutside') ;
    xlabel(hcb, units)
    
    % Equalize colorbars of rainfed/irrigated
    if rem(c,2)==0
        tmp = cat(3,thisMap_YX_prev,thisMap_YX) ;
        new_caxis = minmax_ssr(tmp) ;
        caxis(hsp_prev, new_caxis)
        caxis(hsp, new_caxis)
    end
    hsp_prev = hsp ;
    thisMap_YX_prev = thisMap_YX ;
end
hold on
a = axes() ;
ht = title(a,this_sgtitle) ;
a.Visible = 'off' ;
ht.Visible = 'on' ;
ht.Position(2) = 1.06 ;
ht.FontSize = 16 ;
ht.FontWeight = 'bold' ;
hold off
export_fig(out_file, '-r75')
close

fprintf('\n')

end

