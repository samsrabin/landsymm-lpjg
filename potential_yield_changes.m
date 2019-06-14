%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLUM2LPJG figures: Multiple runs: %%%
%%% Potential yield changes%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisVer = 'remapv6p7' ;

incl_N = {'0'} ;


%% Setup

if strcmp(thisVer, 'remapv6p7')
    dir_list = { ...
        '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_2001-2100_remap6p7_forPotYields_rcp45/LPJGPLUM_2001-2100_remap6p7_forPotYields_rcp45_forED_20190225101539' ;
        '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_2001-2100_remap6p7_forPotYields_rcp60/LPJGPLUM_2001-2100_remap6p7_forPotYields_rcp60_forED_20190225101712' ;
        '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_2001-2100_remap6p7_forPotYields_rcp85/LPJGPLUM_2001-2100_remap6p7_forPotYields_rcp85_forED_20190225121232' ;
        } ;
    rcp_list = {'RCP4.5', 'RCP6.0', 'RCP8.5'} ;
else
    error('thisVer (%s) not recognized', thisVer) ;
end

Nrcp = length(dir_list) ;
if Nrcp ~= length(rcp_list)
    error('Mismatch in length between dir_list and rcp_list')
end

dir_out = sprintf('%s/figures_yieldDiff_%s/', ...
    '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam', ...
    thisVer) ;
if ~exist(dir_out,'dir')
    mkdir(dir_out)
end


%% Import

for r = 1:Nrcp
    
    thisTopDir = dir_list{r} ;
    yearDirs = dir(sprintf('%s/*-*',thisTopDir)) ;
    if r==1
        Npds = length(yearDirs) ;
        yearDirs_orig = yearDirs ;
        tmp = strsplit(yearDirs(1).name,'-') ;
        pdLength = str2double(tmp{2}) - str2double(tmp{1}) + 1 ;
        tmp2 = str2double(strsplit(strjoin({yearDirs_orig.name},'-'),'-')) ;
        yearList_1 = tmp2(1:2:end) ;
        yearList_N = tmp2(2:2:end) ;
        clear tmp*
    else
        if ~isequal({yearDirs.name}, {yearDirs_orig.name})
            error('Mismatch between this directory''s year list and original')
        end
    end
    
    for p = 1:Npds
        fprintf('Importing RCP %d/%d, period %d/%d...\n', r, Nrcp, p, Npds) ;
        file_in = sprintf('%s/%s/yield.out', thisTopDir, yearDirs(p).name) ;
        struct_in = lpjgu_matlab_readTable_then2map(file_in, 'verboseIfNoMat', false) ;
        incl = contains(struct_in.varNames, incl_N) ;        
        if r==1 && p==1
            % Get list of all CFTs (including, e.g., i1000)
            crop_list_all = struct_in.varNames ;
            crop_list_all(~incl) = [] ;
            Ncrops_all = length(crop_list_all) ;
            
            % Get crop names
            crop_list = crop_list_all(cellfun(@isempty,regexp(crop_list_all,sprintf('.*i%s', incl_N{1})))) ;
            crop_list = strrep(crop_list, incl_N{1}, '') ;
            Ncrops = length(crop_list) ;
            
            % Make empty array for yields
            yield_YXcpr = nan([size(struct_in.maps_YXv,1) size(struct_in.maps_YXv,2) Ncrops_all Npds Nrcp]) ;
        elseif ~isequal(crop_list_all, struct_in.varNames(incl))
            error('Mismatch in crop list')
        end
        yield_YXcpr(:,:,:,p,r) = struct_in.maps_YXv(:,:,incl) ;
        
        clear struct_in
    end
    clear yearDirs

end
disp('Done.')


%% Make figures: Potential yield in a given period

% Options %%%%%%
spacing = [0.05 0.05] ; % v h
thisPos = [1    33   638   772] ;
ylims = 69:360 ;
fontSize = 14 ;
this_colormap = 'jet' ;
% yrs = 2041:2050 ; % What PLUM considers for 2046-2055
yrs = 2086:2095 ; % What PLUM considers for 2091-2100
%%%%%%%%%%%%%%%%

% Check for valid year inputs
y1 = yrs(1) ;
yN = yrs(end) ;
if ~any(yearList_1==y1)
    error('Code can''t handle start year (%d) that''s not in yearList_1', y1)
elseif ~any(yearList_N==yN)
    error('Code can''t handle end year (%d) that''s not in yearList_N', yN)
end
yy1 = find(yearList_1==y1) ;
yyN = find(yearList_N==yN) ;

ir_list = {'','i'} ;
for n = 1:length(incl_N)
    thisN = incl_N{n} ;
    for c = 1:Ncrops
        thisCrop = crop_list{c} ;
        
        for i = 1:2
            
            % Get index of this CFT
            thisCFT = [thisCrop ir_list{i} thisN] ;
            ii = find(strcmp(crop_list_all, thisCFT)) ;
            if isempty(ii)
                error('%s not found in crop_list_all', thisCFT) ;
            elseif length(ii) > 1
                error('Multiple matches for %s found in crop_list_all', thisCFT) ;
            end
            
            maps_YXr = 10*squeeze(mean(yield_YXcpr(:,:,ii,yy1:yyN,:),4)) ;
            clims = [0 1] * max(max(max(abs(maps_YXr)))) ;
            figure('Color','w','Position',thisPos) ;
            for r = 1:Nrcp
                thisRCP = rcp_list{r} ;
                subplot_tight(Nrcp,1,r,spacing) ;
                pcolor(maps_YXr(ylims,:,r)) ;
                shading flat; axis equal tight off
                caxis(clims) ;
                colormap(this_colormap)
                hcb = colorbar ;
                title(hcb,'\Delta tons ha^{-1}') ;
                title(sprintf('%s: %s (%d-%d)', thisCFT, thisRCP, y1, yN)) ;
                set(gca,'FontSize',fontSize)
            end
            file_out = sprintf('%s/%s_yield_%d-%d.png', dir_out, thisCFT, y1, yN) ;
            export_fig(file_out, '-r300') ;
            close ;
        end
        
    end
end



%% Make figures: Check differences in baseline yield

% Options %%%%%%
spacing = [0.05 0.05] ; % v h
ylims = 69:360 ;
fontSize = 14 ;
%%%%%%%%%%%%%%%%

figure('Color','w','Position',figurePos) ;
hf1 = gcf ;

for n = 1:length(incl_N)
    thisN = incl_N{n} ;
    for c = 1:Ncrops
        thisCrop = crop_list{c} ;
        
        thisCFT = [thisCrop thisN] ;
        ii = find(strcmp(crop_list_all, thisCFT)) ;
        if isempty(ii)
            error('%s not found in crop_list_all', thisCFT) ;
        elseif length(ii) > 1
            error('Multiple matches for %s found in crop_list_all', thisCFT) ;
        end
        maps_YXr = squeeze(yield_YXcpr(:,:,ii,1,:)) ;
        clims = [0 max(max(max(maps_YXr)))] ;
        
        figure(hf1) ;
        subplot_tight(3,2,c,spacing) ;
        map_YX = max(maps_YXr,[],3) - min(maps_YXr,[],3) ;
        if any(map_YX>0 & min(maps_YXr,[],3)==0)
            error('This doesn''t work when any min yield = 0') ;
        end
        map_YX = 100 * map_YX ./ min(maps_YXr,[],3) ;
        pcolor(map_YX(ylims,:)) ;
        shading flat; axis equal tight off
        caxis(clims) ;
        colorbar
        title(sprintf('%s: max-min yield (%% of min)', thisCFT)) ;
        set(gca,'FontSize',fontSize)
        
%         figure(2, 'Color','w','Position',thisPos) ;
%         for r = 1:Nrcp
%             thisRCP = rcp_list{r} ;
%             subplot_tight(Nrcp,1,r,spacing) ;
%             pcolor(maps_YXr(ylims,:,r)) ;
%             shading flat; axis equal tight off
%             caxis(clims) ;
%             colorbar
%             title(sprintf('%s: %s', thisCFT, thisRCP)) ;
%             set(gca,'FontSize',fontSize)
%         end
%         close(2) ;
        
        
        
    end
end



%% Make figures: Change in yield from baseline

% Options %%%%%%
spacing = [0.05 0.05] ; % v h
thisPos = [1    33   638   772] ;
ylims = 69:360 ;
fontSize = 14 ;
this_colormap = 'brbg_ssr' ;
%%%%%%%%%%%%%%%%

ir_list = {'','i'} ;

for n = 1:length(incl_N)
    thisN = incl_N{n} ;
    for c = 1:Ncrops
        thisCrop = crop_list{c} ;
        
        for i = 1:2
            
            thisCFT = [thisCrop ir_list{i} thisN] ;
            ii = find(strcmp(crop_list_all, thisCFT)) ;
            if isempty(ii)
                error('%s not found in crop_list_all', thisCFT) ;
            elseif length(ii) > 1
                error('Multiple matches for %s found in crop_list_all', thisCFT) ;
            end
            maps_YXr = 10*squeeze(yield_YXcpr(:,:,ii,end,:) - yield_YXcpr(:,:,ii,1,:)) ;
            clims = [-1 1] * max(max(max(abs(maps_YXr)))) ;
            figure('Color','w','Position',thisPos) ;
            for r = 1:Nrcp
                thisRCP = rcp_list{r} ;
                subplot_tight(Nrcp,1,r,spacing) ;
                pcolor(maps_YXr(ylims,:,r)) ;
                shading flat; axis equal tight off
                caxis(clims) ;
                colormap(brewermap(64,this_colormap))
                hcb = colorbar ;
                title(hcb,'\Delta tons ha^{-1}') ;
                title(sprintf('%s: %s', thisCFT, thisRCP)) ;
                set(gca,'FontSize',fontSize)
            end
            file_out = sprintf('%s/%s_yieldDiff_%s_to_%s.png', dir_out, thisCFT, yearDirs(1).name, yearDirs(end).name) ;
            export_fig(file_out, '-r300') ;
            close ;
        end
        
    end
end







