%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make figures showing N sensitivity %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ggcm_list = {'LPJ-GUESS', 'LPJmL', 'pDSSAT', 'EPIC-TAMU'} ;
% ggcm_list = {'EPIC-TAMU', 'LPJmL', 'pDSSAT'} ;
% ggcm_list = {'LPJ-GUESS'} ;
gcm = 'IPSL-CM5A-MR_r1i1p1' ;
rcp = 'rcp45' ;
thisVer = '20200310' ;

% Behaviors
% excl_lowBL_agmerra = true ;
% excl_lowBL_emu = true ;
% interp_infs = true ;
remove_outliers = true ;


%% Setup

lu_file = '/Users/Shared/PLUM/input/remaps_v6p7/LU.remapv6p7.2010.txt.mat' ;
cf_file = '/Users/Shared/PLUM/input/remaps_v6p7/cropfracs.remapv6p7.2010.txt.mat' ;
landarea_file = '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc' ;

Nggcm = length(ggcm_list) ;

inDir = sprintf('/Users/Shared/GGCMI2PLUM_sh/send_to_plum/%s_%s_v%s', gcm, rcp, thisVer) ;
if remove_outliers
    inDir = [inDir '_rmol'] ;
end

which_file = 'yield' ;

getbasename = @(x) regexprep(x,'i?\d\d\d$','') ;
getbasenamei = @(x) regexprep(x,'\d\d\d$','') ;
getN = @(x) x(end-2:end) ;


%% Import and process land use areas

% Land use
bl_lu = lpjgu_matlab_readTable_then2map(lu_file) ;

% Crop fractions
bl_cf = lpjgu_matlab_readTable_then2map(cf_file) ;
bl_cf.maps_YXv(:,:,contains(bl_cf.varNames, {'ExtraCrop', 'Miscanthus'})) = [] ;
bl_cf.varNames(contains(bl_cf.varNames, {'ExtraCrop', 'Miscanthus'})) = [] ;

% Gridcell area
gcelArea_YXqd = 1e6*double(transpose(ncread(landarea_file,'carea'))) ;
land_frac_YXqd = 1 - double(flipud(transpose(ncread(landarea_file,'icwtr')))) ;
landArea_YXqd = gcelArea_YXqd .* land_frac_YXqd ;
tmp = gcelArea_YXqd(:,1:2:1440) + gcelArea_YXqd(:,2:2:1440) ;
gcelArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;

% Crop area
croparea_m2_YXv = bl_cf.maps_YXv ...
    .* repmat(bl_lu.maps_YXv(:,:,strcmp(bl_lu.varNames, 'CROPLAND')), [1 1 size(bl_cf.maps_YXv,3)]) ...
    .* repmat(gcelArea_YX, [1 1 size(bl_cf.maps_YXv,3)]) ;
clear *_YXqd tmp


%% Get results for each model
NvarsI = 0 ;
NcropsI = 0 ;
for m = 1:Nggcm
    
    % Which GGCM?
    ggcm = ggcm_list{m} ;
    disp(ggcm)
    
    % Import
    inDir_ggcm = sprintf('%s/emu_%s', inDir, ggcm) ;
    in_file = sprintf('%s/future_%s.mat', inDir_ggcm, which_file) ;
    if m==1
        tmp = load(in_file) ;
        croparea_m2_xv = lpjgu_matlab_maps2xvy(croparea_m2_YXv, tmp.data_fu_lpj.list2map) ;
        data_fu_lpj = tmp.data_fu_lpj ;
        data_fu_emu = tmp.data_fu_emu ;
        precision = 'double' ;
        data_fu_emu = rmfield(data_fu_emu, 'garr_xvt') ;
        data_fu_emu.garr_xvtm = nan([size(tmp.data_fu_emu.garr_xvt) Nggcm], precision) ;
        data_fu_out = tmp.data_fu_out ;
        data_fu_out = rmfield(data_fu_out, 'garr_xvt') ;
        data_fu_out.garr_xvtm = nan([size(tmp.data_fu_out.garr_xvt) Nggcm], precision) ;
        Ntpers = size(data_fu_out.garr_xvtm,3) ;
        [varNames_lpj, cropList_lpj] = ...
            e2p_get_names([], data_fu_lpj.varNames, ...
            getbasename, getbasenamei, getN) ;
        cropList_lpj_asEmu_mc = cell(Nggcm, length(cropList_lpj)) ;
    else
        tmp = load(in_file, 'data_fu_emu', 'data_fu_out') ;
    end
    data_fu_emu.garr_xvtm(:,1:size(tmp.data_fu_emu.garr_xvt,2),:,m) = tmp.data_fu_emu.garr_xvt ;
    data_fu_out.garr_xvtm(:,:,:,m) = tmp.data_fu_out.garr_xvt ;
    [varNames_emu, cropList_emu] = ...
        e2p_get_names([], tmp.data_fu_emu.varNames, ...
        getbasename, getbasenamei, getN) ;
    NvarsI = max(NvarsI, length(varNames_emu)) ;
    
    % Update cell array of varNames
    if m==1
        varNames_emu_mc = varNames_emu ;
    elseif NvarsI > size(varNames_emu_mc,2)
        tmp = cell(size(varNames_emu_mc,1), NvarsI-size(varNames_emu_mc,2)) ;
        tmp(:) = {'NA'} ;
        varNames_emu_mc = cat(2, varNames_emu_mc, tmp) ;
        clear tmp
    elseif length(varNames_emu) < NvarsI
        tmp = cell(1, NvarsI-length(varNames_emu)) ;
        tmp(:) = {'NA'} ;
        varNames_emu = cat(2, varNames_emu, tmp) ;
    end
    varNames_emu_mc(m,:) = varNames_emu ;
    
    % Update cell array of cropLists
    NcropsI = max(NcropsI, length(cropList_emu)) ;
    if m==1
        cropList_emu_mc = cropList_emu ;
    elseif NcropsI > size(cropList_emu_mc,2)
        tmp = cell(m-1, NcropsI-size(cropList_emu_mc,2)) ;
        tmp(:) = {'NA'} ;
        cropList_emu_mc = cat(2, cropList_emu_mc, tmp) ;
        clear tmp
    elseif length(cropList_emu) < NcropsI
        tmp = cell(1, NcropsI-length(cropList_emu)) ;
        tmp(:) = {'NA'} ;
        cropList_emu = cat(2, cropList_emu, tmp) ;
    end
    cropList_emu_mc(m,:) = cropList_emu ;
    
    % Update cell array of lpj_asEmu
    cropList_lpj_asEmu_mc(m,:) = e2p_translate_crops( ...
        cropList_lpj, cropList_emu) ;
    clear tmp
end


%% Get crop-level data

prod_emu_ntmc = zeros(3, Ntpers, Nggcm, length(cropList_lpj)) ;
prod_out_ntmc = zeros(3, Ntpers, Nggcm, length(cropList_lpj)) ;
prod_lpj_ntc = zeros(3, Ntpers, length(cropList_lpj)) ;
irr_text = {'', 'i'} ;
for c = 1:length(cropList_lpj)
    for ii = 1:2
        thisCrop = [cropList_lpj{c} irr_text{ii}] ;
        disp(thisCrop)
        
        area_x = croparea_m2_xv(:, strcmp(bl_cf.varNames, thisCrop)) ;
        
        % Emulator
        % Cycle through models here because they have different variables
        % included. Not needed with Final Outputs section below because all
        % have the same variable list.
        for m = 1:Nggcm
            thisCrop_emu = [cropList_lpj_asEmu_mc{m,c} irr_text{ii}] ;
            yield_emu_xnt = data_fu_emu.garr_xvtm(:, ...
                strcmp(getbasenamei(varNames_emu_mc(m,:)), thisCrop_emu), :, m) ;
            
            % Save excluded cells (for applying to final outputs)
            excl_x11 = isnan(yield_emu_xnt(:,1,1)) ;
            if m == 1
                excl_xntm = false([size(yield_emu_xnt) Nggcm]) ;
            end
            excl_xntm(:,:,:,m) = repmat(excl_x11, [1 3 Ntpers]) ;
            
            prod_emu_xnt = yield_emu_xnt .* repmat(area_x, [1 3 Ntpers]) ;
            prod_emu_ntmc(:,:,m,c) = prod_emu_ntmc(:,:,m,c) ...
                + squeeze(nansum(prod_emu_xnt,1)) ;
        end
        
        % Final outputs
        yield_out_xntm = data_fu_out.garr_xvtm(:, ...
            strcmp(getbasenamei(data_fu_out.varNames), thisCrop), :, :) ;
        prod_out_xntm = yield_out_xntm .* repmat(area_x, [1 3 Ntpers Nggcm]) ;
        prod_out_xntm(excl_xntm) = NaN ;
        prod_out_ntmc(:,:,:,c) = prod_out_ntmc(:,:,:,c) ...
            + squeeze(nansum(prod_out_xntm,1)) ;
        
        % LPJ-GUESS simulation
        yield_lpj_xnt = data_fu_lpj.garr_xvt(:, ...
            strcmp(getbasenamei(data_fu_lpj.varNames), thisCrop), :) ;
        prod_lpj_xnt = yield_lpj_xnt .* repmat(area_x, [1 3 Ntpers]) ;
        prod_lpj_ntc(:,:,c) = prod_lpj_ntc(:,:,c) ...
            + squeeze(nansum(prod_lpj_xnt,1)) ;
        
    end
end
pctDiff_emu_ntmc = ...
    100*(prod_emu_ntmc(2:3,:,:,:) - prod_emu_ntmc([1 1],:,:,:)) ./ prod_emu_ntmc([1 1],:,:,:) ;
pctDiff_out_ntmc = ...
    100*(prod_out_ntmc(2:3,:,:,:) - prod_out_ntmc([1 1],:,:,:)) ./ prod_out_ntmc([1 1],:,:,:) ;
pctDiff_lpj_ntc = ...
    100*(prod_lpj_ntc(2:3,:,:) - prod_lpj_ntc([1 1],:,:)) ./ prod_lpj_ntc([1 1],:,:) ;

disp('Done')


%% Make plot: Emulators

overall_title = 'N sensitivity: Emulator outputs (dashed 10-60, solid 10-200)' ;

% Options
thisPos = figurePos ;
spacing = [0.1 0.05] ; % v, h
lineStyles = {'--', '-'} ; % Nstep1, Nstep2
lineWidth = 2 ;
fontSize = 14 ;

figure('Color', 'w', 'Position', thisPos) ;
legend_key = [{'LPJG sim.'} ; ggcm_list'] ;

for c = 1:length(cropList_lpj)
    subplot_tight(2, 3, c, spacing) ;
    thisCrop = cropList_lpj{c} ;
    for n = 1:2
        hold on
        thisLineStyle = lineStyles{n} ;
        h = plot(1:Ntpers, pctDiff_lpj_ntc(n,:,c), ...
            [thisLineStyle 'k'], ...
            'LineWidth', lineWidth) ;
        set(gca, 'ColorOrderIndex', 1) ;
        out_tm = squeeze(pctDiff_emu_ntmc(n,:,:,c)) ;
        hp = plot(1:Ntpers, out_tm, ...
            thisLineStyle, ...
            'LineWidth', lineWidth) ;
        h = [h; hp] ;
        hold off
        if strcmp(thisLineStyle, '-')
            hl = legend(h, legend_key, ...
                'Location', 'EastOutside', ...
                'Orientation', 'Vertical') ;
        end
    end
    title(thisCrop)
    set(gca, ...
        'FontSize', fontSize, ...
        'XLim', [1 Ntpers])
    xlabel('Decade')
    ylabel('% difference (N_{lo} to N_{hi})')
end

% Add overall title
hold on; ax2 = axes(gcf); hold off
ax2.Position = [0 0 1 1] ;
ht = text(ax2, 0.5, 0.97, ...
    overall_title, ...
    'FontSize', fontSize*1.5, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center') ;
ht.Units = 'normalized' ;
ht.Position(1) = 0.5 ;
ax2.Visible = 'off' ;



%% Make plot: Final outputs

overall_title = 'N sensitivity: To PLUM (dashed 10-60, solid 10-200)' ;

% Options
thisPos = figurePos ;
spacing = [0.1 0.05] ; % v, h
lineStyles = {'--', '-'} ; % Nstep1, Nstep2
lineWidth = 2 ;
fontSize = 14 ;

figure('Color', 'w', 'Position', thisPos) ;
legend_key = [{'LPJG sim.'} ; ggcm_list'] ;

for c = 1:length(cropList_lpj)
    subplot_tight(2, 3, c, spacing) ;
    thisCrop = cropList_lpj{c} ;
    for n = 1:2
        hold on
        thisLineStyle = lineStyles{n} ;
        h = plot(1:Ntpers, pctDiff_lpj_ntc(n,:,c), ...
            [thisLineStyle 'k'], ...
            'LineWidth', lineWidth) ;
        set(gca, 'ColorOrderIndex', 1) ;
        out_tm = squeeze(pctDiff_out_ntmc(n,:,:,c)) ;
        hp = plot(1:Ntpers, out_tm, ...
            thisLineStyle, ...
            'LineWidth', lineWidth) ;
        h = [h; hp] ;
        hold off
        if strcmp(thisLineStyle, '-')
            hl = legend(h, legend_key, ...
                'Location', 'EastOutside', ...
                'Orientation', 'Vertical') ;
        end
    end
    title(thisCrop)
    set(gca, ...
        'FontSize', fontSize, ...
        'XLim', [1 Ntpers])
    xlabel('Decade')
    ylabel('% difference (N_{lo} to N_{hi})')
end

% Add overall title
hold on; ax2 = axes(gcf); hold off
ax2.Position = [0 0 1 1] ;
ht = text(ax2, 0.5, 0.97, ...
    overall_title, ...
    'FontSize', fontSize*1.5, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center') ;
ht.Units = 'normalized' ;
ht.Position(1) = 0.5 ;
ax2.Visible = 'off' ;


%% Bar graph

overall_title = 'Global mean N sensitivity' ;

% Options
thisPos = figurePos ;
spacing = [0.1 0.05] ; % v, h
lineStyles = {'--', '-'} ; % Nstep1, Nstep2
lineWidth = 2 ;
fontSize = 14 ;

figure('Color', 'w', 'Position', thisPos) ;
% legend_key = [{'LPJG sim.'} ; ggcm_list'] ;
tmp = {'10-60 before', '10-60 after', '10-200 before', '10-200 after'} ;
X = categorical(tmp);
X = reordercats(X,tmp);

for c = 1:length(cropList_lpj)
    subplot_tight(2, 3, c, spacing) ;
    thisCrop = cropList_lpj{c} ;
    
    % Bars
    bar_data = [] ;
    bar_data(:,1) = squeeze(mean(pctDiff_emu_ntmc(1,:,:,c),2)) ;
    bar_data(:,2) = squeeze(mean(pctDiff_out_ntmc(1,:,:,c),2)) ;
    bar_data(:,3) = squeeze(mean(pctDiff_emu_ntmc(2,:,:,c),2)) ;
    bar_data(:,4) = squeeze(mean(pctDiff_out_ntmc(2,:,:,c),2)) ;
    hbars = bar(X, bar_data') ;
    title(thisCrop)
    ylabel('% difference (N_{lo} to N_{hi})')
    
    % Lines for LPJ-GUESS simulation
    hold on
    tmp = {'10-60 before', '10-60 after'} ;
    Xtmp = categorical(tmp);
    Xtmp = reordercats(Xtmp,tmp);
    hl = plot(Xtmp, ...
        [1 1]*mean(pctDiff_lpj_ntc(1,:,c),2), ...
        '-k', ...
        'LineWidth', 2.5) ;
    set(hl, 'XData', [0 2.4])
    tmp = {'10-200 before', '10-200 after'} ;
    Xtmp = categorical(tmp);
    Xtmp = reordercats(Xtmp,tmp);
    hl = plot(Xtmp, ...
        [1 1]*mean(pctDiff_lpj_ntc(2,:,c),2), ...
        '-k', ...
        'LineWidth', 2.5) ;
    set(hl, 'XData', [2.6 4.5])
    hold off
    
    legend(hbars, ggcm_list, ...
        'Location', 'SouthOutside', ...
        'Orientation', 'Horizontal') ;
    set(gca, ...
        'FontSize', fontSize)
end

% Add overall title
hold on; ax2 = axes(gcf); hold off
ax2.Position = [0 0 1 1] ;
ht = text(ax2, 0.5, 0.97, ...
    overall_title, ...
    'FontSize', fontSize*1.5, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center') ;
ht.Units = 'normalized' ;
ht.Position(1) = 0.5 ;
ax2.Visible = 'off' ;
