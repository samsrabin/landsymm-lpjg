%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make fertilizer-yield curves for e2p outputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ggcm_list = {'LPJ-GUESS', 'GEPIC', 'EPIC-TAMU', 'pDSSAT'} ;

% topDir = ['/Volumes/Reacher/GGCMI/AgMIP.output/CMIP_emulated_work/' ...
%     'A1_v2.5_UKESM1-0-LL_ssp126_20210528_ph2bl_intpinfs_rmolend_fake1k'] ;
% topDir = ['/Volumes/Reacher/GGCMI/AgMIP.output/CMIP_emulated_work/' ...
%     'A1_v2.5_UKESM1-0-LL_ssp126_20210531_ph2bl_intpinfs_rmolend_fake1k'] ;
% topDir = ['/Volumes/Reacher/GGCMI/AgMIP.output/CMIP_emulated_work/' ...
%     'h'] ;
% remapVer = '12_g2p' ;
topDir = '/Volumes/Reacher/GGCMI/AgMIP.output/CMIP_emulated_work/A1_v2.5_UKESM1-0-LL_ssp126_20210923_ph2bl_intpinfs_rmolend_fake1k' ;


%% Setup

if ~strcmp(ggcm_list{1}, 'LPJ-GUESS')
    warning('Code relies on ggcm_list{1} being LPJ-GUESS; setting this.')
    if any(strcmp(ggcm_list, 'LPJ-GUESS'))
        ggcm_list(strcmp(ggcm_list, 'LPJ-GUESS')) = [] ;
    end
    ggcm_list = [{'LPJ-GUESS'} ggcm_list] ;
end

Nggcm = length(ggcm_list) ;

irrList = {'', 'i'} ;
Nirr = length(irrList) ;

% Cell area (convert from km2 to ha)
cell_area_YXqd = transpose(ncread( ...
    '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc', ...
    'carea')) ;
cell_area_YX = aggregate_land_area(cell_area_YXqd,0.5,0.5) ;
cell_area_YX = cell_area_YX * 100 ;

% Should really save these into the MAT files
future_y1 = 2005 ;
baseline_y1 = 2001 ;
baseline_yN = 2010 ;
future_ts = 10 ; % Number of years in future time step
future_yN_emu = 2084 ;
Nyears_ts = future_ts ;
tsN_y1 = floor(future_yN_emu/Nyears_ts)*Nyears_ts ;
ts1_list = future_y1:Nyears_ts:tsN_y1 ;
tsN_list = ts1_list + Nyears_ts - 1 ;
Ntpers = length(ts1_list) ;

% Set up directory for output figure files
outDir = sprintf('%s/fert-yield_curves', topDir) ;
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

% Import calibration factors
calib_factors = cell(1, Nggcm) ;
for g = 1:Nggcm
    thisGGCM = ggcm_list{g} ;
    thisDir = dir(sprintf('%s/*%s*', topDir, thisGGCM)) ;
    if length(thisDir) ~= 1
        error('Expected to find 1 match for %s; found %d', ...
            thisGGCM, length(thisDir))
    end
    thisDir = sprintf('%s/%s', topDir, thisDir.name) ;
    cf_file = sprintf('%s/calib_factors.csv', thisDir) ;
    calib_factors{g} = readtable(cf_file) ;
end


%% Import yields
disp('Importing simulated and emulated yields...')

% LPJ-GUESS
matfile = sprintf('%s/sim_LPJ-GUESS/future_yield.mat', topDir) ;
load(matfile, 'data_fu_lpj2') ;
garr_xvt = data_fu_lpj2.garr_xvt ;
data_fu = rmfield(data_fu_lpj2, 'garr_xvt') ;
clear data_fu_lpj2
data_fu.garr_xvtg = nan([size(garr_xvt) Nggcm]) ;
data_fu.garr_xvtg(:,:,:,1) = garr_xvt ;
clear garr_xvt

for g = 2:Nggcm
    thisGGCM = ggcm_list{g} ;
    matfile = sprintf('%s/emu_%s/future_yield.mat', topDir, thisGGCM) ;
    load(matfile, 'data_fu_emu3') ;
    
    % Sanity checks
    if ~isequal(data_fu.list2map, data_fu_emu3.list2map)
        error('list2map mismatch')
    elseif ~isequal(data_fu.varNames, data_fu_emu3.varNames)
        error('varNames mismatch')
    end
    
    data_fu.garr_xvtg(:,:,:,g) = data_fu_emu3.garr_xvt ;
    clear data_fu_emu3
end

% Sort cropIrrN list
[C,I] = sort(data_fu.varNames) ;
data_fu.garr_xvtg = data_fu.garr_xvtg(:,I,:,:) ;
data_fu.varNames = data_fu.varNames(I) ;

% Get crop and N lists
cropList = unique(getbasename(data_fu.varNames)) ;
Ncrops = length(cropList) ;
cropIrrList = unique(getbasenamei(data_fu.varNames)) ;
N_list = unique(getN_num(data_fu.varNames)) ;
Nn = length(N_list) ;

% Get calibration factor array
cfs_cg = nan(Ncrops, Nggcm) ;
for g = 1:Nggcm
    T_cf = calib_factors{g} ;
    [~, ~, IB] = intersect(cropList, T_cf.Crop, 'stable') ;
    cfs_cg(:,g) = T_cf.calibration_factor(IB) ;
    clear T_cf
end

% Embed timestep lists
data_fu.ts1_list = ts1_list ;
data_fu.tsN_list = tsN_list ;

disp('Done.')

% % What is the maximum yield seen?
% [M, I] = max(data_fu.garr_xvtg(:)) ;
% [x, v, t, g] = ind2sub(size(data_fu.garr_xvtg), I) ;
% fprintf('%s %s %d-%d: %g t/ha\n', ...
%     ggcm_list{g}, data_fu.varNames{v}, ts1_list(t), tsN_list(t), ...
%     M*10)


%% Import land use
disp('Importing land use...')

% Determine which maps should be used
if isequal(shiftdim(sort(cropList)), shiftdim(sort({'CerealsC3', 'CerealsC4', ...
        'FruitAndVeg', 'Oilcrops', 'Pulses', 'Rice', 'StarchyRoots', ...
        'Sugar'})))
    LUdir = '/Users/Shared/PLUM/input/remaps_v8c' ;
else
    disp(cropList)
    error('Unknown LU directory for this cropList')
end

% Get file paths
filename_landuse = get_luinput_filepath(LUdir, 'LU') ;
filename_cropfrac = get_luinput_filepath(LUdir, 'cropfracs') ;

% Import land uses
thisYear = 2015 ;
LU = lpjgu_matlab_read2geoArray(filename_landuse) ;
garr_xv = LU.garr_xvy(:,:,LU.yearList==thisYear) ;
LU = rmfield(LU, {'garr_xvy', 'yearList'}) ;
LU.garr_xv = garr_xv ;
clear garr_xv

% Import crop fractions (already only 1 year)
cropfracs = lpjgu_matlab_read2geoArray(filename_cropfrac) ;
toRemove = contains(cropfracs.varNames, {'ExtraCrop', 'Miscanthus'}) ;
cropfracs.garr_xv(:,toRemove) = [] ;
cropfracs.varNames(toRemove) = [] ;

% Combine wheats, if needed
if any(strcmp(cropfracs.varNames, 'CerealsC3s'))
    c3_x = ...
        cropfracs.garr_xv(:,strcmp(cropfracs.varNames, 'CerealsC3s')) ...
        + cropfracs.garr_xv(:,strcmp(cropfracs.varNames, 'CerealsC3w')) ;
    c3i_x = ...
        cropfracs.garr_xv(:,strcmp(cropfracs.varNames, 'CerealsC3si')) ...
        + cropfracs.garr_xv(:,strcmp(cropfracs.varNames, 'CerealsC3wi')) ;
    toRemove = contains(cropfracs.varNames, 'CerealsC3') ;
    cropfracs.garr_xv(:,toRemove) = [] ;
    cropfracs.varNames(toRemove) = [] ;
    cropfracs.garr_xv = cat(2, cropfracs.garr_xv, c3_x, c3i_x) ;
    cropfracs.varNames = [cropfracs.varNames {'CerealsC3', 'CerealsC3i'}] ;
end

% Consistency check
if ~isequal(cropList, unique(getbasename(cropfracs.varNames)))
    error('cropList mismatch')
elseif ~isequal(cropIrrList, unique(getbasenamei(cropfracs.varNames)))
    error('cropIrrList mismatch')
end

% Align gridlists
LU = align_gridlist(LU, data_fu.list2map) ;
cropfracs = align_gridlist(cropfracs, data_fu.list2map) ;

% Convert to crop area
cropareas = cropfracs ;
cropareas.garr_xv = cropfracs.garr_xv ...
    .* LU.garr_xv(:,strcmp(LU.varNames,'CROPLAND')) ...
    .* cell_area_YX(data_fu.list2map) ;
clear cropfracs LU

disp('Done.')


%% Make curves: Each time period, all crops in one figure, area-wtd
warning off export_fig:exportgraphics

% Figure options
thisPos = figurePos ;
fontSize = 12 ;
lineWidth = 1.25 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisCO = colororder ;

output_units = {'prod', 'yield'} ;

for t = 1:Ntpers
    
    for u = 1:length(output_units)
        
        prod_or_yield = output_units{u} ;
        
        make_figure_areaweighted(cropareas, data_fu, t, ggcm_list, cropList, irrList, ...
            thisPos, fontSize, lineWidth, thisCO, cfs_cg, prod_or_yield)
        export_fig(sprintf('%s/fertyield_%s_%d-%d.png', ...
            outDir, prod_or_yield, data_fu.ts1_list(t), data_fu.tsN_list(t)), ...
            '-r200')
        close
        
%         for g = 1:Nggcm
%             thisGGCM = ggcm_list{g} ;
%             fprintf('%s...\n', thisGGCM)
%             data_fu_tmp = data_fu ;
%             data_fu_tmp.garr_xvtg = data_fu.garr_xvtg(:,:,:,g) ;
%             make_figure_areaweighted(cropareas, data_fu_tmp, t, ggcm_list(g), cropList, irrList, ...
%                 thisPos, fontSize, lineWidth, thisCO(g,:), cfs_cg, prod_or_yield)
%             export_fig(sprintf('%s/fertyield_%s_%d-%d_%s.png', ...
%                 outDir, prod_or_yield, data_fu.ts1_list(t), data_fu.tsN_list(t), ...
%                 thisGGCM), '-r200')
%             close
%         end
    end
    
end

disp('Done!')


%% Make curves: Each time period, all crops in one figure, percentile
warning off export_fig:exportgraphics

% What percentile of yield?
this_prctile = 100 ;
% Figure options
thisPos = figurePos ;
fontSize = 12 ;
lineWidth = 1.25 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisCO = colororder ;

for t = 1:Ntpers

    make_figure_allcellsequal(data_fu, t, ggcm_list, cropList, irrList, ...
        thisPos, fontSize, lineWidth, thisCO, cfs_cg, this_prctile)
    export_fig(sprintf('%s/fertyield_%g-pctleYield_%d-%d.png', ...
        outDir, this_prctile, data_fu.ts1_list(t), data_fu.tsN_list(t)), ...
        '-r200')
    close
    
%     for g = 1:Nggcm
%         thisGGCM = ggcm_list{g} ;
%         fprintf('%s...\n', thisGGCM)
%         data_fu_tmp = data_fu ;
%         data_fu_tmp.garr_xvtg = data_fu.garr_xvtg(:,:,:,g) ;
%         make_figure_allcellsequal(data_fu_tmp, t, ggcm_list(g), cropList, irrList, ...
%             thisPos, fontSize, lineWidth, thisCO(g,:), cfs_cg, this_prctile)
%         export_fig(sprintf('%s/fertyield_%g-pctleYield_%d-%d_%s.png', ...
%             outDir, this_prctile, data_fu.ts1_list(t), data_fu.tsN_list(t), ...
%             thisGGCM), '-r200')
%         close
%     end

end


    
disp('Done!')


%% FUNCTIONS

function make_figure_areaweighted(cropareas, data_fu, t, ggcm_list, cropList, irrList, ...
    thisPos, fontSize, lineWidth, thisCO, cfs_cg, prod_or_yield)

Nggcm = length(ggcm_list) ;
Ncrops = length(cropList) ;

theLegend = [ ...
    strcat(ggcm_list, ' (rainfed)') ;
    strcat(ggcm_list, ' (irrigated)') ...
    ] ;
theLegend = theLegend(:) ;

figure('Color', 'w', 'Position', thisPos)
tiledlayout('flow')

for cc = 1:Ncrops+1
    c = min(cc,Ncrops) ;
    is_fakeplot = cc==Ncrops+1 ;    
    
    thisCrop = cropList{c} ;
    is_thisCrop = contains( ...
        data_fu.varNames, thisCrop) ;
    varNames_thisCrop = ...
        data_fu.varNames(is_thisCrop) ;
    data_thisCrop_xvg = ...
        data_fu.garr_xvtg(:,is_thisCrop,t,:) ;
    
    thisCrop_area_x = nansum(...
        cropareas.garr_xv(:,contains(cropareas.varNames, thisCrop)), 2) ;
    
    nexttile
    for g = 1:Nggcm
        
        for ii = 1:length(irrList)
            
            % Rainfed dashed, irrigated solid
            if ii == 1
                lineMarkerStyle = '--o' ;
                markerSize = 6.5 ;
            else
                lineMarkerStyle = '-o' ;
                markerSize = 6.5 ;
            end
            
            % Get and sort data for this crop*irrigation
            thisCropIrr = [thisCrop irrList{ii}] ;
            is_thisCropIrr = ~cellfun(@isempty, regexp( ...
                varNames_thisCrop, ...
                [thisCropIrr '\d+'])) ;
            varNames_thisCropIrr = ...
                varNames_thisCrop(is_thisCropIrr) ;
            N_thisCropIrr = ...
                getN_num(varNames_thisCropIrr) ;
            I_is_thisCropIrr = find(is_thisCropIrr) ;
            [N_thisCropIrr, I_sort] = sort(N_thisCropIrr) ;
            I_is_thisCropIrr = I_is_thisCropIrr(I_sort) ;
%             varNames_thisCropIrr = varNames_thisCropIrr(I_sort) ;
            
            % Get values
            data_thisCropIrr_n = nansum( ...
                data_thisCrop_xvg(:,I_is_thisCropIrr,g) ...
                .* thisCrop_area_x, ...
                1) ;
            if strcmp(prod_or_yield, 'prod')
                % Convert t to Mt
                data_thisCropIrr_n = data_thisCropIrr_n * 1e-6 ;
            elseif strcmp(prod_or_yield, 'yield')
                data_thisCropIrr_n = data_thisCropIrr_n / sum(thisCrop_area_x) ;
            else
                error('prod_or_yield %s not recognized', prod_or_yield)
            end
            
            % Multiply by calibration factors
            data_thisCropIrr_n = data_thisCropIrr_n ...
                * cfs_cg(c,g) ;
            
            % Make fakeplot lines invisible
            if is_fakeplot
                data_thisCropIrr_n = nan(size(data_thisCropIrr_n)) ;
            end
            
            hold on
%             set(gca, 'ColorOrderIndex', g) ;
            thisColor = thisCO(g,:) ;
            hl = plot(N_thisCropIrr, ...
                data_thisCropIrr_n, ...
                lineMarkerStyle, ...
                'MarkerSize', markerSize, ...
                'LineWidth', lineWidth, ...
                'Color', thisColor, ...
                'MarkerEdgeColor', thisColor) ;
            set(gca, ...
                'XTick', N_thisCropIrr, ...
                'FontSize', fontSize)
            if ii == 2
                set(hl, 'MarkerFaceColor', thisColor)
            end
            
            % Add calibration factor text
            if length(ggcm_list) == 1 && ii==1 && ~is_fakeplot
                thisText = sprintf('Calib. factor = %0.3f', ...
                    cfs_cg(c,g)) ;
                text(0.5, 0.05, thisText, ...
                    'FontSize', fontSize, ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'center') ;
            end
            
            hold off
            
        end % Loop through irrigation
    end % Loop through GGCMs
    if ~is_fakeplot
        title(sprintf('%s %d-%d', ...
            thisCrop, data_fu.ts1_list(t), data_fu.tsN_list(t)))
        xlabel('Fertilization (kg/ha)')
        if strcmp(prod_or_yield, 'prod')
            ylabel('Production (Mt)')
        elseif strcmp(prod_or_yield, 'yield')
            ylabel('Yield (t/ha)')
        else
            error('prod_or_yield %s not recognized', prod_or_yield)
        end
    end
end % Loop through crops

axis off
legend(theLegend, ...
    'Location', 'west', ...
    'FontSize', fontSize + 2) ;

end


function make_figure_allcellsequal(data_fu, t, ggcm_list, cropList, irrList, ...
    thisPos, fontSize, lineWidth, thisCO, cfs_cg, this_prctile)

Nggcm = length(ggcm_list) ;
Ncrops = length(cropList) ;

theLegend = [ ...
    strcat(ggcm_list, ' (rainfed)') ;
    strcat(ggcm_list, ' (irrigated)') ...
    ] ;
theLegend = theLegend(:) ;

figure('Color', 'w', 'Position', thisPos)
tiledlayout('flow')

for cc = 1:Ncrops+1
    c = min(cc,Ncrops) ;
    is_fakeplot = cc==Ncrops+1 ;    
    
    thisCrop = cropList{c} ;
    is_thisCrop = contains( ...
        data_fu.varNames, thisCrop) ;
    varNames_thisCrop = ...
        data_fu.varNames(is_thisCrop) ;
    data_thisCrop_xvg = ...
        data_fu.garr_xvtg(:,is_thisCrop,t,:) ;
        
    nexttile
    for g = 1:Nggcm
        
        for ii = 1:length(irrList)
            
            % Rainfed dashed, irrigated solid
            if ii == 1
                lineMarkerStyle = '--o' ;
                markerSize = 6.5 ;
            else
                lineMarkerStyle = '-o' ;
                markerSize = 6.5 ;
            end
            
            % Get and sort data for this crop*irrigation
            thisCropIrr = [thisCrop irrList{ii}] ;
            is_thisCropIrr = ~cellfun(@isempty, regexp( ...
                varNames_thisCrop, ...
                [thisCropIrr '\d+'])) ;
            varNames_thisCropIrr = ...
                varNames_thisCrop(is_thisCropIrr) ;
            N_thisCropIrr = ...
                getN_num(varNames_thisCropIrr) ;
            I_is_thisCropIrr = find(is_thisCropIrr) ;
            [N_thisCropIrr, I_sort] = sort(N_thisCropIrr) ;
            I_is_thisCropIrr = I_is_thisCropIrr(I_sort) ;
%             varNames_thisCropIrr = varNames_thisCropIrr(I_sort) ;
            
            % Get values
            data_thisCropIrr_n = prctile( ...
                data_thisCrop_xvg(:,I_is_thisCropIrr,g), ...
                this_prctile, ...
                1) ;
            
            % Convert to t/ha
            data_thisCropIrr_n = 10 * data_thisCropIrr_n ;
            
            % Multiply by calibration factors
            data_thisCropIrr_n = data_thisCropIrr_n ...
                * cfs_cg(c,g) ;
            
            % Make fakeplot lines invisible
            if is_fakeplot
                data_thisCropIrr_n = nan(size(data_thisCropIrr_n)) ;
            end
            
            hold on
%             set(gca, 'ColorOrderIndex', g) ;
            thisColor = thisCO(g,:) ;
            hl = plot(N_thisCropIrr, ...
                data_thisCropIrr_n, ...
                lineMarkerStyle, ...
                'MarkerSize', markerSize, ...
                'LineWidth', lineWidth, ...
                'Color', thisColor, ...
                'MarkerEdgeColor', thisColor) ;
            set(gca, ...
                'XTick', N_thisCropIrr, ...
                'FontSize', fontSize)
            if ii == 2
                set(hl, 'MarkerFaceColor', thisColor)
            end
            
            % Add calibration factor text
            if length(ggcm_list) == 1 && ii==1 && ~is_fakeplot
                thisText = sprintf('Calib. factor = %0.3f', ...
                    cfs_cg(c,g)) ;
                text(0.5, 0.05, thisText, ...
                    'FontSize', fontSize, ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'center') ;
            end
            
            hold off
            
        end % Loop through irrigation
    end % Loop through GGCMs
    if ~is_fakeplot
        title(sprintf('%s %d-%d: %g percentile yield', ...
            thisCrop, data_fu.ts1_list(t), data_fu.tsN_list(t), this_prctile))
        xlabel('Fertilization (kg/ha)')
        ylabel('Yield (t/ha)')
    end
end % Loop through crops

axis off
legend(theLegend, ...
    'Location', 'west', ...
    'FontSize', fontSize + 2) ;

end


function S = align_gridlist(S, list2map)

[C,~,IB] = intersect(list2map, S.list2map, 'stable') ;

if ~isequal(C, list2map)
    error('Not all cells in target list2map were found in input struct')
end

S.garr_xv = S.garr_xv(IB,:) ;
S.list2map = list2map ;

end


function land_area_YX_out = aggregate_land_area(land_area_YX_in,xres_out,yres_out)

% Get info first
xres_in = 360 / size(land_area_YX_in,2) ;
yres_in = 180 / size(land_area_YX_in,1) ;
if xres_out < xres_in
    error(['xres_out ' num2str(xres_out) ' is too small given xres_in=' num2str(xres_in).'])
    
elseif yres_out < yres_in
    error(['yres_out ' num2str(yres_out) ' is too small given yres_in=' num2str(yres_in) '.'])
    
elseif xres_in~=xres_out || yres_in~=yres_out
    
    xratio = xres_out / xres_in ;
    yratio = yres_out / yres_in ;
    if ~isint(xratio) || ~isint(yratio)
        error(['xratio (' num2str(xratio) ') and yratio (' num2str(yratio) ') must both be integers!'])
    end
    
    % Aggregate X
    land_area_YX_tmp = zeros(size(land_area_YX_in,1),size(land_area_YX_in,2)/xratio) ;
    for i = 1:xratio
        land_area_YX_tmp = land_area_YX_tmp + land_area_YX_in(:,i:xratio:end) ;
    end
    
    % Check
    globalland_area_in = sum(sum(land_area_YX_in)) ;
    if abs(sum(sum(land_area_YX_tmp)) - globalland_area_in)>1e-4
        error('Error in X aggregation.')
    end
    
    % Aggregate Y
    land_area_YX_out = zeros(size(land_area_YX_in,1)/yratio,size(land_area_YX_tmp,2)) ;
    for j = 1:yratio
        land_area_YX_out = land_area_YX_out + land_area_YX_tmp(j:yratio:end,:) ;
    end
    
    % Check
    if abs(sum(sum(land_area_YX_out)) - globalland_area_in)>1e-4
        error('Error in Y aggregation.')
    end
    
else
    land_area_YX_out = land_area_YX_in ;
    
end


end


function filepath = get_luinput_filepath(LUdir, file_prefix)

filelist = dir(sprintf('%s/%s.*.txt*', LUdir, file_prefix)) ;
filepath = sprintf('%s/%s', ...
    LUdir, regexprep(filelist(1).name, '.txt.*', '.txt'));

end
