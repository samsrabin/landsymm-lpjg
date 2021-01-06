%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make figures from PLUM outputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topDir = '/Volumes/Reacher/G2P/outputs_PLUM/ggcmi/test_20210101' ;
runList = { ...
    'EPIC-TAMU (emu)', 'epictamu_emu/s1_26' ;
    'GEPIC (emu)', 'gepic_emu/s1_26' ;
    'LPJ-GUESS (sim)', 'lpjguess_sim/s1_26' ;
    } ;


%% Setup

Nruns = size(runList, 1) ;

fileList = {'landuse', 'landuse_etc', 'cropfracs', 'fert', 'irrig', 'yield'} ;
% fileList = {'landuse', 'landuse_etc', 'cropfracs'} ;
Nfiles = length(fileList) ;

reqfields = {'lonlats', 'list2map', 'varNames', 'garr_xvy', 'yearList'} ;

models_legend = runList(:,1) ;

outDir = sprintf('%s/figs', topDir) ;
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

pngres = '-r100' ;


%% Import

for f = 1:Nfiles
    thisFile = fileList{f} ;
    fprintf('Importing %s:\n', thisFile)
    
    for r = 1:Nruns
        thisRun = runList{r,1} ;
        fprintf('   %s (%d/%d)...\n', thisRun, r, Nruns)
        
        inFile = sprintf('%s/%s-matlab/%s', topDir, runList{r,2}, thisFile) ;
        in_tmp = lpjgu_matlab_read2geoArray(inFile) ;
        
        % Remove pasture from everything except land use fractions
        if ~strcmp(thisFile, 'LandCoverFract')
            in_tmp.garr_xvy(:,strcmp(in_tmp.varNames, 'Pasture'),:) = [] ;
            in_tmp.varNames(strcmp(in_tmp.varNames, 'Pasture')) = [] ;
        end
        
        % Make sure all fields are present
        fieldok = isfield(in_tmp, reqfields) ;
        if ~all(fieldok)
            error('Imported struct missing field(s): %s', strjoin(reqfields(~fieldok)))
        end
        
%         % Move setaside to end
%         is_setaside = contains(in_tmp.varNames, 'setaside') ;
%         if any(is_setaside)
%             in_tmp.garr_xvy(:,is_setaside,:) = [] ;
%             in_tmp.varNames(is_setaside) = [] ;
%         end
%         clear is_setaside
        
        % Convert Mha to km2, if needed
        [~, IA] = intersect(in_tmp.varNames, ...
            {'area', 'suitable', 'protected', ...
            'natural', 'cropland', 'pasture', 'barren', 'urban'}) ;
        if ~isempty(IA)
            in_tmp.garr_xvy(:,IA,:) = 1e4 * in_tmp.garr_xvy(:,IA,:) ;
        end
        clear IA
        
        % Set up or check structures
        if r==1
            out_struct = in_tmp ;
            out_struct = rmfield(out_struct, 'garr_xvy') ;
            out_struct.garr_xvyr = nan([size(in_tmp.garr_xvy) Nruns]) ;
            if ~exist('yearList', 'var')
                yearList = in_tmp.yearList ;
            else
                if ~isequal(yearList, in_tmp.yearList)
                    error('yearList mismatch')
                end
            end
        else
            
            if ~isequal(out_struct.lonlats, in_tmp.lonlats)
                error('Gridlist mismatch in %s: %d vs. %d', ...
                    inFile, y1, thisYear)
            elseif ~isequal(out_struct.varNames, in_tmp.varNames)
                error('varNames mismatch in %s: %d vs. %d', ...
                    inFile, y1, thisYear)
            elseif ~isequal(out_struct.yearList, in_tmp.yearList)
                error('yearList mismatch in %s: %d vs. %d', ...
                    inFile, y1, thisYear)
            end
        end
        out_struct.garr_xvyr(:,:,:,r) = in_tmp.garr_xvy ;
        clear in_tmp
    end
    
    eval([thisFile ' = out_struct ;']) ;
    clear out_struct
end

if exist('cropfracs', 'var')
    cropList = cropfracs.varNames ;
elseif exist('fert', 'var')
    cropList = fert.varNames ;
elseif exist('irrig', 'var')
    cropList = irrig.varNames ;
elseif exist('yield', 'var')
    cropList = yield.varNames ;
end

% Get crop areas, converting km2 to ha
croparea_ha_xvyr = 100 * repmat(landuse.garr_xvyr(:,strcmp(landuse.varNames, 'cropland'),:,:), ...
    [1 length(cropfracs.varNames) 1 1]) .* cropfracs.garr_xvyr ;

disp('Done.')


%% Make overview time series

% Options %%%%%%%%%%%%%%%%%%%%
thisPos = figurePos ;
spacing = [0.1 0.05] ; % y x
lineWidth = 3 ;
fontSize = 14 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'w', 'Position', thisPos) ;

% Land use area
plot_lu_list = {'cropland', 'pasture', 'natural'} ;
for v = 1:length(plot_lu_list)
    thisVar = plot_lu_list{v} ;
    vv = find(strcmp(landuse.varNames, thisVar)) ;
    subplot_tight(2, 3, v, spacing)
    
    area_yr = squeeze(sum(landuse.garr_xvyr(:,vv,:,:),1)) ;
    % Convert km2 to Mha
    area_yr = area_yr * 1e-4 ;
    units = 'Mha' ;
    
    plot(yearList, area_yr, ...
        'LineWidth', lineWidth) ;
    
    legend(models_legend, 'Location', 'Best')
    title(thisVar) ;
    xlabel('Year')
    ylabel(units)
    set(gca, ...
        'XLim', minmax_ssr(yearList), ...
        'FontSize', fontSize)
end

% Crop production
subplot_tight(2, 3, 4, spacing) ;
data_yr = squeeze(sum(sum((1e-9*yield.garr_xvyr) .* croparea_ha_xvyr,1),2)) ;
units = 'Mt dry matter' ;
plot(yearList, data_yr, ...
    'LineWidth', lineWidth) ;
legend(models_legend, 'Location', 'Best')
title('crop production') ;
xlabel('Year')
ylabel(units)
set(gca, ...
    'XLim', minmax_ssr(yearList), ...
    'FontSize', fontSize)

% Fertilizer
subplot_tight(2, 3, 5, spacing) ;
data_yr = squeeze(sum(sum((1e-9*fert.garr_xvyr) .* croparea_ha_xvyr,1),2)) ;
units = 'Mt N' ;
plot(yearList, data_yr, ...
    'LineWidth', lineWidth) ;
legend(models_legend, 'Location', 'Best')
title('fertilizer use') ;
xlabel('Year')
ylabel(units)
set(gca, ...
    'XLim', minmax_ssr(yearList), ...
    'FontSize', fontSize) 

% Irrigation
subplot_tight(2, 3, 6, spacing) ;
data_yr = squeeze(sum(sum(irrig.garr_xvyr .* (1e-15*1e6*croparea_ha_xvyr),1),2)) ;
units = 'Kkm^3' ;
plot(yearList, data_yr, ...
    'LineWidth', lineWidth) ;
legend(models_legend, 'Location', 'Best')
title('irrigation use') ;
xlabel('Year')
ylabel(units)
set(gca, ...
    'XLim', minmax_ssr(yearList), ...
    'FontSize', fontSize)

% Save figure
export_fig(sprintf('%s/ts_overview.pdf', outDir), gcf)
close


%% Per-crop time series

% thisTitle = 'Crop area' ;
thisTitle = 'Fertilizer' ;
% thisTitle = 'Irrigation' ;

if strcmp(thisTitle, 'Crop area')
    data_vyr = squeeze(sum(croparea_ha_xvyr, 1)) ;
    % Convert ha to Mha
    data_vyr = data_vyr * 1e-6 ;
    units = 'Mha' ;
elseif strcmp(thisTitle, 'Fertilizer')
    % Convert kg to Mt
    data_vyr = squeeze(sum((1e-9*fert.garr_xvyr) .* croparea_ha_xvyr, 1)) ;
    units = 'Mt N' ;
elseif strcmp(thisTitle, 'Irrigation')
    % Convert Mt to Kkm3
    data_vyr = squeeze(sum(irrig.garr_xvyr .* (1e-15*1e6*croparea_ha_xvyr), 1)) ;
    units = 'Kkm3' ;
else
    error('thisTitle %s not recognized', thisTitle)
end

make_percrop_fig(data_vyr, cropList, models_legend, yearList, units, thisTitle)
clear data_vyr units

% Save figure
export_fig(sprintf('%s/ts_crops_%s.pdf', outDir, strrep(lower(thisTitle), ' ', '')), gcf)
% close


%% FUNCTIONS

function make_percrop_fig(data_vyr, cropList, models_legend, yearList, ...
    units, thisTitle)

% Options %%%%%%%%%%%%%%%%%%%%
thisPos = figurePos ;
spacing = [0.08 0.05] ; % y x
lineWidth = 3 ;
fontSize = 14 ;
ny = 3 ;
nx = 4 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'w', 'Position', thisPos) ;

Ncrops = length(cropList) ;

for c = 1:Ncrops
    thisCrop = cropList{c} ;
    subplot_tight(ny, nx, c, spacing) ;
    
    data_yr = squeeze(data_vyr(c,:,:)) ;
    
    
    plot(yearList, data_yr, ...
        'LineWidth', lineWidth) ;
    
    title(thisCrop) ;
%     xlabel('Year')
    ylabel(units)
    set(gca, ...
        'XLim', minmax_ssr(yearList), ...
        'FontSize', fontSize)

end

legendflex(models_legend, ...
    'ref', gca, ...
    'anchor', {'w','w'}, ...
    'buffer', [1/nx 0], 'bufferunit', 'normalized')

% Add main title
hold on
hat = axes ;
hat.Position = [0 0 1 1] ;
ht = text(hat, 100, 100, ...
    thisTitle, ...
    'FontSize', fontSize+4, 'FontWeight', 'bold') ;
ht.Units = 'normalized' ;
ht.HorizontalAlignment = 'center' ;
ht.Position(1) = 0.5 ;
ht.Position(2) = 0.97 ;
hat.Visible = 'off' ;
hold off

end















