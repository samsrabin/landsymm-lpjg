function data_fu_out = e2p_apply_deltas( ...
    data_bl_thisBL, data_bl_emu, data_fu_emu, deltas_emu_xvt, ...
    cropList_thisBL, cropList_thisBL_asEmu, varNames_thisBL, ...
    list2map, which_file, figure_visibility)

% warning('TROUBLESHOOTING')
% deltas_emu_xvt = ones(size(deltas_emu_xvt)) ;

if any(any(any(isinf(deltas_emu_xvt))))
    error('Infinite deltas were supposed to have been fixed by now!')
end

Ntpers = size(deltas_emu_xvt,3) ;

% Set up output structure
data_fu_out.list2map = list2map ;
data_fu_out.lonlats = data_bl_thisBL.lonlats ;
data_fu_out.varNames = varNames_thisBL ;
data_fu_out.y1s = data_fu_emu.y1s ;
data_fu_out.yNs = data_fu_emu.yNs ;
if isfield(data_bl_thisBL, 'actually_emu_char')
    data_fu_out.actually_emu_char = data_bl_thisBL.actually_emu_char ;
end

% Get corresponding array indices
tmp_i_thisBL = nan(size(data_fu_out.varNames)) ;
tmp_i_emu = nan(size(data_fu_out.varNames)) ;
for v = 1:length(tmp_i_emu)
    
% % %     % TROUBLESHOOTING
% % %     tmp_x = data_bl_thisBL.garr_xv(:,v) ;
% % %     tmp_x(tmp_x>10 & ~isinf(tmp_x)) = 10 ;
% % %     figure; hist(tmp_x(tmp_x>0 & ~isinf(tmp_x)));
% % %     pause(3)
% % %     close
% % %     continue
    
    thisVar_thisBL = data_fu_out.varNames{v} ;
    thisCrop_thisBL = getbasename(thisVar_thisBL) ;
    thisN = thisVar_thisBL(regexp(thisVar_thisBL,'\d\d\d$'):end) ; % regexp: Get last 3 digits (must also be last 3 characters)
    thisI = strrep(strrep(thisVar_thisBL,thisCrop_thisBL,''), thisN, '') ;
%     fprintf('%s = %s + %s + %s\n', thisVar_thisBL, thisCrop_thisBL, thisI, thisN)
    thisCrop_emu = cropList_thisBL_asEmu{strcmp(cropList_thisBL, thisCrop_thisBL)} ;
    thisVar_emu = sprintf('%s%s%s', thisCrop_emu, thisI, thisN) ;
    
    ii = find(strcmp(data_bl_thisBL.varNames, thisVar_thisBL)) ;
    if length(ii) ~= 1
        error('Error finding %s in data_bl_thisBL.varNames (%d found)', thisVar_thisBL, length(ii)) ;
    end
    tmp_i_thisBL(v) = ii ;
    
    ii = find(strcmp(data_fu_emu.varNames, thisVar_emu)) ;
    if length(ii) ~= 1
        error('Error finding %s in data_fu_emu.varNames (%d found)', thisVar_emu, length(ii)) ;
    end
    tmp_i_emu(v) = ii ;
    
%     % TROUBLESHOOTING
%     fprintf('%s <-- %s (ind emu %d)\n', thisVar_thisBL, thisVar_emu, tmp_i_emu(v)) ;
end
if length(tmp_i_thisBL) ~= length(unique(tmp_i_thisBL))
    error('length(tmp_i_thisBL) ~= length(unique(tmp_i_thisBL))')
end

% Actually apply the deltas
deltas_tmp_xvt = deltas_emu_xvt(:,tmp_i_emu,:) ;
byield_tmp_xvt = repmat(data_bl_thisBL.garr_xv(:,tmp_i_thisBL), [1 1 Ntpers]) ;
positive_delta_on_zero_baseline_xv = any(byield_tmp_xvt==0 & deltas_tmp_xvt>0, 3) ;
if any(any(positive_delta_on_zero_baseline_xv))
    warning('%d cell-variables have positive deltas applied to zero baseline', ...
        length(find(positive_delta_on_zero_baseline_xv)))
end
data_fu_out.garr_xvt = deltas_tmp_xvt .* byield_tmp_xvt ;

% Sort variable names
[data_fu_out.varNames, I] = sort(data_fu_out.varNames) ;
data_fu_out.garr_xvt = data_fu_out.garr_xvt(:,I,:) ;
clear I

% Test for consistent deltas
a = data_fu_out.garr_xvt(:,strcmp(data_fu_out.varNames,'CerealsC4200'),4) ./ ...
    data_bl_thisBL.garr_xv(:,strcmp(data_bl_thisBL.varNames,'CerealsC4200'),:) ;
b = data_fu_emu.garr_xvt(:,strcmp(data_fu_emu.varNames,'maize'),4) ./ ...
    data_bl_emu.garr_xv(:,strcmp(data_bl_emu.varNames,'maize'),:) ;
a2 = a(~isnan(a) & ~isnan(b) & ~isinf(a) & ~isinf(b)) ;
b2 = b(~isnan(a) & ~isnan(b) & ~isinf(a) & ~isinf(b)) ;
if ~isequaln(round(a2, 12), round(b2, 12))
    error('Inconsistent deltas in test!')
end
clear a b a2 b2


end



function make_fig(y_thisBL_bl, y_thisBL_fu, y_out_fu, ...
    varName, y1, yN)

figure('Color', 'w', 'Position', [1 34 720 771], 'Visible', figure_visibility) ;
xlims = [0 0] ;

h1 = subplot(3,1,1) ;
hist(y_thisBL_bl);
xlims = [min(xlims(1), h1.XLim(1)) max(xlims(2), h1.XLim(2))] ;
title(sprintf('LPJ BL: %s (max %0.3f)', varName, max(y_thisBL_bl)))

% h2 = subplot(3,1,2) ;
% hist(y_thisBL_fu);
% xlims = [min(xlims(1), h2.XLim(1)) max(xlims(2), h2.XLim(2))] ;
% title(sprintf('LPJ %d-%d: %s (max %0.3f)', y1, yN, varName, max(y_thisBL_fu)))

h2 = subplot(3,1,2) ;
ygt0 = y_thisBL_fu(y_thisBL_fu>0) ;
scaled_MAD = (-1/(sqrt(2)*erfcinv(3/2)))*median(abs(ygt0-median(ygt0))) ;
upper_outl_thresh = median(ygt0) + 3*scaled_MAD ;
is_upper_outlier = ygt0 > upper_outl_thresh ;
hist(y_thisBL_fu);
xlims = [min(xlims(1), h2.XLim(1)) max(xlims(2), h2.XLim(2))] ;
title(sprintf('Output %d-%d: %s (%d outliers>%0.3f)', ...
    y1, yN, ...
    varName, ...
    length(find(is_upper_outlier)), ...
    upper_outl_thresh))

h3 = subplot(3,1,3) ;
ygt0 = y_out_fu(y_out_fu>0) ;
scaled_MAD = (-1/(sqrt(2)*erfcinv(3/2)))*median(abs(ygt0-median(ygt0))) ;
upper_outl_thresh = median(ygt0) + 3*scaled_MAD ;
is_upper_outlier = ygt0 > upper_outl_thresh ;
hist(y_out_fu);
xlims = [min(xlims(1), h3.XLim(1)) max(xlims(2), h3.XLim(2))] ;
title(sprintf('Output %d-%d: %s (%d outliers>%0.3f)', ...
    y1, yN, ...
    varName, ...
    length(find(is_upper_outlier)), ...
    upper_outl_thresh))

h1.XLim = xlims ;
h2.XLim = xlims ;
h3.XLim = xlims ;


end










