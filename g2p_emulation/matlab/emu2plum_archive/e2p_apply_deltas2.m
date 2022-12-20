function [data_fu_out, data_bl_thisBL, data_bl_emu, data_fu_emu, deltas_emu_xvt] ...
    = e2p_apply_deltas2( ...
    data_bl_thisBL, data_bl_emu, data_fu_emu, deltas_emu_xvt, ...
    cropList_thisBL, cropList_thisBL_asEmu, varNames_thisBL, ...
    list2map, which_file, figure_visibility, interp_infs)

% warning('TROUBLESHOOTING')
% deltas_emu_xvt = ones(size(deltas_emu_xvt)) ;

Ntpers = size(deltas_emu_xvt,3) ;

% Get corresponding array indices
Nvars = length(data_bl_thisBL.varNames) ;
tmp_i_thisBL = nan(Nvars, 1) ;
tmp_i_emu = nan(Nvars, 1) ;
for v = 1:length(tmp_i_emu)
    
% % %     % TROUBLESHOOTING
% % %     tmp_x = data_bl_thisBL.garr_xv(:,v) ;
% % %     tmp_x(tmp_x>10 & ~isinf(tmp_x)) = 10 ;
% % %     figure; hist(tmp_x(tmp_x>0 & ~isinf(tmp_x)));
% % %     pause(3)
% % %     close
% % %     continue
    
    thisVar_thisBL = data_bl_thisBL.varNames{v} ;
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

% Rearrange
data_bl_thisBL.garr_xv = data_bl_thisBL.garr_xv(:,tmp_i_thisBL) ;
data_bl_thisBL.varNames = data_bl_thisBL.varNames(tmp_i_thisBL) ;
data_bl_thisBL.actually_emu_char = data_bl_thisBL.actually_emu_char(tmp_i_thisBL) ;
data_bl_emu.garr_xv = data_bl_emu.garr_xv(:,tmp_i_emu) ;
data_bl_emu.varNames = data_bl_emu.varNames(tmp_i_emu) ;
data_fu_emu.garr_xvt = data_fu_emu.garr_xvt(:,tmp_i_emu,:) ;
data_fu_emu.varNames = data_fu_emu.varNames(tmp_i_emu) ;
deltas_emu_xvt = deltas_emu_xvt(:,tmp_i_emu,:) ;

% Set up output structure
data_fu_out.list2map = list2map ;
data_fu_out.lonlats = data_bl_thisBL.lonlats ;
data_fu_out.varNames = varNames_thisBL ;
data_fu_out.y1s = data_fu_emu.y1s ;
data_fu_out.yNs = data_fu_emu.yNs ;
if isfield(data_bl_thisBL, 'actually_emu_char')
    data_fu_out.actually_emu_char = data_bl_thisBL.actually_emu_char ;
end

% Actually apply the deltas
data_fu_out.garr_xvt = deltas_emu_xvt .* repmat(data_bl_thisBL.garr_xv, [1 1 Ntpers]) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Deal with 0 baseline --> positive future (results in delta=Inf) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Might want to instead limit to (e.g.) 99.9th percentile of deltas

% Find cells that need fixing
data_fu_emu_xvMax = max(data_fu_emu.garr_xvt, [], 3) ;
[isbad_deltas_xv, isbad_thisBL_xv] = getbad(deltas_emu_xvt, data_bl_thisBL, data_fu_emu_xvMax) ;

% Warn for cell-variables where a positive delta will be applied to a zero
% baseline: Not fixing these! (At least, not yet.)
isbad_thisBL_notDeltas_xv = isbad_thisBL_xv & ~isbad_deltas_xv ;
Nbad_thisBL_notDeltas = length(find(isbad_thisBL_notDeltas_xv)) ;
if Nbad_thisBL_notDeltas > 0
    warning('Not fixing 0 used baseline --> positive deltas for %d cell-variables', Nbad_thisBL_notDeltas)
end
isbad_thisBL_andDeltas_xv = isbad_thisBL_xv & isbad_deltas_xv ;


Nbad_deltas = length(find(isbad_deltas_xv)) ;
if Nbad_deltas > 0 && ~interp_infs
    warning('%d cell-variables in data_fu_emu.garr_xvt have future positive but baseline emulator zero, resulting in delta=Inf. Will NOT fix.', ...
        Nbad_deltas)
elseif Nbad_deltas > 0 && interp_infs
    warning('%d cell-variables in data_fu_emu.garr_xvt have future positive but baseline emulator zero, resulting in delta=Inf. Fixing...', ...
        Nbad_deltas)
    t = 0 ;
    alreadyFixed_xv = false(size(isbad_deltas_xv)) ;
    while any(any(isbad_deltas_xv))
        t = t + 1 ;
        if t > Ntpers
            error('???')
        end
        
        for v = 1:Nvars
                        
            data_fu_emu_x1fut = data_fu_emu.garr_xvt(:,v,t:end) ;
            data_fu_emu_x1this = data_fu_emu_x1fut(:,:,1) ;
            fixnow_x = isbad_deltas_xv(:,v) & data_fu_emu_x1this > 0;
            if any(fixnow_x)
                fprintf('t %d, v %d: %d bad\n', t, v, length(find(fixnow_x)))
                data_bl_thisBL_x = data_bl_thisBL.garr_xv(:,v) ;
                actual_NaN_x = isnan(data_bl_thisBL_x) ;
                
                %             fixnow_x(badcell)
                
                % Sanity check
                if any(alreadyFixed_xv(fixnow_x,v))
                    tmp = find(fixnow_x) ;
                    I = find(alreadyFixed_xv(fixnow_x,v), 1) ;
                    error('This cell (%d) should already have been fixed!', tmp(I))
                end
                
                if any(any(isinf(data_fu_emu_x1fut(fixnow_x,:))))
                    error('any(any(isinf(data_fu_emu_x1fut(fixnow_x,:))))')
                end
                
                % Set deltas relative to this year
                deltas_tmp_x1fut = data_fu_emu_x1fut ./ data_fu_emu_x1this ;
                deltas_emu_xvt(fixnow_x,v,t:end) = deltas_tmp_x1fut(fixnow_x,:,:) ;
                
                % Get original map
                orig_YX = nan(360,720) ;
                data_bl_thisBL_x(fixnow_x) = NaN ;
                orig_YX(list2map) = data_bl_thisBL_x ;
                
                % Set baseline as interpolated value for this timestep
                intp_YX = inpaint_nans(orig_YX, 4) ;
                data_bl_thisBL_x = intp_YX(list2map) ;
                data_bl_thisBL_x(actual_NaN_x) = NaN ;
                data_bl_thisBL.garr_xv(:,v) = data_bl_thisBL_x ;
                
                % Update isbad_deltas_xv
                isbad_deltas_xv = ...
                    getbad(deltas_emu_xvt, data_bl_thisBL, data_fu_emu_xvMax) ;
                alreadyFixed_xv(fixnow_x,v) = true ;
            end
            
        end
    end
    
    % Make sure you've eliminated all infinite deltas
    if any(any(any(isinf(deltas_emu_xvt))))
        error('Inf remaining in deltas_emu_xvt')
    end
    
    % Warn if any de-Inf'ed deltas are applied to zero baselines
    deltas_tmp_xvMax = max(deltas_emu_xvt, [], 3) ;
    Nweird = length(find(isbad_thisBL_andDeltas_xv & deltas_tmp_xvMax > 0 & data_bl_thisBL.garr_xv == 0)) ;
    if Nweird > 0
        warning('%d de-Inf''ed cell-variables have positive deltas applied to zero baselines. This could be an acceptable result (interpolated baseline value == 0), or it could indicate a problem.', ...
            Nweird)
    end    
    
    data_fu_out.garr_xvt = deltas_emu_xvt .* repmat(data_bl_thisBL.garr_xv, [1 1 Ntpers]) ;
end

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


function [isbad_deltas_xv, isbad_thisBL_xv] = getbad(deltas_emu_xvt, data_bl_thisBL, data_fu_emu_xvMax)

isbad_deltas_xv = any(isinf(deltas_emu_xvt), 3) ;
isbad_thisBL_xv = data_bl_thisBL.garr_xv==0 & data_fu_emu_xvMax > 0 ;

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










