function deltas_emu_xvt = e2p_get_deltas(...
    data_bl_emu, data_fu_emu, interp_infs, cropList_emu, ...
    which_file, ...
    used_emuCrops, ...
    save_interp_figs, outDir_interp_figs, ggcm, figure_visibility, ...
    when_remove_outliers, outDir_ggcm, renderer)

verbose = false ;

list2map = data_bl_emu.list2map ;

% Get deltas. Produces NaNs where excluded, because exclusions were applied
% above.
Ntpers = size(data_fu_emu.garr_xvt,3) ;
emu_bl_xvt = repmat(data_bl_emu.garr_xv, [1 1 Ntpers]) ;
deltas0_emu_xvt = data_fu_emu.garr_xvt ./ emu_bl_xvt ;

% Set delta=0 where both emulated baseline and future had 0. Note, this
% only applies for ZERO, not NaN.
deltas0_emu_xvt(emu_bl_xvt==0 & data_fu_emu.garr_xvt==0) = 0 ;

% Deal with 0 baseline --> positive future (results in delta=Inf)
% Might want to instead limit to (e.g.) 99.9th percentile of deltas
isbad_xvt = isinf(deltas0_emu_xvt) ;
isbad_vt = squeeze(any(isbad_xvt, 1)) ;
isbad = find(isbad_xvt) ;
Nbad = length(isbad) ;

% Sanity check
if ~isequal(size(shiftdim(cropList_emu)), size(shiftdim(used_emuCrops)))
    error('Size mismatch between cropList_emu and used_emuCrops')
end

if Nbad > 0 && ~interp_infs
    warning('%d elements of data_fu_emu.garr_xvt are positive but were 0 in baseline, resulting in delta=Inf. Will NOT fix.', ...
        length(find(isbad)))
    deltas_emu_xvt = deltas0_emu_xvt ;
elseif Nbad > 0 && interp_infs
    fprintf('Interpolating %d elements of deltas_emu_xvt where baseline was 0 but future is positive (%d interpolations)...\n', ...
        length(find(isbad)), length(find(isbad_vt)))
    
%     if Nbad > 100
%         warning('Large number of interpolations. Could you make this more efficient by using each cell''s first interpolated value as its "baseline"?')
%     end
    
    if strcmp(when_remove_outliers, 'before_interp')
        disp('    Removing outliers...')
        
        data_tmp.garr_xvt = deltas0_emu_xvt ;
        is_inf_xvt = isinf(data_tmp.garr_xvt) ;
        data_tmp.garr_xvt(is_inf_xvt) = 0 ;
        data_tmp.varNames = data_fu_emu.varNames ;
        [data_tmp, outlier_info] = e2p_remove_outliers(data_tmp, which_file) ;
        data_tmp.garr_xvt(is_inf_xvt) = Inf ;
        deltas0_emu_xvt = data_tmp.garr_xvt ;
        clear data_tmp
        e2p_check_correct_zeros(deltas0_emu_xvt, ...
            which_file, data_fu_emu.varNames, ...
            'Future', @getbasenamei)
                
        % Save outlier info
        e2p_save_outlier_info(outlier_info, outDir_ggcm, which_file, data_fu_emu.y1s, data_fu_emu.yNs)
    end
    
    
    
    deltas_emu_xvt = deltas0_emu_xvt ;
    for v = 1:length(data_fu_emu.varNames)
        thisVar_emu = data_fu_emu.varNames{v} ;
        if ~used_emuCrops(strcmp(cropList_emu,getbasename(thisVar_emu)))
            if verbose
                fprintf('    %s: Skipping (not needed)\n', thisVar_emu) ;
            end
            continue
        elseif ~any(isbad_vt(v,:))
            if verbose
                fprintf('    %s: Skipping (no Inf)\n', thisVar_emu) ;
            end
            continue
        end
        if verbose
            fprintf('    %s:\n', thisVar_emu) ;
        end
                
        for t = 1:Ntpers
            
            fixnow_x = isbad_xvt(:,v,t) ;
            if ~any(fixnow_x)
                continue
            end
            
            deltas0_emu_x = deltas0_emu_xvt(:,v,t) ;
            
            if ~any(deltas0_emu_x>0)
                error('No positive values!')
            end
            
% % %             % TROUBLESHOOTING
% % %             deltas0_emu_x(deltas0_emu_x>10 & ~isinf(deltas0_emu_x)) = 10 ;
% % %             figure; hist(deltas0_emu_x(deltas0_emu_x>0 & ~isinf(deltas0_emu_x)));
% % %             pause(3)
% % %             close
            
% % %             % TROUBLESHOOTING
% % %             fprintf('%0.2f, %0.2f\n', ...
% % %                 prctile(deltas0_emu_x(deltas0_emu_x>0 & ~isinf(deltas0_emu_x)),99), ...
% % %                 prctile(deltas0_emu_x(deltas0_emu_x>0 & ~isinf(deltas0_emu_x)),99.9)) ;
            
            % There should be infs here!
            if ~any(isinf(deltas0_emu_x))
                error('Why are there no Inf values here?')
            end
            if verbose
                fprintf('        %d/%d (%d bad)...', t, Ntpers, length(find(isinf(deltas0_emu_x)))) ; %#ok<UNRCH>
            end
            
            % Record cells that were already NaN
            already_nan = isnan(deltas0_emu_x) ;
            
            % Interpolate over Infs
            deltas0_emu_x(already_nan) = 0 ;
            deltas0_emu_x(isinf(deltas0_emu_x)) = NaN ;
            tmp_YX = nan(360,720) ;
            tmp_YX(list2map) = deltas0_emu_x ;
            intp_YX = inpaint_nans(tmp_YX, 4) ;
            deltas_emu_x = intp_YX(list2map) ;
            if ~any(deltas0_emu_x>0)
                error('No positive values!')
            end
            
%             % TROUBLESHOOTING
%             new_YX = nan(360,720) ;
%             new_YX(list2map) = deltas0_emu_x ;
%             shademap(orig_YX) ; title('original')
%             shademap(new_YX) ; title('interpolated')
            
            % Fill back in NaNs
            deltas_emu_x(already_nan) = NaN ;
            if ~any(~isnan(deltas_emu_x))
                error('No non-NaN values!')
            end
            
            % Put into output array
            deltas_emu_xvt(:,v,t) = deltas_emu_x ;
            
            % Fix any subsequent infinite deltas for cells we just fixed
            if t < Ntpers
                deltas0_emu_x1Fut = deltas0_emu_xvt(:,v,t+1:end) ;
                Nbad_fut_toFix = length(find(isinf(deltas0_emu_x1Fut(fixnow_x,:,:)))) ;
                if Nbad_fut_toFix > 0
                    if verbose
                        fprintf(' + %d in future timesteps', Nbad_fut_toFix)
                    end
                    deltas0_emu_x1Fut(fixnow_x,:,:) = ...
                        deltas_emu_x(fixnow_x) ...
                        .* data_fu_emu.garr_xvt(fixnow_x,v,t+1:end) ...
                        ./ data_fu_emu.garr_xvt(fixnow_x,v,t) ;
                    deltas_emu_xvt(fixnow_x,v,t+1:end) = deltas0_emu_x1Fut(fixnow_x,:,:) ;
                end
                if verbose
                    fprintf('\n')
                end
            end
            
            % Update check for infinite deltas
            isbad_xvt = isinf(deltas_emu_xvt) ;
            
            clear intp_YX
            if ~any(any(any(isbad_xvt)))
                break
            end
        end
        
        if save_interp_figs
            
            % Make maps
            orig_YX = nan(360,720) ;
            orig_YX(list2map) = max(deltas0_emu_xvt(:,v,:), [], 3) ;
            intp_YX = nan(360,720) ;
            intp_YX(list2map) = max(deltas_emu_xvt(:,v,:), [], 3) ;
            
            figure('Color', 'w', 'Position', figurePos, 'Visible', figure_visibility) ;
            spacing = [0.05 0.01] ; % v h
            fontSize = 14 ;
            ybnds = 65:360 ;
            
            % Map highlighting Inf cells
            map_infs_YX = double(isinf(orig_YX)) ;
            map_infs_YX(isnan(orig_YX)) = NaN ;
            subplot_tight(1,3,1,spacing)
            pcolor(map_infs_YX(ybnds,:)); shading flat; axis equal tight off
            colormap(gca, flipud(colormap('parula')))
            title('Inf cells (blue)')
            set(gca, 'FontSize', fontSize)
            
            % Map before interpolation
            map_before_YX = orig_YX ;
            map_before_YX(isinf(orig_YX)) = NaN ;
            h1 = subplot_tight(2,3,2:3,spacing) ;
            pcolor(map_before_YX(ybnds,:)); shading flat; axis equal tight off
            colormap(gca, 'jet')
            hcb = colorbar ;
            ylabel(hcb, '\Delta (future/baseline)')
            title('Before interp. (Infs are white)')
            set(gca, 'FontSize', fontSize)
            
            % Map after interpolation
            map_after_YX = intp_YX ;
            map_after_YX(isnan(orig_YX)) = NaN ;
            h2 = subplot_tight(2,3,5:6,spacing) ;
            pcolor(map_after_YX(ybnds,:)); shading flat; axis equal tight off
            colormap(gca, 'jet')
            hcb = colorbar ;
            ylabel(hcb, '\Delta (future/baseline)')
            title('After interp.')
            set(gca, 'FontSize', fontSize)
            
            % Standardize caxis
            tmp_beforeafter = cat(3, map_before_YX, map_after_YX) ;
            new_caxis = [ ...
                min(min(min(tmp_beforeafter))), ...
                max(max(max(tmp_beforeafter)))] ;
            caxis(h1, new_caxis) ;
            caxis(h2, new_caxis) ;
            
            % Add info
            hold on; ax2 = axes(); hold off
            ax2.Position = [0 0 1 1] ;
            ax2.Visible = 'off' ;
            info_text = sprintf( ...
                '%s:\n%s\n(maximum)', ...
                ggcm, ...
                strrep(thisVar_emu, '_', '\_')) ;
            text(ax2, 0.1, 0.9, ...
                info_text, ...
                'FontSize', fontSize*2, ...
                'FontWeight', 'bold')
            
            % Save figure
            if ~exist(outDir_interp_figs, 'dir')
                mkdir(outDir_interp_figs)
            end
            filename = sprintf('%s/intp_%s_%s.png', ...
                outDir_interp_figs, ggcm, ...
                strrep(thisVar_emu, '_', '')) ;
            export_fig(filename, '-r150', renderer) ;
            close
        end
        
        clear orig_YX
    end
else
    if verbose
        disp('No values of data_fu_emu.garr_xvt are positive but were 0 in baseline.')
    end
    deltas_emu_xvt = deltas0_emu_xvt ;
end
clear emu_bl_xvt


end
