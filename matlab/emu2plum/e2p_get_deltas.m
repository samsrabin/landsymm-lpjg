function deltas_emu_xvt = e2p_get_deltas(...
    data_bl_emu, data_fu_emu, interp_infs, cropList_emu, ...
    getbasename, getbasenamei, which_file, ...
    used_emuCrops, list2map, ...
    save_interp_figs, outDir_interp_figs, ggcm, figure_visibility)

verbose = false ;

% Get deltas. Produces NaNs where excluded, because exclusions were applied
% above.
Ntpers = size(data_fu_emu.garr_xvt,3) ;
emu_bl_xvt = repmat(data_bl_emu.garr_xv, [1 1 Ntpers]) ;
deltas0_emu_xvt = data_fu_emu.garr_xvt ./ emu_bl_xvt ;

e2p_check_correct_zeros(deltas0_emu_xvt, which_file, getbasenamei(data_fu_emu.varNames))

% Set delta=0 where both emulated baseline and future had 0. Note, this
% only applies for ZERO, not NaN.
deltas0_emu_xvt(emu_bl_xvt==0 & data_fu_emu.garr_xvt==0) = 0 ;

e2p_check_correct_zeros(deltas0_emu_xvt, which_file, getbasenamei(data_fu_emu.varNames))

% Deal with 0 baseline --> positive future (results in delta=Inf)
% Might want to instead limit to (e.g.) 99.9th percentile of deltas
isbad = find(data_fu_emu.garr_xvt>0 & emu_bl_xvt==0) ;

if ~isempty(isbad) && ~interp_infs
    warning('%d elements of data_fu_emu.garr_xvt are positive but were 0 in baseline, resulting in delta=Inf. Will NOT fix.', ...
        length(find(isbad)))
    deltas_emu_xvt = deltas0_emu_xvt ;
elseif ~isempty(isbad) && interp_infs
    fprintf('Interpolating %d elements of deltas_emu_xvt where baseline was 0 but future is positive...\n', ...
        length(find(isbad)))
    deltas_emu_xvt = deltas0_emu_xvt ;
    for v = 1:length(data_fu_emu.varNames)
        thisVar_emu = data_fu_emu.varNames{v} ;
        if ~used_emuCrops(strcmp(cropList_emu,getbasename(thisVar_emu)))
            if verbose
                fprintf('    %s: Skipping (not needed)\n', thisVar_emu) ;
            end
            continue
        elseif ~any(any(isinf(deltas0_emu_xvt(:,v,:))))
            if verbose
                fprintf('    %s: Skipping (no Inf)\n', thisVar_emu) ;
            end
            continue
        end
        if verbose
            fprintf('    %s:\n', thisVar_emu) ;
        end
        for t = 1:Ntpers
            
            tmp_x = deltas0_emu_xvt(:,v,t) ;
            
            if ~any(tmp_x>0)
                error('No positive values!')
            end
            
% % %             % TROUBLESHOOTING
% % %             tmp_x(tmp_x>10 & ~isinf(tmp_x)) = 10 ;
% % %             figure; hist(tmp_x(tmp_x>0 & ~isinf(tmp_x)));
% % %             pause(3)
% % %             close
            
% % %             % TROUBLESHOOTING
% % %             fprintf('%0.2f, %0.2f\n', ...
% % %                 prctile(tmp_x(tmp_x>0 & ~isinf(tmp_x)),99), ...
% % %                 prctile(tmp_x(tmp_x>0 & ~isinf(tmp_x)),99.9)) ;
            
            % Skip if no Infs
            if ~any(isinf(tmp_x))
                if verbose
                    fprintf('        %d/%d skipped\n', t, Ntpers) ;
                end
                continue
            else
                if verbose
                    fprintf('        %d/%d (%d bad)...\n', t, Ntpers, length(find(isinf(tmp_x)))) ;
                end
            end
            
            % Record cells that were already NaN
            already_nan = isnan(tmp_x) ;
            
            % Make map of original
            orig_YX = nan(360,720) ;
            orig_YX(list2map) = tmp_x ;
            
            % Interpolate over Infs
            tmp_x(already_nan) = 0 ;
            tmp_x(isinf(tmp_x)) = NaN ;
            tmp_YX = nan(360,720) ;
            tmp_YX(list2map) = tmp_x ;
            intp_YX = inpaint_nans(tmp_YX, 4) ;
            
%             % TROUBLESHOOTING
%             intp_YX(isnan(orig_YX)) = NaN ;
%             shademap(orig_YX) ;
%             shademap(intp_YX) ;
%             stop
            
            tmp_x = intp_YX(list2map) ;
            if ~any(tmp_x>0)
                error('No positive values!')
            end
            
            tmp_x(already_nan) = NaN ;
            if ~any(~isnan(tmp_x))
                error('No non-NaN values!')
            end
            
            deltas_emu_xvt(:,v,t) = tmp_x ;
            e2p_check_correct_zeros(deltas_emu_xvt(:,:,t), which_file, getbasenamei(data_fu_emu.varNames))
            
            if save_interp_figs && t==Ntpers
                
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
                    '%s:\n%s\n(last timestep)', ...
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
                export_fig(filename, '-r150') ;
                close
            end
            
            clear orig_YX intp_YX
        end
    end
elseif isempty(isbad)
    disp('No values of data_fu_emu.garr_xvt are positive but were 0 in baseline.')
    deltas_emu_xvt = deltas0_emu_xvt ;
end
clear emu_bl_xvt


end