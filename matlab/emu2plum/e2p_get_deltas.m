function deltas_emu_xvt = e2p_get_deltas(...
    data_bl_emu, data_fu_emu, interp_infs, cropList_emu, ...
    getbasename, getbasenamei, which_file, ...
    used_emuCrops, list2map)

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
            clear orig_YX intp_YX
        end
    end
elseif isempty(isbad)
    disp('No values of data_fu_emu.garr_xvt are positive but were 0 in baseline.')
    deltas_emu_xvt = deltas0_emu_xvt ;
end
clear emu_bl_xvt


end