function deltas_emu_xvt = e2p_get_deltas2(...
    data_bl_emu, data_fu_emu, interp_infs, cropList_emu, ...
    which_file, ...
    used_emuCrops, list2map, ...
    save_interp_figs, outDir_interp_figs, ggcm, figure_visibility, ...
    when_remove_outliers, outDir_ggcm, renderer)

verbose = false ;

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
isbad_xvt = data_fu_emu.garr_xvt>0 & emu_bl_xvt==0 ;
isbad = find(isbad_xvt) ;
Nbad = length(isbad) ;

% Sanity check
if any(any(any(isbad_xvt & ~isinf(deltas0_emu_xvt)))) ...
|| any(any(any(~isbad_xvt & isinf(deltas0_emu_xvt))))
    error('???')
end

% Sanity check
if ~isequal(size(shiftdim(cropList_emu)), size(shiftdim(used_emuCrops)))
    error('Size mismatch between cropList_emu and used_emuCrops')
end

Ntpers = size(deltas0_emu_xvt, 3) ;

deltas_emu_xvt = deltas0_emu_xvt ;

% % Have to reshape so you can work with booleans across what would otherwise
% % be multiple dimensions
% deltas_emu_wt = reshape(deltas0_emu_xvt, [prod(size(deltas0_emu_xvt,1:2)) Ntpers]) ;
% emu_fu_wt = reshape(data_fu_emu.garr_xvt, [prod(size(deltas0_emu_xvt,1:2)) Ntpers]) ;
% isbad_wt = reshape(isbad_xvt, [prod(size(deltas0_emu_xvt,1:2)) Ntpers]) ;
% 
% if Nbad > 0
%     fprintf('%d cell-variables have baseline 0 but future positive. Deltas will be relative to first non-zero timestep.\n', ...
%         length(find(isbad)))
%     
%     t = 0 ;
%     while any(isinf(deltas_emu_wt))
%         t = t + 1
%         if t > Ntpers
%             error('Inf(s) remaining???')
%         end
%         
%         isbad_w = isbad_wt(:,t) ;
%         
%         % Set subsequent deltas relative to this time period
%         % (and set this time period's delta to 1)
%         deltas_emu_wt(isbad_w,t:end) = emu_fu_wt(isbad_w,t:end) ...
%             ./ repmat(deltas_emu_wt(isbad_w,t), [1 Ntpers-t+1]) ;
%         
%         % Make sure no infinite deltas remain for this time period or
%         % previous
%         if any(any(isinf(deltas_emu_wt(:,1:t))))
%             error('How do you still have infinite deltas for this timestep and/or previous?')
%         end
%         
%     end
% elseif verbose && isempty(isbad)
%     disp('No values of data_fu_emu.garr_xvt are positive but were 0 in baseline.')
% end
% 
% deltas_emu_xvt = reshape(deltas_emu_wt, size(deltas0_emu_xvt)) ;
% 
% % Make sure no infinite deltas remain
% if any(any(any(isinf(deltas_emu_xvt))))
%     error('How do you still have infinite deltas? Try getting rid of this code in this function and moving it to apply function')
% end

end
