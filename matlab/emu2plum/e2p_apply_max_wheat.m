function Sout = e2p_apply_max_wheat(Sin, outDir)

Sout = Sin ;

actually_emu_char = {} ;
if isfield(Sin, 'actually_emu')
    actually_emu_char = cell(size(Sin.actually_emu)) ;
    actually_emu_char(~Sin.actually_emu) = {'sim'} ;
    actually_emu_char(Sin.actually_emu) = {'*EMU*'} ;
end

% Load existing info
is_ww_max_file = sprintf('%s/is_ww_max.mat', outDir) ;
load(is_ww_max_file, 'is_ww_max_bl_gW', 'is_ww_max_fu_gWt', 'winter_wheats') ;

Nww = length(winter_wheats) ;
if isfield(Sin,'garr_xv')
    out_xvt = Sin.garr_xv ;
elseif isfield(Sin,'garr_xvt')
    out_xvt = Sin.garr_xvt ;
else
    error('???')
end
for w = 1:Nww
    
    % Get indices
    thisWW = winter_wheats{w} ;
    [i_thisWW, i_thisSW, i_thisMW, Sout.varNames] = ...
        e2p_wheatInds(thisWW, Sout.varNames) ;
    
    % Get winter and spring wheat values
    if isfield(Sin,'garr_xv')
        tmpWW_x1t_in = Sin.garr_xv(:,i_thisWW) ;
        tmpSW_x1t_in = Sin.garr_xv(:,i_thisSW) ;
        is_ww_max_gWt = is_ww_max_bl_gW ;
    elseif isfield(Sin,'garr_xvt')
        tmpWW_x1t_in = Sin.garr_xvt(:,i_thisWW,:) ;
        tmpSW_x1t_in = Sin.garr_xvt(:,i_thisSW,:) ;
        is_ww_max_gWt = is_ww_max_fu_gWt ;
    else
        error('???')
    end
    
    % Assign values from max wheat
    tmp_x1t_out = nan(size(out_xvt,1), 1, size(out_xvt,3)) ;
    tmp_x1t_out(is_ww_max_gWt(:,w,:)) = tmpWW_x1t_in(is_ww_max_gWt(:,w,:)) ;
    tmp_x1t_out(~is_ww_max_gWt(:,w,:)) = tmpSW_x1t_in(~is_ww_max_gWt(:,w,:)) ;
    out_xvt(:,i_thisMW,:) = tmp_x1t_out ;
    
    % Assign character-based designator for whether it was simulated or
    % emulated
    if isfield(Sin, 'actually_emu')
        if Sin.actually_emu(i_thisWW) == Sin.actually_emu(i_thisSW)
            actually_emu_char{i_thisMW} = actually_emu_char{i_thisWW} ; %#ok<AGROW>
        else
            actually_emu_char{i_thisMW} = ...
                [actually_emu_char{i_thisWW} actually_emu_char{i_thisSW}] ; %#ok<AGROW>
        end
    end
end

if isfield(Sin, 'actually_emu')
    Sout.actually_emu_char = actually_emu_char ;
end

if isfield(Sin,'garr_xv')
    Sout.garr_xv = out_xvt ;
elseif isfield(Sin,'garr_xvt')
    Sout.garr_xvt = out_xvt ;
else
    error('???')
end

end