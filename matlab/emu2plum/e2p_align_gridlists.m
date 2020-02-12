function [data_bl_lpj, data_fu_lpj, data_bl_emu, data_fu_emu, list2map] = ...
    e2p_align_gridlists(data_bl_lpj, data_fu_lpj, data_bl_emu, data_fu_emu)

if ~isequal(data_bl_lpj.list2map, data_fu_lpj.list2map)
    error('Mismatch between baseline and future gridlists: LPJ-GUESS')
elseif ~isequal(data_bl_emu.list2map, data_fu_emu.list2map)
    error('Mismatch between baseline and future gridlists: Emulator')
end

% Skip all this if gridlists already align
if isequal(data_bl_lpj.list2map, data_bl_emu.list2map)
    list2map = data_bl_lpj.list2map ;

else
    
    [list2map, gridlist_IA, gridlist_IB] = intersect(data_bl_lpj.list2map, data_bl_emu.list2map) ;
    
    data_bl_lpj.garr_xv = data_bl_lpj.garr_xv(gridlist_IA,:) ;
    data_bl_lpj.list2map = list2map ;
    data_bl_lpj.lonlats = data_bl_lpj.lonlats(gridlist_IA,:) ;
    
    data_fu_lpj.garr_xvt = data_fu_lpj.garr_xvt(gridlist_IA,:,:) ;
    data_fu_lpj.list2map = list2map ;
    data_fu_lpj.lonlats = data_fu_lpj.lonlats(gridlist_IA,:) ;
    
    data_bl_emu.garr_xv = data_bl_emu.garr_xv(gridlist_IB,:) ;
    data_bl_emu.list2map = list2map ;
    data_bl_emu.lonlats = data_bl_emu.lonlats(gridlist_IB,:) ;
    
    data_fu_emu.garr_xvt = data_fu_emu.garr_xvt(gridlist_IB,:,:) ;
    data_fu_emu.list2map = list2map ;
    data_fu_emu.lonlats = data_fu_emu.lonlats(gridlist_IB,:) ;
end


end