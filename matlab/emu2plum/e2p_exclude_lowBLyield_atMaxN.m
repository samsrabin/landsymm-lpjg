function exclude_xc = e2p_exclude_lowBLyield_atMaxN( ...
    varNames_emu, cropList_emu_basei, Nlist, ...
    yieldBL_xv, thresh_kgm2)

% Get info
maxN = max(str2double(Nlist)) ;
Ncropsi_emu = length(cropList_emu_basei) ;

% Setup
exclude_xc = false(size(yieldBL_xv,1), Ncropsi_emu) ;

% Exclude a given crop-irrig if it had low yield at maximum N
for c = 1:Ncropsi_emu
    
    % Get index of thisCropNNN
    thisCrop = cropList_emu_basei{c} ;
    v = find(strcmp(varNames_emu, sprintf('%s%d', thisCrop, maxN))) ;
    if length(v) ~= 1
        error('Error finding v (%d found)', length(v)) ;
    end
    
    % Exclude
    exclude_xc(:,c) = yieldBL_xv(:,v) < thresh_kgm2 ;
end

end