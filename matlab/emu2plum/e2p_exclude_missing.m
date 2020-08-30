function missing_xc = e2p_exclude_missing( ...
    varNames_emu_basei, cropList_emu_basei, Nlist, ...
    data_xv)
% data_xv can also have a third dimension

% Get info
Ncropsi_emu = length(cropList_emu_basei) ;

% Setup
missing_xc = false(size(data_xv,1), Ncropsi_emu) ;

% Exclude a given crop-irrig if it had low yield at maximum N
for c = 1:Ncropsi_emu
    
    thisCrop = cropList_emu_basei{c} ;
    
    % Get indices of thisCrop
    thisCrop_i = find(strcmp(varNames_emu_basei,thisCrop)) ;
    if length(thisCrop_i)~=length(Nlist)
        error('Error finding isThisCrop (%d found)', length(thisCrop_i))
    end
    
    % Find where it's missing
    missing_xc(:,c) = isnan(mean(mean(data_xv(:,thisCrop_i,:),2),3)) ;
    
end

end