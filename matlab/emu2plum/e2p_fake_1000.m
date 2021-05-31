function data_fu_out = e2p_fake_1000( ...
    data_fu_out, data_fu_lpj, Nlist_lpj, Nlist_emu, ...
    which_file)

% Get crop lists
cropList_out = unique(getbasename(data_fu_out.varNames)) ;
cropList_lpj = unique(getbasename(data_fu_lpj.varNames)) ;

% Set up array of fake 1000s
Ncells = length(data_fu_out.list2map) ;
Ncrops = length(cropList_lpj) ;
irrList = {'', 'i'} ;
Nirr = length(irrList) ;
Nts = size(data_fu_out.garr_xvt,3) ;
Ncrops_fake = Ncrops*Nirr ;
% fake1000_xvt = nan(Ncells, Ncrops_fake, Nts) ;
% fake1000_varNames = cell(1, Ncrops_fake) ;
if isfield(data_fu_out, 'actually_emu_char')
    actually_emu_char_1000 = cell(1, Ncrops_fake) ;
end

% Get strings corresponding to "200" and "1000" for LPJ-GUESS and emulator
Nlist_num_lpj = str2double(Nlist_lpj) ;
N200_str_lpj = Nlist_lpj{Nlist_num_lpj==200} ;
N1000_str_lpj = Nlist_lpj{Nlist_num_lpj==1000} ;
Nlist_num_emu = str2double(Nlist_emu) ;
N200_str_emu = Nlist_emu{Nlist_num_emu==200} ;

% Process
for c = 1:Ncrops
    thisCrop = cropList_lpj{c} ;
    for ii = 1:Nirr
        thisIrr = irrList{ii} ;
        thisCropIrr = [thisCrop thisIrr] ;
        thisCrop_200_lpj = [thisCropIrr N200_str_lpj] ;
        thisCrop_200_emu = [thisCropIrr N200_str_emu] ;
        thisCrop_1000 = [thisCropIrr N1000_str_lpj] ;
        
        % Find stands
        ind_lpj_200 = find_stand(data_fu_lpj.varNames, thisCropIrr, 200) ;
        ind_lpj_1000 = find_stand(data_fu_lpj.varNames, thisCropIrr, 1000) ;
        ind_out_200 = find_stand(data_fu_out.varNames, thisCropIrr, 200) ;
        ind_out_1000 = find_stand(data_fu_out.varNames, thisCropIrr, 1000) ;
        
        % Calculate deltas
        delta_x1t = calculate_deltas( ...
            data_fu_lpj, ind_lpj_200, ind_lpj_1000, which_file) ;
        
        % Multiply onto existing N200 to get fake N1000
        thisInd = (c-1)*Nirr + ii ;
        data_fu_out.garr_xvt(:,ind_out_1000,:) = ...
            data_fu_out.garr_xvt(:,ind_out_200,:) .* delta_x1t ;
        
        if isfield(data_fu_out, 'actually_emu_char')
            data_fu_out.actually_emu_char{ind_out_1000} = ...
                data_fu_out.actually_emu_char{ind_out_200} ;
        end
    end
end

e2p_check_correct_zeros(data_fu_out.garr_xvt, ...
    which_file, data_fu_out.varNames, ...
    'Future', @getbasenamei)
    
end


function delta_x1t = calculate_deltas(data_fu, ind_loN, ind_hiN, which_file)

delta_x1t = data_fu.garr_xvt(:,ind_hiN,:) ./ data_fu.garr_xvt(:,ind_loN,:) ;
delta_x1t(data_fu.garr_xvt(:,ind_loN,:)==0 & data_fu.garr_xvt(:,ind_hiN,:)==0) = 1 ;

% If any infinite, first try excluding very small values
if any(isinf(delta_x1t))
    if strcmp(which_file, 'yield')
        small_thresh = 0.001 ; % 0.001 kgC/m2 yield
    elseif strcmp(which_file, 'gsirrigation')
        small_thresh = 20 ; % 20 mm irrigation
    else
        error('which_file (%s) not recognized for e2p_fake_1000().', which_file)
    end
    is_small = isinf(delta_x1t) & ...
        data_fu.garr_xvt(:,ind_hiN,:) <= small_thresh ;
    delta_x1t(is_small) = 1 ;
end
if any(isinf(delta_x1t))
    error('Infs remain in delta_x1t. Deal with these.')
end

end


function ind = find_stand(varNames, thisCropIrr, thisN)

ind = find( ...
    strcmp(getbasenamei(varNames), thisCropIrr) ...
    & getN_num(varNames)==thisN) ;
if length(ind) ~= 1
    error('Error finding %s: %d found', ...
        thisCropIrr, length(ind))
end

end
