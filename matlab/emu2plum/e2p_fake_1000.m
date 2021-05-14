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
        
        % Find stand with 200 kgN/ha in LPJ-GUESS
        ind_lpj_200 = find( ...
            strcmp(getbasenamei(data_fu_lpj.varNames), thisCropIrr) ...
            & getN_num(data_fu_lpj.varNames)==200) ;
        if length(ind_lpj_200) ~= 1
            error('Error finding thisCrop_200: %d found', length(ind_lpj_200))
        end
        
        % Find stand with 1000 kgN/ha in LPJ-GUESS
        ind_lpj_1000 = find( ...
            strcmp(getbasenamei(data_fu_lpj.varNames), thisCropIrr) ...
            & getN_num(data_fu_lpj.varNames)==1000) ;
        if length(ind_lpj_1000) ~= 1
            error('Error finding ind_lpj_1000: %d found', length(ind_lpj_1000))
        end
        
        % Find stand with 200 kgN/ha in output
        ind_out_200 = find( ...
            strcmp(getbasenamei(data_fu_out.varNames), thisCropIrr) ...
            & getN_num(data_fu_out.varNames)==200) ;
        if length(ind_out_200) ~= 1
            error('Error finding ind_out_200: %d found', length(ind_out_200))
        end
        
        % Find stand with 1000 kgN/ha in output
        ind_out_1000 = find( ...
            strcmp(getbasenamei(data_fu_out.varNames), thisCropIrr) ...
            & getN_num(data_fu_out.varNames)==1000) ;
        if length(ind_out_1000) ~= 1
            error('Error finding ind_out_1000: %d found', length(ind_out_1000))
        end
                
        % Calculate deltas
        delta_x1t = data_fu_lpj.garr_xvt(:,ind_lpj_1000,:) ./ data_fu_lpj.garr_xvt(:,ind_lpj_200,:) ;
        delta_x1t(data_fu_lpj.garr_xvt(:,ind_lpj_200,:)==0 & data_fu_lpj.garr_xvt(:,ind_lpj_1000,:)==0) = 1 ;
        
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
                data_fu_lpj.garr_xvt(:,ind_lpj_1000,:) <= small_thresh ;
            delta_x1t(is_small) = 1 ;
        end
        if any(isinf(delta_x1t))
            error('Infs remain in delta_x1t. Deal with these.')
        end
        
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
    which_file, unique(getbasenamei(data_fu_out.varNames)), ...
    'Future', @getbasenamei)

end