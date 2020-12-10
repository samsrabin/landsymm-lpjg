function data_fu_out = e2p_fake_1000( ...
    data_fu_out, data_bl_lpj0, cropList_lpj, Nlist_lpj0, Nlist_emu, ...
    which_file)

% Set up array of fake 1000s
Ncells = length(data_fu_out.list2map) ;
Ncrops = length(cropList_lpj) ;
irrList = {'', 'i'} ;
Nirr = length(irrList) ;
Nts = size(data_fu_out.garr_xvt,3) ;
Ncrops_fake = Ncrops*Nirr ;
fake1000_xvt = nan(Ncells, Ncrops_fake, Nts) ;
fake1000_varNames = cell(1, Ncrops_fake) ;

% Get strings corresponding to "200" and "1000" for LPJ-GUESS and emulator
Nlist_num_lpj = str2double(Nlist_lpj0) ;
N200_str_lpj = Nlist_lpj0{Nlist_num_lpj==200} ;
N1000_str_lpj = Nlist_lpj0{Nlist_num_lpj==1000} ;
Nlist_num_emu = str2double(Nlist_emu) ;
N200_str_emu = Nlist_emu{Nlist_num_emu==200} ;

% Process
for c = 1:Ncrops
    thisCrop = cropList_lpj{c} ;
    for ii = 1:Nirr
        thisIrr = irrList{ii} ;
        thisCrop_200_lpj = [thisCrop thisIrr N200_str_lpj] ;
        thisCrop_200_emu = [thisCrop thisIrr N200_str_emu] ;
        thisCrop_1000 = [thisCrop thisIrr N1000_str_lpj] ;
        
        % Find stand with 200 kgN/ha in LPJ-GUESS
        ind_lpj_200 = find(strcmp(data_bl_lpj0.varNames, thisCrop_200_lpj)) ;
        if length(ind_lpj_200) ~= 1
            error('Error finding thisCrop_200: %d found', length(ind_lpj_200))
        end
        
        % Find stand with 1000 kgN/ha in LPJ-GUESS
        ind_lpj_1000 = find(strcmp(data_bl_lpj0.varNames, thisCrop_1000)) ;
        if length(ind_lpj_1000) ~= 1
            error('Error finding ind_lpj_1000: %d found', length(ind_lpj_1000))
        end
        
        % Find stand with 200 kgN/ha in output
        ind_out_200 = find(strcmp(data_fu_out.varNames, thisCrop_200_emu)) ;
        if length(ind_out_200) ~= 1
            error('Error finding ind_out_200: %d found', length(ind_out_200))
        end
                
        % Calculate deltas
        delta_x = data_bl_lpj0.garr_xv(:,ind_lpj_1000) ./ data_bl_lpj0.garr_xv(:,ind_lpj_200) ;
        delta_x(data_bl_lpj0.garr_xv(:,ind_lpj_200)==0 & data_bl_lpj0.garr_xv(:,ind_lpj_1000)==0) = 1 ;
        
        % If any infinite, first try excluding very small values
        if any(isinf(delta_x))
            if strcmp(which_file, 'yield')
                small_thresh = 0.001 ; % 0.001 kgC/m2 yield
            elseif strcmp(which_file, 'gsirrigation')
                small_thresh = 20 ; % 20 mm irrigation
            else
                error('which_file (%s) not recognized for e2p_fake_1000().', which_file)
            end
            is_small = isinf(delta_x) & ...
                data_bl_lpj0.garr_xv(:,ind_lpj_1000) <= small_thresh ;
            delta_x(is_small) = 1 ;
        end
        if any(isinf(delta_x))
            error('Infs remain in delta_x. Deal with these.')
        end
        
        % Multiply onto existing N200 to get fake N1000
        delta_x1t = repmat(delta_x, [1 1 Nts]) ;
        thisInd = (c-1)*Nirr + ii ;
        fake1000_xvt(:,thisInd,:) = ...
            data_fu_out.garr_xvt(:,ind_out_200,:) .* delta_x1t ;
        fake1000_varNames{thisInd} = thisCrop_1000 ;
    end
end

data_fu_out.varNames = [data_fu_out.varNames fake1000_varNames] ;
data_fu_out.garr_xvt = cat(2, data_fu_out.garr_xvt, fake1000_xvt) ;

end