function data_fu_out = e2p_fake_1000( ...
    data_fu_out, data_fu_lpj, Nlist_lpj, Nlist_emu, ...
    which_file, scale_200to1000)

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
        delta_200_lpjg_x1t = get_deltas_betweenNs( ...
            data_fu_lpj, ind_lpj_200, ind_lpj_1000, which_file, ...
            thisCropIrr, 'LPJ-GUESS', 200, 1000) ;
        
        % If needed, change 200-1000 deltas according to ratio of
        % 60-200 deltas
        if scale_200to1000
            ind_lpj_60 = find_stand(data_fu_lpj.varNames, thisCropIrr, 60) ;
            ind_out_60 = find_stand(data_fu_out.varNames, thisCropIrr, 60) ;
            delta_60_lpj_x1t = get_deltas_betweenNs( ...
                data_fu_lpj, ind_lpj_60, ind_lpj_200, which_file, ...
                thisCropIrr, 'LPJ-GUESS', 60, 200) ;
            [delta_60_out_x1t, use_lpj_x1t] = get_deltas_betweenNs( ...
                data_fu_out, ind_out_60, ind_out_200, which_file, ...
                thisCropIrr, 'emulator', 60, 200) ;
            delta_60_ratio_x1t = delta_60_out_x1t ./ delta_60_lpj_x1t ;
            
            % Handle cells with zero delta in LPJ-GUESS
            delta_60_ratio_x1t(delta_60_lpj_x1t==0) = 1 ;
            
            % Adjust deltas, assuming that there's no way increasing N from
            % 200 to 100 should *decrease* yield or irrigation requirement
            delta_x1t = max(1, delta_200_lpjg_x1t .* delta_60_ratio_x1t) ;
            
            % Handle cells with infinite emulated N60-N200 where N200 was
            % too high to be ignored
            delta_x1t(use_lpj_x1t) = ...
                delta_200_lpjg_x1t(use_lpj_x1t) ;
        else
            delta_x1t = delta_200_lpjg_x1t ;
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
    which_file, data_fu_out.varNames, ...
    'Future', @getbasenamei)
    
end


function [delta_x1t, use_lpj_x1t] = get_deltas_betweenNs(data_fu, ind_loN, ind_hiN, ...
    which_file, thisCropIrr, whichModel, loN, hiN)

% Calculate deltas
data_lo_x1t = data_fu.garr_xvt(:,ind_loN,:) ;
data_hi_x1t = data_fu.garr_xvt(:,ind_hiN,:) ;
delta_x1t = data_hi_x1t ./ data_lo_x1t ;

% MATLAB sometimes gives -Inf when it should be positive Inf????????
delta_x1t(isinf(delta_x1t) & delta_x1t<0 ...
    & data_hi_x1t>=0 & data_lo_x1t>=0) = Inf ;

delta_x1t(data_fu.garr_xvt(:,ind_loN,:)==0 & data_fu.garr_xvt(:,ind_hiN,:)==0) = 1 ;

% If any infinite, first try excluding very small values
if any(any(isinf(delta_x1t)))
    if strcmp(which_file, 'yield')
        small_thresh = 0.01 ; % kgC/m2 yield
        small_thresh_units = 'kgC/m2 yield' ;
    elseif strcmp(which_file, 'gsirrigation')
        small_thresh = 20 ; % mm irrigation
        small_thresh_units = 'mm irrigation' ;
    else
        error('which_file (%s) not recognized for e2p_fake_1000().', which_file)
    end
    is_small = isinf(delta_x1t) & ...
        data_fu.garr_xvt(:,ind_hiN,:) <= small_thresh ;
    delta_x1t(is_small) = 1 ;
end

% Make sure no values are negative
if any(any(delta_x1t < 0 ))
    error('Trouble getting %s %s deltas between N%d and N%d: Negative value(s)', ...
        whichModel, thisCropIrr, loN, hiN)
end

% If infinite values remain, either throw error or mark for processing in
% next step
use_lpj_x1t = false(size(delta_x1t)) ;
if any(any(isinf(delta_x1t)))
    hiN_x1t = data_fu.garr_xvt(:,ind_hiN,:) ;
    bad_hiN = hiN_x1t(isinf(delta_x1t)) ;
    Nbad = length(find(isinf(delta_x1t))) ;
    if strcmp(whichModel, 'emulator') && loN==60 && hiN==200
        warning(['Trouble getting %s %s deltas between N%d and N%d: %d cell-' ...
            'timesteps have zero at low N but at high N exceed "small value" ' ...
            'threshold of %g %s (median %g, mean %g, max %g). ' ...
            'These will use unmodified LPJ-GUESS deltas.'], ...
            whichModel, thisCropIrr, loN, hiN, Nbad, ...
            small_thresh, small_thresh_units, median(bad_hiN), ...
            mean(bad_hiN), max(bad_hiN))
        use_lpj_x1t(isinf(delta_x1t)) = true ;
    else
        error(['Trouble getting %s %s deltas between N%d and N%d: %d cell-' ...
            'timesteps have zero at low N but at high N exceed "small value" ' ...
            'threshold of %g %s (median %g, mean %g, max %g).'], ...
            whichModel, thisCropIrr, loN, hiN, Nbad, ...
            small_thresh, small_thresh_units, median(bad_hiN), ...
            mean(bad_hiN), max(bad_hiN))
    end
end
if any(any(isinf(delta_x1t) & ~use_lpj_x1t))
    error('Trouble getting %s %s deltas between N%d and N%d: How do you STILL have infinite values?', ...
        whichModel, thisCropIrr, loN, hiN)
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
