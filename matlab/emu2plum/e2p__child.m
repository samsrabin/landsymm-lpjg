warning('on','all')

for g = 1:length(gcm_list)
    gcm = gcm_list{g} ;
    
    switch gcm
        case 'UKESM1-0-LL'
            gcm_prefix = 'ukesm' ;
        otherwise
            error('What prefix should be used for GCM %s?', gcm)
    end
    
    for s = 1:length(ssp_list)
        ssp = ssp_list{s} ;

        outDir = sprintf('%s_work/A%d_%s_%s_%s_%s', ...
            topDir_emu, adaptation, emuVer, gcm, ssp, thisVer) ;
        outDir = [outDir '_ph2bl'] ; %#ok<AGROW>
        if interp_infs
            outDir = [outDir '_intpinfs'] ;
        end
        if remove_outliers
            outDir = [outDir '_rmol' when_remove_outliers] ; %#ok<AGROW>
        end
        if fake1k
            outDir = [outDir '_fake1k'] ; %#ok<AGROW>
        end
        if ~excl_lowBL_agmerra
            outDir = [outDir '_ignLoPh2'] ; %#ok<AGROW>
        end
        if ~excl_lowBL_emu
            outDir = [outDir '_ignLoEm'] ; %#ok<AGROW>
        end
        if ~exist(outDir, 'dir')
            mkdir(outDir) ;
        end
        
        % Copy log to file
        diary('off')
        diaryfile = sprintf('%s/matlab_log.txt', outDir) ;
        if exist(diaryfile, 'file')
            delete(diaryfile) ;
        end
        diary(diaryfile)
        diary('on')
            
        outDir_lpj = sprintf('%s/sim_LPJ-GUESS', outDir) ;
        outDir_excl_figs_inCrops = sprintf('%s/excl_figs_GGCMIcrops', outDir) ;
        outDir_excl_figs_outCrops = sprintf('%s/excl_figs_PLUMcrops', outDir) ;
        outDir_interp_figs = sprintf('%s/interp_figs', outDir) ;
        outDir_yield_figs = sprintf('%s/yield_figs', outDir) ;
        outDir_irrig_figs = sprintf('%s/irrig_figs', outDir) ;
        
        if ~exist(outDir_lpj, 'dir')
            mkdir(outDir_lpj) ;
        end
        
        % Import LPJ-GUESS yield and irrigation
        if strcmp(which_system, 'mymac')
            tmp = sprintf('/Volumes/Reacher/G2P/outputs_LPJG/remap12_2016/%s_actual_2015soc_default/outputs', ...
                gcm_prefix) ;
        elseif strcmp(which_system, 'keal')
            tmp = sprintf('/pd/data/lpj/sam/ggcmi2plum/lpj-guess_runs/remap12_2016/%s_actual_2015soc_default/outputs', ...
                gcm_prefix) ;
        else
            error('Error parsing which_system "%s" for topDir_lpj', which_system)
        end
        if ~exist(tmp, 'dir')
            error('%s not found', tmp)
        end
        outdirs = dir(sprintf('%s/outForPLUM*', tmp)) ;
        topDir_lpj = sprintf('%s/%s', ...
            outdirs(end).folder, outdirs(end).name) ;
        clear outdirs tmp
        disp('Importing LPJ-GUESS yield...')
        which_file = 'yield' ;
        data_fu_lpj_yield = e2p_import_fu_lpj(future_y1-1, future_ts, future_yN_lpj, topDir_lpj, ...
            which_file, get_unneeded) ;
        e2p_check_correct_zeros(data_fu_lpj_yield.garr_xvt, ...
            which_file, data_fu_lpj_yield.varNames, ...
            'Future', @getbasenamei)
        [varNames_lpj, cropList_lpj, ...
            varNames_lpj_basei, cropList_lpj_basei, ...
            Nlist_lpj_char, ~] = ...
            e2p_get_names({}, data_fu_lpj_yield.varNames, ...
            get_unneeded) ;
        cropList_out_basei = cropIrrList_out ;
        gridlist.list_to_map = data_fu_lpj_yield.list2map ;
        gridlist.lonlats = data_fu_lpj_yield.lonlats ;
        gridlist.mask_YX = false(360, 720) ;
        gridlist.mask_YX(gridlist.list_to_map) = true ;
        Ncells = length(gridlist.list_to_map) ;
        gridlist_target = {gridlist.lonlats gridlist.list_to_map} ;
        disp('Importing LPJ-GUESS irrigation...')
        which_file = 'gsirrigation' ;
        data_fu_lpj_irrig = e2p_import_fu_lpj(future_y1-1, future_ts, future_yN_lpj, topDir_lpj, ...
            which_file, get_unneeded, gridlist_target) ;
        e2p_check_correct_zeros(data_fu_lpj_irrig.garr_xvt, ...
            which_file, data_fu_lpj_irrig.varNames, ...
            'Future', @getbasenamei)
        [~, cropList_lpj2, ...
            varNames_lpj_basei2, cropList_lpj_basei2, ...
            Nlist_lpj_char2, ~] = ...
            e2p_get_names({}, data_fu_lpj_irrig.varNames, ...
            get_unneeded) ;
        if ~isequal(cropList_lpj,cropList_lpj2)
            error('Mismatch between cropList_lpj for yield vs. gsirrigation')
        elseif ~isequal(varNames_lpj_basei,varNames_lpj_basei2)
            error('Mismatch between varNames_lpj_basei for yield vs. gsirrigation')
        elseif ~isequal(cropList_lpj_basei,cropList_lpj_basei2)
            error('Mismatch between cropList_lpj_basei for yield vs. gsirrigation')
        elseif ~isequal(Nlist_lpj_char,Nlist_lpj_char2)
            error('Mismatch between Nlist_lpj_char2 for yield vs. gsirrigation')
        end
        clear cropList_lpj2 varNames_lpj_basei2 cropList_lpj_basei2 Nlist_lpj_char2
        
        % Get max wheat and irrigation, if needed
        if any(strcmp(cropList_out, 'CerealsC3')) && ~any(strcmp(cropList_lpj, 'CerealsC3'))
            disp('Getting max wheats...')
            [data_fu_lpj_yield.garr_xvt, data_fu_lpj_yield.varNames, is_ww_max_fu_gWt, winter_wheats] = ...
                e2p_get_max_wheat(data_fu_lpj_yield.garr_xvt, data_fu_lpj_yield.varNames) ;
            out_file = sprintf('%s/is_ww_max.mat', topDir_lpj) ;
            save(out_file, 'is_ww_max_fu_gWt', 'winter_wheats') ;
            clear is_ww_max_fu_gWt winter_wheats
            data_fu_lpj_irrig = e2p_apply_max_wheat(data_fu_lpj_irrig, topDir_lpj) ;
            [varNames_lpj, cropList_lpj, ...
                varNames_lpj_basei, cropList_lpj_basei] = ...
                e2p_get_names(data_fu_lpj_yield.varNames, data_fu_lpj_irrig.varNames, ...
                get_unneeded) ;
        end
        
        disp('Done importing LPJ-GUESS yield and irrigation.')
        
        % Loop through crop emulators
        for ggcm_counter = 1:length(ggcm_list)
            ggcm = ggcm_list{ggcm_counter} ;
            topDir_phase2 = sprintf('%s/%s/phase2', topDir_agmipout, ggcm) ;
            outDir_ggcm = sprintf('%s/emu_%s', outDir, ggcm) ;
            if ~exist(outDir_ggcm, 'dir')
                mkdir(outDir_ggcm) ;
            end
            
            for w = 1:length(whichfile_list)
                which_file = whichfile_list{w} ;
                fprintf('%s %s %s %s\n', gcm, ssp, ggcm, which_file)
                
                % Get LPJ-GUESS simulation data for this file
                if strcmp(which_file, 'yield')
                    data_fu_lpj = data_fu_lpj_yield ;
                elseif strcmp(which_file, 'gsirrigation')
                    data_fu_lpj = data_fu_lpj_irrig ;
                else
                    error('which_file (%s) not recognized', which_file)
                end

                % If allowed to import existing file, then do so
                out_file = sprintf('%s/future_%s.mat', outDir_ggcm, which_file) ;
                did_load_existing = load_existing_file && exist(out_file, 'file') ;
                if did_load_existing
                    error('load_existing_file only tested with use_lpjg_baseline')
                    
                    fprintf('    Importing future_%s.mat...\n', which_file) ;

                    load(out_file, 'data_fu_lpj', 'data_fu_emu', 'data_fu_out')

                    % Translate crop names
                    [varNames_emu, cropList_emu, ...
                            varNames_emu_basei, cropList_emu_basei, ...
                            Nlist_emu_char, ~] = ...
                            e2p_get_names([], data_fu_emu.varNames, ...
                            get_unneeded) ;
                    [cropList_out_asEmu, used_emuCrops_out] = e2p_translate_crops_2emu( ...
                        cropList_out, cropList_emu, 'PLUM') ;
                    [cropList_out_asLpj, used_lpjCrops_out] = e2p_translate_crops_2lpj( ...
                        cropList_out, cropList_lpj, 'PLUM') ;
                else

                    %% Import emulator outputs

                    if contains(which_file, {'yield', 'gsirrigation'})

                        disp('    Importing emulator outputs...')

                        try
                            [data_bl_emu, data_fu_emu] = e2p_import_emu( ...
                                topDir_emu, gcm, ggcm, ssp, which_file, ...
                                cropList_in, gridlist, ts1_list, tsN_list, Nlist_emu, ...
                                baseline_yN, future_yN_emu, irrList_in, irrList_out, ...
                                emuVer, adaptation) ;
                        catch ME
                            if strcmp(ME.identifier, 'e2p:e2p_import_emu:noFilesFound')
                                warning('No files found for %s %s %s; skipping', ...
                                    gcm, ssp, ggcm)
                                continue
                            else
                                rethrow(ME)
                            end
                        end

                        e2p_check_correct_zeros(data_bl_emu.garr_xv, ...
                            which_file, data_bl_emu.varNames, ...
                            'Baseline', @getbasenamei, ...
                            true)
                        e2p_check_correct_zeros(data_fu_emu.garr_xvt, ...
                            which_file, data_fu_emu.varNames, ...
                            'Future', @getbasenamei, ...
                            true)

                        [varNames_emu, cropList_emu, ...
                            varNames_emu_basei, cropList_emu_basei, ...
                            Nlist_emu_char, ~] = ...
                            e2p_get_names(data_bl_emu.varNames, data_fu_emu.varNames, ...
                            get_unneeded) ;

                        disp('    Done.')

                    elseif ~strcmp(which_file, 'anpp')
                        error('which_file (%s) not recognized', which_file)
                    end
                    
                    false_xc = false(length(data_bl_emu.list2map),length(cropList_emu_basei)) ;
                    fprintf('    Importing AgMERRA %s...\n', which_file)
                    if ~isequal(data_fu_lpj.list2map, data_bl_emu.list2map)
                        error('gridlist mismatch')
                    end
                    data_bl_agm = e2p_get_agmerra_baseline(...
                        varNames_emu, topDir_phase2, topDir_emu, ggcm, ...
                        data_bl_emu.list2map, data_fu_lpj.lonlats, ...
                        adaptation, which_file) ;
                    if ~any(any(~isnan(data_bl_agm.garr_xv)))
                        error('data_bl_agm.garr_xv is all NaN')
                    end
                    % Translate to output crops
                    [varNames_agm, cropList_agm, ...
                        varNames_agm_basei, cropList_agm_basei, ...
                        Nlist_agm_char, ~] = ...
                        e2p_get_names(data_bl_agm.varNames, [], ...
                        get_unneeded) ;
                    [cropList_agm_asEmu, used_emuCrops_agm] = e2p_translate_crops_2emu( ...
                        cropList_agm, cropList_emu, 'Phase 2 sims') ;
                    % Sanity check
                    if ~isequal(size(shiftdim(cropList_emu)), size(shiftdim(used_emuCrops_agm)))
                        error('Size mismatch between cropList_emu and used_emuCrops_agm')
                    end
                    % Find missing crops
                    missing_agmerra_xc = e2p_exclude_missing( ...
                        varNames_emu_basei, cropList_emu_basei, Nlist_emu_char, ...
                        data_bl_agm.garr_xv) ;

                    
                    %% Get and apply exclusions, if doing so
                    
                    % Where do we exclude based on missing
                    % emulation?
                    [missing_emu_bl_xc, all0_emu_bl_c] = e2p_exclude_missing( ...
                        varNames_emu_basei, cropList_emu_basei, Nlist_emu_char, ...
                        data_bl_emu.garr_xv) ;
                    [missing_emu_fu_xc, all0_emu_fu_c] = e2p_exclude_missing( ...
                        varNames_emu_basei, cropList_emu_basei, Nlist_emu_char, ...
                        data_fu_emu.garr_xvt) ;
                    missing_emu_xc = missing_emu_bl_xc | missing_emu_fu_xc ;
                    
                    if strcmp(which_file, 'yield')
                        
                        isexcl_lowBL_emu_xc = false_xc ;
                        
                        % Where do we exclude based on low AgMERRA yield at max N
                        % (or existing exclusions)?
                        disp('    Excluding based on low AgMERRA yield at max N...')
                        isexcl_lowBL_agmerra_xc = e2p_exclude_lowBLyield_atMaxN( ...
                            varNames_emu, cropList_emu_basei, Nlist_emu_char, ...
                            data_bl_agm.garr_xv, low_yield_threshold_kgm2) ;
                        if ~any(any(~isexcl_lowBL_agmerra_xc))
                            error('All cells excluded because of low AgMERRA yield')
                        end
                        
                        % Where do we exclude based on NaN or low baseline-year emulated yield at
                        % max N (or existing exclusions)?
                        if excl_lowBL_emu
                            disp('    Excluding based on NaN/low baseline-year emulated yield at max N...')
                            isexcl_lowBL_emu_xc = e2p_exclude_lowBLyield_atMaxN( ...
                                varNames_emu, cropList_emu_basei, Nlist_emu_char, ...
                                data_bl_emu.garr_xv, low_yield_threshold_kgm2) ;
                        end
                        
                        if ~any(any(~isexcl_lowBL_emu_xc))
                            error('All cells excluded because of low AgMERRA yield and/or low baseline-year emulated yield')
                        end
                        
                        exclude_lowBLyield_xc = isexcl_lowBL_agmerra_xc | isexcl_lowBL_emu_xc ;
                        missing_yield_xc = missing_emu_xc | missing_agmerra_xc ;
                        
                        % Set up for figures
                        excl_vecs{1} = missing_emu_xc ;
                        excl_vecs{2} = missing_agmerra_xc ;
                        excl_vecs{3} = isexcl_lowBL_emu_xc ;
                        excl_vecs{4} = isexcl_lowBL_agmerra_xc ;
                        
                    elseif strcmp(which_file, 'gsirrigation')
                        disp('    Applying yield-based exclusions...')
                        exclude_xc_file = sprintf('%s/missing_or_excluded.mat', outDir_ggcm) ;
                        load(exclude_xc_file) ;
                        
                        % Set up for figures
                        excl_vecs{1} = missing_emu_xc ;
                        excl_vecs{2} = missing_yield_xc ;
                        excl_vecs{3} = exclude_lowBLyield_xc ;
                        excl_vecs{4} = all0_emu_bl_c | all0_emu_fu_c ;
                    else
                        error('which_file (%s) not recognized', which_file)
                    end
                    
                    % Save exclusion figures, if doing so
                    if save_excl_figs
                        disp('    Saving exclusion figures...')
                        e2p_save_excl_figs_outCrops( ...
                            ggcm, which_file, gridlist, excl_vecs, ...
                            cropList_out, cropList_out_basei, ...
                            cropList_emu, cropList_emu_basei, ...
                            figure_visibility, figure_extension, ...
                            outDir_excl_figs_outCrops, ...
                            overwrite_existing_figs, renderer)
                        e2p_save_excl_figs_inCrops(...
                            ggcm, which_file, gridlist, excl_vecs, ...
                            cropList_emu_basei, figure_visibility, ...
                            figure_extension, outDir_excl_figs_inCrops, ...
                            overwrite_existing_figs, renderer)
                    end
                    
                    % Apply exclusions
                    exclude_xc = missing_yield_xc | missing_emu_xc ...
                        | exclude_lowBLyield_xc  ;
                    for c = 1:length(cropList_emu_basei)
                        thisCrop = cropList_emu_basei{c} ;
                        thisCrop_i = find(strcmp(varNames_emu_basei,thisCrop)) ;
                        if length(thisCrop_i)~=length(Nlist_emu_char)
                            error('Error finding isThisCrop (%d found)', length(thisCrop_i))
                        end
                        
                        % If processing irrigation files, we don't care
                        % about rainfed crops. We know they're all zero
                        % anyway because of checks in
                        % e2p_check_correct_zeros() above.
                        if strcmp(which_file,'gsirrigation') && ...
                                ~strcmp(thisCrop(end), 'i')
                            continue
                        end
                        
                        % If irrigation, check for any missing cells
                        % that weren't missing in yield even after
                        % exclusions
                        if strcmp(which_file,'gsirrigation')
                            isbad = isnan(data_bl_emu.garr_xv(:,thisCrop_i)) & ~exclude_xc(:,c) ;
                            if any(any(isbad))
                                warning('%s: %d cells that were included in yield are NaN in irrig baseline emulation', ...
                                    thisCrop, length(find(isbad)))
                            end
                            isbad = isnan(data_fu_emu.garr_xvt(:,thisCrop_i,:)) & repmat(~exclude_xc(:,c),[1 length(thisCrop_i) Ntpers]) ;
                            if any(any(any(isbad)))
                                warning('%s: %d cells that were included in yield are NaN in irrig future emulation', ...
                                    thisCrop, length(find(isbad)))
                            end
                        end
                        
                        % Apply exclusions to emulated baseline
                        tmp = data_bl_emu.garr_xv(:,thisCrop_i) ;
                        tmp(exclude_xc(:,c),:) = NaN ;
                        data_bl_emu.garr_xv(:,thisCrop_i) = tmp ;
                        clear tmp
                        
                        % Apply exclusions to emulated future
                        tmp = data_fu_emu.garr_xvt(:,thisCrop_i,:) ;
                        tmp(exclude_xc(:,c),:,:) = NaN ;
                        data_fu_emu.garr_xvt(:,thisCrop_i,:) = tmp ;
                        clear tmp
                    end
                    
                    disp('    Done.')
                    
                    % Sanity checks
                    if ~any(any(~isnan(data_bl_emu.garr_xv)))
                        error('data_bl_emu.garr_xv is all NaN!')
                    end
                    if ~any(any(any(~isnan(data_fu_emu.garr_xvt))))
                        error('data_fu_emu.garr_xvt is all NaN!')
                    end


                    %% Get max wheats

                    disp('    Getting max wheats...')

                    refresh_vars = true ;
                    if strcmp(which_file, 'yield')
                        [data_bl_emu.garr_xv, data_bl_emu.varNames, is_ww_max_bl_gW, winter_wheats] = ...
                            e2p_get_max_wheat(data_bl_emu.garr_xv, data_bl_emu.varNames) ;
                        [data_fu_emu.garr_xvt, data_fu_emu.varNames, is_ww_max_fu_gWt, winter_wheats_test] = ...
                            e2p_get_max_wheat(data_fu_emu.garr_xvt, data_fu_emu.varNames) ;
                        if ~isequal(winter_wheats, winter_wheats_test)
                            error('Mismatch between winter wheat lists from data_bl_emu vs. data_fu_emu')
                        end
                        [data_bl_agm.garr_xv, data_bl_agm.varNames, ...
                            is_ww_max_bl_gW, winter_wheats_test, ...
                            data_bl_agm.actually_emu_char] = ...
                            e2p_get_max_wheat(data_bl_agm.garr_xv, data_bl_agm.varNames, ...
                            data_bl_agm.actually_emu) ;
                        if ~isequal(winter_wheats, winter_wheats_test)
                            error('Mismatch between winter wheat lists from data_bl_emu vs. data_bl_agm')
                        end
                    elseif strcmp(which_file, 'gsirrigation')
                        data_bl_emu = e2p_apply_max_wheat(data_bl_emu, outDir_ggcm) ;
                        data_fu_emu = e2p_apply_max_wheat(data_fu_emu, outDir_ggcm) ;
                        data_bl_agm = e2p_apply_max_wheat(data_bl_agm, outDir_ggcm) ;
                    elseif ~strcmp(which_file, 'anpp')
                        error('which_file (%s) not recognized', which_file)
                    else
                        refresh_vars = false ;
                    end

                    e2p_check_correct_zeros(data_bl_emu.garr_xv, ...
                        which_file, data_bl_emu.varNames, ...
                        'Baseline', @getbasenamei)
                    e2p_check_correct_zeros(data_fu_emu.garr_xvt, ...
                        which_file, data_fu_emu.varNames, ...
                        'Future', @getbasenamei)
                    e2p_check_correct_zeros(data_bl_agm.garr_xv, ...
                        which_file, data_bl_agm.varNames, ...
                        'Baseline', @getbasenamei)
                    
                    % Refresh variable and crop lists
                    if refresh_vars
                        [varNames_emu, cropList_emu, ...
                            varNames_emu_basei, cropList_emu_basei, ...
                            Nlist_emu_char, ~] = ...
                            e2p_get_names(data_bl_emu.varNames, data_fu_emu.varNames, ...
                            get_unneeded) ;
                        [varNames_agm, cropList_agm, ...
                            varNames_agm_basei, cropList_agm_basei, ...
                            Nlist_agm_char, ~] = ...
                            e2p_get_names(data_bl_agm.varNames, [], ...
                            get_unneeded) ;
                    end
                    
                    % Rearrange and rename AgMERRA baseline to match
                    % all possible variable names, expanding with
                    % NaN as needed where N values are missing.
                    data_bl_out.list2map = data_bl_agm.list2map ;
                    data_bl_out.lonlats = data_bl_agm.lonlats ;
                    [I_agm, I_out] = e2p_translate_crops_agm2out(...
                        varNames_agm, varNames_out_allN) ;
                    data_bl_out.varNames = varNames_out_allN ;
                    data_bl_out.actually_emu_char = cell(1,length(data_bl_out.varNames)) ;
                    data_bl_out.actually_emu_char(I_out) = ...
                        data_bl_agm.actually_emu_char(I_agm) ;
                    tmp = find(cellfun(@isempty, data_bl_out.actually_emu_char)) ;
                    data_bl_out.actually_emu_char(tmp) = ...
                        strcat(data_bl_out.actually_emu_char(tmp-1), 'xLPJG200-1000') ;
                    clear tmp
                    data_bl_out.garr_xv = nan(Ncells, ...
                        length(data_bl_out.varNames)) ;
                    data_bl_out.garr_xv(:,I_out) = data_bl_agm.garr_xv(:,I_agm) ;
                    clear tmp_garr_xv
                    
                    % Do the same thing for emulated baseline.
                    tmp_garr_xv = data_bl_emu.garr_xv ;
                    data_bl_emu = rmfield(data_bl_emu, 'garr_xv') ;
                    data_bl_emu.list2map = data_bl_emu.list2map ;
                    data_bl_emu.lonlats = gridlist.lonlats ;
                    [I_emu, I_out] = e2p_translate_crops_agm2out(...
                        varNames_emu, varNames_out_allN) ;
                    data_bl_emu.varNames = varNames_out_allN ;
                    data_bl_emu.garr_xv = nan(Ncells, length(data_bl_emu.varNames)) ;
                    data_bl_emu.garr_xv(:,I_out) = tmp_garr_xv(:,I_emu,:) ;
                    clear tmp_garr_xv
                    
                    % Do the same thing for emulated future.
                    tmp_garr_xvt = data_fu_emu.garr_xvt ;
                    data_fu_emu = rmfield(data_fu_emu, 'garr_xvt') ;
                    data_fu_emu.list2map = data_fu_emu.list2map ;
                    data_fu_emu.lonlats = gridlist.lonlats ;
                    [I_emu, I_out] = e2p_translate_crops_agm2out(...
                        varNames_emu, varNames_out_allN) ;
                    data_fu_emu.varNames = varNames_out_allN ;
                    data_fu_emu.garr_xvt = nan(Ncells, length(data_fu_emu.varNames), Ntpers) ;
                    data_fu_emu.garr_xvt(:,I_out,:) = tmp_garr_xvt(:,I_emu,:) ;
                    clear tmp_garr_xvt
                    
                    % Do the same thing for LPJ-GUESS simulated future.
                    tmp_garr_xvt = data_fu_lpj.garr_xvt ;
                    data_fu_lpj = rmfield(data_fu_lpj, 'garr_xvt') ;
                    data_fu_lpj.list2map = data_fu_lpj.list2map ;
                    data_fu_lpj.lonlats = gridlist.lonlats ;
                    [I_lpj, I_out] = e2p_translate_crops_lpj2out(...
                        varNames_lpj, varNames_out_allN) ;
                    data_fu_lpj.varNames = varNames_out_allN ;
                    data_fu_lpj.garr_xvt = nan(Ncells, length(data_fu_lpj.varNames), Ntpers) ;
                    data_fu_lpj.garr_xvt(:,I_out,:) = tmp_garr_xvt(:,I_lpj,:) ;
                    clear tmp_garr_xvt
                    
                    % Now that you've expanded to match PLUM (output) 
                    % variable names, these are the same.
                    Nlist_agm_char = Nlist_lpj_char ;
                    Nlist_emu_char = Nlist_lpj_char ;
                    
                    % Get indices of variables that should NOT be missing
                    ok_missing_v = get_ok_missing(varNames_out_allN, Nlist_emu) ;
                    
                    disp('    Done.')


                    % Harmonize LPJ-GUESS and emulator data

                    disp('    Harmonizing LPJ-GUESS and emulator data...')

                    % Align gridlists
                    [data_fu_lpj, data_bl_emu, data_fu_emu, list2map] = ...
                        e2p_align_gridlists(data_fu_lpj, data_bl_emu, data_fu_emu) ;

                    e2p_check_correct_zeros(data_fu_lpj.garr_xvt, ...
                        which_file, data_fu_lpj.varNames, ...
                        'Future', @getbasenamei)
                    tmp_garr_xv = data_bl_emu.garr_xv(:,~ok_missing_v) ;
                    varNames_NOTok_missing = data_fu_emu.varNames(~ok_missing_v) ;
                    e2p_check_correct_zeros(tmp_garr_xv, ...
                        which_file, varNames_NOTok_missing, ...
                        'Baseline', @getbasenamei)
                    tmp_garr_xvt = data_fu_emu.garr_xvt(:,~ok_missing_v,:) ;
                    e2p_check_correct_zeros(tmp_garr_xvt, ...
                        which_file, varNames_NOTok_missing, ...
                        'Future', @getbasenamei)
                    clear tmp_garr_xv tmp_garr_xvt

                    % Make sure that N lists match
                    if ~isequal(Nlist_agm_char,Nlist_emu_char)
                        error('Mismatch in N levels between Phase 2 baseline and emulator future')
                    end

                    % Translate crop names
                    [cropList_out_asEmu, used_emuCrops_out] = e2p_translate_crops_2emu( ...
                        cropList_out, cropList_emu, 'PLUM') ;
                    [cropList_agm_asEmu, used_emuCrops_agm] = e2p_translate_crops_2emu( ...
                        cropList_agm, cropList_emu, 'Phase2 sims') ;
                    [cropList_out_asLpj, used_lpjCrops_out] = e2p_translate_crops_2lpj( ...
                        cropList_out, cropList_lpj, 'PLUM') ;
                    % Sanity check
                    if ~isequal(size(shiftdim(cropList_emu)), size(shiftdim(used_emuCrops_agm)))
                        error('Size mismatch between cropList_emu and used_emuCrops_agm')
                    end

                    disp('    Done.')
                    

                    % Get deltas
                    if contains(which_file, {'yield','gsirrigation'})
                        disp('    Getting deltas...')

                        deltas_emu_xvt = e2p_get_deltas(...
                            data_bl_emu, data_fu_emu, interp_infs, cropList_emu, ...
                            which_file, ...
                            used_emuCrops_agm, list2map, ...
                            save_interp_figs, outDir_interp_figs, ggcm, figure_visibility, ...
                            when_remove_outliers, outDir_ggcm, renderer) ;
                        
                        % Check
                        tmp_garr_xvt = deltas_emu_xvt(:,~ok_missing_v,:) ;
                        e2p_check_correct_zeros(tmp_garr_xvt, ...
                            which_file, varNames_NOTok_missing, ...
                            'Future', @getbasenamei)
                        clear tmp_garr_xvt
                        
                        disp('    Done.')
                    elseif ~strcmp(which_file, 'anpp')
                        error('which_file (%s) not recognized', which_file)
                    end
                    
                    % Apply emulator deltas to chosen baseline.
                    disp('    Applying deltas...')
                    data_fu_out = e2p_apply_deltas( ...
                        data_bl_out, data_bl_emu, data_fu_emu, deltas_emu_xvt, ...
                        cropList_out, cropList_out, varNames_out_allN, ...
                        list2map, which_file, figure_visibility) ;
                    % Check
                    [ok_missing_v, varNames_NOTok_missing] = ...
                        get_ok_missing(data_fu_out.varNames, Nlist_emu) ;
                    tmp_garr_xvt = data_fu_out.garr_xvt(:,~ok_missing_v,:) ;
                    e2p_check_correct_zeros(tmp_garr_xvt, ...
                        which_file, varNames_NOTok_missing, ...
                        'Future', @getbasenamei)
                    
                    % Sort variable names in data_fu_lpj (was previously
                    % done in e2p_apply_deltas()
                    if ~isequal(data_fu_lpj.varNames, sort(data_fu_lpj.varNames))
                        [data_fu_lpj.varNames, I] = sort(data_fu_lpj.varNames) ;
                        data_fu_lpj.garr_xvt = data_fu_lpj.garr_xvt(:,I,:) ;
                        clear I
                        e2p_check_correct_zeros(data_fu_lpj.garr_xvt, ...
                            which_file, data_fu_lpj.varNames, ...
                            'Future', @getbasenamei)
                    end
                    
                    % Get fake N1000 values; trim unneeded N200
                    data_fu_out = e2p_fake_1000(data_fu_out, data_fu_lpj, ...
                        Nlist_lpj_char, Nlist_emu_char, which_file) ;
                    unneededN200 = find(getN_num(data_fu_out.varNames)==200) ;
                    data_fu_out.garr_xvt(:,unneededN200,:) = [] ;
                    data_fu_out.varNames(unneededN200) = [] ;
                    data_fu_out.actually_emu_char(unneededN200) = [] ;
                    e2p_check_correct_zeros(data_fu_out.garr_xvt, ...
                        which_file, data_fu_out.varNames, ...
                        'Future', @getbasenamei)
                    
                    % Remove outliers
                    if strcmp(when_remove_outliers, 'end')
                        disp('    Removing outliers...')

                        [data_fu_lpj, outlier_info_lpj] = e2p_remove_outliers(data_fu_lpj, which_file) ;
                        e2p_check_correct_zeros(data_fu_lpj.garr_xvt, ...
                            which_file, data_fu_lpj.varNames, ...
                            'Future', @getbasenamei)

                        [data_fu_out, outlier_info_out] = e2p_remove_outliers(data_fu_out, which_file) ;
                        e2p_check_correct_zeros(data_fu_out.garr_xvt, ...
                            which_file, data_fu_out.varNames, ...
                            'Future', @getbasenamei)

                        % Save info
                        if strcmp(ggcm, ggcm_list{1})
                            e2p_save_outlier_info(outlier_info_lpj, outDir_lpj, which_file, data_fu_out.y1s, data_fu_out.yNs)
                        end
                        e2p_save_outlier_info(outlier_info_out, outDir_ggcm, which_file, data_fu_out.y1s, data_fu_out.yNs)
                    end
                    
                    % Consistency check
                    if isfield(data_fu_out, 'actually_emu_char') ...
                    && ~isequal(size(shiftdim(data_fu_out.varNames)), size(shiftdim(data_fu_out.actually_emu_char)))
                        error('data_fu_out.actually_emu_char must be the same size as data_fu_out.varNames')
                    end

                end
                
                % Sort. Technically unnecessary, but maybe PLUM isn't
                % robust to non-sorted crops.
                data_fu_lpj = do_sort(data_fu_lpj) ;
                data_fu_out = do_sort(data_fu_out) ;

                disp('    Done.')


                %% Save outputs
                
                if ~did_load_existing
                    disp('    Saving MAT file(s)...')

                    % Save outputs (one file)
                    out_file = sprintf('%s/future_%s.mat', outDir_ggcm, which_file) ;
                    save(out_file, 'data_fu_lpj', 'data_fu_emu', 'data_fu_out')

                    % Save exclusion info
                    if strcmp(which_file, 'yield')
                        out_file = sprintf('%s/missing_or_excluded.mat', outDir_ggcm) ;
                        save(out_file, 'exclude_lowBLyield_xc', 'missing_yield_xc')
                    end
                    clear exclude_xc exclude_lowBLyield_xc missing_yield_xc
                    
                    % For yield, save is_ww_max_bl_gW*
                    if strcmp(which_file, 'yield')
                        out_file = sprintf('%s/is_ww_max.mat', outDir_ggcm) ;
                        save(out_file, 'is_ww_max_bl_gW', 'is_ww_max_fu_gWt', 'winter_wheats') ;
                    elseif ~contains(which_file, {'anpp', 'gsirrigation'})
                        error('which_file (%s) not recognized', which_file)
                    end
                end

                % Save outputs (for PLUM)
                if save_txt_files_emu || save_txt_files_lpjg
                    
                    % Import ANPP and runoff, if needed
                    save_anpp_runoff = save_txt_files_lpjg ...
                        & strcmp(ggcm, ggcm_list{1}) ...
                        & strcmp(which_file, 'yield') ;
                    if save_anpp_runoff
                        disp('Importing LPJ-GUESS ANPP and runoff...')
                        data_fu_lpj_anpp = e2p_import_fu_lpj(...
                            future_y1-1, future_ts, future_yN_lpj, ...
                            topDir_lpj, 'anpp', [], ...
                            gridlist_target) ;
                        data_fu_lpj_runoff = e2p_import_fu_lpj(...
                            future_y1-1, future_ts, future_yN_lpj, ...
                            topDir_lpj, 'tot_runoff', [], ...
                            gridlist_target) ;
                        out_header_cell_anpp = [{'Lon', 'Lat'} data_fu_lpj_anpp.varNames] ;
                        out_header_cell_runoff = [{'Lon', 'Lat'} data_fu_lpj_runoff.varNames] ;
                    end
                    
                    disp('    Saving txt files for PLUM...')
                    out_header_cell = [{'Lon', 'Lat'} data_fu_out.varNames] ;
                    for t = 1:Ntpers
                        y1 = data_fu_out.y1s(t) ;
                        yN = data_fu_out.yNs(t) ;

                        if save_txt_files_lpjg && strcmp(ggcm, ggcm_list{1})
                            % Trim unneeded N200 from LPJ-GUESS
                            unneededN200 = getN_num(data_fu_lpj.varNames)==200 ;
                            tmp_xv = data_fu_lpj.garr_xvt(:,~unneededN200,t) ;
                            fprintf('    %d/%d lpj', t, Ntpers)
                            e2p_save(outDir_lpj, y1, yN, out_header_cell, ...
                                data_fu_lpj.lonlats, tmp_xv, which_file, ...
                                false)
                            
                            % Save ANPP and runoff
                            if save_anpp_runoff
                                fprintf('    %d/%d lpj ANPP', t, Ntpers)
                                e2p_save(outDir_lpj, y1, yN, out_header_cell_anpp, ...
                                    data_fu_lpj.lonlats, data_fu_lpj_anpp.garr_xvt(:,:,t), ...
                                    'anpp', false)
                                fprintf('    %d/%d lpj runoff', t, Ntpers)
                                e2p_save(outDir_lpj, y1, yN, out_header_cell_runoff, ...
                                    data_fu_lpj.lonlats, data_fu_lpj_runoff.garr_xvt(:,:,t), ...
                                    'tot_runoff', false)
                            end
                        end
                        if save_txt_files_emu
                            fprintf('    %d/%d emu', t, Ntpers)
                            e2p_save(outDir_ggcm, y1, yN, out_header_cell, ...
                                data_fu_out.lonlats, data_fu_out.garr_xvt(:,:,t), which_file, ...
                                overwrite_existing_txt)
                        end

                    end
                    clear data_fu_lpj_anpp data_fu_lpj_runoff
                end

                %% Save yield diagnostic figures, if doing so
                if save_out_figs
                    disp('    Saving output figures...')
                    if strcmp(which_file, 'yield')
                        e2p_save_out_figs(data_fu_lpj, ...
                            data_fu_emu, data_fu_out, ...
                            ggcm, outDir_yield_figs, ...
                            which_file, cropList_out_asEmu, figure_visibility, figure_extension, ...
                            which_out_figs, overwrite_existing_figs, renderer)
                    elseif strcmp(which_file, 'gsirrigation')
                        e2p_save_out_figs(data_fu_lpj, ...
                            data_fu_emu, data_fu_out, ...
                            ggcm, outDir_irrig_figs, ...
                            which_file, cropList_out_asEmu, figure_visibility, figure_extension, ...
                            which_out_figs, overwrite_existing_figs, renderer)
                    else
                        warning('which_file (%s) not recognized for save_yield_figs; skipping.', which_file)
                    end
                    clear tmp
                end
                fprintf('Done with %s %s %s %s.\n', gcm, ssp, ggcm, which_file)
                
            end
            
            fprintf('Done with %s %s %s.\n', gcm, ssp, ggcm)

        end % Loop through crop models

        fprintf('Done with %s %s.\n', gcm, ssp)
        diary('off')
        
        clear data_fu_lpj_yield data_fu_lpj_irrig

    end % Loop through SSPs
    
    fprintf('Done with %s.\n', gcm)
end % Loop through climate models

disp(' ')
disp('All done!')


%% FUNCTIONS

function [ok_missing_v, varNames_NOTok_missing] = get_ok_missing(varNames, Nlist_emu)

tmp_Ns_out = cellfun(@str2double, getN_char(varNames)) ;
C = setdiff(tmp_Ns_out, Nlist_emu) ;
ok_missing_v = false(size(varNames)) ;
for n = 1:length(C)
    ok_missing_v(tmp_Ns_out==C(n)) = true ;
end

varNames_NOTok_missing = varNames(~ok_missing_v) ;

end

function S = do_sort(S)

% Sort N strings
Nlist = unique(getN_char(S.varNames)) ;
Nlist_num = str2double(Nlist) ;
[Nlist_num, I] = sort(Nlist_num) ; %#ok<ASGLU>
Nlist = Nlist(I) ;

% Sort crop names (with irrig indicator)
cropList = unique(getbasenamei(S.varNames)) ;
cropList = sort(cropList) ;

% Get new index order
Nn = length(Nlist) ;
Ncrops = length(cropList) ;
Nvars = Nn*Ncrops ;
new_order = nan(1,Nvars) ;
for c = 1:Ncrops
    thisCrop = cropList{c} ;
    for n = 1:Nn
        thisCropN = [thisCrop Nlist{n}] ;
        thisInd = find(strcmp(S.varNames, thisCropN)) ;
        if length(thisInd) ~= 1
            error('Error finding thisCropN: %d found', length(thisInd))
        end
        new_order((c-1)*Nn+n) = thisInd ;
    end
end

% Rearrange
if ~isequal(1:Nvars, new_order)
    S.garr_xvt = S.garr_xvt(:,new_order,:) ;
    S.varNames = S.varNames(new_order) ;
end

end





