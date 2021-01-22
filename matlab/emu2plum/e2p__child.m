warning('on','all')

for g = 1:length(gcm_list)
    gcm = gcm_list{g} ;
    
    for s = 1:length(ssp_list)
        ssp = ssp_list{s} ;

        outDir = sprintf('%s_work/A%d_%s_%s_%s_%s', ...
            topDir_emu, adaptation, emuVer, gcm, ssp, thisVer) ;
        if use_ph2_baseline
            outDir = [outDir '_ph2bl'] ; %#ok<AGROW>
        elseif ~use_lpjg_baseline
            error('What suffix do I use for this baseline?')
        end
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
            
        outDir_lpj = sprintf('%s/sim_LPJ-GUESS', outDir) ;
        outDir_excl_figs_inCrops = sprintf('%s/excl_figs_GGCMIcrops', outDir) ;
        outDir_excl_figs_outCrops = sprintf('%s/excl_figs_PLUMcrops', outDir) ;
        outDir_interp_figs = sprintf('%s/interp_figs', outDir) ;
        outDir_yield_figs = sprintf('%s/yield_figs', outDir) ;
        outDir_irrig_figs = sprintf('%s/irrig_figs', outDir) ;

        if ~exist(outDir, 'dir')
            mkdir(outDir) ;
        end
        if ~exist(outDir_lpj, 'dir')
            mkdir(outDir_lpj) ;
        end

        for ggcm_counter = 1:length(ggcm_list)

            ggcm = ggcm_list{ggcm_counter} ;
            topDir_phase2 = sprintf('%s/AgMIP.output/%s/phase2', topdir_sh, ggcm) ;
            outDir_ggcm = sprintf('%s/emu_%s', outDir, ggcm) ;
            if ~exist(outDir_ggcm, 'dir')
                mkdir(outDir_ggcm) ;
            end
            
            for w = 1:length(whichfile_list)
                which_file = whichfile_list{w} ;
                fprintf('%s %s %s %s\n', gcm, ssp, ggcm, which_file)
                
                % Get LPJ-GUESS simulation data for this file
                if strcmp(which_file, 'yield')
                    data_bl_lpj = data_bl_lpj_yield ;
                    data_fu_lpj = data_fu_lpj_yield ;
                    data_bl_lpj0 = data_bl_lpj0_yield ;
                    data_fu_lpj0 = data_fu_lpj0_yield ;
                elseif strcmp(which_file, 'gsirrigation')
                    data_bl_lpj = data_bl_lpj_irrig ;
                    data_fu_lpj = data_fu_lpj_irrig ;
                    data_bl_lpj0 = data_bl_lpj0_irrig ;
                    data_fu_lpj0 = data_fu_lpj0_irrig ;
                else
                    error('which_file (%s) not recognized', which_file)
                end

                % If allowed to import existing file, then do so
                out_file = sprintf('%s/future_%s.mat', outDir_ggcm, which_file) ;
                did_load_existing = load_existing_file && exist(out_file, 'file') ;
                if did_load_existing
                    if ~use_lpjg_baseline
                        error('load_existing_file only tested with use_lpjg_baseline')
                    end
                    
                    fprintf('    Importing future_%s.mat...\n', which_file) ;

                    load(out_file, 'data_fu_lpj', 'data_fu_emu', 'data_fu_out')

                    % Translate crop names
                    [varNames_emu, cropList_emu, ...
                            varNames_emu_basei, cropList_emu_basei, ...
                            Nlist_emu, ~] = ...
                            e2p_get_names([], data_fu_emu.varNames, ...
                            getN, get_unneeded) ;
                    [cropList_lpj_asEmu, used_emuCrops_lpj] = e2p_translate_crops_lpj2emu( ...
                        cropList_lpj, cropList_emu) ;
                else

                    %% Import emulator outputs

                    if contains(which_file, {'yield', 'gsirrigation'})

                        disp('    Importing emulator outputs...')

                        try
                            [data_bl_emu, data_fu_emu] = e2p_import_emu( ...
                                topDir_emu, gcm, ggcm, ssp, which_file, ...
                                cropList_in, gridlist, ts1_list, tsN_list, Nlist, ...
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
                            Nlist_emu, ~] = ...
                            e2p_get_names(data_bl_emu.varNames, data_fu_emu.varNames, ...
                            getN, get_unneeded) ;

                        disp('    Done.')

                    elseif ~strcmp(which_file, 'anpp')
                        error('which_file (%s) not recognized', which_file)
                    end
                    
                    false_xc = false(length(data_bl_emu.list2map),length(cropList_emu_basei)) ;
                    missing_agmerra_xc = false_xc ;
                    if (excl_lowBL_agmerra && strcmp(which_file, 'yield')) ...
                            || use_ph2_baseline
                        fprintf('    Importing AgMERRA %s...', which_file)
                        if ~isequal(data_fu_lpj.list2map, data_bl_emu.list2map)
                            error('gridlist mismatch')
                        end
                        data_bl_agm = e2p_get_agmerra_yield(...
                            varNames_emu, topDir_phase2, ggcm, ...
                            data_bl_emu.list2map, data_fu_lpj.lonlats, getN, ...
                            adaptation, which_file, use_ph2_baseline) ;
                        if ~any(any(~isnan(data_bl_agm.garr_xv)))
                            error('data_bl_agm.garr_xv is all NaN')
                        end
                        %% Translate to output crops
                        [varNames_agm, cropList_agm, ...
                            varNames_agm_basei, cropList_agm_basei, ...
                            Nlist_agm, ~] = ...
                            e2p_get_names(data_bl_agm.varNames, [], ...
                            getN, get_unneeded) ;
                        [cropList_agm_asEmu, used_emuCrops_agm] = e2p_translate_crops_lpj2emu( ...
                            cropList_agm, cropList_emu) ;
                        % Sanity check
                        if ~isequal(size(shiftdim(cropList_emu)), size(shiftdim(used_emuCrops_agm)))
                            error('Size mismatch between cropList_emu and used_emuCrops_agm')
                        end
                        % Find missing crops
                        missing_agmerra_xc = e2p_exclude_missing( ...
                            varNames_emu_basei, cropList_emu_basei, Nlist, ...
                            data_bl_agm.garr_xv) ;
                    end


                    %% Get and apply exclusions, if doing so
                    
                    % Where do we exclude based on missing
                    % emulation?
                    [missing_emu_bl_xc, all0_emu_bl_c] = e2p_exclude_missing( ...
                        varNames_emu_basei, cropList_emu_basei, Nlist, ...
                        data_bl_emu.garr_xv) ;
                    [missing_emu_fu_xc, all0_emu_fu_c] = e2p_exclude_missing( ...
                        varNames_emu_basei, cropList_emu_basei, Nlist, ...
                        data_fu_emu.garr_xvt) ;
                    missing_emu_xc = missing_emu_bl_xc | missing_emu_fu_xc ;
                    
                    if strcmp(which_file, 'yield')
                        
                        isexcl_lowBL_agmerra_xc = false_xc ;
                        isexcl_lowBL_emu_xc = false_xc ;
                        
                        % Where do we exclude based on low AgMERRA yield at max N
                        % (or existing exclusions)?
                        if excl_lowBL_agmerra || use_ph2_baseline
                            disp('    Excluding based on low AgMERRA yield at max N...')
                            isexcl_lowBL_agmerra_xc = e2p_exclude_lowBLyield_atMaxN( ...
                                varNames_emu, cropList_emu_basei, Nlist_emu, ...
                                data_bl_agm.garr_xv, low_yield_threshold_kgm2) ;
                        end
                        if ~any(any(~isexcl_lowBL_agmerra_xc))
                            error('All cells excluded because of low AgMERRA yield')
                        end
                        
                        % Where do we exclude based on NaN or low baseline-year emulated yield at
                        % max N (or existing exclusions)?
                        if excl_lowBL_emu
                            disp('    Excluding based on NaN/low baseline-year emulated yield at max N...')
                            isexcl_lowBL_emu_xc = e2p_exclude_lowBLyield_atMaxN( ...
                                varNames_emu, cropList_emu_basei, Nlist_emu, ...
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
                    
                    %% Save exclusion figures, if doing so
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
                    
                    exclude_xc = missing_yield_xc | missing_emu_xc ...
                        | exclude_lowBLyield_xc  ;
                    
                    for c = 1:length(cropList_emu_basei)
                        thisCrop = cropList_emu_basei{c} ;
                        thisCrop_i = find(strcmp(varNames_emu_basei,thisCrop)) ;
                        if length(thisCrop_i)~=length(Nlist_emu)
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
                        
                        %% If irrigation, check for any missing cells
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
                        
                        % Apply to emulated baseline
                        tmp = data_bl_emu.garr_xv(:,thisCrop_i) ;
                        tmp(exclude_xc(:,c),:) = NaN ;
                        data_bl_emu.garr_xv(:,thisCrop_i) = tmp ;
                        clear tmp
                        
                        % Apply to emulated future
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
                        if use_ph2_baseline
                            [data_bl_agm.garr_xv, data_bl_agm.varNames, is_ww_max_bl_gW, winter_wheats_test] = ...
                                e2p_get_max_wheat(data_bl_agm.garr_xv, data_bl_agm.varNames) ;
                            if ~isequal(winter_wheats, winter_wheats_test)
                                error('Mismatch between winter wheat lists from data_bl_emu vs. data_bl_agm')
                            end
                        end
                    elseif strcmp(which_file, 'gsirrigation')
                        data_bl_emu = e2p_apply_max_wheat(data_bl_emu, outDir_ggcm) ;
                        data_fu_emu = e2p_apply_max_wheat(data_fu_emu, outDir_ggcm) ;
                        if use_ph2_baseline
                            data_bl_agm = e2p_apply_max_wheat(data_bl_agm, outDir_ggcm) ;
                        end
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
                    if use_ph2_baseline
                        e2p_check_correct_zeros(data_bl_agm.garr_xv, ...
                            which_file, data_bl_agm.varNames, ...
                            'Baseline', @getbasenamei)
                    end
                    
                    % Refresh variable and crop lists
                    if refresh_vars
                        [varNames_emu, cropList_emu, ...
                            varNames_emu_basei, cropList_emu_basei, ...
                            Nlist_emu, ~] = ...
                            e2p_get_names(data_bl_emu.varNames, data_fu_emu.varNames, ...
                            getN, get_unneeded) ;
                        if use_ph2_baseline
                            [varNames_agm, cropList_agm, ...
                                varNames_agm_basei, cropList_agm_basei, ...
                                Nlist_agm, ~] = ...
                                e2p_get_names(data_bl_agm.varNames, [], ...
                                getN, get_unneeded) ;
                        end
                    end
                    
                    %% Set up output baseline structure
                    if use_lpjg_baseline
                        data_bl_out = data_bl_lpj ;
                    elseif use_ph2_baseline
                        % Rearrange and rename AgMERRA baseline to match
                        % LPJ-GUESS (output) variable names
                        data_bl_out.list2map = data_bl_agm.list2map ;
                        data_bl_out.lonlats = data_bl_agm.lonlats ;
                        I = e2p_translate_crops_agm2out(...
                            varNames_agm, varNames_out) ;
                        data_bl_out.varNames = varNames_out ;
                        data_bl_out.garr_xv = data_bl_agm.garr_xv(:,I) ;
                    elseif ~(use_lpjg_baseline || use_ph2_baseline)
                        error('Which baseline are you using? I can''t set up data_bl_out.')
                    end

                    disp('    Done.')


                    %% Harmonize LPJ-GUESS and emulator data

                    disp('    Harmonizing LPJ-GUESS and emulator data...')

                    % Align gridlists
                    [data_bl_lpj, data_fu_lpj, data_bl_emu, data_fu_emu, list2map] = ...
                        e2p_align_gridlists(data_bl_lpj, data_fu_lpj, data_bl_emu, data_fu_emu) ;

                    e2p_check_correct_zeros(data_bl_lpj.garr_xv, ...
                        which_file, data_bl_lpj.varNames, ...
                        'Baseline', @getbasenamei)
                    e2p_check_correct_zeros(data_fu_lpj.garr_xvt, ...
                        which_file, data_fu_lpj.varNames, ...
                        'Future', @getbasenamei)
                    e2p_check_correct_zeros(data_bl_emu.garr_xv, ...
                        which_file, data_bl_emu.varNames, ...
                        'Baseline', @getbasenamei)
                    e2p_check_correct_zeros(data_fu_emu.garr_xvt, ...
                        which_file, data_fu_emu.varNames, ...
                        'Future', @getbasenamei)

                    % Make sure that N lists match
                    if use_lpjg_baseline && ~isequal(Nlist_lpj,Nlist_emu)
                        error('Mismatch in N levels between LPJ-GUESS baseline and emulator future')
                    elseif use_ph2_baseline && ~isequal(Nlist_agm,Nlist_emu)
                        error('Mismatch in N levels between Phase 2 baseline and emulator future')
                    elseif ~(use_lpjg_baseline || use_ph2_baseline)
                        error('Which baseline are you using? I can''t check for matching N levels.')
                    end

                    % Translate crop names
                    [cropList_lpj_asEmu, used_emuCrops_lpj] = e2p_translate_crops_lpj2emu( ...
                        cropList_lpj, cropList_emu) ;
                    cropList_out_asEmu = cropList_lpj_asEmu ;
                    if use_lpjg_baseline
                        used_emuCrops = used_emuCrops_lpj ;
                    elseif use_ph2_baseline
                        [cropList_agm_asEmu, used_emuCrops_agm] = e2p_translate_crops_lpj2emu( ...
                            cropList_agm, cropList_emu) ;
                        used_emuCrops = used_emuCrops_agm ;
                    else
                        error('What used_emuCrops should I use?')
                    end
                    % Sanity check
                    if ~isequal(size(shiftdim(cropList_emu)), size(shiftdim(used_emuCrops)))
                        error('Size mismatch between cropList_emu and used_emuCrops')
                    end

                    disp('    Done.')


                    %% Get deltas
                    if contains(which_file, {'yield','gsirrigation'})
                        disp('    Getting deltas...')

                        deltas_emu_xvt = e2p_get_deltas(...
                            data_bl_emu, data_fu_emu, interp_infs, cropList_emu, ...
                            which_file, ...
                            used_emuCrops, list2map, ...
                            save_interp_figs, outDir_interp_figs, ggcm, figure_visibility, ...
                            when_remove_outliers, outDir_ggcm, renderer) ;
                        e2p_check_correct_zeros(deltas_emu_xvt, ...
                            which_file, data_fu_emu.varNames, ...
                            'Future', @getbasenamei)
                        disp('    Done.')

                    elseif ~strcmp(which_file, 'anpp')
                        error('which_file (%s) not recognized', which_file)
                    end


                    %% Apply emulator deltas to chosen baseline.
                    disp('    Applying deltas...')
                    data_fu_out = e2p_apply_deltas( ...
                        data_bl_out, data_bl_emu, data_fu_emu, deltas_emu_xvt, ...
                        cropList_out, cropList_out_asEmu, varNames_out, ...
                        list2map, which_file, figure_visibility) ;
                    
                    % Sort variable names in data_fu_lpj (was previously
                    % done in e2p_apply_deltas()
                    if ~isequal(data_fu_lpj.varNames, sort(data_fu_lpj.varNames))
                        [data_fu_lpj.varNames, I] = sort(data_fu_lpj.varNames) ;
                        data_fu_lpj.garr_xvt = data_fu_lpj.garr_xvt(:,I,:) ;
                        clear I
                        e2p_check_correct_zeros(data_fu_lpj.garr_xvt, ...
                            which_file, data_fu_lpj.varNames, ...
                            'Baseline', @getbasenamei)
                    end
                    
                    
                    %% Get fake N1000 values
                    if fake1k
                        data_fu_lpj = e2p_fake_1000(data_fu_lpj, data_bl_lpj0, ...
                            cropList_lpj, Nlist_lpj0, Nlist_emu, which_file) ;
                        data_fu_out = e2p_fake_1000(data_fu_out, data_bl_lpj0, ...
                            cropList_lpj, Nlist_lpj0, Nlist_emu, which_file) ;
                    end
                    
                    %% Remove outliers
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
                        e2p_save_outlier_info(outlier_info_lpj, outDir_lpj, which_file, data_fu_out.y1s, data_fu_out.yNs)
                        e2p_save_outlier_info(outlier_info_out, outDir_ggcm, which_file, data_fu_out.y1s, data_fu_out.yNs)
                    end

                end
                
                % Sort. Technically unnecessary, but maybe PLUM isn't
                % robust to non-sorted crops.
                data_fu_lpj = do_sort(data_fu_lpj, getN) ;
                data_fu_out = do_sort(data_fu_out, getN) ;

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
                if save_txt_files
                    disp('    Saving txt files for PLUM...')
                    out_header_cell = [{'Lon', 'Lat'} data_fu_out.varNames] ;
                    for t = 1:Ntpers
                        y1 = data_fu_out.y1s(t) ;
                        yN = data_fu_out.yNs(t) ;

                        if strcmp(ggcm, ggcm_list{1})
                            fprintf('    %d/%d lpj', t, Ntpers)
                            e2p_save(outDir_lpj, y1, yN, out_header_cell, ...
                                data_fu_lpj.lonlats, data_fu_lpj.garr_xvt(:,:,t), which_file, ...
                                false)
                        end
                        fprintf('    %d/%d emu', t, Ntpers)
                        e2p_save(outDir_ggcm, y1, yN, out_header_cell, ...
                            data_fu_out.lonlats, data_fu_out.garr_xvt(:,:,t), which_file, ...
                            overwrite_existing_txt)

                    end
                end

                % Save yield diagnostic figures, if doing so
                if save_out_figs
                    disp('    Saving output figures...')
                    if ~isempty(data_fu_lpj0) && ~save_out_figs_Nth0
                        tmp = [] ;
                    else
                        tmp = data_fu_lpj0 ;
                    end
                    if strcmp(which_file, 'yield')
                        e2p_save_out_figs(data_fu_lpj, tmp, ...
                            data_fu_emu, data_fu_out, ...
                            ggcm, getN, outDir_yield_figs, ...
                            which_file, cropList_lpj_asEmu, figure_visibility, figure_extension, ...
                            which_out_figs, overwrite_existing_figs, renderer, ...
                            use_lpjg_baseline, use_ph2_baseline)
                    elseif strcmp(which_file, 'gsirrigation')
                        e2p_save_out_figs(data_fu_lpj, tmp, ...
                            data_fu_emu, data_fu_out, ...
                            ggcm, getN, outDir_irrig_figs, ...
                            which_file, cropList_lpj_asEmu, figure_visibility, figure_extension, ...
                            which_out_figs, overwrite_existing_figs, renderer, ...
                            use_lpjg_baseline, use_ph2_baseline)
                    else
                        warning('which_file (%s) not recognized for save_yield_figs; skipping.', which_file)
                    end
                    clear tmp
                end

                fprintf('Done with %s %s %s %s.\n', gcm, ssp, ggcm, which_file)
                
            end
            
            fprintf('Done with %s %s %s.\n', gcm, ssp, ggcm)

        end

        fprintf('Done with %s %s.\n', gcm, ssp)

    end
    
    fprintf('Done with %s.\n', gcm)
end

disp(' ')
disp('All done!')


%% FUNCTIONS

function S = do_sort(S, getN)

% Sort N strings
Nlist = unique(getN(S.varNames)) ;
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





