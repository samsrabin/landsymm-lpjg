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
        
        outDir_lpj = sprintf('%s/sim_LPJ-GUESS', outDir) ;
        outDir_excl_figs_inCrops = sprintf('%s/excl_figs_GGCMIcrops', outDir) ;
        outDir_excl_figs_outCrops = sprintf('%s/excl_figs_PLUMcrops', outDir) ;
        outDir_interp_figs = sprintf('%s/interp_figs', outDir) ;
        outDir_combineCrops_figs = sprintf('%s/combineCrops_figs', outDir) ;
        outDir_yield_figs = sprintf('%s/yield_figs', outDir) ;
        outDir_irrig_figs = sprintf('%s/irrig_figs', outDir) ;
        
        if ~exist(outDir_lpj, 'dir')
            mkdir(outDir_lpj) ;
        end
        
        % Import LPJ-GUESS yield and irrigation
        if ~exist('topDir_lpj', 'var') || ~exist(topDir_lpj, 'dir')
            if exist('topDir_lpj', 'var') && ~exist(topDir_lpj, 'dir')
                warning('topDir_lpj (%s) not found; using latest directory instead', ...
                    topDir_lpj)
            end
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
        end
        fprintf('topDir_lpj: %s\n', topDir_lpj)
        clear outdirs tmp
        disp('Importing LPJ-GUESS yield...')
        which_file = 'yield' ;
        data_fu_lpj0_yield = e2p_import_fu_lpj(future_y1-1, future_ts, future_yN_lpj, topDir_lpj, ...
            which_file, ssp) ;
        e2p_check_correct_zeros(data_fu_lpj0_yield.garr_xvt, ...
            which_file, data_fu_lpj0_yield.varNames, ...
            'Future', @getbasenamei)
        [varNames_lpj0, cropList_lpj0, ...
            varNames_lpj0_basei, cropList_lpj0_basei, ...
            Nlist_lpj0_char, ~] = ...
            e2p_get_names({}, data_fu_lpj0_yield.varNames) ;
        cropList_cf_basei = cropIrrList_cf ;
        cropList_out_basei = cropIrrList_out ;
        gridlist.list_to_map = data_fu_lpj0_yield.list2map ;
        gridlist.lonlats = data_fu_lpj0_yield.lonlats ;
        gridlist.mask_YX = false(360, 720) ;
        gridlist.mask_YX(gridlist.list_to_map) = true ;
        Ncells = length(gridlist.list_to_map) ;
        gridlist_target = {gridlist.lonlats gridlist.list_to_map} ;
        disp('Importing LPJ-GUESS irrigation...')
        which_file = 'gsirrigation' ;
        data_fu_lpj0_irrig = e2p_import_fu_lpj(future_y1-1, future_ts, future_yN_lpj, topDir_lpj, ...
            which_file, ssp, gridlist_target) ;
        e2p_check_correct_zeros(data_fu_lpj0_irrig.garr_xvt, ...
            which_file, data_fu_lpj0_irrig.varNames, ...
            'Future', @getbasenamei)
        [~, cropList_lpj0i, ...
            varNames_lpj0_basei2, cropList_lpj0_basei2, ...
            Nlist_lpj0_char2, ~] = ...
            e2p_get_names({}, data_fu_lpj0_irrig.varNames) ;
        if ~isequal(cropList_lpj0, cropList_lpj0i)
            error('Mismatch between cropList_lpj0 for yield vs. gsirrigation')
        elseif ~isequal(varNames_lpj0_basei, varNames_lpj0_basei2)
            error('Mismatch between varNames_lpj0_basei for yield vs. gsirrigation')
        elseif ~isequal(cropList_lpj0_basei, cropList_lpj0_basei2)
            error('Mismatch between cropList_lpj0_basei for yield vs. gsirrigation')
        elseif ~isequal(Nlist_lpj0_char, Nlist_lpj0_char2)
            error('Mismatch between Nlist_lpj0_char2 for yield vs. gsirrigation')
        end
        clear cropList_lpj02 varNames_lpj0_basei2 cropList_lpj0_basei2 Nlist_lpj0_char2
        
        % Rename crops to match GGCMI names
        if any(strcmp(cropList_lpj0, 'CerealsC3'))
            error('This version assumes that you simulated spring and winter wheat separately')
        end
        cropList_lpj1 = e2p_translate_crops_2emu(cropList_lpj0, cropList_in, 'LPJ-GUESS') ;
        if length(cropList_lpj1) ~= length(cropList_lpj0)
            error('Error renaming LPJ-GUESS crops to match GGCMI names: Length mismatch???')
        elseif length(cropList_lpj1) ~= length(unique(cropList_lpj1))
            error('Error renaming LPJ-GUESS crops to match GGCMI names: Non-unique names in result')
        end
        varNames_lpj1 = varNames_lpj0 ;
        varNames_lpj1_basei = varNames_lpj0_basei ;
        cropList_lpj1_basei = cropList_lpj0_basei ;
        for c = 1:length(cropList_lpj1)
            thisCrop0 = cropList_lpj0{c} ;
            thisCrop1 = cropList_lpj1{c} ;
            varNames_lpj1 = strrep(varNames_lpj1, thisCrop0, thisCrop1) ;
            varNames_lpj1_basei = strrep(varNames_lpj1_basei, thisCrop0, thisCrop1) ;
            cropList_lpj1_basei = strrep(cropList_lpj1_basei, thisCrop0, thisCrop1) ;
        end
        data_fu_lpj1_yield = data_fu_lpj0_yield;
        data_fu_lpj1_yield.varNames = varNames_lpj1 ;
        data_fu_lpj1_irrig = data_fu_lpj0_irrig;
        data_fu_lpj1_irrig.varNames = varNames_lpj1 ;
        if save_txt_files_lpjg
            disp('Importing LPJ-GUESS ANPP...')
            data_fu_lpj_anpp = e2p_import_fu_lpj(...
                future_y1-1, future_ts, future_yN_lpj, ...
                topDir_lpj, 'anpp', ssp, gridlist_target) ;
            out_header_cell_anpp = [{'Lon', 'Lat'} data_fu_lpj_anpp.varNames] ;
            lpj_anpp_precision = get_precision(data_fu_lpj_anpp.garr_xvt) ;
            
            disp('Importing LPJ-GUESS runoff...')
            data_fu_lpj_runoff = e2p_import_fu_lpj(...
                future_y1-1, future_ts, future_yN_lpj, ...
                topDir_lpj, 'tot_runoff', ssp, gridlist_target) ;
            out_header_cell_runoff = [{'Lon', 'Lat'} data_fu_lpj_runoff.varNames] ;
            lpj_runoff_precision = get_precision(data_fu_lpj_runoff.garr_xvt) ;
        end
                
        %% Get LPJ-GUESS calibration factors
        cf_lpj = e2p_get_CFs(cropList_cf, 'LPJ-GUESS', cfDir, combineCrops, true) ;
        
        % Split crops up to match calibration factor list; combine as
        % necessary to match PLUM-expected list, burning in calibration
        % factors of combined crops.
        disp('Splitting and combining...')
        [data_fu_lpj2_yield, combineCrops_lpj, ~, cf_lpj] = ...
            e2p_split_combine_burn(data_fu_lpj1_yield, ...
            combineCrops, cropList_out, cf_lpj) ;
        save_cf_csv(cf_lpj, outDir_lpj)
        out_file = sprintf('%s/combineCrops_lpj_%s.mat', topDir_lpj, ssp) ;
        save_combineCrops(combineCrops_lpj, 'yield', out_file)
        data_fu_lpj2_irrig = e2p_apply_max(data_fu_lpj1_irrig, out_file) ;
        if save_combineCrops_figs
            e2p_save_combineCrops_figs(combineCrops_lpj, gridlist, ...
                'LPJ-GUESS', figure_visibility, ...
                figure_extension, outDir_combineCrops_figs, ...
                overwrite_existing_figs, renderer)
        end
        
        % Rename variables to match desired outputs
        data_fu_lpj2_yield = rename_and_rearrange(data_fu_lpj2_yield, ...
            varNames_out_allN) ;
        data_fu_lpj2_irrig = rename_and_rearrange(data_fu_lpj2_irrig, ...
            varNames_out_allN) ;
        
        % Sort. Technically unnecessary, but maybe PLUM isn't
        % robust to non-sorted crops.
        data_fu_lpj2_yield = do_sort(data_fu_lpj2_yield) ;
        data_fu_lpj2_irrig = do_sort(data_fu_lpj2_irrig) ;
        
        [varNames_lpj2, cropList_lpj2, ...
            varNames_lpj2_basei, cropList_lpj2_basei] = ...
            e2p_get_names(data_fu_lpj2_yield.varNames, data_fu_lpj2_irrig.varNames) ;
                
        disp('Done importing LPJ-GUESS yield and irrigation.')
                        
        % Loop through crop emulators
        for ggcm_counter = 1:length(ggcm_list)
            ggcm = ggcm_list{ggcm_counter} ;
            topDir_phase2 = sprintf('%s/%s/phase2', topDir_agmipout, ggcm) ;
            outDir_ggcm = sprintf('%s/emu_%s', outDir, ggcm) ;
            if ~exist(outDir_ggcm, 'dir')
                mkdir(outDir_ggcm) ;
            end
            
            % Copy log to file
            diary('off')
            diaryfile = sprintf('%s/matlab_log_%s.txt', ...
                outDir_ggcm, datetime('now', 'Format', 'yyyyMMdd_HHmmss')) ;
            if exist(diaryfile, 'file')
                delete(diaryfile) ;
            end
            diary(diaryfile)
            diary('on')
            
            % Get this emulator's calibration factors
            cf_emu0 = e2p_get_CFs(cropList_cf, ggcm, cfDir, combineCrops, true) ;
            
            for w = 1:length(whichfile_list)
                which_file = whichfile_list{w} ;
                fprintf('%s %s %s %s\n', gcm, ssp, ggcm, which_file)
                
                % Get LPJ-GUESS simulation data for this file
                if strcmp(which_file, 'yield')
                    data_fu_lpj1 = data_fu_lpj1_yield ;
                    data_fu_lpj2 = data_fu_lpj2_yield ;
                elseif strcmp(which_file, 'gsirrigation')
                    data_fu_lpj1 = data_fu_lpj1_irrig ;
                    data_fu_lpj2 = data_fu_lpj2_irrig ;
                else
                    error('which_file (%s) not recognized', which_file)
                end

                % If allowed to import existing file, then do so
                out_file = sprintf('%s/future_%s.mat', outDir_ggcm, which_file) ;
                did_load_existing = load_existing_file && exist(out_file, 'file') ;
                if did_load_existing
                    error('load_existing_file only tested with use_lpjg_baseline')
                    
                    fprintf('    Importing future_%s.mat...\n', which_file) ;

                    load(out_file, 'data_fu_lpj1', 'data_fu_emu', 'data_fu_out')
                    
                    % Translate crop names
                    [varNames_emu, cropList_emu, ...
                        varNames_emu_basei, cropList_emu_basei, ...
                        Nlist_emu0_char, ~] = ...
                        e2p_get_names([], data_fu_emu.varNames) ;
                else
                    
                    % Import emulator outputs
                    
                    if contains(which_file, {'yield', 'gsirrigation'})
                        
                        disp('    Importing emulator outputs...')
                        try
                            [data_bl_emu0, data_fu_emu0] = e2p_import_emu( ...
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
                        
                        e2p_check_correct_zeros(data_bl_emu0.garr_xv, ...
                            which_file, data_bl_emu0.varNames, ...
                            'Baseline', @getbasenamei, ...
                            true)
                        e2p_check_correct_zeros(data_fu_emu0.garr_xvt, ...
                            which_file, data_fu_emu0.varNames, ...
                            'Future', @getbasenamei, ...
                            true)

                        [varNames_emu0, cropList_emu0, ...
                            varNames_emu0_basei, cropList_emu0_basei, ...
                            Nlist_emu0_char, ~] = ...
                            e2p_get_names(data_bl_emu0.varNames, data_fu_emu0.varNames) ;

                        disp('    Done.')

                    elseif ~strcmp(which_file, 'anpp')
                        error('which_file (%s) not recognized', which_file)
                    end
                    
                    false_xc = false(length(data_bl_emu0.list2map),length(cropList_emu0_basei)) ;
                    if ~isequal(data_fu_lpj1.list2map, data_bl_emu0.list2map)
                        error('gridlist mismatch')
                    end
                    if emulated_baseline
                        fprintf('    Using emulated %s as "true" baseline...\n', which_file)
                        data_bl_agm = data_bl_emu0 ;
                        data_bl_agm.list2map = gridlist.list_to_map ;
                        data_bl_agm.lonlats = gridlist.lonlats ;
                        data_bl_agm.actually_emuBL = true(size(data_bl_agm.varNames)) ;
                    else
                        fprintf('    Importing AgMERRA %s...\n', which_file)
                        data_bl_agm = e2p_import_bl_ggcmi(...
                            varNames_emu, topDir_phase2, topDir_emu, ggcm, ...
                            data_bl_emu0.list2map, data_fu_lpj1.lonlats, ...
                            adaptation, which_file, force_consistent_baseline) ;
                    end
                    
                    % Sanity check
                    if ~any(any(~isnan(data_bl_agm.garr_xv)))
                        error('data_bl_agm.garr_xv is all NaN')
                    end
                    
                    % Convert booleans to chars
                    actually_emuBL_char = cell(size(data_bl_agm.actually_emuBL)) ;
                    actually_emuBL_char(~data_bl_agm.actually_emuBL) = {'sim'} ;
                    actually_emuBL_char(data_bl_agm.actually_emuBL) = {'*EMU*'} ;
                    data_bl_agm.actually_emuBL_char = actually_emuBL_char ;
                    clear actually_emuBL_char
                    
                    % Get info about Phase 2 variables
                    [varNames_agm, cropList_agm, ...
                        varNames_agm_basei, cropList_agm_basei, ...
                        ~, ~] = ...
                        e2p_get_names(data_bl_agm.varNames, []) ;
                    [~, used_emuCrops_agm] = e2p_translate_crops_2emu( ...
                        cropList_agm, cropList_emu0, 'Phase 2 sims') ;
                    % Sanity check
                    if ~isequal(size(shiftdim(cropList_emu0)), size(shiftdim(used_emuCrops_agm)))
                        error('Size mismatch between cropList_emu0 and used_emuCrops_agm')
                    end
                    

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% Emulator and/or phase 2 exclusions %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % Where do we exclude based on missing
                    % emulation (or Phase 2, if using)?
                    missing_agmerra_xc = e2p_exclude_missing( ...
                        varNames_emu0_basei, cropList_emu0_basei, Nlist_emu0_char, ...
                        data_bl_agm.garr_xv) ;
                    [missing_emu_bl_xc, all0_emu_bl_c] = e2p_exclude_missing( ...
                        varNames_emu0_basei, cropList_emu0_basei, Nlist_emu0_char, ...
                        data_bl_emu0.garr_xv) ;
                    [missing_emu_fu_xc, all0_emu_fu_c] = e2p_exclude_missing( ...
                        varNames_emu0_basei, cropList_emu0_basei, Nlist_emu0_char, ...
                        data_fu_emu0.garr_xvt) ;
                    missing_emu_xc = missing_emu_bl_xc | missing_emu_fu_xc ;
                    
                    exclude_xc_file = sprintf('%s/missing_or_excluded.mat', outDir_ggcm) ;
                    if strcmp(which_file, 'yield')
                        
                        isexcl_lowBL_emu_xc = false_xc ;
                        
                        % Where do we exclude based on low AgMERRA yield at max N
                        % (or existing exclusions)?
                        disp('    Excluding based on low AgMERRA yield at max N...')
                        isexcl_lowBL_agmerra_xc = e2p_exclude_lowBLyield_atMaxN( ...
                            varNames_emu0, cropList_emu0_basei, Nlist_emu0_char, ...
                            data_bl_agm.garr_xv, low_yield_threshold_kgm2) ;
                        if ~any(any(~isexcl_lowBL_agmerra_xc))
                            error('All cells excluded because of low AgMERRA yield')
                        end
                        
                        % Where do we exclude based on NaN or low baseline-year emulated yield at
                        % max N (or existing exclusions)?
                        if excl_lowBL_emu
                            disp('    Excluding based on NaN/low baseline-year emulated yield at max N...')
                            isexcl_lowBL_emu_xc = e2p_exclude_lowBLyield_atMaxN( ...
                                varNames_emu0, cropList_emu0_basei, Nlist_emu0_char, ...
                                data_bl_emu0.garr_xv, low_yield_threshold_kgm2) ;
                        end
                        
                        if ~any(any(~isexcl_lowBL_emu_xc))
                            error('All cells excluded because of low AgMERRA yield and/or low baseline-year emulated yield')
                        end
                        
                        exclude_lowBLyield_xc = isexcl_lowBL_agmerra_xc | isexcl_lowBL_emu_xc ;
                        missing_yield_xc = missing_emu_xc | missing_agmerra_xc ;
                        
                        % Save exclusion info
                        save(exclude_xc_file, ...
                            'exclude_lowBLyield_xc', 'missing_yield_xc')
                        
                        % Set up for figures
                        excl_vecs_emu0{1} = missing_emu_xc ;
                        excl_vecs_emu0{2} = missing_agmerra_xc ;
                        excl_vecs_emu0{3} = isexcl_lowBL_emu_xc ;
                        excl_vecs_emu0{4} = isexcl_lowBL_agmerra_xc ;
                        
                    elseif strcmp(which_file, 'gsirrigation')
                        disp('    Applying yield-based exclusions...')
                        load(exclude_xc_file) ;
                        
                        % Set up for figures
                        excl_vecs_emu0{1} = missing_emu_xc ;
                        excl_vecs_emu0{2} = missing_yield_xc ;
                        excl_vecs_emu0{3} = exclude_lowBLyield_xc ;
                        excl_vecs_emu0{4} = all0_emu_bl_c | all0_emu_fu_c ;
                        if size(excl_vecs_emu0{4}, 2) == 1
                            excl_vecs_emu0{4} = excl_vecs_emu0{4}' ;
                        end
                    else
                        error('which_file (%s) not recognized', which_file)
                    end
                    
                    % Save exclusion figures, if doing so
                    if save_excl_figs
                        disp('    Saving exclusion figures (input crops)...')
                        e2p_save_excl_figs_inCrops(...
                            ggcm, which_file, gridlist, excl_vecs_emu0, ...
                            cropList_emu0_basei, figure_visibility, ...
                            figure_extension, outDir_excl_figs_inCrops, ...
                            overwrite_existing_figs, renderer)
                    end
                    
                    % Apply exclusions
                    exclude_xc = missing_yield_xc | missing_emu_xc ...
                        | exclude_lowBLyield_xc  ;
                    for c = 1:length(cropList_emu0_basei)
                        thisCrop = cropList_emu0_basei{c} ;
                        thisCrop_i = find(strcmp(varNames_emu0_basei,thisCrop)) ;
                        if length(thisCrop_i)~=length(Nlist_emu0_char)
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
                            isbad = isnan(data_bl_emu0.garr_xv(:,thisCrop_i)) & ~exclude_xc(:,c) ;
                            if any(any(isbad))
                                warning('%s: %d cells that were included in yield are NaN in irrig baseline emulation', ...
                                    thisCrop, length(find(isbad)))
                            end
                            isbad = isnan(data_fu_emu0.garr_xvt(:,thisCrop_i,:)) & repmat(~exclude_xc(:,c),[1 length(thisCrop_i) Ntpers]) ;
                            if any(any(any(isbad)))
                                warning('%s: %d cells that were included in yield are NaN in irrig future emulation', ...
                                    thisCrop, length(find(isbad)))
                            end
                        end
                        
                        % Apply exclusions to emulated baseline
                        tmp = data_bl_emu0.garr_xv(:,thisCrop_i) ;
                        tmp(exclude_xc(:,c),:) = NaN ;
                        data_bl_emu0.garr_xv(:,thisCrop_i) = tmp ;
                        clear tmp
                        
                        % Apply exclusions to emulated future
                        tmp = data_fu_emu0.garr_xvt(:,thisCrop_i,:) ;
                        tmp(exclude_xc(:,c),:,:) = NaN ;
                        data_fu_emu0.garr_xvt(:,thisCrop_i,:) = tmp ;
                        clear tmp
                    end
                    
                    disp('    Done.')
                    
                    % Sanity checks
                    if ~any(any(~isnan(data_bl_emu0.garr_xv)))
                        error('data_bl_emu0.garr_xv is all NaN!')
                    end
                    if ~any(any(any(~isnan(data_fu_emu0.garr_xvt))))
                        error('data_fu_emu0.garr_xvt is all NaN!')
                    end
                    
                    % Align gridlists
                    [data_fu_lpj1, data_bl_emu0, data_fu_emu0] = ...
                        e2p_align_gridlists(data_fu_lpj1, data_bl_emu0, data_fu_emu0) ;
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% Apply emulator deltas to Phase 2 baseline %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    [ok_missing_v, varNames_NOTok_missing] = ...
                        get_ok_missing(data_fu_emu0.varNames, Nlist_emu) ;
                    
                    % Get deltas
                    if contains(which_file, {'yield','gsirrigation'})
                        disp('    Getting deltas...')

                        deltas_emu_xvt = e2p_get_deltas(...
                            data_bl_emu0, data_fu_emu0, interp_infs, cropList_emu0, ...
                            which_file, ...
                            used_emuCrops_agm, ...
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
                    data_fu_emu1 = e2p_apply_deltas( ...
                        data_bl_agm, data_bl_emu0, data_fu_emu0, deltas_emu_xvt) ;
                    % Check
                    tmp_garr_xvt = data_fu_emu1.garr_xvt(:,~ok_missing_v,:) ;
                    e2p_check_correct_zeros(tmp_garr_xvt, ...
                        which_file, varNames_NOTok_missing, ...
                        'Future', @getbasenamei)
                    
                    % Update names (v1 should be the same as v0)
                    [varNames_emu1, cropList_emu1, ...
                        varNames_emu1_basei, cropList_emu1_basei, ...
                        Nlist_emu1_char] = ...
                        e2p_get_names(data_fu_emu0.varNames, ...
                        data_fu_emu1.varNames) ;
                    
                    % Copy this over
                    data_fu_emu1.actually_emuBL_char = ...
                        data_bl_agm.actually_emuBL_char ;
                    
                    disp('    Done.')
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% Get N1000; remove outliers %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % Rearrange and emulated future to match data_fu_lpj1,
                    % expanding with NaN as needed where N values are missing.
                    tmp_garr_xvt = data_fu_emu1.garr_xvt ;
                    data_fu_emu2 = rmfield(data_fu_emu1, 'garr_xvt') ;
                    data_fu_emu2.list2map = data_fu_emu1.list2map ;
                    data_fu_emu2.lonlats = gridlist.lonlats ;
                    % Rename Ns, if needed
                    C = intersect(Nlist_lpj0_char, Nlist_emu1_char) ;
                    if ~isequal(sort(C), sort(Nlist_emu1_char))
                        tmp_varNames = data_fu_emu1.varNames ;
                        for n = 1:length(Nlist_emu)
                            thisN_emu = Nlist_emu(n) ;
                            I_thisN_lpj = find(Nlist_lpj == thisN_emu) ;
                            if length(I_thisN_lpj) ~= 1
                                error('Expected 1 match; found %d', I_thisN_lpj)
                            end
                            thisN_emu_char = Nlist_emu1_char{n} ;
                            thisN_lpj_char = Nlist_lpj0_char{I_thisN_lpj} ;
                            if ~strcmp(thisN_emu_char, thisN_lpj_char)
                                tmp_varNames = strrep(tmp_varNames, ...
                                    thisN_emu_char, thisN_lpj_char) ;
                            end
                        end
                        data_fu_emu1.varNames = tmp_varNames ;
                        [~, ~, ~, ~, Nlist_emu1_char] = ...
                            e2p_get_names(data_fu_emu1.varNames, []) ;
                        if ~isequal(varNames_emu0, varNames_emu1)
                            error('~isequal(varNames_emu0, varNames_emu1)')
                        end
                        data_fu_emu0.varNames = data_fu_emu1.varNames ;
                        [varNames_emu0, cropList_emu0, ...
                            varNames_emu0_basei, cropList_emu0_basei, ...
                            Nlist_emu0_char] = ...
                            e2p_get_names(data_fu_emu0.varNames, []) ;
                    end
                    % (Continue)
                    [~, I_lpj, I_emu] = intersect(data_fu_lpj1.varNames, ...
                        data_fu_emu1.varNames, 'stable') ;
                    if length(I_emu) ~= length(data_fu_emu1.varNames)
                        error('Not all data_fu_emu1.varNames found in data_fu_lpj1.varNames?')
                    end
                    data_fu_emu2.varNames = data_fu_lpj1.varNames ;
                    data_fu_emu2.garr_xvt = nan(Ncells, length(data_fu_lpj1.varNames), Ntpers) ;
                    data_fu_emu2.garr_xvt(:,I_lpj,:) = tmp_garr_xvt(:,I_emu,:) ;
                    clear tmp_garr_xvt
                    data_fu_emu2 = rmfield(data_fu_emu2, 'actually_emuBL') ;
                    data_fu_emu2.actually_emuBL_char = cell(size(data_fu_emu2.varNames)) ;
                    [C, I1, I2] = intersect(data_fu_emu1.varNames, ...
                        data_fu_emu2.varNames, 'stable') ;
                    if ~isequal(shiftdim(C), shiftdim(data_fu_emu1.varNames))
                        error('Expected each member of data_fu_emu1.varNames to be found exactly once in data_fu_emu2.varNames')
                    end
                    data_fu_emu2.actually_emuBL_char(I2) = ...
                        data_fu_emu1.actually_emuBL_char(I1) ;

                    % Check
                    e2p_check_correct_zeros(data_fu_lpj1.garr_xvt, ...
                        which_file, data_fu_lpj1.varNames, ...
                        'Future', @getbasenamei)
                    ok_missing_v = get_ok_missing(data_fu_emu2.varNames, Nlist_emu) ;
                    varNames_NOTok_missing = data_fu_emu2.varNames(~ok_missing_v) ;
                    tmp_garr_xvt = data_fu_emu2.garr_xvt(:,~ok_missing_v,:) ;
                    e2p_check_correct_zeros(tmp_garr_xvt, ...
                        which_file, varNames_NOTok_missing, ...
                        'Future', @getbasenamei)
                    clear tmp_garr_xv tmp_garr_xvt
                    
                    % Get fake N1000 values
                    data_fu_emu2 = e2p_fake_1000(data_fu_emu2, data_fu_lpj1, ...
                        Nlist_lpj0_char, Nlist_emu1_char, which_file, ...
                        scale_200to1000) ;
                    
                    % Update names
                    [varNames_emu2, cropList_emu2, ...
                        varNames_emu2_basei, cropList_emu2_basei, ...
                        Nlist_emu2_char] = ...
                        e2p_get_names(data_fu_emu2.varNames, ...
                        []) ;
                                        
                    % Remove outliers
                    if strcmp(when_remove_outliers, 'end')
                        disp('    Removing outliers...')

                        [data_fu_lpj1, outlier_info_lpj] = e2p_remove_outliers(data_fu_lpj1, which_file) ;
                        e2p_check_correct_zeros(data_fu_lpj1.garr_xvt, ...
                            which_file, data_fu_lpj1.varNames, ...
                            'Future', @getbasenamei)

                        [data_fu_emu2, outlier_info_out] = e2p_remove_outliers(data_fu_emu2, which_file) ;
                        e2p_check_correct_zeros(data_fu_emu2.garr_xvt, ...
                            which_file, data_fu_emu2.varNames, ...
                            'Future', @getbasenamei)

                        % Save info
                        if strcmp(ggcm, ggcm_list{1})
                            e2p_save_outlier_info(outlier_info_lpj, outDir_lpj, which_file, data_fu_emu2.y1s, data_fu_emu2.yNs)
                        end
                        e2p_save_outlier_info(outlier_info_out, outDir_ggcm, which_file, data_fu_emu2.y1s, data_fu_emu2.yNs)
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% Split, combine, and burn %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    disp('Splitting and combining...')

                    refresh_vars = true ;
                    emu_combineCrops_file = ...
                        sprintf('%s/combineCrops_fu_emu.mat', outDir_ggcm) ;
                    
                    if strcmp(which_file, 'yield')
                        [data_fu_emu3, combineCrops_emu3, ~, cf_emu3, ...
                            excl_vecs_emu3] = ...
                            e2p_split_combine_burn(data_fu_emu2, ...
                            combineCrops, cropList_out, cf_emu0, ...
                            excl_vecs_emu0) ;
                        save_cf_csv(cf_emu3, outDir_ggcm)
                        % For yield, save result combineCrops
                        save_combineCrops(combineCrops_emu3, which_file, ...
                            emu_combineCrops_file)
                        if save_combineCrops_figs
                            e2p_save_combineCrops_figs(combineCrops_emu3, gridlist, ...
                                ggcm, figure_visibility, ...
                                figure_extension, outDir_combineCrops_figs, ...
                                overwrite_existing_figs, renderer)
                        end
                    elseif strcmp(which_file, 'gsirrigation')
                        [data_fu_emu3, excl_vecs_emu3] = e2p_apply_max( ...
                            data_fu_emu2, emu_combineCrops_file, ...
                            excl_vecs_emu0) ;
                    elseif ~strcmp(which_file, 'anpp')
                        error('which_file (%s) not recognized', which_file)
                    else
                        refresh_vars = false ;
                    end
                    
                    e2p_check_correct_zeros(data_fu_emu3.garr_xvt, ...
                        which_file, data_fu_emu3.varNames, ...
                        'Future', @getbasenamei)
                    
                    % Rename variables to match desired outputs
                    data_fu_emu3 = rename_and_rearrange(data_fu_emu3, ...
                        varNames_out_allN) ;
                    
                    % Update names
                    [varNames_emu3, cropList_emu3, ...
                        varNames_emu3_basei, cropList_emu3_basei, ...
                        Nlist_emu3_char] = ...
                        e2p_get_names(data_fu_emu3.varNames, ...
                        []) ;
                    
                    % Save exclusion figures, if doing so
                    if save_excl_figs
                        disp('    Saving exclusion figures (output crops)...')
                        e2p_save_excl_figs_outCrops( ...
                            ggcm, which_file, gridlist, excl_vecs_emu3, ...
                            cropList_emu3, cropList_emu3_basei, ...
                            figure_visibility, figure_extension, ...
                            outDir_excl_figs_outCrops, ...
                            overwrite_existing_figs, renderer)
                    end

                    % Sort. Technically unnecessary, but maybe PLUM isn't
                    % robust to non-sorted crops.
                    data_fu_emu3 = do_sort(data_fu_emu3) ;
                                        
                    % Save MAT files
                    disp('    Saving MAT file(s)...')
                    % Save outputs (one file)
                    out_file = sprintf('%s/future_%s.mat', outDir_ggcm, which_file) ;
                    save(out_file, 'data_fu_emu3')
                    if ggcm_counter == 1
                        out_file = sprintf('%s/future_%s.mat', outDir_lpj, which_file) ;
                        save(out_file, 'data_fu_lpj2')
                    end
                    
                    % Consistency check
                    if isfield(data_fu_emu3, 'actually_emuBL_char') ...
                    && ~isequal(size(shiftdim(data_fu_emu3.varNames)), size(shiftdim(data_fu_emu3.actually_emuBL_char)))
                        error('data_fu_emu3.actually_emuBL_char must be the same size as data_fu_emu3.varNames')
                    end

                end % if did_load_existing
                
                disp('    Done.')

                % Save yield diagnostic figures, if doing so
                if save_out_figs
                    disp('    Preparing for output figures...')
                    
                    emu0s_combineCrops_file = strrep(emu_combineCrops_file, ...
                        '.mat', '.0s.mat') ;
                    if strcmp(which_file, 'yield')
                        [data_fu_emu0s, combineCrops_emu0s, ~, cf_emu0s] = ...
                            e2p_split_combine_burn(data_fu_emu0, ...
                            combineCrops, cropList_out, cf_emu0) ;
                        if ~isequal(cf_emu0s, cf_emu3)
                            error('isequal(cf_emu0s, cf_emu3)')
                        end
                        % For yield, save result combineCrops
                        save_combineCrops(combineCrops_emu0s, which_file, ...
                            emu0s_combineCrops_file)
                    elseif strcmp(which_file, 'gsirrigation')
                        data_fu_emu0s = e2p_apply_max(data_fu_emu0, ...
                            emu0s_combineCrops_file) ;
                    elseif ~strcmp(which_file, 'anpp')
                        error('which_file (%s) not recognized', which_file)
                    end
                    
                    disp('    Saving output figures...')
                    if strcmp(which_file, 'yield')
                        e2p_save_out_figs(data_fu_lpj2, ...
                            data_fu_emu0s, data_fu_emu3, ...
                            ggcm, outDir_yield_figs, ...
                            which_file, figure_visibility, figure_extension, ...
                            which_out_figs, overwrite_existing_figs, renderer, ...
                            cf_lpj, cf_emu3)
                    elseif strcmp(which_file, 'gsirrigation')
                        e2p_save_out_figs(data_fu_lpj2, ...
                            data_fu_emu0s, data_fu_emu3, ...
                            ggcm, outDir_irrig_figs, ...
                            which_file, figure_visibility, figure_extension, ...
                            which_out_figs, overwrite_existing_figs, renderer, ...
                            cf_lpj, cf_emu3)
                    else
                        warning('which_file (%s) not recognized for save_yield_figs; skipping.', which_file)
                    end
                    clear tmp
                end
                
                %% Save text file outputs (for PLUM)
                
                data_fu_emu3 = trim_unneeded_N200(data_fu_emu3, which_file) ;
                if save_txt_files_emu || save_txt_files_lpjg
                    
                    % Will we be saving LPJ-GUESS ANPP and runoff?
                    save_anpp_runoff = save_txt_files_lpjg ...
                        & strcmp(ggcm, ggcm_list{1}) ...
                        & strcmp(which_file, 'yield') ;
                    
                    disp('    Saving txt files for PLUM...')
                    out_header_cell = [{'Lon', 'Lat'} data_fu_emu3.varNames] ;
                    for t = 1:Ntpers
                        y1 = data_fu_emu3.y1s(t) ;
                        yN = data_fu_emu3.yNs(t) ;

                        if save_txt_files_lpjg && strcmp(ggcm, ggcm_list{1})
                            % Trim unneeded N200 from LPJ-GUESS
                            unneededN200 = getN_num(data_fu_lpj2.varNames)==200 ;
                            tmp_xv = data_fu_lpj2.garr_xvt(:,~unneededN200,t) ;
                            fprintf('    %d/%d lpj', t, Ntpers)
                            e2p_save(outDir_lpj, y1, yN, out_header_cell, ...
                                data_fu_lpj2.lonlats, tmp_xv, which_file, ...
                                overwrite_existing_txt)
                            
                            % Save ANPP and runoff
                            if save_anpp_runoff
                                fprintf('    %d/%d lpj ANPP', t, Ntpers)
                                e2p_save(outDir_lpj, y1, yN, out_header_cell_anpp, ...
                                    data_fu_lpj2.lonlats, data_fu_lpj_anpp.garr_xvt(:,:,t), ...
                                    'anpp', overwrite_existing_txt, lpj_anpp_precision)
                                fprintf('    %d/%d lpj runoff', t, Ntpers)
                                e2p_save(outDir_lpj, y1, yN, out_header_cell_runoff, ...
                                    data_fu_lpj2.lonlats, data_fu_lpj_runoff.garr_xvt(:,:,t), ...
                                    'tot_runoff', overwrite_existing_txt, lpj_runoff_precision)
                            end
                        end
                        if save_txt_files_emu
                            fprintf('    %d/%d emu', t, Ntpers)
                            e2p_save(outDir_ggcm, y1, yN, out_header_cell, ...
                                data_fu_emu3.lonlats, data_fu_emu3.garr_xvt(:,:,t), which_file, ...
                                overwrite_existing_txt)
                        end

                    end
                    clear data_fu_lpj1 data_fu_lpj_anpp data_fu_lpj_runoff
                end
                
                
                clear data_fu_lpj2
                fprintf('Done with %s %s %s %s.\n', gcm, ssp, ggcm, which_file)
                
            end % loop through whichfile_list
            
            fprintf('Done with %s %s %s.\n', gcm, ssp, ggcm)
            diary('off')
            
        end % Loop through crop models

        fprintf('Done with %s %s.\n', gcm, ssp)
        
        clear data_fu_lpj*_yield data_fu_lpj*_irrig

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


function save_combineCrops(combineCrops, which_file, out_file)

if strcmp(which_file, 'yield')
    save(out_file, 'combineCrops') ;
elseif ~contains(which_file, {'anpp', 'gsirrigation'})
    error('which_file (%s) not recognized', which_file)
end

end


function S = rename_and_rearrange(S, varNames_out)

has_time = ~isfield(S, 'garr_xv') ;
if has_time
    tmp_garr = S.garr_xvt ;
    S = rmfield(S, 'garr_xvt') ;
else
    tmp_garr = S.garr_xv ;
    S = rmfield(S, 'garr_xv') ;
end

[I_in, I_out] = e2p_translate_crops_agm2out(...
    S.varNames, varNames_out) ;
S.varNames = varNames_out ;
if isfield(S, 'actually_emuBL_char')
    S.actually_emuBL_char = S.actually_emuBL_char(I_in) ;
end

Ncells = size(tmp_garr, 1) ;
Nvar_out = length(S.varNames) ;
if has_time
    S.garr_xvt = nan(Ncells, Nvar_out, size(tmp_garr, 3)) ;
    S.garr_xvt(:,I_out,:) = tmp_garr(:,I_in,:) ;
else
    S.garr_xv = nan(Ncells, Nvar_out) ;
    S.garr_xv(:,I_out) = tmp_garr(:,I_in) ;
end
clear tmp_garr_xv

end


function S = trim_unneeded_N200(S, which_file)

unneededN200 = find(getN_num(S.varNames)==200) ;
S.garr_xvt(:,unneededN200,:) = [] ;
S.varNames(unneededN200) = [] ;
if isfield(S, 'actually_emuBL_char')
    S.actually_emuBL_char(unneededN200) = [] ;
end
e2p_check_correct_zeros(S.garr_xvt, ...
    which_file, S.varNames, ...
    'Future', @getbasenamei)

end


function save_cf_csv(cf_table, outDir_thisGGCM)

writetable(cf_table, ...
    sprintf('%s/calib_factors.csv', outDir_thisGGCM))

end


function precision = get_precision(A)

precision = 1 ;
while max(abs(minmax_ssr(A - round(A, precision)))) > 1e-11
    precision = precision + 1 ;
    if precision == 11
        error('Possible problem finding precision')
    end
end

end
