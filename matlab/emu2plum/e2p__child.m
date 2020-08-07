for g = 1:length(gcm_list)
    gcm = gcm_list{g} ;
    
    for s = 1:length(ssp_list)
        ssp = ssp_list{s} ;

        outDir = sprintf('%s_forPLUM/%s_%s_v%s', ...
            topDir_emu, gcm, ssp, thisVer) ;
        if remove_outliers
            outDir = [outDir '_rmol'] ; %#ok<AGROW>
        end
        outDir_lpj = sprintf('%s/sim_LPJ-GUESS', outDir) ;
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

            %% Setup
            ggcm = ggcm_list{ggcm_counter} ;
            fprintf('%s %s %s\n', gcm, ssp, ggcm)
            topDir_phase2 = sprintf('%s/AgMIP.output/%s/phase2', topdir_sh, ggcm) ;
            outDir_ggcm = sprintf('%s/emu_%s', outDir, ggcm) ;
            if ~exist(outDir_ggcm, 'dir')
                mkdir(outDir_ggcm) ;
            end

            % If allowed to import existing file, then do so
            out_file = sprintf('%s/future_%s.mat', outDir_ggcm, which_file) ;
            did_load_existing = load_existing_file && exist(out_file, 'file') ;
            if did_load_existing
                fprintf('    Importing future_%s.mat...\n', which_file) ;

                load(out_file, 'data_fu_lpj', 'data_fu_emu', 'data_fu_out')

                % Translate crop names
                [varNames_emu, cropList_emu, ...
                        varNames_emu_basei, cropList_emu_basei, ...
                        Nlist_emu, ~] = ...
                        e2p_get_names([], data_fu_emu.varNames, ...
                        getbasename, getbasenamei, getN) ;
                [cropList_lpj_asEmu, used_emuCrops] = e2p_translate_crops( ...
                    cropList_lpj, cropList_emu) ;
            else

                %% Import emulator outputs

                if contains(which_file, {'yield', 'gsirrigation'})

                    disp('    Importing emulator outputs...')

                    try
                        [data_bl_emu, data_fu_emu] = e2p_import_emu( ...
                            topDir_emu, gcm, ggcm, ssp, which_file, ...
                            cropList_in, gridlist, ts1_list, tsN_list, Nlist, ...
                            baseline_yN, future_yN_emu, irrList_in, irrList_out) ;
                    catch ME
                        if strcmp(ME.identifier, 'e2p:e2p_import_emu:noFilesFound')
                            warning('No files found for %s %s %s; skipping', ...
                                gcm, ssp, ggcm)
                            continue
                        else
                            rethrow(ME)
                        end
                    end

                    e2p_check_correct_zeros(data_bl_emu.garr_xv, which_file, getbasenamei(data_bl_emu.varNames))
                    e2p_check_correct_zeros(data_fu_emu.garr_xvt, which_file, getbasenamei(data_fu_emu.varNames))

                    [varNames_emu, cropList_emu, ...
                        varNames_emu_basei, cropList_emu_basei, ...
                        Nlist_emu, ~] = ...
                        e2p_get_names(data_bl_emu.varNames, data_fu_emu.varNames, ...
                        getbasename, getbasenamei, getN) ;

                    disp('    Done.')

                elseif ~strcmp(which_file, 'anpp')
                    error('which_file (%s) not recognized', which_file)
                end


                %% Get and apply exclusions, if doing so

                if contains(which_file, {'yield','gsirrigation'}) && (excl_lowBL_agmerra || excl_lowBL_emu)
                    if strcmp(which_file, 'yield')

                        % Set up exclusion array
                        exclude_xc = false(length(data_bl_emu.list2map),length(cropList_emu_basei)) ;

                        % Where do we exclude based on low AgMERRA yield at max N
                        % (or existing exclusions)?
                        if excl_lowBL_agmerra
                            disp('    Excluding based on low AgMERRA yield...')
                            %% Get AgMERRA yield
                            yield_agmerraBL_xv = e2p_get_agmerra_yield(...
                                varNames_emu, topDir_phase2, ggcm, data_bl_emu.list2map, getbasenamei, getN) ;
                            if ~any(any(~isnan(yield_agmerraBL_xv)))
                                error('yield_agmerraBL_xv is all NaN')
                            end
                            % Exclude
                            exclude_xc = e2p_exclude_lowBLyield_atMaxN( ...
                                varNames_emu, cropList_emu_basei, Nlist_emu, ...
                                yield_agmerraBL_xv, 0.01, exclude_xc) ;
                        end

                        if ~any(any(~exclude_xc))
                            error('All cells excluded because of low AgMERRA yield')
                        end

                        % Where do we exclude based on low baseline-year emulated yield at
                        % max N (or existing exclusions)?
                        if excl_lowBL_emu
                            disp('    Excluding based on low baseline-year emulated yield...')
                            exclude_xc = e2p_exclude_lowBLyield_atMaxN( ...
                                varNames_emu, cropList_emu_basei, Nlist_emu, ...
                                data_bl_emu.garr_xv, 0.01, exclude_xc) ;
                        end

                        if ~any(any(~exclude_xc))
                            error('All cells excluded because of low AgMERRA yield and/or low baseline-year emulated yield')
                        end

                        % Error checks
                    elseif strcmp(which_file, 'gsirrigation')
                        disp('    Applying yield-based exclusions...')
                        exclude_xc_file = sprintf('%s/exclude_xc.mat', outDir_ggcm) ;
                        load(exclude_xc_file) ;
                    else
                        error('which_file (%s) not recognized', which_file)
                    end

                    for c = 1:length(cropList_emu_basei)
                        thisCrop = cropList_emu_basei{c} ;
                        thisCrop_i = find(strcmp(varNames_emu_basei,thisCrop)) ;
                        if length(thisCrop_i)~=length(Nlist_emu)
                            error('Error finding isThisCrop (%d found)', length(thisCrop_i))
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

                elseif ~strcmp(which_file,'anpp')
                    error('which_file (%s) not recognized', which_file)
                end

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
                        error('Mismatch between winter wheat lists')
                    end
                elseif strcmp(which_file, 'gsirrigation')
                    data_bl_emu = e2p_apply_max_wheat(data_bl_emu, outDir) ;
                    data_fu_emu = e2p_apply_max_wheat(data_fu_emu, outDir) ;
                elseif ~strcmp(which_file, 'anpp')
                    error('which_file (%s) not recognized', which_file)
                else
                    refresh_vars = false ;
                end

                e2p_check_correct_zeros(data_bl_emu.garr_xv, which_file, getbasenamei(data_bl_emu.varNames))
                e2p_check_correct_zeros(data_fu_emu.garr_xvt, which_file, getbasenamei(data_fu_emu.varNames))

                % Refresh variable and crop lists
                if refresh_vars
                    [varNames_emu, cropList_emu, ...
                        varNames_emu_basei, cropList_emu_basei, ...
                        Nlist_emu, ~] = ...
                        e2p_get_names(data_bl_emu.varNames, data_fu_emu.varNames, ...
                        getbasename, getbasenamei, getN) ; %#ok<ASGLU>
                end

                disp('    Done.')


                %% Harmonize LPJ-GUESS and emulator data

                disp('    Harmonizing LPJ-GUESS and emulator data...')

                % Align gridlists
                [data_bl_lpj, data_fu_lpj, data_bl_emu, data_fu_emu, list2map] = ...
                    e2p_align_gridlists(data_bl_lpj, data_fu_lpj, data_bl_emu, data_fu_emu) ;

                e2p_check_correct_zeros(data_bl_lpj.garr_xv, which_file, getbasenamei(data_bl_lpj.varNames))
                e2p_check_correct_zeros(data_fu_lpj.garr_xvt, which_file, getbasenamei(data_fu_lpj.varNames))
                e2p_check_correct_zeros(data_bl_emu.garr_xv, which_file, getbasenamei(data_bl_emu.varNames))
                e2p_check_correct_zeros(data_fu_emu.garr_xvt, which_file, getbasenamei(data_fu_emu.varNames))

                % Make sure that N lists match
                if ~isequal(Nlist_lpj,Nlist_emu)
                    error('Mismatch in N levels between LPJ-GUESS and emulator')
                end

                % Translate crop names
                [cropList_lpj_asEmu, used_emuCrops] = e2p_translate_crops( ...
                    cropList_lpj, cropList_emu) ;

                disp('    Done.')


                %% Get and apply deltas (and remove outliers, if doing so)

                if contains(which_file, {'yield','gsirrigation'})
                    disp('    Getting deltas...')

                    deltas_emu_xvt = e2p_get_deltas(...
                        data_bl_emu, data_fu_emu, interp_infs, cropList_emu, ...
                        getbasename, getbasenamei, which_file, ...
                        used_emuCrops, list2map, ...
                        save_interp_figs, outDir_interp_figs, ggcm, figure_visibility) ;
                    e2p_check_correct_zeros(deltas_emu_xvt, which_file, getbasenamei(data_fu_emu.varNames))
                    disp('    Done.')

                elseif ~strcmp(which_file, 'anpp')
                    error('which_file (%s) not recognized', which_file)
                end

                disp('    Applying deltas...')

                % Applies emulator deltas to LPJ-GUESS baseline.
                % data_fu_lpj is not affected, except to have its variables sorted for
                % consistency.
                [data_fu_lpj, data_fu_out] = e2p_apply_deltas( ...
                    data_bl_lpj, data_fu_lpj, data_bl_emu, data_fu_emu, deltas_emu_xvt, ...
                    cropList_lpj, cropList_lpj_asEmu, varNames_lpj, ...
                    list2map, getbasename, getbasenamei, which_file, figure_visibility) ;
                e2p_check_correct_zeros(data_fu_out.garr_xvt, which_file, getbasenamei(data_fu_out.varNames))

                if remove_outliers
                    disp('    Removing outliers...')

                    if strcmp(which_file,'yield')
                        smad_mult = 3 ;
                    elseif strcmp(which_file,'gsirrigation')
                        warning('Might have to change gsirrigation smad_mult once properly pre-thresholding')
                        smad_mult = 3 ;
                    else
                        error('which_file (%s) not recognized', which_file)
                    end

                    [data_fu_lpj, outlier_info_lpj] = e2p_remove_outliers(data_fu_lpj, smad_mult) ;
                    e2p_check_correct_zeros(data_fu_lpj.garr_xvt, which_file, getbasenamei(data_fu_lpj.varNames))

                    [data_fu_out, outlier_info_out] = e2p_remove_outliers(data_fu_out, smad_mult) ;
                    e2p_check_correct_zeros(data_fu_out.garr_xvt, which_file, getbasenamei(data_fu_out.varNames))

                    %% Save info
                    outlier_info_cols = string([ ...
                        repmat('y',[Ntpers 1]) ...
                        num2str(shiftdim(data_fu_out.y1s)) ...
                        repmat('_',[Ntpers 1]) ...
                        num2str(shiftdim(data_fu_out.yNs))]) ;
                    e2p_save_outlier_info(outlier_info_lpj, outDir_lpj, which_file, outlier_info_cols)
                    e2p_save_outlier_info(outlier_info_out, outDir_ggcm, which_file, outlier_info_cols)
                    clear outlier_info_cols
                end

            end

            disp('    Done.')


            %% Save outputs

            if ~did_load_existing
                disp('    Saving...')

                % Save outputs (one file)
                out_file = sprintf('%s/future_%s.mat', outDir_ggcm, which_file) ;
                save(out_file, 'data_fu_lpj', 'data_fu_emu', 'data_fu_out')

                % Save exclusion info
                out_file = sprintf('%s/exclude_xc.mat', outDir_ggcm) ;
                save(out_file, 'exclude_xc')

                % Save outputs (for PLUM)
                out_header_cell = [{'Lon', 'Lat'} data_fu_out.varNames] ;
                for t = 1:Ntpers
                    fprintf('    %d/%d...\n', t, Ntpers)
                    y1 = data_fu_out.y1s(t) ;
                    yN = data_fu_out.yNs(t) ;

                    if strcmp(ggcm, ggcm_list{1})
                        e2p_save(outDir_lpj, y1, yN, out_header_cell, ...
                            data_fu_lpj.lonlats, data_fu_lpj.garr_xvt(:,:,t), which_file, ...
                            interp_infs, remove_outliers, false)
                    end
                    e2p_save(outDir_ggcm, y1, yN, out_header_cell, ...
                        data_fu_out.lonlats, data_fu_out.garr_xvt(:,:,t), which_file, ...
                        interp_infs, remove_outliers, true)

                end

                % For yield, save is_ww_max_bl_gW*
                if strcmp(which_file, 'yield')
                    out_file = sprintf('%s/is_ww_max.mat', outDir_ggcm) ;
                    save(out_file, 'is_ww_max_bl_gW', 'is_ww_max_fu_gWt', 'winter_wheats') ;
                elseif ~contains(which_file, {'anpp', 'gsirrigation'})
                    error('which_file (%s) not recognized', which_file)
                end

            end

            % Save yield diagnostic figures, if doing so
            if save_out_figs
                disp('    Saving output figures...')
                if strcmp(which_file, 'yield')
                    e2p_save_out_figs(data_fu_lpj, data_fu_emu, data_fu_out, ...
                        ggcm, getbasename, getbasenamei, getN, outDir_yield_figs, ...
                        which_file, cropList_lpj_asEmu, figure_visibility)
                elseif strcmp(which_file, 'gsirrigation')
                    e2p_save_out_figs(data_fu_lpj, data_fu_emu, data_fu_out, ...
                        ggcm, getbasename, getbasenamei, getN, outDir_irrig_figs, ...
                        which_file, cropList_lpj_asEmu, figure_visibility)
                else
                    warning('which_file (%s) not recognized for save_yield_figs; skipping.', which_file)
                end
            end



            fprintf('Done with %s %s %s.\n', gcm, ssp, ggcm)

        end

        fprintf('Done with %s %s.\n', gcm, ssp)

    end
    
    fprintf('Done with %s.\n', gcm)
end

disp(' ')
disp('All done!')


