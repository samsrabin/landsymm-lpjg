function [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
    inDir,base_year,yearList, ...
    landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
    is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, thisVer, ...
    is_orig, fruitveg_sugar_2oil, allow_unveg)

combineCrops = isempty(PLUMtoLPJG) ;

S_nfert_out = [] ;
S_irrig_out = [] ;

% Does MAT-file already exist?
MATfile = [inDir '.processed.' thisVer 'mat'] ;
MATfile = sprintf('%s.processed.%d-%d.%smat', inDir, yearList(1), yearList(end), thisVer) ;
disp(MATfile)
if exist(MATfile,'file')
    MATfile_info = dir(MATfile) ;
    [~, TXTfile] = unix(['ls -t ' inDir '/*/*.' thisVer 'txt | head -n 1  | tr -d ''\n''']) ;
    TXTfile_info = dir(TXTfile) ;
    
    if isempty(TXTfile_info) || MATfile_info.datenum > TXTfile_info.datenum
        disp('   Loading MAT file...')
        if combineCrops
            load(MATfile,'S_out') ;
        else
            load(MATfile,'S_out','S_nfert_out','S_irrig_out') ;
        end
    else
        warning('.mat file exists but is older than latest .txt file. Regenerating.')
        if combineCrops
            [S_out, ~, ~] = ...
                generate_struct( ...
                inDir,base_year,yearList, ...
                landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
                is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, ...
                fruitveg_sugar_2oil) ;
        else
            [S_out, S_nfert_out, S_irrig_out] = ...
            generate_struct( ...
                inDir,base_year,yearList, ...
                landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
                is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, ...
                fruitveg_sugar_2oil) ;
        end
        % Save, catching exceptions to avoid having to re-read everything
        disp('   Saving MAT file...')
        try
            if combineCrops
                save(MATfile,'S_out') ;
            elseif ~combineCrops
                save(MATfile,'S_out', 'S_nfert_out', 'S_irrig_out') ;
            end
        catch ME
            if batchStartupOptionUsed
                rethrow(ME)
            else 
                keyboard
            end
        end
    end
elseif combineCrops
    [S_out, ~, ~] = ...
        generate_struct( ...
            inDir,base_year,yearList, ...
            landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
            is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, ...
            is_orig, fruitveg_sugar_2oil) ;
    % Save, catching exceptions to avoid having to re-read everything
    disp('   Saving MAT file...')
    try
        save(MATfile,'S_out') ;
    catch ME
        if batchStartupOptionUsed
            rethrow(ME)
        else 
            keyboard
        end
    end
else
    [S_out, S_nfert_out, S_irrig_out] = ...
        generate_struct( ...
            inDir,base_year,yearList, ...
            landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
            is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, ...
            is_orig, fruitveg_sugar_2oil, allow_unveg) ;
    % Save, catching exceptions to avoid having to re-read everything
    disp('   Saving MAT file...')
    try
        save(MATfile,'S_out', 'S_nfert_out', 'S_irrig_out') ;
    catch ME
        if batchStartupOptionUsed
            rethrow(ME)
        else 
            keyboard
        end
    end
end

disp('Done.')


end



function [S, S_nfert, S_irrig] = generate_struct(...
    inDir, base_year, yearList, ...
    landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
    is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, ...
    is_orig, fruitveg_sugar_2oil, allow_unveg)

combineCrops = isempty(PLUMtoLPJG) ;
S_nfert = [] ;
S_irrig = [] ;

S.varNames = LUnames ;
if ~combineCrops
    S_nfert.varNames = LPJGcrops ;
    S_irrig.varNames = LPJGcrops ;
end

S.yearList = yearList ;
if ~combineCrops
    S_nfert.yearList = yearList ;
    S_irrig.yearList = yearList ;
end
Nyears = length(yearList) ;

if ~is2deg
    landArea_2deg_YX = PLUMharm_aggregate(landArea_YX,0.5,2) ;
else
    landArea_2deg_YX = landArea_YX ;
end

cf_kgNha_kgNm2 = 1e-4 ;

for y = 1:Nyears
    
    thisYear = yearList(y) ;
    clear *_exists
    txt_exists = false ;
    gz_exists = false ;
    
    file_in_lu = fullfile(inDir, num2str(thisYear), 'LandCoverFract.mat') ;
    if ~is_orig
        file_in_lu = strrep(file_in_lu, '.mat', sprintf('.base%d.mat', base_year)) ;
    end
    mat_exists = exist(file_in_lu, 'file') ;
    if ~mat_exists
        file_in_lu = strrep(file_in_lu, '.mat', '.txt') ;
        txt_exists = exist(file_in_lu, 'file') ;
        if ~txt_exists
            gz_exists = exist([file_in_lu '.gz'], 'file') ;
        end
    end
    
    if ~(mat_exists || txt_exists || gz_exists)
        warning('\n%d (%s) does not exist! Aborting.\n', thisYear, ...
            strrep(file_in_lu, '.txt.gz', '[.mat,.txt(.gz)]'))
        break
    end
    
    if y>1 && rem(y,5)==1
        fprintf('\n   %d... ',thisYear)
    elseif y==1
        fprintf('   %d... ',thisYear)
    elseif y==Nyears
        fprintf('%d...\n',thisYear)
    else
        fprintf('%d... ',thisYear)
    end
    if mat_exists % Each loads as structure out_y1
        
        % Import land use areas
        tmp = load(file_in_lu) ;
        LU_tmp_in = tmp.out_y1 ; clear tmp
        LUarea_tmp_in = LU_tmp_in ;
        if is2deg
            error('You need to pass in half-degree land area for next step...')
        end
        LUarea_tmp_in.maps_YXv = LU_tmp_in.maps_YXv ...
            .* repmat(landArea_YX,[1 1 size(LU_tmp_in.maps_YXv,3)]) ;
        if is2deg
            LUarea_tmp_in.maps_YXv = PLUMharm_aggregate( ...
                LUarea_tmp_in.maps_YXv, 0.5, 2) ;
        end
        
        % Import crop areas
        if combineCrops
            cropFrac_tmp_in.maps_YXv = ...
                ones(size(LUarea_tmp_in.maps_YXv, 1:2)) ;
%                 LUarea_tmp_in.maps_YXv(:,:,strcmp(LUarea_tmp_in.varNames,'CROPLAND')) ...
%                 ./ repmat(thisLandArea_YX,[1 1 size(LU_tmp_in.maps_YXv,3)]) ;
            cropFrac_tmp_in.varNames = LPJGcrops ;
        else
            file_in_cf = strrep(file_in_lu,'LandCoverFract','CropFract') ;
            tmp = load(file_in_cf) ;
            cropFrac_tmp_in = tmp.out_y1 ; clear tmp
        end
        cropArea_tmp_in = cropFrac_tmp_in ;
        Ncrops = size(cropFrac_tmp_in.maps_YXv,3) ;
        cropArea_tmp_in.maps_YXv = cropFrac_tmp_in.maps_YXv .* repmat(LUarea_tmp_in.maps_YXv(:,:,strcmp(LUarea_tmp_in.varNames,'CROPLAND')), [1 1 Ncrops]) ;
        
        % Combine areas
        PLUM_in.varNames = [cropArea_tmp_in.varNames LUarea_tmp_in.varNames(~strcmp(LUarea_tmp_in.varNames,'CROPLAND'))] ;
        PLUM_in.maps_YXv = cat(3, cropArea_tmp_in.maps_YXv, LUarea_tmp_in.maps_YXv(:,:,~strcmp(LUarea_tmp_in.varNames,'CROPLAND'))) ;
        
        if ~combineCrops
            % Import managements
            file_in_nfert = strrep(file_in_lu,'LandCoverFract','Fert') ;
            tmp = load(file_in_nfert) ;
            nfert_in = tmp.out_y1 ; clear tmp
            file_in_irrig = strrep(file_in_lu,'LandCoverFract','Irrig') ;
            tmp = load(file_in_irrig) ;
            irrig_in = tmp.out_y1 ; clear tmp
            
            % Convert kgN/ha to kgN/m2
            nfert_in.maps_YXv = nfert_in.maps_YXv * cf_kgNha_kgNm2 ;
        end
    else
        if is2deg
            try
                [~, ~, ~, PLUM_in, nfert_in, irrig_in] = PLUMharm_processPLUMin_areaCrops(...
                    file_in_lu,landArea_YX, landArea_2deg_YX, LUnames, bareFrac_y0_YX, ...
                    [], [], PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method, fruitveg_sugar_2oil, ...
                    allow_unveg) ;
            catch ME
                if batchStartupOptionUsed
                    rethrow(ME)
                else 
                    keyboard
                end
            end
        else
            [PLUM_in, nfert_in, irrig_in, ~, ~, ~] = PLUMharm_processPLUMin_areaCrops(...
                file_in_lu,landArea_YX, landArea_2deg_YX, LUnames, bareFrac_y0_YX, ...
                [], [], PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method, fruitveg_sugar_2oil, ...
                allow_unveg) ;
        end
    end
    if y==1
        S.maps_YXvy = nan(size(PLUM_in.maps_YXv,1), ...
            size(PLUM_in.maps_YXv,2), ...
            length(LUnames), Nyears, 'single') ;
        if ~combineCrops
            S_nfert.maps_YXvy = nan(size(nfert_in.maps_YXv,1), ...
                size(nfert_in.maps_YXv,2), ...
                length(LPJGcrops), Nyears, 'single') ;
            S_irrig.maps_YXvy = nan(size(nfert_in.maps_YXv,1), ...
                size(nfert_in.maps_YXv,2), ...
                length(LPJGcrops), Nyears, 'single') ;
        end
    end
    S.maps_YXvy(:,:,:,y) = PLUM_in.maps_YXv ;
    if ~combineCrops
        S_nfert.maps_YXvy(:,:,:,y) = nfert_in.maps_YXv ;
        S_irrig.maps_YXvy(:,:,:,y) = irrig_in.maps_YXv ;
    end

end

if ~isfield(S,'maps_YXvy')
    error('No files read!')
end


end
