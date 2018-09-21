function [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
    inDir,base_year,yearList, ...
    landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
    is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, thisVer, ...
    read_MATs, is_orig)

% Does MAT-file already exist?
MATfile = [inDir '.processed.' thisVer 'mat'] ;
if false%exist(MATfile,'file')
    MATfile_info = dir(MATfile) ;
    [~, TXTfile] = unix(['ls -t ' inDir '/*/*.' thisVer 'txt | head -n 1  | tr -d ''\n''']) ;
    TXTfile_info = dir(TXTfile) ;
    if MATfile_info.datenum > TXTfile_info.datenum
        S_out = [] ;
        S_nfert_out = [] ;
        S_irrig_out = [] ;
        disp('   Loading MAT file...')
        load(MATfile,'S_out','S_nfert_out','S_irrig_out') ;
    else
        warning('.mat file exists but is older than latest .txt file. Regenerating.')
        [S_out, S_nfert_out, S_irrig_out] = ...
            generate_struct( ...
                inDir,base_year,yearList, ...
                landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
                is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, read_MATs) ;
        % Save, catching exceptions to avoid having to re-read everything
        disp('   Saving MAT file...')
        try
            save(MATfile,'S_out', 'S_nfert_out', 'S_irrig_out') ;
        catch ME
            keyboard
        end
    end
else
    [S_out, S_nfert_out, S_irrig_out] = ...
        generate_struct( ...
            inDir,base_year,yearList, ...
            landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
            is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, ...
            read_MATs, is_orig) ;
%     % Save, catching exceptions to avoid having to re-read everything
%     disp('   Saving MAT file...')
%     try
%         save(MATfile,'S_out', 'S_nfert_out', 'S_irrig_out') ;
%     catch ME
%         keyboard
%     end
end

disp('Done.')


end



function [S, S_nfert, S_irrig] = generate_struct(...
    inDir, base_year, yearList, ...
    landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
    is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, ...
    read_MATs, is_orig)

S.varNames = LUnames ;
S_nfert.varNames = LPJGcrops ;
S_irrig.varNames = LPJGcrops ;

S.yearList = yearList ;
S_nfert.yearList = yearList ;
S_irrig.yearList = yearList ;
Nyears = length(yearList) ;

landArea_2deg_YX = PLUMharm_aggregate(landArea_YX,0.5,2) ;

cf_kgNha_kgNm2 = 1e-4 ;

for y = 1:Nyears
    
    thisYear = yearList(y) ;
    if read_MATs
        file_in_lu = [inDir '/' num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.mat'] ;
        if is2deg
            file_in_lu = strrep(file_in_lu,'.mat','.2deg.mat') ;
        end
    else
        if is_orig
            file_in_lu = [inDir '/' num2str(thisYear) '/LandCoverFract.txt'] ;
        else
            file_in_lu = [inDir '/' num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.txt'] ;
        end
    end
    
    if exist(file_in_lu,'file')
        if y>1 && rem(y,5)==1
            fprintf('\n   %d... ',thisYear)
        elseif y==1
            fprintf('   %d... ',thisYear)
        elseif y==Nyears
            fprintf('%d...\n',thisYear)
        else
            fprintf('%d... ',thisYear)
        end
        if read_MATs % Each loads as structure out_y1
            if is2deg
                thisLandArea_YX = landArea_2deg_YX ;
            else
                thisLandArea_YX = landArea_YX ;
            end
            
            % Import land use areas
            tmp = load(file_in_lu) ;
            LU_tmp_in = tmp.out_y1 ; clear tmp
            LUarea_tmp_in = LU_tmp_in ;
            LUarea_tmp_in.maps_YXv = LU_tmp_in.maps_YXv ...
                .* repmat(thisLandArea_YX,[1 1 size(LU_tmp_in.maps_YXv,3)]) ;
            
            % Import crop areas
            file_in_cf = strrep(file_in_lu,'LandCoverFract','CropFract') ;
            tmp = load(file_in_cf) ;
            cropFrac_tmp_in = tmp.out_y1 ; clear tmp
            cropArea_tmp_in = cropFrac_tmp_in ;
            Ncrops = size(cropFrac_tmp_in.maps_YXv,3) ;
            cropArea_tmp_in.maps_YXv = cropFrac_tmp_in.maps_YXv .* repmat(LUarea_tmp_in.maps_YXv(:,:,strcmp(LUarea_tmp_in.varNames,'CROPLAND')), [1 1 Ncrops]) ;
            
            % Combine areas
            PLUM_in.varNames = [cropArea_tmp_in.varNames LUarea_tmp_in.varNames(~strcmp(LUarea_tmp_in.varNames,'CROPLAND'))] ;
            PLUM_in.maps_YXv = cat(3, cropArea_tmp_in.maps_YXv, LUarea_tmp_in.maps_YXv(:,:,~strcmp(LUarea_tmp_in.varNames,'CROPLAND'))) ;
            
            % Import managements
            file_in_nfert = strrep(file_in_lu,'LandCoverFract','Fert') ;
            tmp = load(file_in_nfert) ;
            nfert_in = tmp.out_y1 ; clear tmp
            file_in_irrig = strrep(file_in_lu,'LandCoverFract','Irrig') ;
            tmp = load(file_in_irrig) ;
            irrig_in = tmp.out_y1 ; clear tmp
            
            % Convert kgN/ha to kgN/m2
            nfert_in.maps_YXv = nfert_in.maps_YXv * cf_kgNha_kgNm2 ;
        else
            if is2deg
                [~, ~, ~, PLUM_in, nfert_in, irrig_in] = PLUMharm_processPLUMin_areaCrops(...
                    file_in_lu,landArea_YX, landArea_2deg_YX, LUnames, bareFrac_y0_YX, ...
                    [], [], PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method) ;
            else
                [PLUM_in, nfert_in, irrig_in, ~, ~, ~] = PLUMharm_processPLUMin_areaCrops(...
                    file_in_lu,landArea_YX, landArea_2deg_YX, LUnames, bareFrac_y0_YX, ...
                    [], [], PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method) ;
            end
        end
        if y==1
            S.maps_YXvy = nan(size(PLUM_in.maps_YXv,1), ...
                              size(PLUM_in.maps_YXv,2), ...
                              length(LUnames), Nyears, 'single') ;
            S_nfert.maps_YXvy = nan(size(nfert_in.maps_YXv,1), ...
                                    size(nfert_in.maps_YXv,2), ...
                                    length(LPJGcrops), Nyears, 'single') ;
            S_irrig.maps_YXvy = nan(size(nfert_in.maps_YXv,1), ...
                                    size(nfert_in.maps_YXv,2), ...
                                    length(LPJGcrops), Nyears, 'single') ;
        end
        S.maps_YXvy(:,:,:,y) = PLUM_in.maps_YXv ;
        S_nfert.maps_YXvy(:,:,:,y) = nfert_in.maps_YXv ;
        S_irrig.maps_YXvy(:,:,:,y) = irrig_in.maps_YXv ;
    else
        warning('\n%d does not exist! Aborting.\n', thisYear)
        break
    end

end

if ~isfield(S,'maps_YXvy')
    error('No files read!')
end


end