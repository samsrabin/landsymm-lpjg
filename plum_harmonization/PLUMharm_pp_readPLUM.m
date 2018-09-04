function [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
    inDir,base_year,yearList, ...
    landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
    is2deg, bareFrac_y0_YX, norm2extra, inpaint_method, thisVer)

% Does MAT-file already exist?
MATfile = [inDir '.processed.' thisVer 'mat'] ;
if exist(MATfile,'file')
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
                is2deg, bareFrac_y0_YX, norm2extra, inpaint_method) ;
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
            is2deg, bareFrac_y0_YX, norm2extra, inpaint_method) ;
    % Save, catching exceptions to avoid having to re-read everything
    disp('   Saving MAT file...')
    try
        save(MATfile,'S_out', 'S_nfert_out', 'S_irrig_out') ;
    catch ME
        keyboard
    end
end

disp('Done.')


end



function [S, S_nfert, S_irrig] = generate_struct(...
    inDir, base_year, yearList, ...
    landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
    is2deg, bareFrac_y0_YX, norm2extra, inpaint_method)

S.varNames = LUnames ;
S_nfert.varNames = LPJGcrops ;
S_irrig.varNames = LPJGcrops ;
Nyears = length(yearList) ;

is_orig = ~isempty(bareFrac_y0_YX) ;
landArea_2deg_YX = PLUMharm_aggregate(landArea_YX,0.5,2) ;

for y = 1:Nyears
    
    thisYear = yearList(y) ;
    if is_orig
        file_in = [inDir '/' num2str(thisYear) '/LandCoverFract.txt'] ;
    else
        file_in = [inDir '/' num2str(thisYear) '/LandCoverFract.base' num2str(base_year) '.txt'] ;
    end
    
    if exist(file_in,'file')
        if y>1 && rem(y,5)==1
            fprintf('\n   %d... ',thisYear)
        elseif y==1
            fprintf('   %d... ',thisYear)
        elseif y==Nyears
            fprintf('%d...\n',thisYear)
        else
            fprintf('%d... ',thisYear)
        end
        if is2deg
            [~, ~, ~, PLUM_in, nfert_in, irrig_in] = PLUMharm_processPLUMin_areaCrops(...
                file_in,landArea_YX, landArea_2deg_YX, LUnames, bareFrac_y0_YX, ...
                [], [], PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method) ;
        else
            [PLUM_in, nfert_in, irrig_in, ~, ~, ~] = PLUMharm_processPLUMin_areaCrops(...
                file_in,landArea_YX, landArea_2deg_YX, LUnames, bareFrac_y0_YX, ...
                [], [], PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method) ;
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