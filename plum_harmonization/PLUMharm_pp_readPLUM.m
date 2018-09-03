function S_out = PLUMharm_pp_readPLUM(inDir,base_year,yearList, ...
    landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
    is2deg, bareFrac_y0_YX, norm2extra)

is_orig = ~isempty(bareFrac_y0_YX) ;

% Does MAT-file already exist?
MATfile = [inDir '.processed.mat'] ;
if exist(MATfile,'file')
    MATfile_info = dir(MATfile) ;
    [~, TXTfile] = unix(['ls -t ' inDir '/*/*.txt | head -n 1  | tr -d ''\n''']) ;
    TXTfile_info = dir(TXTfile) ;
    if MATfile_info.datenum > TXTfile_info.datenum
        S_out = [] ;
        disp('Loading MAT file...')
        load(MATfile,'S_out') ;
    else
        warning('.mat file exists but is older than latest .txt file. Regenerating.')
        S_out = generate_struct(inDir,base_year,yearList, ...
            landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
            is2deg, bareFrac_y0_YX, norm2extra) ;
        disp('Saving MAT file...')
        try
            save(MATfile,'S_out') ;
        catch ME
            keyboard
        end
    end
else
    S_out = generate_struct(inDir,base_year,yearList, ...
        landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
        is2deg, bareFrac_y0_YX, norm2extra) ;
    disp('Saving MAT file...')
    try
        save(MATfile,'S_out') ;
    catch ME
        keyboard
    end
end

disp('Done.')


end



function S = generate_struct(inDir, base_year, yearList, ...
    landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
    is2deg, bareFrac_y0_YX, norm2extra)

S.varNames = LUnames ;
Nyears = length(yearList) ;

is_orig = ~isempty(bareFrac_y0_YX) ;

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
            [~,PLUM_in] = PLUMharm_processPLUMin_areaCrops(...
                file_in,landArea_YX, LUnames, bareFrac_y0_YX, ...
                PLUMtoLPJG, LPJGcrops, norm2extra) ;
        else
            [PLUM_in,~] = PLUMharm_processPLUMin_areaCrops(...
                file_in,landArea_YX, LUnames, bareFrac_y0_YX, ...
                PLUMtoLPJG, LPJGcrops, norm2extra) ;
        end
        if y==1
            S.maps_YXvy = nan(size(PLUM_in.maps_YXv,1), ...
                              size(PLUM_in.maps_YXv,2), ...
                              length(LUnames), Nyears, 'single') ;
        end
        S.maps_YXvy(:,:,:,y) = PLUM_in.maps_YXv ;
    else
        warning('\n%d does not exist! Aborting.\n', thisYear)
        break
    end

end

if ~isfield(S,'maps_YXvy')
    error('No files read!')
end


end