function hyde_this_yx = get_HYDE_LUs( ...
    hyde_years, yearList, inDir_hyde, hyde_files, ...
    thisName, i_ok)

Ncells = length(i_ok) ;
Nyears = length(yearList) ;
hyde_this_yx = nan(Nyears,Ncells) ;

for y = 1:Nyears
    
    tic
    thisYear = yearList(y) ;
    fprintf('%d:\n',thisYear) ;
    
    % Get HYDE land use maps (km2), linearly interpolating
    % between available years
    if y==1
        i = find(hyde_years==thisYear) ;
        if isempty(i)
            error('Deal with exact start year not being in HYDE')
        end
        year1 = yearList(1) ;
        fprintf('   Reading HYDE %d...\n',year1) ;
        filename_zip = sprintf('%s/%s',inDir_hyde,hyde_files(i).name) ;
        
        hyde_this_1_YX = read_HYDE_file(thisName,year1,filename_zip) ;
        
        hyde_this_1_x = hyde_this_1_YX(i_ok) ;
        if any(hyde_this_1_x<0)
            error('Failed sanity test for hyde_this_1_x!')
        end
        
    elseif any(hyde_years==thisYear)
        hyde_this_1_x = hyde_this_N_x ;
        year1 = yearN ;
    end
    if y==1 || any(hyde_years==thisYear)
        i = find(hyde_years==thisYear) + 1 ;
        yearN = hyde_years(i) ;
        fprintf('   Reading HYDE %d...\n',yearN) ;
        filename_zip = sprintf('%s/%s',inDir_hyde,hyde_files(i).name) ;
        
        hyde_this_N_YX = read_HYDE_file(thisName,yearN,filename_zip) ;
        
        hyde_this_N_x = hyde_this_N_YX(i_ok) ;
        if any(hyde_this_N_x<0)
            error('Failed sanity test for hyde_this_N_x!')
        end
        
    end
    hyde_this_x = (yearN-thisYear)/(yearN-year1) * hyde_this_1_x ...
                + (thisYear-year1)/(yearN-year1) * hyde_this_N_x ;
    hyde_this_yx(y,:) = hyde_this_x' ;
end


end


function out_YX = read_HYDE_file(thisName,yearN,filename_zip)

filename_this = sprintf('%s%dAD.asc',thisName,yearN) ;
unix(sprintf('unzip -qq %s %s',filename_zip,filename_this)) ;
out_YX = flipud(dlmread(filename_this,' ',6,0)) ;
unix(sprintf('rm %s',filename_this)) ;

end