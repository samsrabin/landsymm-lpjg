function get_LU_timeseries(thisSSPdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get land use time series for all PLUM outputs in a given SSP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%
%%% Setup %%%
%%%%%%%%%%%%%

cd(thisSSPdir) ;

PUTURBANHERE = 'BARREN' ;

% Get LUH2 land area (m2)
tmp = pwd ;
if strcmp(tmp(1:5),'/User')
    onMac = true ;
elseif strcmp(tmp(1:5),'/pfs/')
    onMac = false ;
else
    error('What system are you on?')
end
clear tmp
if onMac
    landarea_file = '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc' ;
else
    landarea_file = '/home/fh1-project-lpjgpi/lr8247/PLUM/input/LUH2/supporting/staticData_quarterdeg.nc' ;
end
gcel_area_YXqd = 1e6*double(transpose(ncread(landarea_file,'carea'))) ;
land_frac_YXqd = 1 - double(flipud(transpose(ncread(landarea_file,'icwtr')))) ;
landArea_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear *_YXqd tmp


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Import PLUM outputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirList = dir('s*') ;
Ndirs = length(dirList) ;

for d = 1:Ndirs
    
    thisDir = dirList(d).name ;
    
    % Get list of years
    tmp = dir(sprintf('%s/2*',thisDir)) ;
    yearList_tmp = cellfun(@str2num,{tmp(:).name}') ;
    if d==1
        yearList = yearList_tmp ;
        Nyears = length(yearList) ;
    elseif ~isequal(yearList, yearList_tmp)
        error(['Not all PLUM replicates have same number of years!\n' ...
               '(%s %d-%d, previously %s %d-%d)'], ...
            thisDir, min(yearList_tmp), max(yearList_tmp), dirList(1).name, min(yearList), max(yearList))
    end
    clear tmp yearList_tmp
    
    % Import each year
    for y = 1:Nyears
        thisYear = yearList(y) ;
        fprintf('%s/%d... ', thisDir, thisYear) ;
        
        % Import
        thisFile = sprintf('%s/%d/LandCoverFract.txt', thisDir, thisYear) ;
        S_lcf = lpjgu_matlab_readTable_then2map(thisFile,'verboseIfNoMat',false,'force_mat_save',true) ;
        
        % Move URBAN into PUTURBANHERE; remove URBAN
        if any(strcmp(S_lcf.varNames,'URBAN'))
            S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,PUTURBANHERE)) = ...
                sum(S_lcf.maps_YXv(:,:,contains(S_lcf.varNames,{PUTURBANHERE,'URBAN'})),3) ;
            S_lcf.maps_YXv(:,:,strcmp(S_lcf.varNames,'URBAN')) = [] ;
            S_lcf.varNames(strcmp(S_lcf.varNames,'URBAN')) = [] ;
        end
        
        % Get this year's totals
        if d==1 && y==1
            LUnames = S_lcf.varNames ;
            Nlu = length(LUnames) ;
            LUts_yvd = nan(Nyears, Nlu, Ndirs) ;
            landArea_YXv = repmat(landArea_YX,[1 1 Nlu]) ;
        elseif ~isequal(S_lcf.varNames,LUnames)
            fprintf('\n')
            error('Not all PLUMreplicate-years have the same LU names! (%s %d)\n', ...
                thisDir, thisYear) ;
        end
        LUts_yvd(y,:,d) = permute(nansum(nansum(S_lcf.maps_YXv.*landArea_YXv, 2), 1), [1 3 2]) ;
        
        
        fprintf('Done.\n')
    end % for each year
end % for each dir

save('LU_timeseries.mat', 'LUts_yvd', 'yearList', 'LUnames', 'dirList')


end


