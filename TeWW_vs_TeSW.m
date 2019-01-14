%%%%%%%%%%%%%%%%%%%%%
%%% TeWW or TeSW? %%%
%%%%%%%%%%%%%%%%%%%%%

inDir = '/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/justEvenWheats_1850-2005/output-2018-01-13-153758' ;

 inFile_cropfracs = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/LU/crop_fractions_fromMIRCA_PLUM_1yr_20180105b.MP.out' ;
outFile_cropfracs = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/LU/crop_fractions_fromMIRCA_PLUM_XXXX-XXXX_20180105b.MP.assignWWorSW.out' ;
outPrec_LC = 3 ;

 inFile_nfert = '/Users/Shared/PLUM/input/Nfert/LUH2/Nfert_fromLUH2.min5kgha.1700-2015.64493.20161216.txt' ;
outFile_nfert = '/Users/Shared/PLUM/input/Nfert/LUH2/Nfert_fromLUH2.min5kgha.XXXX-XXXX.YYYYY.20161216.assignWWorSW.txt' ;
outPrec_mgmtInputs = 6 ;


%% Setup

cd '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work'
addpath(genpath(pwd))

inDir = find_PLUM2LPJG_run(inDir) ;

% Anonymous function for checking whether any version of a file exists
anyfileexist = @(in_file) ...
    exist(in_file,'file') ...
    || exist([in_file '.maps.mat'],'file') ...
    || exist([in_file '.mat'],'file') ...
    || exist([in_file '.gz'],'file') ;

% Get CFTs of interest
CFTnames = {'TeWW','TeSW','TeWWi','TeSWi'} ;
Ncfts = length(CFTnames) ;

% For output
force_overwrite = true ;
fclose_every = 1000 ;
pct_progress = 5 ;
lons_all = -179.75:0.5:179.75 ;
lats_all = -89.75:0.5:89.75 ;
lons_map = repmat(lons_all,[length(lats_all) 1]) ;
lats_map = repmat(lats_all',[1 length(lons_all)]) ;



%% Import

disp('Importing...')
yield = lpjgu_matlab_readTable_then2map([inDir 'yield.out'],'force_mat_save',true,'verboseIfNoMat',true) ;
hdate1 = lpjgu_matlab_readTable_then2map([inDir 'hdate1.out'],'force_mat_save',true,'verboseIfNoMat',true) ;
hdate2 = lpjgu_matlab_readTable_then2map([inDir 'hdate2.out'],'force_mat_save',true,'verboseIfNoMat',true) ;

% Trim unnecessary CFTs
disp('Finishing...')
[~,IA] = setdiff(yield.varNames,CFTnames) ;
if length(yield.varNames)-length(IA) ~= Ncfts
    error('length(yield.varNames)-length(IA) ~= Ncfts')
elseif length(hdate1.varNames)-length(IA) ~= Ncfts
    error('length(hdate1.varNames)-length(IA) ~= Ncfts')
elseif length(hdate2.varNames)-length(IA) ~= Ncfts
    error('length(hdate2.varNames)-length(IA) ~= Ncfts')
end
yield.maps_YXvy(:,:,IA,:) = [] ;
yield.varNames(IA) = [] ;
hdate1.maps_YXvy(:,:,IA,:) = [] ;
hdate1.varNames(IA) = [] ;
hdate2.maps_YXvy(:,:,IA,:) = [] ;
hdate2.varNames(IA) = [] ;

% Get info
yearList_in = yield.yearList ;
Nyears = length(yearList_in) ;
isOK_YX = ~isnan(yield.maps_YXvy(:,:,1,1)) ;
if ~isequal(size(lons_map),size(isOK_YX))
    error('~isequal(size(lons_map),size(isOK_YX))')
end
disp('Done.')


% %% Find where TeWW(i) is better (simple)
% 
% TeWW_better = squeeze(yield.maps_YXvy(:,:,strcmp(yield.varNames,'TeWW'),:) >= yield.maps_YXvy(:,:,strcmp(yield.varNames,'TeSW'),:)) ;
% TeWWi_better = squeeze(yield.maps_YXvy(:,:,strcmp(yield.varNames,'TeWWi'),:) >= yield.maps_YXvy(:,:,strcmp(yield.varNames,'TeSWi'),:)) ;


%% Find where TeWW(i) is better (complicated)

disp('Finding where TeWW(i) is better...')

y1s_in = ceil((min(yearList_in)+1)/5)*5+1:5:floor(max(yearList_in)/5)*5-4 ;
yNs_in = y1s_in+4 ;

y1s_out = y1s_in + 5 ;
yNs_out = yNs_in + 5 ;

Npds = length(y1s_in) ;

TeWW_better_YXp = false(size(yield.maps_YXvy,1),size(yield.maps_YXvy,2),Npds) ;
TeWWi_better_YXp = false(size(yield.maps_YXvy,1),size(yield.maps_YXvy,2),Npds) ;
for p = 1:Npds
    % Get indices
    thisY1 = y1s_in(p) ;
    thisYN = yNs_in(p) ;
    I = yield.yearList>=thisY1 & yield.yearList<=thisYN ;
    % Disp info
    thisY1_out = y1s_out(p) ;
    thisYN_out = yNs_out(p) ;
    disp(['   ' num2str(thisY1) '-' num2str(thisYN) ' --> ' num2str(thisY1_out) '-' num2str(thisYN_out)])
    % Get data
    tmpYield_YXv = sum(yield.maps_YXvy(:,:,:,I),4) ;
    tmpHdate1_YXvy = hdate1.maps_YXvy(:,:,:,I) ;
    tmpHdate2_YXvy = hdate2.maps_YXvy(:,:,:,I) ;
    Nharvs_YXv = sum(tmpHdate1_YXvy>0,4) + sum(tmpHdate2_YXvy>0,4) ;
    % Calculate mean per-season yield
    tmpYield_TeWW = tmpYield_YXv(:,:,strcmp(yield.varNames,'TeWW')) ./ Nharvs_YXv(:,:,strcmp(yield.varNames,'TeWW')) ;
    tmpYield_TeSW = tmpYield_YXv(:,:,strcmp(yield.varNames,'TeSW')) ./ Nharvs_YXv(:,:,strcmp(yield.varNames,'TeSW')) ;
    tmpYield_TeWWi = tmpYield_YXv(:,:,strcmp(yield.varNames,'TeWWi')) ./ Nharvs_YXv(:,:,strcmp(yield.varNames,'TeWWi')) ;
    tmpYield_TeSWi = tmpYield_YXv(:,:,strcmp(yield.varNames,'TeSWi')) ./ Nharvs_YXv(:,:,strcmp(yield.varNames,'TeSWi')) ;
    % Get which is better
    TeWW_better_YXp(:,:,p) = tmpYield_TeWW >= tmpYield_TeSW ;
    TeWWi_better_YXp(:,:,p) = tmpYield_TeWWi >= tmpYield_TeSWi ;
end

disp('   Finishing...')

% Expand
TeWW_better_YXy = repmat(TeWW_better_YXp,[1 1 1 5]) ;
TeWW_better_YXy = permute(TeWW_better_YXy,[1 2 4 3]) ;
TeWW_better_YXy = reshape(TeWW_better_YXy,[size(TeWW_better_YXp,1) size(TeWW_better_YXp,2) Npds*5]) ;
TeWWi_better_YXy = repmat(TeWWi_better_YXp,[1 1 1 5]) ;
TeWWi_better_YXy = permute(TeWWi_better_YXy,[1 2 4 3]) ;
TeWWi_better_YXy = reshape(TeWWi_better_YXy,[size(TeWWi_better_YXp,1) size(TeWWi_better_YXp,2) Npds*5]) ;
yearList_tmp = y1s_out(1):yNs_out(end) ;

% Prepend first "out" year to line up beginnings of in and out
Nmissing = length(find(yearList_in<y1s_out(1))) ;
TeWW_better_YXy = cat(3,repmat(TeWW_better_YXy(:,:,1),[1 1 Nmissing]),TeWW_better_YXy) ;
TeWWi_better_YXy = cat(3,repmat(TeWWi_better_YXy(:,:,1),[1 1 Nmissing]),TeWWi_better_YXy) ;
yearList_out = [yearList_in ; yearList_tmp(end-4:end)'] ;
Nyears_out = length(yearList_out) ;
if Nyears_out ~= size(TeWW_better_YXy,3) || Nyears_out ~= size(TeWWi_better_YXy,3)
    error('Nyears_out ~= size(TeWW_better_YXy,3) || Nyears_out ~= size(TeWWi_better_YXy,3)')
end
disp('Done.')


%% Import cropfracs

cropfracs_in = lpjgu_matlab_readTable_then2map(inFile_cropfracs,'force_mat_save',true,'verboseIfNoMat',true) ;

% CFTnames = {'TeWWi','TeSWi','TeCoi','TrRii','TeWW','TeSW','TeCo','TrRi'} ;
% if ~isequal(CFTnames,cropfracs_in.varNames)
%     error('Rework code to work with this exact set of cropfracs_in.varNames!')
% end
CFTnames = cropfracs_in.varNames ;
Ncfts = length(CFTnames) ;

% Get cropfracs_in.maps_YXvy, if needed
if ~isfield(cropfracs_in,'maps_YXvy')
    cropfracs_in.maps_YXvy = repmat(cropfracs_in.maps_YXv,[1 1 1 Nyears_out]) ;
    cropfracs_in.yearList = yearList_out ;
end

% Make sure yearLists match!
if ~isequal(cropfracs_in.yearList,yearList_out)
    error('~isequal(cropfracs_in.yearList,yearList_out)')
end

% Set up for output
cropfracs_out = cropfracs_in ;
if isfield(cropfracs_out,'maps_YXv')
    cropfracs_out = rmfield(cropfracs_out,'maps_YXv') ;
end


%% Import Nfert

nfert_in = lpjgu_matlab_readTable_then2map(inFile_nfert,'force_mat_save',true,'verboseIfNoMat',true) ;

% Set up for output
nfert_out = nfert_in ;
nfert_out.maps_YXvy(:,:,:,nfert_in.yearList<yearList_out(1) | nfert_in.yearList>yearList_out(end)) = [] ;
nfert_out.yearList(nfert_in.yearList<yearList_out(1) | nfert_in.yearList>yearList_out(end)) = [] ;
if ~isequal(nfert_out.yearList,yearList_out)
    error('~isequal(nfert_out.yearList,yearList_out)')
end



%% Move original TeSW(i) to TeXX(i)

%%% To keep wheat separate from oilcrops/starchy roots/pulses

% Cropfracs
cropfracs_out.maps_YXvy(:,:,end+1,:) = cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeSW'),:) ;
cropfracs_out.varNames{end+1} = 'TeXX' ;
cropfracs_out.maps_YXvy(:,:,strcmp(CFTnames,'TeSW'),:) = 0 ;
cropfracs_out.maps_YXvy(:,:,end+1,:) = cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeSWi'),:) ;
cropfracs_out.varNames{end+1} = 'TeXXi' ;
cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeSWi'),:) = 0 ;

% Nfert
nfert_out.maps_YXvy(:,:,end+1,:) = nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeSW'),:) ;
nfert_out.varNames{end+1} = 'TeXX' ;
nfert_out.maps_YXvy(:,:,strcmp(CFTnames,'TeSW'),:) = 0 ;
nfert_out.maps_YXvy(:,:,end+1,:) = nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeSWi'),:) ;
nfert_out.varNames{end+1} = 'TeXXi' ;
nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeSWi'),:) = 0 ;


%% Assign original TeWW(i) to either TeWW(i) or TeSW(i)

% Do it: cropfracs, rainfed
tmpWW_YXy = squeeze(cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeWW'),:)) ;
tmpSW_YXy = squeeze(cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeSW'),:)) ;
tmpSW_YXy(~TeWW_better_YXy) = tmpWW_YXy(~TeWW_better_YXy) ;
tmpWW_YXy(~TeWW_better_YXy) = 0 ;
cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeWW'),:) = reshape(tmpWW_YXy,[size(tmpWW_YXy,1) size(tmpWW_YXy,2) 1 size(tmpWW_YXy,3)]) ;
cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeSW'),:) = reshape(tmpSW_YXy,[size(tmpSW_YXy,1) size(tmpSW_YXy,2) 1 size(tmpSW_YXy,3)]) ;
clear tmp*

% Do it: cropfracs, irrigated
tmpWWi_YXy = squeeze(cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeWWi'),:)) ;
tmpSWi_YXy = squeeze(cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeSWi'),:)) ;
tmpSWi_YXy(~TeWWi_better_YXy) = tmpWWi_YXy(~TeWWi_better_YXy) ;
tmpWWi_YXy(~TeWWi_better_YXy) = 0 ;
cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeWWi'),:) = reshape(tmpWWi_YXy,[size(tmpWWi_YXy,1) size(tmpWWi_YXy,2) 1 size(tmpWWi_YXy,3)]) ;
cropfracs_out.maps_YXvy(:,:,strcmp(cropfracs_out.varNames,'TeSWi'),:) = reshape(tmpSWi_YXy,[size(tmpSWi_YXy,1) size(tmpSWi_YXy,2) 1 size(tmpSWi_YXy,3)]) ;
clear tmp*

cropfracs_out.maps_YXvy(~repmat(isOK_YX,[1 1 size(cropfracs_out,3) Nyears_out])) = NaN ;


% Do it: nfert, rainfed
tmpWW_YXy = squeeze(nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeWW'),:)) ;
tmpSW_YXy = squeeze(nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeSW'),:)) ;
tmpSW_YXy(~TeWW_better_YXy) = tmpWW_YXy(~TeWW_better_YXy) ;
tmpWW_YXy(~TeWW_better_YXy) = 0 ;
nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeWW'),:) = reshape(tmpWW_YXy,[size(tmpWW_YXy,1) size(tmpWW_YXy,2) 1 size(tmpWW_YXy,3)]) ;
nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeSW'),:) = reshape(tmpSW_YXy,[size(tmpSW_YXy,1) size(tmpSW_YXy,2) 1 size(tmpSW_YXy,3)]) ;
clear tmp*

% Do it: nfert, irrigated
tmpWWi_YXy = squeeze(nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeWWi'),:)) ;
tmpSWi_YXy = squeeze(nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeSWi'),:)) ;
tmpSWi_YXy(~TeWWi_better_YXy) = tmpWWi_YXy(~TeWWi_better_YXy) ;
tmpWWi_YXy(~TeWWi_better_YXy) = 0 ;
nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeWWi'),:) = reshape(tmpWWi_YXy,[size(tmpWWi_YXy,1) size(tmpWWi_YXy,2) 1 size(tmpWWi_YXy,3)]) ;
nfert_out.maps_YXvy(:,:,strcmp(nfert_out.varNames,'TeSWi'),:) = reshape(tmpSWi_YXy,[size(tmpSWi_YXy,1) size(tmpSWi_YXy,2) 1 size(tmpSWi_YXy,3)]) ;
clear tmp*

nfert_out.maps_YXvy(~repmat(isOK_YX,[1 1 size(nfert_out,3) Nyears_out])) = NaN ;


%% Save cropfracs

disp('Preparing to save cropfracs...')

% Get header and formatSpec
out_header_cropfracs = 'Lon Lat Year' ;
for c = 1:length(cropfracs_out.varNames)
    out_header_cropfracs = [out_header_cropfracs ' ' cropfracs_out.varNames{c}] ;
end
out_formatSpec_cropfracs = ['%4.2f %4.2f %4.0f' repmat([' %5.' num2str(outPrec_LC) 'f'],[1 length(cropfracs_out.varNames)])] ;
out_formatSpec_cropfracs = [out_formatSpec_cropfracs '\n'] ;

% Set up output files
outFile_cropfracs = strrep(outFile_cropfracs,'XXXX-XXXX',[num2str(yearList_out(1)) '-' num2str(yearList_out(end))]) ;
outFile_cropfracs = strrep(outFile_cropfracs,'YYYYY',num2str(Ncells_out)) ;
if exist(outFile_cropfracs,'file') && ~force_overwrite
    error(['Outfile already exists! (' outFile_cropfracs ')']) ;
end
fid1_cropfracs = fopen(outFile_cropfracs, 'w') ;
fprintf(fid1_cropfracs,'%s \n',out_header_cropfracs);

cropfracs_maps_yvYX = permute(cropfracs_out.maps_YXvy,[4 3 1 2]) ;

Ncells_out = length(yield.list_to_map) ;
disp('Saving...')
for c = 1:Ncells_out
    
    i = yield.list_to_map(c) ;
    [I,J] = ind2sub([size(cropfracs_out.maps_YXvy,1) size(cropfracs_out.maps_YXvy,2)],i) ;
    thisMat = cropfracs_maps_yvYX(:,:,I,J) ;
    thisLon = lons_map(i) ;
    thisLat = lats_map(i) ;
    
    thisOut = [thisLon*ones(Nyears_out,1) ...
               thisLat*ones(Nyears_out,1) ...
               yearList_out ...
               thisMat] ;
           
    if c>1 && rem(c-1,fclose_every)==0
        fid1_cropfracs = fopen(outFile_cropfracs, 'a+') ;
    end
    
    fprintf(fid1_cropfracs,out_formatSpec_cropfracs,thisOut') ;
    
    if rem(c,fclose_every)==0 || c==Ncells_out
        fclose(fid1_cropfracs) ;
    end
    
    if rem(c,pct_progress*floor(Ncells_out/100))==0
        compltn = 100*(c)/Ncells_out ;
        if compltn>0
            disp(['   ' num2str(round(compltn)) '% complete...'])
        end
    end
end
disp('Done.')


%% Save nfert

disp('Preparing to save nfert...')

% Get header and formatSpec
out_header_nfert = 'Lon Lat Year' ;
for c = 1:length(nfert_out.varNames)
    out_header_nfert = [out_header_nfert ' ' nfert_out.varNames{c}] ;
end
out_formatSpec_nfert = ['%4.2f %4.2f %4.0f' repmat([' %5.' num2str(outPrec_mgmtInputs) 'f'],[1 length(nfert_out.varNames)])] ;
out_formatSpec_nfert = [out_formatSpec_nfert '\n'] ;

% Set up output files
outFile_nfert = strrep(outFile_nfert,'XXXX-XXXX',[num2str(yearList_out(1)) '-' num2str(yearList_out(end))]) ;
outFile_nfert = strrep(outFile_nfert,'YYYYY',num2str(Ncells_out)) ;
if exist(outFile_nfert,'file') && ~force_overwrite
    error(['Outfile already exists! (' outFile_nfert ')']) ;
end
fid1_nfert = fopen(outFile_nfert, 'w') ;
fprintf(fid1_nfert,'%s \n',out_header_nfert);

nfert_maps_yvYX = permute(nfert_out.maps_YXvy,[4 3 1 2]) ;

Ncells_out = length(yield.list_to_map) ;
disp('Saving...')
for c = 1:Ncells_out
    
    i = yield.list_to_map(c) ;
    [I,J] = ind2sub([size(nfert_out.maps_YXvy,1) size(nfert_out.maps_YXvy,2)],i) ;
    thisMat = nfert_maps_yvYX(:,:,I,J) ;
    thisLon = lons_map(i) ;
    thisLat = lats_map(i) ;
    
    thisOut = [thisLon*ones(Nyears_out,1) ...
               thisLat*ones(Nyears_out,1) ...
               yearList_out ...
               thisMat] ;
           
    if c>1 && rem(c-1,fclose_every)==0
        fid1_nfert = fopen(outFile_nfert, 'a+') ;
    end
    
    fprintf(fid1_nfert,out_formatSpec_nfert,thisOut') ;
    
    if rem(c,fclose_every)==0 || c==Ncells_out
        fclose(fid1_nfert) ;
    end
    
    if rem(c,pct_progress*floor(Ncells_out/100))==0
        compltn = 100*(c)/Ncells_out ;
        if compltn>0
            disp(['   ' num2str(round(compltn)) '% complete...'])
        end
    end
           
end
disp('Done.')
