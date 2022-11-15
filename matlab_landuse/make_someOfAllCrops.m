%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make version with some of all crops %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inDir = '/Users/Shared/SAI-LandSyMM/input/LU/remaps_v10_f09_g17' ;
remove_miscanthus = true ;
notSomeOfAllIrrig = false ;
firstSomeOfAllYear = -Inf ;
% firstSomeOfAllYear = 1995 ;


%% Setup

% Which remap version?
tmp = strsplit(inDir, '/') ;
thisVer = strrep(tmp{end}, 'remaps_', '') ;
clear tmp

% Get input files
inFile_lu = find_input_file(sprintf('%s/LU.remap%s*.txt*', inDir, thisVer)) ;
inFile_cf = find_input_file(sprintf('%s/cropfracs.remap%s*.txt*', inDir, thisVer)) ;

% Get output files
outFile_lu = get_outFile(inFile_lu, firstSomeOfAllYear) ;
outFile_cf = get_outFile(inFile_cf, firstSomeOfAllYear) ;
if notSomeOfAllIrrig
    outFile_cf = strrep(outFile_cf, 'EachCrop', 'EachRainfedCrop') ;
end

outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
donation_order = {'PASTURE','NATURAL','BARREN'} ;
mincropfrac = 2*10^-outPrec ;


%% Import

% Land uses
disp('Importing land uses...')
lu_in = lpjgu_matlab_readTable(inFile_lu,'do_save_MAT',true) ;
lu_in_header = lu_in.Properties.VariableNames ;

% Crop fractions
disp('Importing crop fractions...')
cf_in = lpjgu_matlab_readTable(inFile_cf,'do_save_MAT',true) ;
if remove_miscanthus
    [~,IA] = find(contains(cf_in.Properties.VariableNames, 'Miscanthus')) ;
    cf_in(:,IA) = [] ;
end
cf_in_header = cf_in.Properties.VariableNames ;

% Land area (km2)
disp('Importing land area...')
luMap_in = lpjgu_matlab_readTable_then2map(inFile_lu,'force_mat_save',true) ;
list2map = luMap_in.list_to_map ;
clear luMap_in
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
landArea_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear *_YXqd
%%%%% Get as vectors
landArea_x = landArea_YX(list2map) ;
Nyears = length(unique(lu_in.Year)) ;
tmp = transpose(repmat(landArea_x,[1 Nyears])) ;
landArea_x_allYrs = tmp(:) ;

disp('Done importing.')


%% Convert

lu_hasyears = any(strcmp(lu_in.Properties.VariableNames,'Year')) ;
cf_hasyears = any(strcmp(cf_in.Properties.VariableNames,'Year')) ;
if ~isequal(lu_in.Lon,cf_in.Lon) || ~isequal(lu_in.Lat,cf_in.Lat)
    lonlat_lu = unique([lu_in.Lon lu_in.Lat], 'rows', 'stable') ;
    lonlat_cf = unique([cf_in.Lon cf_in.Lat], 'rows', 'stable') ;
    lu_hasyears_cf_not = lu_hasyears && ~cf_hasyears ...
        &&  isequal(lonlat_lu, lonlat_cf) ;
    if lu_hasyears_cf_not
        warning('lu_in has years but cf_in doesn''t. Shouldn''t be a problem.')
    else
        error('This code assumes equal gridlists for lu_in and cf_in!')
    end
end

% Get Years column(s)
if lu_hasyears
    incl_years_lu = lu_in.Year >= firstSomeOfAllYear ;
else
    incl_years_lu = true(size(lu_in.Lon)) ;
end
if cf_hasyears
    incl_years_cf = cf_in.Year >= firstSomeOfAllYear ;
else
    incl_years_cf = true(size(cf_in.Lon)) ;
end

lu_out = table2array(lu_in) ;
cf_out = table2array(cf_in) ;

moved_area_lu = zeros(size(landArea_x_allYrs)) ;
moved_area_cf = zeros(size(cf_in.Lon)) ;

if mincropfrac>0
    
    % Set up
    iCrop = find(strcmp(lu_in_header(4:end),'CROPLAND')) ;
    lu_tmp = lu_out(incl_years_lu,4:end) ;
    lu_tmp = round(lu_tmp,outPrec) ;
    [~,IA] = setdiff(cf_in_header,{'Lon','Lat','Year'},'stable') ;
    list_crops = cf_in_header(IA) ;
    Ncrops = length(IA) ;
    cf_tmp = cf_out(incl_years_cf,IA) ;
    if any(sum(cf_tmp,2)==0)
        all_zero = sum(cf_tmp,2)==0 ;
        cf_tmp(all_zero,:) = 1/size(cf_tmp, 2) ;
    end
    cf_tmp = cf_tmp ./ repmat(sum(cf_tmp,2),[1 size(cf_tmp,2)]) ;
    
    % Some cropland everywhere
    no_cropland = lu_tmp(:,iCrop)<mincropfrac ;
    no_cropland_orig = no_cropland ;
    
    i = 0 ;
    while(any(no_cropland))
        i = i+1 ;
        if i > length(donation_order)
            error('GET CROPLAND FROM SOMEWHERE')
        end
        this_donor = donation_order{i} ;
        iThis = find(strcmp(lu_in_header(4:end),this_donor)) ;
        involved = lu_tmp(:,iThis)>=2*mincropfrac & no_cropland ;
        if any(involved)
            transfer_amt = mincropfrac-lu_tmp(involved,iCrop) ;
            transfer_amt_area = transfer_amt .* landArea_x_allYrs(involved) ;
            warning('Giving some from %s to CROPLAND (%d cells, %0.1f km2).', this_donor, length(find(involved)), sum(transfer_amt_area)); pause(0.1)
            moved_area_lu(involved) = moved_area_lu(involved) + transfer_amt_area ;
            lu_tmp(involved,iThis) = lu_tmp(involved,iThis) - transfer_amt ;
            lu_tmp(involved,iCrop) = lu_tmp(involved,iCrop) + transfer_amt ;
            no_cropland = lu_tmp(:,iCrop)==0 ;
        end
    end
    lu_out(incl_years_lu,4:end) = lu_tmp ;
    clear lu_tmp
    
    % NOW JUST USING "EVERYWHERE ELSE" ALGORITHM HERE TO ALLOW FOR WHEN LU
    % HAS YEARS BUT CF DOESN'T
%     % Some of each type: Where there was no cropland
%     cf_tmp(repmat(no_cropland_orig,[1 Ncrops])) = 1/Ncrops ;
%     if max(sum(cf_tmp,2))>=1+2*10^(-outPrec)
%         error('max(round(sum_cf))>1')
%     end
    
    % Some of each type: Everywhere else
    for c = 1:Ncrops
        thisCrop_name = list_crops{c} ;
        isIrr = strcmp(thisCrop_name(end), 'i') ;
        % Find rows with 0 for this crop
        thiscrop = cf_tmp(:,c) ;
        iszerothiscrop = thiscrop<mincropfrac ;
        if any(iszerothiscrop) && (~isIrr || ~notSomeOfAllIrrig)
            % Find the crop that currently has the greatest area
            maxcropfrac_Xv = repmat(max(cf_tmp,[],2),[1 Ncrops]) ;
            ismaxcropfrac_Xv = cf_tmp==maxcropfrac_Xv ;
            % Make sure there's only one ismaxcropfrac in each row
            for i = fliplr(2:Ncrops)
                tmp = ismaxcropfrac_Xv(:,i) ;
                sumtoleft = sum(ismaxcropfrac_Xv(:,1:(i-1)),2) ;
                tmp(sumtoleft>0) = false ;
                ismaxcropfrac_Xv(:,i) = tmp ;
            end
            cf_tmp(ismaxcropfrac_Xv & iszerothiscrop) = cf_tmp(ismaxcropfrac_Xv & iszerothiscrop) - mincropfrac ;
            cf_tmp(iszerothiscrop,c) = cf_tmp(iszerothiscrop,c) + mincropfrac ;
            moved_tmp = mincropfrac*lu_out(iszerothiscrop,iCrop).*landArea_x_allYrs(iszerothiscrop) ;
            moved_area_cf(iszerothiscrop) = moved_area_cf(iszerothiscrop) + moved_tmp ;
            warning('Stealing some for %s (%0.1e km2)',thisCrop_name,sum(moved_tmp)); pause(0.1)
        end
    end
    
    % Save
    cf_out(incl_years_cf,IA) = cf_tmp ;
    clear cf_tmp
    
end

% Fake a temporal split if needed
if ~lu_hasyears
    [lu_out, lu_in_header] = fake_yearsplit(lu_in, lu_out, lu_in_header, firstSomeOfAllYear) ;
end
if ~cf_hasyears
    [cf_out, cf_in_header] = fake_yearsplit(cf_in, cf_out, cf_in_header, firstSomeOfAllYear) ;
end


disp('Done')



%% Save

overwrite = true ;

disp('Saving LU...')
lpjgu_matlab_saveTable(lu_in_header, lu_out, outFile_lu,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy) ;

disp('Saving cropfrac...')
lpjgu_matlab_saveTable(cf_in_header, cf_out, outFile_cf,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy) ;


%% FUNCTIONS

function outFile = get_outFile(inFile, firstSomeOfAllYear)

outFile = strrep(strrep(inFile, '.gz', ''), '.mat', '') ;
outFile = strrep(outFile, '.txt', '.someOfEachCrop.txt') ;

if firstSomeOfAllYear > -Inf
    outFile = strrep(inFile, '.txt', ...
        sprintf('.someOfEachCrop.from%d.txt', firstSomeOfAllYear)) ;
end

end


function A = add_yearcol(A, fillyear)

A = cat(2, A, nan(size(A(:,1)))) ;
A(:,4:end) = A(:,3:end-1) ;
A(:,3) = fillyear ;

end


function [A_out, header_cell] = fake_yearsplit(T_in, A_out, header_cell, firstSomeOfAllYear)

A_in = add_yearcol(table2array(T_in), firstSomeOfAllYear-1) ;
A_out = add_yearcol(A_out, firstSomeOfAllYear) ;
A_out = cat(1, A_in, A_out) ;
header_cell = [header_cell(1:2) {'Year'} header_cell(3:end)] ;

end


function inFile = find_input_file(inPattern)

filelist = dir(inPattern) ;

if isempty(filelist)
    error('No file found matching pattern %s', inFile_lu)
end

if length(filelist) > 1
    basenames = cell(length(filelist), 1) ;
    for f = 1:length(filelist)
        basenames{f} = strrep(strrep(strrep(strrep(filelist(f).name, '.txt', ''), '.gz', ''), '.mat', ''), '.garr', '') ;
    end
    basename = unique(basenames) ;
    if length(basename) > 1
        basename
        error('Multiple matches found for %s', inPattern)
    end
    filelist = filelist(1) ;
end

inFile = sprintf('%s/%s', filelist(1).folder, filelist(1).name) ;

end


