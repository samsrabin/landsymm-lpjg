%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make version with some of all crops %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inFile_lu = '/Users/Shared/PLUM/input/remaps_v5/LU.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% inFile_cf = '/Users/Shared/PLUM/input/remaps_v5/cropfracs.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt' ;

% inFile_lu = '/Users/Shared/PLUM/input/remaps_v6/LU.remapv6.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% inFile_cf = '/Users/Shared/PLUM/input/remaps_v6/cropfracs.remapv6.20180214.ecFertIrr0.setaside0103.m4.txt' ;

% inFile_lu = '/Users/Shared/PLUM/input/remaps_v6p3/LU.remapv6p3.txt' ;
% inFile_cf = '/Users/Shared/PLUM/input/remaps_v6p3/cropfracs.remapv6p3.txt' ;

% inFile_lu = '/Users/Shared/PLUM/input/remaps_v6p6/LU.remapv6p6.txt' ;
% inFile_cf = '/Users/Shared/PLUM/input/remaps_v6p6/cropfracs.remapv6p6.txt' ;
% outFile_lu = '/Users/Shared/PLUM/input/remaps_v6p6/LU.remapv6p6.someOfEachCrop.txt' ;
% outFile_cf = '/Users/Shared/PLUM/input/remaps_v6p6/cropfracs.remapv6p6.someOfEachCrop.txt' ;

% inFile_lu = '/Users/Shared/PLUM/input/remaps_v6p7/LU.remapv6p7.txt' ;
% inFile_cf = '/Users/Shared/PLUM/input/remaps_v6p7/cropfracs.remapv6p7.txt' ;
% outFile_lu = '/Users/Shared/PLUM/input/remaps_v6p7/LU.remapv6p7.someOfEachCrop.txt' ;
% outFile_cf = '/Users/Shared/PLUM/input/remaps_v6p7/cropfracs.remapv6p7.someOfEachCrop.txt' ;

% inFile_lu = '/Users/Shared/PLUM/input/remaps_v7a/LU.remapv7a.txt' ;
% inFile_cf = '/Users/Shared/PLUM/input/remaps_v7a/cropfracs.remapv7a.txt' ;
% outFile_lu = '/Users/Shared/PLUM/input/remaps_v7a/LU.remapv7a.someOfEachCrop.txt' ;
% outFile_cf = '/Users/Shared/PLUM/input/remaps_v7a/cropfracs.remapv7a.someOfEachCrop.txt' ;

inFile_lu = '/Users/Shared/PLUM/input/remaps_v7b/LU.remapv7b.txt' ;
inFile_cf = '/Users/Shared/PLUM/input/remaps_v7b/cropfracs.remapv7b.txt' ;
outFile_lu = '/Users/Shared/PLUM/input/remaps_v7b/LU.remapv7b.someOfEachCrop.txt' ;
outFile_cf = '/Users/Shared/PLUM/input/remaps_v7b/cropfracs.remapv7b.someOfEachCrop.txt' ;


%% Setup

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

if ~isequal(lu_in.Lon,cf_in.Lon) || ~isequal(lu_in.Lat,cf_in.Lat)
    error('This code assumes equal gridlists for lu_in and cf_in!')
end

lu_out = table2array(lu_in) ;
cf_out = table2array(cf_in) ;

moved_area = zeros(size(landArea_x_allYrs)) ;

if mincropfrac>0
    
    % Set up
    iCrop = find(strcmp(lu_in_header(4:end),'CROPLAND')) ;
    lu_tmp = lu_out(:,4:end) ;
    lu_tmp = round(lu_tmp,outPrec) ;
    [~,IA] = setdiff(cf_in_header,{'Lon','Lat','Year'},'stable') ;
    cf_tmp = cf_out(:,IA) ;
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
            moved_area(involved) = moved_area(involved) + transfer_amt_area ;
            lu_tmp(involved,iThis) = lu_tmp(involved,iThis) - transfer_amt ;
            lu_tmp(involved,iCrop) = lu_tmp(involved,iCrop) + transfer_amt ;
            no_cropland = lu_tmp(:,iCrop)==0 ;
        end
    end
    lu_out(:,4:end) = lu_tmp ;
    
    % Some of each type: Where there was no cropland
    Ncrops = length(IA) ;
    cf_tmp(repmat(no_cropland_orig,[1 Ncrops])) = 1/Ncrops ;
    if max(sum(cf_tmp,2))>=1+2*10^(-outPrec)
        error('max(round(sum_cf))>1')
    end
    
    % Some of each type: Everywhere else
    for c = 1:Ncrops
        % Find rows with 0 for this crop
        thiscrop = cf_tmp(:,c) ;
        iszerothiscrop = thiscrop<mincropfrac ;
        if any(iszerothiscrop)
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
            moved_tmp = mincropfrac*lu_tmp(iszerothiscrop,iCrop).*landArea_x_allYrs(iszerothiscrop) ;
            moved_area(iszerothiscrop) = moved_area(iszerothiscrop) + moved_tmp ;
            warning('Stealing some for %s (%0.1e km2)',cf_in_header{3+c},sum(moved_tmp)); pause(0.1)
        end
    end
    
    % Save
    cf_out(:,IA) = cf_tmp ;
    %         clear cf_tmp Ncrops
    
end
disp('Done')



%% Save

lpjgu_matlab_saveTable(lu_in_header, lu_out, outFile_lu,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', false, ...
    'fancy', fancy) ;

lpjgu_matlab_saveTable(cf_in_header, cf_out, outFile_cf,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', false, ...
    'fancy', fancy) ;




