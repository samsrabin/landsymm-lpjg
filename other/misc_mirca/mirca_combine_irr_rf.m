%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Manipulating MIRCA data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd '/Users/sam/Geodata/MIRCA/harvested_area_grids_26crops_30mn'
addpath(genpath(pwd))

list_crops = {'Wheat' ; 'Maize' ; 'Rice' ; 'Barley' ; 'Rye' ; 'Millet' ;
              'Sorghum' ; 'Soybeans' ; 'Sunflower' ; 'Potatoes' ; 'Cassava' ;
              'Sugarcane' ; 'Sugarbeet' ; 'Oilpalm' ; 'RapeseedCanola' ;
              'GroundnutsPeanuts' ; 'Pulses' ; 'Citrus' ; 'Datepalm' ;
              'GrapesVine' ; 'Cotton' ; 'Cocoa' ; 'Coffee' ;
              'OtherPerennials' ; 'FodderGrasses' ; 'OtherAnnuals'} ;          
          

%% Import crops

file_list_ir = dir('originals/*irc*.asc') ;
file_list_rf = dir('originals/*rfc*.asc') ;
if length(file_list_ir) ~= length(file_list_rf)
    error('length(file_list_ir) ~= length(file_list_rf)')
end
Ncrops = length(file_list_ir) ;
if Ncrops ~= length(list_crops)
    error('Ncrops ~= length(list_crops)')
end

maps_ir_YXc = nan(360,720,Ncrops) ;
maps_rf_YXc = nan(360,720,Ncrops) ;
for c = 1:Ncrops
    disp(['Importing: ' list_crops{c} ' (' num2str(c) ' of ' num2str(Ncrops) ')'])
    maps_ir_YXc(:,:,c) = import_mirca(['originals/' file_list_ir(c).name]) ;
    maps_rf_YXc(:,:,c) = import_mirca(['originals/' file_list_rf(c).name]) ;
end


%% Combine

maps_area_YXc = maps_ir_YXc + maps_rf_YXc ;
maps_prop_YXc = maps_area_YXc ./ repmat(sum(maps_area_YXc,3),[1 1 Ncrops]) ;

list_TomIgn = {'Sugarcane';'Oilpalm';'Citrus';'Datepalm';'GrapesVine';
               'Cotton';'Cocoa';'Coffee';'OtherAnnuals';'OtherPerennials'} ;

maps_TomIgn_area_YX = zeros(360,720) ;
maps_TomIgn_prop_YX = zeros(360,720) ;
for c = 1:length(list_TomIgn)
    maps_TomIgn_area_YX = maps_TomIgn_area_YX + maps_area_YXc(:,:,find_string_in_cell_exact(list_crops,list_TomIgn{c})) ;
    maps_TomIgn_prop_YX = maps_TomIgn_prop_YX + maps_prop_YXc(:,:,find_string_in_cell_exact(list_crops,list_TomIgn{c})) ;
end

maps_TomInc_prop_YX = 1 - maps_TomIgn_prop_YX ;



%% Save

out_header = {'NCOLS 720' ;
              'NROWS 360' ;
              'XLLCORNER -180' ;
              'YLLCORNER -90' ;
              'CELLSIZE 0.5' ;
              'NODATA_VALUE -9999'} ;

for c = 1:Ncrops
    disp(['Saving: ' list_crops{c} ' (' num2str(c) ' of ' num2str(Ncrops) ')'])
    save_ascii_map(['combined_eachCrop_area/annCropArea_' list_crops{c} '.asc'],...
                    maps_area_YXc(:,:,c),out_header) ;
    save_ascii_map(['combined_eachCrop_prop/annCropProp_' list_crops{c} '.asc'],...
                    maps_prop_YXc(:,:,c),out_header) ;
end ; clear c

disp('Saving: Other...')
save_ascii_map('combined_eachCrop_area/annCropArea_TomIgnored.asc',maps_TomIgn_area_YX,out_header)
save_ascii_map('combined_eachCrop_prop/annCropProp_TomIgnored.asc',maps_TomIgn_prop_YX,out_header)
save_ascii_map('combined_eachCrop_prop/annCropProp_TomIncluded.asc',maps_TomInc_prop_YX,out_header)
tmp = nan(size(maps_TomInc_prop_YX)) ;
tmp(maps_TomInc_prop_YX==0) = 1 ;
save_ascii_map('combined_eachCrop_prop/annCropProp_TomIncluded_isZero.asc',tmp,out_header)
clear tmp


%% Check against Hurtt data

% Import Hurtt gridlist
which_hurtt_gridlist = 'expanded' ;
if strcmp(which_hurtt_gridlist,'original')
    scripts_dir = '/Users/sam/Documents/Dropbox/lpj-guess-plum/input/unifying_gridlist/extrapd_to_hurttList/handy_MATLAB_scripts/' ;
    addpath(genpath(scripts_dir))
    thisDir = '/Users/sam/Documents/Dropbox/lpj-guess-plum/input/unifying_gridlist/' ;
    filename_hurttRNDM = 'lists_for_MATLAB/gridlist_hurtt_RNDM.csv' ;
    import_hurtt_RNDM
elseif strcmp(which_hurtt_gridlist,'expanded')
    filename_hurttRNDM = '/Users/sam/Documents/Dropbox/lpj-guess-plum/input/unifying_gridlist/lists_for_MATLAB/gridlist_hurtt_grossfrac.csv' ;
    gridlist_hurttRNDM = readtable(filename_hurttRNDM) ;
    gridlist_hurttRNDM = sortrows(gridlist_hurttRNDM,{'Lon','Lat'},{'ascend','ascend'}) ;
    hurtt_lons = gridlist_hurttRNDM.Lon ;
    hurtt_lats = gridlist_hurttRNDM.Lat ;
    Ncells_hurtt = length(hurtt_lons) ;
else
    error('Invalid "which_hurtt_gridlist".')
end
% Import
in_dir = '/Users/Shared/PLUM/input/LU/' ;
in_filename = 'netfrac_hurtt_hist_1901_2014.out' ;
[in_path,in_filename_noext,in_ext] = fileparts([in_dir in_filename]) ;
tmp_filename = [in_path '/' in_filename_noext '.mat'] ;
if ~exist(tmp_filename,'file')
    in_data = dlmread([in_dir in_filename],'',1,0) ;
    in_lons = in_data(:,1) ;
    in_lats = in_data(:,2) ;
    in_years = in_data(:,3) ;
    in_data(:,1:3) = [] ;
    in_total = sum(in_data,2) ;
    save(tmp_filename,'in_lons','in_lats','in_years','in_data','in_total')
else
    load(tmp_filename)
end
clear tmp_filename
lu_lons = in_lons ;
lu_lats = in_lats ;
lu_years = in_years ;
lu_data = in_data ; 
lu_total = in_total ;
clear in_*
years_list = unique(lu_years) ;
Nyears = length(years_list) ;
Ncells_lu = length(lu_years)/Nyears ;
lu_lats_1yr = lu_lats(lu_years==1901) ;
lu_lons_1yr = lu_lons(lu_years==1901) ;
lu_data_u = lu_data(:,1) ;
lu_data_c = lu_data(:,2) ;
lu_data_p = lu_data(:,3) ;
lu_data_f = lu_data(:,4) ;
lu_data_n = lu_data(:,5) ;
lu_data_b = lu_data(:,6) ;

% Make mask
lons = -180:0.5:179.5 ;
lats = -90:0.5:89.5 ;
lons_map = repmat(lons,[length(lats) 1]) ;
lats_map = repmat(lats',[1 length(lons)]) ;
hurtt_mask = false(360,720) ;
for y = 1:length(lats)
     thislat = lats(y) ;
     hurtt_lons_in_thislat = hurtt_lons(hurtt_lats==thislat) ;
     if isempty(hurtt_lons_in_thislat)
         continue
     end
     for x = 1:length(lons)
         thislon = lons(x) ;
         if ~isempty(find(hurtt_lons_in_thislat==thislon,1))
             hurtt_mask(y,x) = true ;
         end
     end ; clear x
end ; clear y

% Make maps
lons = -180:0.5:179.5 ;
lats = -90:0.5:89.5 ;
lons_map = repmat(lons,[length(lats) 1]) ;
lats_map = repmat(lats',[1 length(lons)]) ;
list_to_map_lu = nan(Ncells_lu,1) ;
disp('Getting indices to convert lists to maps...')
progress = 0 ;
progress_step_pct = 5 ;
tic ;
for c = 1:Ncells_lu
    list_to_map_lu(c) = find(lats_map==lu_lats_1yr(c) & lons_map==lu_lons_1yr(c)) ;
    % Update progress
    if rem(c,ceil(Ncells_lu*progress_step_pct/100))==0
        progress = progress + progress_step_pct ;
        disp(['   ' num2str(progress) '% complete (' toc_hms(toc) ')'])
    end
end
disp('Converting list to maps...')
maps_lu_YXcy = nan(length(lats),length(lons),6,Nyears) ;
for y = 1:Nyears
    thisYear = years_list(y) ;
    disp(num2str(thisYear))
    lu_data_inThisYear = lu_data(lu_years==thisYear,:) ;
    thisYear_YXc = nan(length(lats),length(lons),6) ;
    for c = 1:6
        tmp = nan(length(lats),length(lons)) ;
        tmp(list_to_map_lu) = lu_data_inThisYear(:,c) ;
        thisYear_YXc(:,:,c) = tmp ;
        clear tmp
    end
    clear lu_data_inThisYear
    maps_lu_YXcy(:,:,:,y) = thisYear_YXc ;
    clear thisYear_YXc
end
disp('Done.')


%% Where does the MIRCA dataset have 0 where LU cropland > 0?
maps_lu_c_YXy = squeeze(maps_lu_YXcy(:,:,2,:)) ;

mirca_crop_area_YX = sum(maps_area_YXc,3) ;
mirca_crop_area_YXy = repmat(mirca_crop_area_YX,[1 1 Nyears]) ;

disp(['There are ' num2str(length(find(maps_lu_c_YXy>0 & mirca_crop_area_YXy==0))) ' cells where crop fractions dataset has 0 but'])
disp(['land use dataset has >0 (' num2str(roundn(100*length(find(maps_lu_c_YXy>0 & mirca_crop_area_YXy==0))/length(lu_data_c),-1)) '%).'])

disp(['There are ' num2str(length(find(maps_lu_c_YXy>0 & isnan(repmat(sum(maps_area_YXc,3),[1 1 Nyears]))))) ' cells where crop fractions dataset has NaN but'])
disp(['land use dataset has >0 (' num2str(roundn(100*length(find(maps_lu_c_YXy>0 & isnan(mirca_crop_area_YXy)))/length(lu_data_c),-1)) '%).'])

bad_because_zero = mean(maps_lu_c_YXy>0 & mirca_crop_area_YXy==0,3) ;
bad_because_nans = mean(maps_lu_c_YXy>0 & isnan(mirca_crop_area_YXy),3) ;
good_because_zero = mean(maps_lu_c_YXy==0 & mirca_crop_area_YXy==0,3) ;
bad_because_zero(~hurtt_mask) = NaN ;
bad_because_nans(~hurtt_mask) = NaN ;
good_because_zero(~hurtt_mask) = NaN ;

figure ;
colormap('jet')
subplot_tight(3,1,1,0.05)
pcolor(bad_because_zero) ; shading flat ; axis equal tight
caxis([0 1]) ; colorbar ;
title(['MIRCA==0 but CROP > 0 (' num2str(roundn(100*length(find(maps_lu_c_YXy>0 & mirca_crop_area_YXy==0))/length(lu_data_c),-1)) '%).'])
set(gca,'XTick',[],'YTick',[])
subplot_tight(3,1,2,0.05)
pcolor(bad_because_nans) ; shading flat ; axis equal tight
caxis([0 1]) ; colorbar ;
title(['MIRCA missing but CROP > 0 (' num2str(roundn(100*length(find(maps_lu_c_YXy>0 & isnan(mirca_crop_area_YXy)))/length(lu_data_c),-1)) '%).'])
set(gca,'XTick',[],'YTick',[])
subplot_tight(3,1,3,0.05)
pcolor(good_because_zero) ; shading flat ; axis equal tight
caxis([0 1]) ; colorbar ;
title(['MIRCA==0 and CROP==0 (' num2str(roundn(100*length(find(maps_lu_c_YXy==0 & mirca_crop_area_YXy==0))/length(lu_data_c),-1)) '%).'])
set(gca,'XTick',[],'YTick',[])

figure ;
hist(maps_lu_c_YXy(maps_lu_c_YXy>0 & mirca_crop_area_YXy==0),20)
histc(maps_lu_c_YXy(maps_lu_c_YXy>0 & mirca_crop_area_YXy==0),0:0.05:1)

figure ;
tmp = mean(maps_lu_c_YXy,3) ;
tmp(~hurtt_mask) = NaN ;
tmp(bad_because_zero==0) = NaN ;
colormap('jet')
pcolor(tmp) ; shading flat ; axis equal tight
caxis([0 max(caxis)])
colorbar ;
title('Mean fraction CROP')
set(gca,'XTick',[],'YTick',[])

figure ;
tmp = maps_lu_c_YXy ;
tmp(~(maps_lu_c_YXy>0 & mirca_crop_area_YXy==0)) = NaN ;
tmp = nanmean(tmp,3) ;
tmp(~hurtt_mask) = NaN ;
tmp(bad_because_zero==0) = NaN ;
colormap('jet')
pcolor(tmp) ; shading flat ; axis equal tight
caxis([0 max(caxis)])
colorbar ;
title('Mean fraction CROP')
set(gca,'XTick',[],'YTick',[])