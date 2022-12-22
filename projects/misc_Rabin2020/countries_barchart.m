%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure after Delzeit et al. (2018) Figure 5 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Get country map and key

% Import
countries_YX = flipud(dlmread('/Users/Shared/PLUM/crop_calib_data/countries/PLUM-specific/country_boundaries62892.noNeg99.extrapd.asc','',6,0)) ;
countries_YX(countries_YX<=0) = NaN ;
countries_key = readtable('/Users/Shared/PLUM/crop_calib_data/countries/PLUM-specific/country_boundaries_codes4.csv') ;
region_list = unique(countries_key.Region) ;
Nregions = length(region_list) ;


%% Import LU data

lu_file = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v12.s1.harm.forLPJG/landcover.txt.gz' ;
lu = lpjgu_matlab_readTable_then2map(lu_file) ;

% Trim extra years
lu.maps_YXvy(:,:,:,lu.yearList<2011) = [] ;
lu.yearList(lu.yearList<2011) = [] ;

% Ignore BARREN
lu.maps_YXvy(:,:,strcmp(lu.varNames,'BARREN'),:) = [] ;
lu.varNames(strcmp(lu.varNames,'BARREN')) = [] ;

% Get area
land_area_YXqd = transpose(ncread('/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc','carea')) ;
land_area_YXqd = double(land_area_YXqd) ;
xres = 360/size(countries_YX,2) ;
yres = 180/size(countries_YX,1) ;
land_area_YX = aggregate_land_area(land_area_YXqd,xres,yres) ;
clear land_area_YXqd
lu.maps_YXvy = lu.maps_YXvy .* repmat(land_area_YX, [1 1 size(lu.maps_YXvy,3) size(lu.maps_YXvy,4)]) ;


%% Get classified data

yearList = [2011 2100] ;
Nyears = length(yearList) ;
Nlu = length(lu.varNames) ;
area_vry = nan(Nlu, Nregions, Nyears) ;

for r = 1:Nregions
    thisRegion = region_list{r} ;
    numCodes_thisReg = countries_key.numCode(strcmp(countries_key.Region,thisRegion)) ;
    [~,IA] = intersect(countries_YX,numCodes_thisReg) ;
    isthisReg = false(size(countries_YX)) ;
    for i = 1:length(numCodes_thisReg)
        isthisReg(countries_YX==numCodes_thisReg(i)) = true ;
    end
    for y = 1:Nyears
        isThisYear = lu.yearList==yearList(y) ;
        for v = 1:Nlu
            tmp = lu.maps_YXvy(:,:,v,isThisYear) ;
            area_vry(v,r,y) = nansum(tmp(isthisReg)) ;
        end
    end
end

if Nyears>2
    error('This won''t work with >2 years.')
end
diff_pct_vr = 100 * (area_vry(:,:,2) - area_vry(:,:,1)) ./ area_vry(:,:,1) ;
diff_pct_vr(area_vry(:,:,1)==0) = NaN ;


%% Make figure

figure('Color','w','Position',[721 34 720 771]) ;
subplot_tight(1,1,1,[0.1 0.3])
barh(categorical(region_list), diff_pct_vr', 'grouped');
xlabel('Area change (%)')
ht = title('Land-use change by region') ;
set(gca,'FontSize',14)
set(ht, 'FontSize', 24)
legend(lu.varNames)




