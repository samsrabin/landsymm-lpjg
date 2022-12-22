%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Diagnosing problems Bart highlighted 2022-06-14 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

potYr = 1850 ;

% thisDir = '/Users/Shared/PLUM/forestry_tests_sam/outForPLUM-2022-06-07-012812/1850pot-1949' ;
% thisDir = sprintf('/Users/Shared/PLUM/forestry_tests_sam/outForPLUM-2022-06-21-001756/%dpot-2100_ssp126', potYr) ;
% thisDir = sprintf('/Users/Shared/PLUM/forestry_tests_sam/outForPLUM-2022-06-21-060203/%dpot-2100_ssp126', potYr) ;
% thisDir = sprintf('/Users/Shared/PLUM/forestry_tests_sam/outForPLUM-2022-06-30-040356/%dpot_1850-2014', potYr) ;
thisDir = sprintf('/Users/Shared/PLUM/forestry_tests_sam/outForPLUM-2022-06-30-040356/%dpot_2015-2100_ssp126', potYr) ;

% fluxList = {'pcutC', 'pcutW', 'plutC', 'plutW'} ;
donor = 'forC' ;


%% Setup

landarea_file = '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc' ;
% Use gridcell area instead of land area
land_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
addpath '/Users/Shared/PLUM/crop_calib_code'
land_area_YX = aggregate_land_area(land_area_YXqd,0.5,0.5) ;
clear land_area_YXqd


%% Potential land-use transitions

fluxList = {'plutC', 'plutW', 'Total'} ;
Nflux = length(fluxList) ;

for f = 1:Nflux
    thisFlux = fluxList{f} ;

    if strcmp(thisFlux, 'Total')
        garr_xvyf(:,:,:,f) = sum(garr_xvyf(:,:,:,1:(f-1)), 4) ;
    else
        data = lpjgu_matlab_read2geoArray(sprintf('%s/landsymm_%s_from_%s.out', thisDir, thisFlux, donor), ...
            'force_mat_save', false, 'force_mat_nosave', true) ;
    
        if f == 1
            yearList = data.yearList ;
            Nyears = length(data.yearList) ;
            Ncells = length(data.list2map) ;
            land_area_x = land_area_YX(data.list2map) ;
            land_area = sum(land_area_x) ;
            list2map = data.list2map ;
            lonlats = data.lonlats ;
            receptorList = data.varNames ;
            Nreceptors = length(receptorList) ;
            garr_xvyf = nan([Ncells Nreceptors Nyears Nflux]) ;
        end
        garr_xvyf(:,:,:,f) = data.garr_xvy ;
    end

end

spacing = [0.075 0.05] ; % v h
fontSize = 14 ;

figure('Color', 'w', 'Position', figurePos) ;
for v = 1:Nreceptors
    subplot_tight(2,2,v,spacing)
    receptor = strrep(receptorList{v}, 'to_', '') ;

    for f = 1:Nflux
        hold on
        data_xy = squeeze(garr_xvyf(:,v,:,f)) ;
        data_xy = data_xy .* repmat(land_area_x / land_area, [1 Nyears]) ;
        plot(yearList, sum(data_xy,1))
        hold off
    end

    this_legend = strrep(receptorList, 'to_', sprintf('%s to ', donor)) ;

    legend(fluxList, 'Location', 'best')

    title(sprintf('%dpot: %s to %s', potYr, donor, receptor))
    ylabel('kgC m^{-2}')
    set(gca, 'FontSize', fontSize)
end


%% Potential wood harvest

fluxList = {'pcutC', 'pcutW', 'Total'} ;
Nflux = length(fluxList) ;

for f = 1:Nflux
    thisFlux = fluxList{f} ;

    if strcmp(thisFlux, 'Total')
        garr_xvyf(:,:,:,f) = sum(garr_xvyf(:,:,:,1:(f-1)), 4) ;
    else
        data = lpjgu_matlab_read2geoArray(sprintf('%s/landsymm_%s_sts.out', thisDir, thisFlux)) ;
    
        if f == 1
            yearList = data.yearList ;
            Nyears = length(data.yearList) ;
            Ncells = length(data.list2map) ;
            land_area_x = land_area_YX(data.list2map) ;
            land_area = sum(land_area_x) ;
            list2map = data.list2map ;
            lonlats = data.lonlats ;
            sourceList = data.varNames ;
            Nsources = length(sourceList) ;
            garr_xvyf = nan([Ncells Nsources Nyears Nflux]) ;
        end
        garr_xvyf(:,:,:,f) = data.garr_xvy ;
    end

end

spacing = [0.075 0.05] ; % v h
fontSize = 14 ;

figure('Color', 'w', 'Position', [1    41   791   764]) ;
for v = 1:Nsources
    subplot_tight(2,1,v,spacing)
    thisSource = strrep(sourceList{v}, 'to_', '') ;

    for f = 1:Nflux
        hold on
        data_xy = squeeze(garr_xvyf(:,v,:,f)) ;
        data_xy = data_xy .* repmat(land_area_x / land_area, [1 Nyears]) ;
        plot(yearList, sum(data_xy,1))
        hold off
    end

    this_legend = strrep(sourceList, 'to_', sprintf('%s to ', donor)) ;

    legend(fluxList, 'Location', 'best')

    title(sprintf('%dpot: %s', potYr, thisSource))
    ylabel('kgC m^{-2}')
    set(gca, 'FontSize', fontSize)
end


