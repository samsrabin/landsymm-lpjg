if ~strcmp(thisVer, 'harm3')
    error('You''ve only done this for harm3 so far. Current thisVer = %s.', thisVer)
end

if ~exist('garr_NTRLfrac_luh1_d9', 'var')
    luh1_file_list = { ...
        '/Volumes/WDMPP_Storage/ExternalGeodata/LUH1/future/states/minicam.v1.mat', ... % RCP4.5
        '/Volumes/WDMPP_Storage/ExternalGeodata/LUH1/future/states/aim.v1.1.mat', ... % RCP6.0
        '/Volumes/WDMPP_Storage/ExternalGeodata/LUH1/future/states/aim.v1.1.mat', ... % RCP6.0
        '/Volumes/WDMPP_Storage/ExternalGeodata/LUH1/future/states/message.v1.mat' ... % RCP8.5
        } ;
    
    % Import
    disp('Importing LUH1...')
    garr_NTRLfrac_luh1_d9 = struct ;
    for r = 1:Nruns
        thisFile = luh1_file_list{r} ;
        if r>1 && strcmp(thisFile, luh1_file_list{r-1})
            garr_NTRLfrac_luh1_d9.garr_xrB(:,r) = garr_NTRLfrac_luh1_d9.garr_xrB(:,r-1) ;
            garr_NTRLfrac_luh1_d9.garr_xrF(:,r) = garr_NTRLfrac_luh1_d9.garr_xrF(:,r-1) ;
            continue
        end
        if r==1
            garr_NTRLfrac_luh1_d9.garr_xrB = nan(length(garr_LU_d9.list2map), Nruns, 'single') ;
            garr_NTRLfrac_luh1_d9.garr_xrF = nan(length(garr_LU_d9.list2map), Nruns, 'single') ;
        end
        
        % Import
        load(thisFile) ; % --> LU_luh1
        
        % Remove unnecessary years
        bad_years = ~(LU_luh1.yearList>=2006 & LU_luh1.yearList<=2010) & ~(LU_luh1.yearList>=2096 & LU_luh1.yearList<=2100) ;
        LU_luh1.maps_YXvy(:,:,:,bad_years) = [] ;
        LU_luh1.yearList(bad_years) = [] ;
        
        % Remove ice/water
        LU_luh1.maps_YXvy = LU_luh1.maps_YXvy ./ sum(LU_luh1.maps_YXvy(:,:,~strcmp(LU_luh1.varNames,'icwr'),:),3) ;
        LU_luh1.maps_YXvy(:,:,strcmp(LU_luh1.varNames,'icwr'),:) = [] ;
        LU_luh1.varNames(strcmp(LU_luh1.varNames,'icwr')) = [] ;
        
        % Extract
        is_othr = contains(LU_luh1.varNames, {'othr','secd'}) ;
        tmp_YX = squeeze(mean(sum(LU_luh1.maps_YXvy(:,:,is_othr,LU_luh1.yearList>=2006 & LU_luh1.yearList<=2010),3),4)) ;
        garr_NTRLfrac_luh1_d9.garr_xrB(:,r) = tmp_YX(garr_LU_d9.list2map) ;
        tmp_YX = squeeze(mean(sum(LU_luh1.maps_YXvy(:,:,is_othr,LU_luh1.yearList>=2096 & LU_luh1.yearList<=2100),3),4)) ;
        garr_NTRLfrac_luh1_d9.garr_xrF(:,r) = tmp_YX(garr_LU_d9.list2map) ;
        
        clear LU_luh1
    end
end

if ~exist('garr_NTRLfrac_plum_d9', 'var')
    is_othr = contains(garr_LU_d9.varNames, {'NATURAL'  'BARREN'}) ;
    garr_NTRLfrac_plum_d9.garr_xrB = repmat(mean(sum(garr_LU_d9.garr_xvyB(:,is_othr,6:end),2),3), [1 Nruns]) ;
    garr_NTRLfrac_plum_d9.garr_xrF = squeeze(mean(sum(garr_LU_d9.garr_xvyr(:,is_othr,6:end,:),2),3)) ;
end
