function [maps1_YXB, maps1_YXr, maps2_YXB, maps2_YXr]= ...
    ESscatter_getMeanOverPeriod(...
        mapstructX, thisVarX, ...
        mapstructY, thisVarY)

maps1_YXB = mean(mapstructX.maps_YXvyB(:,:,strcmp(mapstructX.varNames,thisVarX),:),4) ;
maps1_YXr = squeeze(mean(mapstructX.maps_YXvyr(:,:,strcmp(mapstructX.varNames,thisVarX),:,:),4)) ;
maps2_YXB = mean(mapstructY.maps_YXvyB(:,:,strcmp(mapstructY.varNames,thisVarY),:),4) ;
maps2_YXr = squeeze(mean(mapstructY.maps_YXvyr(:,:,strcmp(mapstructY.varNames,thisVarY),:,:),4)) ;


end