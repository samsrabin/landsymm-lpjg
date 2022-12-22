p95_q20c = prctile(maps_pk_runoff_d9.maps_YXvyB,95,4) ;
p95_q21c = prctile(maps_pk_runoff_d9.maps_YXvyr(:,:,:,:,1),95,4) ;
p05_q20c = prctile(maps_pk_runoff_d9.maps_YXvyB,5,4) ;
p05_q21c = prctile(maps_pk_runoff_d9.maps_YXvyr(:,:,:,:,1),5,4) ;
below_thresh = mean(maps_awater_d1.maps_YXvyB(:,:,strcmp(maps_awater_d1.varNames,'Runoff'),:),4)/365 < 0.01 ;
p95_q20c(below_thresh) = NaN ;
p95_q21c(below_thresh) = NaN ;
p05_q20c(below_thresh) = NaN ;
p05_q21c(below_thresh) = NaN ;
% p95_dq = (p95_q21c - p95_q20c) ./ (p95_q21c + p95_q20c) ;
% p05_dq = (p05_q21c - p05_q20c) ./ (p05_q21c + p05_q20c) ;
tmp95 = p95_q21c-p95_q20c ;
mean_incWhereInc = 100*(nansum(nansum(p95_q21c(tmp95>0))) - nansum(nansum(p95_q20c(tmp95>0)))) / nansum(nansum(p95_q20c(tmp95>0))) ;
mean_decWhereDec = 100*(nansum(nansum(p95_q21c(tmp95<0))) - nansum(nansum(p95_q20c(tmp95<0)))) / nansum(nansum(p95_q20c(tmp95<0))) ;
shademap(p95_q21c-p95_q20c);
title(sprintf('p95: %d pct inc (mean %d), %d pct dec (mean %d)', ...
    round(length(find(~isnan(tmp95) & tmp95>0))/ length(find(~isnan(tmp95)))*100), ...
    round(mean_incWhereInc), ...
    round(length(find(~isnan(tmp95) & tmp95<0))/ length(find(~isnan(tmp95)))*100), ...
    round(mean_decWhereDec))) ;
caxis([-max(abs(caxis)) max(abs(caxis))])
tmp05 = p05_q21c-p05_q20c ;
mean_incWhereInc = 100*(nansum(nansum(p05_q21c(tmp05>0))) - nansum(nansum(p05_q20c(tmp05>0)))) / nansum(nansum(p05_q20c(tmp05>0))) ;
mean_decWhereDec = 100*(nansum(nansum(p05_q21c(tmp05<0))) - nansum(nansum(p05_q20c(tmp05<0)))) / nansum(nansum(p05_q20c(tmp05<0))) ;
shademap(tmp05);
title(sprintf('p05: %d pct dec (mean %d), %d pct inc (mean %d)', ...
    round(length(find(~isnan(tmp05) & tmp05<0))/ length(find(~isnan(tmp05)))*100), ...
    round(mean_decWhereDec), ...
    round(length(find(~isnan(tmp05) & tmp05>0))/ length(find(~isnan(tmp05)))*100), ...
    round(mean_incWhereInc))) ;
caxis([-max(abs(caxis)) max(abs(caxis))])
both = tmp95>0 & tmp05<0 ;
notnan = ~isnan(tmp95.*tmp05) ;
fprintf('p95 inc + p05 dec: %d pct\n', round(length(find(both & notnan))/length(find(notnan))*100))
