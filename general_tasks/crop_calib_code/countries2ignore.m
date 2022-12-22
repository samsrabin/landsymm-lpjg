function ignore_Cc = countries2ignore(croparea_Ccy)

croparea_Cc = nansum(croparea_Ccy,3) ;
ignore_Cc = croparea_Cc==0 ;


end