function save_geotiffs(maps_YXr, filename, runList, R, gtif_missing)

for r = 1:length(runList)
    thisYX = maps_YXr(:,:,r) ;
    thisFile = strrep(filename,'.png',['.' runList{r} '.tif']) ;
    thisFile = strrep(thisFile, '/maps/', '/gtif/') ;
    geotiffwrite_ssr(thisFile, thisYX, R, gtif_missing) ;
end


end