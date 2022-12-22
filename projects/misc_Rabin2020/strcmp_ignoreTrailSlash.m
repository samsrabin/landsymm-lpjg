function tf_out = strcmp_ignoreTrailSlash(dir1,dir2)

dir1 = removeslashifneeded(dir1) ;
dir2 = removeslashifneeded(dir2) ;

tf_out = strcmp(dir1, dir2) ;

end