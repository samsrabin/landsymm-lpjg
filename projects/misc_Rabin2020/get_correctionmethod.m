function correctionmethod = get_correctionmethod(in_dir)

[status, result] = unix(sprintf( ...
    'grep -he "^param \\"%s" %s/main_plum_lpjg.ins', ...
    'correctionmethod', in_dir)) ;
if status~=0
    error('Error in unix() call')
end
tmp = strsplit(result,'"') ;
correctionmethod = tmp{4} ;


end