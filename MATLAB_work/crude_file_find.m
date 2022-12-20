function path_out = crude_file_find(str_in)

path_out = strrep(str_in,'param "file_lu" (str "','') ;
path_out = strrep(path_out,'param "file_lucrop" (str "','') ;
path_out = strrep(path_out,'param "file_nfert" (str "','') ;
path_out = strrep(path_out,'param "file_irrigIntens" (str "','') ;

if ~check_exist(path_out)
    path_out = strrep(path_out,'/project/fh1-project-lpjgpi/lr8247','/Users/Shared') ;
    path_out = strrep(path_out,'/home/fh1-project-lpjgpi/lr8247','/Users/Shared') ;
    path_out = strrep(path_out,'/pfs/data5/home/kit/imk-ifu/lr8247','/Users/Shared') ;
    path_out = strrep(path_out,'/home/kit/imk-ifu/lr8247','/Users/Shared') ;
    if ~check_exist(path_out)
        path_out = strrep(path_out,'/Users/Shared/PLUM/input/PLUMouts_2011-2100','/Users/Shared/PLUM/PLUM_outputs_for_LPJG') ;
        if ~check_exist(path_out)
            error('That didn''t work: file not found! (%s)', str_in)
        end
    end
end



end


function tf_out = check_exist(file_in)

tf_out = ~isempty(dir(file_in)) ...
      || ~isempty(dir([file_in '.gz'])) ...
      || ~isempty(dir([file_in '.mat'])) ;

end
