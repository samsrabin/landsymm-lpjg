function [ts_yr_out, title_suffix, file_suffix] = ...
    rebase_future2baseline(rebase, Nsmth, ts_bl, ts_yr, ignYrs, yrList_future)

title_suffix = '' ;
file_suffix = '' ;

if ~rebase && ignYrs==0 && Nsmth==1
    ts_yr_out = ts_yr ;
else
    
    if ignYrs > 0
        title_suffix = [' (ign. ' num2str(yrList_future(1)) '-' num2str(yrList_future(ignYrs)) ')'] ;
        file_suffix = ['_ign' num2str(ignYrs)] ;
        ts_yr(1:ignYrs,:) = NaN ;
    end
    
    if Nsmth>1
        title_suffix = [title_suffix ' (movmean' num2str(Nsmth) ')'] ;
        file_suffix = [file_suffix '_mm' num2str(Nsmth)] ;
        if rebase
%             warning('Figure out how to do smoothing with rebase! Turning off smoothing.')
        else
            ts_yr_tmp = ts_yr(ignYrs+1:end,:) ;
            ts_yr_tmp = movmean(ts_yr_tmp,Nsmth,1) ;
            ts_yr = [nan(ignYrs,size(ts_yr_tmp,2));ts_yr_tmp] ;
        end
    end
    
    if rebase
        title_suffix = [title_suffix ' (rebased)'] ;
        file_suffix = [file_suffix '_rb'] ;
        ts_yr_out = nan(size(ts_yr)) ;
        for r = 1:size(ts_yr,2)
            yrDiff = ts_bl(end) - ts_yr(ignYrs+1,r) ;
            ts_yr_out(:,r) = ts_yr(:,r) + yrDiff ;
            if Nsmth>1
                ts_yr_out(:,r) = movmean(ts_yr_out(:,r),Nsmth,1) ;
            end
        end
    else
        ts_yr_out = ts_yr ;
    end
    
end