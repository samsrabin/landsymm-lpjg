%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make plot of LU time series for all PLUM outputs in a given SSP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Import data 
file_in = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/ssp12/SSP3/LU_timeseries.mat' ;
load(file_in)

Nlu = length(LUnames) ;
Nyears = length(yearList) ;

[badYr, badDir] = find(squeeze(sum(LUts_yvd,2))<1.2616e14) ;
if ~isempty(badYr)
    warnMsg = '' ;
    for i = 1:length(badYr)
        tmp = sprintf('\n s%d %d\n', badDir(i), yearList(badYr(i))) ;
        warnMsg = [warnMsg tmp] ;
        LUts_yvd(badYr(i),:,badDir(i)) = NaN ;
    end
    warnMsg = ['Something went wrong in the following (setting to NaN):' warnMsg] ;
    warning(warnMsg)
end


%% Get info on minimum and maximum of each LU at end of run

for v = 1:Nlu
    
    thisLU = LUnames{v} ;
    if strcmp(thisLU,'BARREN')
        thisLU = 'CROP+PAST' ;
        eor_thisLU = squeeze(sum(LUts_yvd(end,contains(LUnames,{'CROPLAND','PASTURE'}),:),2)) ;
    else
        eor_thisLU = squeeze(LUts_yvd(end,v,:)) ;
    end
    
    
    thisMinRun = find(eor_thisLU==min(eor_thisLU)) ;
    thisMaxRun = find(eor_thisLU==max(eor_thisLU)) ;
    
    fprintf('Minimum %s:\ts%d\n', thisLU, thisMinRun) ;
    fprintf('Maximum %s:\ts%d\n', thisLU, thisMaxRun) ;
    disp(' ')
end


%% Make figure

figure('Position',figurePos,'Color','w') ;

for v = 1:Nlu
    subplot_tight(2,2,v,0.04)
    thisLU = LUnames{v} ;
    
    theseTS_yd = squeeze(LUts_yvd(:,v,:)) ;
    theseTS_yd = theseTS_yd - repmat(theseTS_yd(1,:),[Nyears 1]) ;
    plot(yearList, theseTS_yd) ;
    title(thisLU)
end


%% Save figure

export_fig(strrep(file_in,'.mat','.pdf'))
close


%%