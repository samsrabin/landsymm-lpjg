function [table1_out, table2_out, table3_out, table4_out, table5_out, table6_out] = make_ES_table(...
    rowInfo, years_endh, years_begf, years_endf, ...
    Nruns, yearList_baseline, yearList_future, struct_in, runColNames, ...
    combine_meanAndErr)

Nvars = size(rowInfo,1) ;
mean_endh = nan(Nvars,1) ;
mean_begf = nan(Nvars,Nruns) ;
mean_endf = nan(Nvars,Nruns) ;
sem_endh = nan(Nvars,1) ;
sem_begf = nan(Nvars,Nruns) ;
sem_endf = nan(Nvars,Nruns) ;
if combine_meanAndErr
    string_endh = cell(Nvars,1) ;
    string_begf = cell(Nvars,Nruns) ;
    string_endf = cell(Nvars,Nruns) ;
end
for c = 1:Nvars
    
    % Get values
    thisVar = rowInfo{c,2} ;
    thisConv = rowInfo{c,3} ;
    mean_endh(c) = thisConv*eval(['mean(struct_in.ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
    sem_endh(c) = thisConv*eval(['std(struct_in.ts_' thisVar '_bl(yearList_baseline>=min(years_endh) & yearList_baseline<=max(years_endh)))']) ;
    mean_begf(c,:) = thisConv*eval(['mean(struct_in.ts_' thisVar '_yr(yearList_future>=min(years_begf) & yearList_future<=max(years_begf),:))']) ;
    sem_begf(c,:) = thisConv*eval(['std(struct_in.ts_' thisVar '_yr(yearList_future>=min(years_begf) & yearList_future<=max(years_begf),:))']) ;
    mean_endf(c,:) = thisConv*eval(['mean(struct_in.ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
    sem_endf(c,:) = thisConv*eval(['std(struct_in.ts_' thisVar '_yr(yearList_future>=min(years_endf) & yearList_future<=max(years_endf),:))']) ;
    
    if combine_meanAndErr
        % Turn into strings
        if strcmp(rowInfo{c,4},'%d')
            thisMean = round(mean_endh(c)) ;
        else
            thisMean = mean_endh(c) ;
        end
        if strcmp(rowInfo{c,4},'%d')
            thisSD = round(sem_endh(c)) ;
        else
            thisSD = sem_endh(c) ;
        end
        string_endh{c} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean thisSD]) ;
        for r = 1:Nruns
            if strcmp(rowInfo{c,4},'%d')
                thisMean_begf = round(mean_begf(c,r)) ;
                thisMean_endf = round(mean_endf(c,r)) ;
            else
                thisMean_begf = mean_begf(c,r) ;
                thisMean_endf = mean_endf(c,r) ;
            end
            if strcmp(rowInfo{c,4},'%d')
                thisSD_begf = round(sem_begf(c,r)) ;
                thisSD_endf = round(sem_endf(c,r)) ;
            else
                thisSD_begf = sem_begf(c,r) ;
                thisSD_endf = sem_endf(c,r) ;
            end
            string_begf{c,r} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean_begf thisSD_begf]) ;
            string_endf{c,r} = sprintf([rowInfo{c,4} ' ± ' rowInfo{c,5}],[thisMean_endf thisSD_endf]) ;
        end
    end

end

table2_out = [] ;

if combine_meanAndErr
    table1_out = table(collate_empties(rowInfo(:,1)),...
        collate_empties(string_endh)) ;
    for r = 1:Nruns
        table1_out = [table1_out collate_twocells(string_begf(:,r),string_endf(:,r))] ;
    end
else
    table1_out = [] ;
    for r = 1:Nruns
        error('figure this out')
        table1_out = [table1_out ; mean_endh] ;
    end
end

table1_out.Properties.VariableNames = [{'Ecosystem_function','Baseline'} runColNames] ;
if ~combine_meanAndErr
    table2_out.Properties.VariableNames = [{'Ecosystem_function','Baseline'} runColNames] ;
end


end