function plot_scen_diff_timeseries(...
    thisVar, ts1, ts2, yearList, conv_fact, ...
    Nsmth, lineWidth, runList, legend_loc, fontSize, units, ...
    title_prefix, do_caps, subplot_letter, dashed_zero_line, do_pct, ...
    expt_titles)

if iscellstr(thisVar)
    data1 = 0 ;
    data2 = 0 ;
    for v = 1:length(thisVar)
        eval(['data1 = data1 + ts1.' thisVar{v} ' ;']) ;
        eval(['data2 = data2 + ts2.' thisVar{v} ' ;']) ;
    end
else
    eval(['data1 = ts1.' thisVar ' ;']) ;
    eval(['data2 = ts2.' thisVar ' ;']) ;
end
if do_pct
    ydata = movmean((data1-data2)./data2*100,Nsmth) ;
    units = ['% (rel. to "' expt_titles{2} '")'] ;
else
    ydata = conv_fact*movmean(data1-data2,Nsmth) ;
end

plot(yearList,ydata,'LineWidth',lineWidth)
legend(runList,'Location',legend_loc,'AutoUpdate','off') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel(units)

title_suffix = '' ;
if Nsmth>1
    title_suffix = [title_suffix ' (Nsmth=' num2str(Nsmth)] ;
end

ht = title([title_prefix ': "' expt_titles{1} '" - "' expt_titles{2} '"' title_suffix]) ;
if ~isempty(subplot_letter)
    letterlabel_align0(subplot_letter,ht,do_caps) ;
end

if dashed_zero_line
    ylims = get(gca,'YLim') ;
    if min(ylims)<0 && max(ylims)>0
        hold on
        xlims = get(gca,'XLim') ;
        plot(xlims,[0 0],'--k') ;
        hold off
    end
end


end