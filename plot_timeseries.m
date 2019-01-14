function ht = plot_timeseries(...
    yearList_baseline, yearList_fao, yearList_future, ...
    ts_in_bl, ts_in_fao, ts_in_yr, ...
    conv_fact, Nsmth, ...
    ylab, thisLegend, title_prefix, title_suffix, lineWidth, fontSize, ...
    skip3rdColor)

plot(yearList_baseline,conv_fact*movmean(ts_in_bl,Nsmth),'-k','LineWidth',lineWidth)
hold on

if ~isempty(yearList_fao) && ~isempty(ts_in_fao)
    plot(yearList_fao,conv_fact*ts_in_fao,'--k','LineWidth',lineWidth)
end

if skip3rdColor
    plot(yearList_future,conv_fact*ts_in_yr(:,1),'LineWidth',lineWidth)
    h = plot(yearList_future,conv_fact*ts_in_yr(:,2),'LineWidth',lineWidth) ;
    set(h,'Color',[255 159 56]/255) ;
    set(gca,'ColorOrderIndex',4) ;
    plot(yearList_future,conv_fact*ts_in_yr(:,3),'LineWidth',lineWidth)
else
    plot(yearList_future,conv_fact*ts_in_yr,'LineWidth',lineWidth)
end

hold off
legend(thisLegend, ...
    'Location','NorthWest') ;
set(gca,'FontSize',fontSize)
xlabel('Year')
ylabel(ylab)
ht = title([title_prefix title_suffix]) ;



end