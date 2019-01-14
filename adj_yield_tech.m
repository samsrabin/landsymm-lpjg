function data_in = adj_yield_tech(data_in, yearList)

tech_chg_rate = 0.002 ;
base_year = 2010 ;

if iscell(data_in)
    for c = 1:length(data_in)
        data_in{c} = doit(data_in{c}, tech_chg_rate, yearList, base_year) ;
    end
else
    data_in = doit(data_in, tech_chg_rate, yearList, base_year) ;
end


end


function data_in = doit(data_in, tech_chg_rate, yearList, base_year)

data_in = data_in .* repmat(...
    1 + tech_chg_rate*(shiftdim(yearList)-base_year),...
    [1 size(data_in,2)]) ;

end