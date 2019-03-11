function data_in = adj_yield_tech(data_in, yearList, varargin)

tech_chg_rate = 0.002 ;
base_year = 2010 ;

if ~isempty(varargin)
    yearList_baseline = yearList ;
    yearList_future = varargin{1} ;
    if length(varargin)>1
        error('At most 1 optional argument allowed! (yearList_future)')
    end
end
    
if iscell(data_in)
    chg_mult_eachYear = get_change_mult(tech_chg_rate, yearList, base_year) ;
    for c = 1:length(data_in)
        data_in{c} = doit(data_in{c}, chg_mult_eachYear) ;
    end
elseif isstruct(data_in)
    data_in.maps_YXvyB = doit_maps_YXvyQ(data_in.maps_YXvyB, get_change_mult(tech_chg_rate, yearList_baseline, base_year)) ;
    data_in.maps_YXvyr = doit_maps_YXvyQ(data_in.maps_YXvyr, get_change_mult(tech_chg_rate, yearList_future, base_year)) ;
else
    chg_mult_eachYear = get_change_mult(tech_chg_rate, yearList, base_year) ;
    data_in = doit(data_in, chg_mult_eachYear) ;
end


end


function chg_mult_eachYear = get_change_mult(tech_chg_rate, yearList, base_year)

chg_mult_eachYear = 1 + tech_chg_rate*(shiftdim(yearList)-base_year) ;

end



function data_in = doit(data_in, chg_mult_eachYear)

data_in = data_in .* repmat(...
    chg_mult_eachYear, ...
    [1 size(data_in,2)]) ;

end



function data_in_YXvyQ = doit_maps_YXvyQ(data_in_YXvyQ, chg_mult_eachYear)

if ndims(data_in_YXvyQ)<4
    error('doit_maps_YXvyQ() assumes at least 4 dimensions of array! (ndims=%d', ndims(data_in_YXvyQ))
end

array_size = size(data_in_YXvyQ) ;
if array_size(4)~=length(chg_mult_eachYear)
    error('doit_maps_YXvyQ() assumes "year" is fourth dimension of array!')
elseif length(find(array_size==length(chg_mult_eachYear)))>1
    warning('doit_maps_YXvyQ() assumes "year" is fourth dimension of array, but there are %d length-matching dimensions.', length(find(array_size==length(chg_mult_eachYear))))
end
chg_rate_vector = permute(chg_mult_eachYear,[4 3 2 1]) ;

repmat_dims = [array_size(1:3) 1] ;
if ndims(data_in_YXvyQ)>4
    repmat_dims = [repmat_dims array_size(5:end)] ;
end
    
data_in_YXvyQ = data_in_YXvyQ .* ...
    repmat(chg_rate_vector, repmat_dims) ;

end