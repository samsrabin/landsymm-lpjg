function data_out = e2p_trim_unneeded(data_in)

data_out = data_in ;

is_unneeded = get_unneeded(data_in.varNames) ;

if ~any(is_unneeded)
    keyboard
end

data_out.varNames(is_unneeded) = [] ;
if isfield(data_in, 'garr_xv')
    data_out.garr_xv(:,is_unneeded) = [] ;
elseif isfield(data_in, 'garr_xvt')
    data_out.garr_xvt(:,is_unneeded,:) = [] ;
else
    error('Neither garr_xv nor garr_xvt found')
end


end