function e2p_check_correct_zeros(data_gvt, which_file, basenamesi)
% If which_file is yield, returns error if no positive cells found for one
% or more crops.
% If which_file is gsirrigation, returns error if no positive cells found
% for one or more irrigated crops, or if any positive cells are found for
% one or more rainfed crops.

if strcmp(which_file, 'yield')
    test_gt0_is_good(data_gvt)
elseif strcmp(which_file, 'gsirrigation')
    test_isirrig = @(x) strcmp(x(end),'i') ;
    isirrig = cellfun(test_isirrig, basenamesi) ;
    data_gvt_rf = data_gvt(:,~isirrig,:) ;
    data_gvt_ir = data_gvt(:,isirrig,:) ;
    test_gt0_is_bad(data_gvt_rf)
    test_gt0_is_good(data_gvt_ir)
else
    error('which_file (%s) not recognized', which_file)
end


end


function test_gt0_is_good(data_gvt)

for t = 1:size(data_gvt, 3)
    if any(~any(data_gvt(:,:,t)>0,1))
        Ifound = find(~any(data_gvt(:,:,t)>0,1)) ;
        Nfound = length(Ifound) ;
        error('No positive cells found for %d/%d variables (first %d) at timestep %d', ...
            Nfound, size(data_gvt, 2), Ifound(1), t)
    end
end

end


function test_gt0_is_bad(data_gvt)

for t = 1:size(data_gvt, 3)
    if any(any(data_gvt(:,:,t)>0,1))
        Ifound = find(~any(data_gvt(:,:,t)>0,1)) ;
        Nfound = length(Ifound) ;
        error('Positive cells found for %d/%d variables (first %d) at timestep %d', ...
            Nfound, size(data_gvt, 2), Ifound(1), t)
    end
end


end