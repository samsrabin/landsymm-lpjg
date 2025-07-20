function PLUMharm_check_LU_map_sum_1(map_YXv)
% Ensure that LU fractions sum to 1 within some tolerance

if ndims(map_YXv) ~= 3
    error('Expected 3d map_YXv; got %d dims', ndims(map_YXv))
end


sum_YX = sum(map_YXv, 3) ;
abs_diff_from_1_YX = abs(1 - sum_YX) ;
max_diff_from_1 = nanmax(nanmax(abs_diff_from_1_YX)) ;
if max_diff_from_1 > 1e-6
    warning('Sum of LU fractions differs from 1 by up to %0.4g', max_diff_from_1)

    % Arbitrarily deciding to allow this if the total divergence is small.
    % If total divergence is large, then you'll start to notice problems
    % like https://github.com/samsrabin/landsymm-lpjg/issues/1
    sum_diff_from_1 = nansum(nansum(abs_diff_from_1_YX)) ;
    if sum_diff_from_1 > 0.1

        % Plot map
        minmax_ssr(sum_YX)
        shademap(sum_YX) ;
        colormap(parula)
        title('Sum of all LU. Expected 1 everywhere; stopping.')

        % Throw error
        error('Stopping because total divergence is >0.1. Try commenting this out to see if you get apparent large gaps after harmonization.')
    end
end

end