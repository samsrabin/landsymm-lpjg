function out_lu = remap_import_lu_note_noveg(out_lu, fill_unveg)
% Note cells with no vegetated land according to LU dataset

[~, Ibare] = intersect(out_lu.varNames, {'BARREN'}) ;
if length(Ibare) ~= 1
    error('Expected 1 match of ''BARREN'' in out_lu.varNames; found %d', ...
        length(Ibare))
end
sum_x1y = sum(out_lu.garr_xvy(:,Ibare,:), 2) ;
bad_x = max(sum_x1y, [], 3)==1 ;
Nbad = length(find(bad_x)) ;
if Nbad
    x = find(bad_x, 1) ;
    lonlat = out_lu.lonlats(x, :) ;
    msg = sprintf('%d cells not vegetated; e.g. lon %0.2f lat %0.2f', ...
        Nbad, lonlat(1), lonlat(2)) ;

    if fill_unveg < 0
        max(sum_x1y(x,:,:), [], 3)
        squeeze(out_lu.garr_xvy(x,:,:))
        error(msg) %#ok<SPERR>
    elseif fill_unveg > 1
        error('Maximum fill_unveg is 1; got %g', fill_unveg)
    end
    
    disp(msg)
    if fill_unveg > 0
        [~, Intrl] = intersect(out_lu.varNames, {'NATURAL'}) ;
        if length(Intrl) ~= 1
            error('Expected 1 match of ''NATURAL'' in out_lu.varNames; found %d', ...
                length(Intrl))
        end
        ntrl_x1y = out_lu.garr_xvy(:,Intrl,:) ;
        bare_x1y = out_lu.garr_xvy(:,Ibare,:) ;
        ntrl_x1y(bare_x1y==1) = fill_unveg ;
        bare_x1y(bare_x1y==1) = 1 - fill_unveg ;
        out_lu.garr_xvy(:,Intrl,:) = ntrl_x1y ;
        out_lu.garr_xvy(:,Ibare,:) = bare_x1y ;

        % Check that this worked
        remap_import_lu_note_noveg(out_lu, -1) ;
    elseif ~strcmp(handle_unveg, 'ignore')
        error('Unrecognized handle_unveg: %s', handle_unveg)
    end

end

end