function out_YXv = PLUMharm_fixTinyNegs(in_YXv, landArea_YXv, ...
    LUnames, outPrec, fixTinyNegs_tol_m2, conserv_tol_area, debugIJ_2deg)
% Rounding errors can result in small negative values. Fix.

Nlu = size(in_YXv, 3) ;
do_debug = ~isempty(debugIJ_2deg) ;

fixTinyNegs_tol_m2 = abs(fixTinyNegs_tol_m2) ;
conserv_tol_area = abs(conserv_tol_area) ;

if fixTinyNegs_tol_m2 > conserv_tol_area
    error('fixTinyNegs_tol_m2 must be ≤ conserv_tol_area')
end

if min(in_YXv(:)) < -fixTinyNegs_tol_m2
    msg = sprintf('Negative value(s) of area exceeding tolerance %g m2 (fixTinyNegs_tol):', ...
        fixTinyNegs_tol_m2);
    non_ntrl_bad = false ;
    for v = 1:Nlu
        this_YX = in_YXv(:,:,v) ;
        thisMin = min(this_YX(:)) ;
        if thisMin >= -fixTinyNegs_tol_m2
            continue
        end
        if length(LUnames) ~= Nlu
            error('length(LUnames) ~= Nlu')
        end
        thisLU = LUnames{v} ;
        if ~strcmp(thisLU, 'NATURAL')
            non_ntrl_bad = true ;
        end
        msg = sprintf('%s\n    %s min\t%g', ...
            msg, thisLU, thisMin) ;
    end
    if min(in_YXv(:)) < -conserv_tol_area
        error(msg)
    else
        msg = sprintf('%s\nStill within conserv_tol_area (%g)', ...
            msg, conserv_tol_area) ;
        if non_ntrl_bad
            msg = sprintf('%s, %s', ...
                msg, ...
                'but this is only expected for NATURAL. Figure out why this happened!') ;
            error(msg) %#ok<*SPERR> 
        else
             msg = sprintf('%s, %s', ...
                msg, ...
                'and it''s only for NATURAL, so zeroing.') ;
            warning(msg) %#ok<*SPWRN> 
        end
    end

end

out_YXv = in_YXv ;

isBad_YX = any(out_YXv<0,3) ;

if do_debug
    dbI = debugIJ_2deg(1) ;
    dbJ = debugIJ_2deg(2) ;
    isdb_YX = false(size(in_YXv, 1:2)) ;
    isdb_YX(dbI,dbJ) = true ;
    
    fprintf('PLUMharm_fixTinyNegs() for (%d,%d):\n', dbI, dbJ)
end

if any(any(isBad_YX))
    Nbad = length(find(isBad_YX)) ;
    isBad_YXv = repmat(isBad_YX,[1 1 Nlu]) ;
    
    % Get garrays
    tmp = out_YXv(isBad_YXv) ;
    if do_debug
        isdb = isdb_YX(isBad_YX) ;
        I_db = find(isdb) ;
    end
    tmp_landArea = landArea_YXv(isBad_YXv) ;
    tmp = tmp ./ tmp_landArea ;
    if any(isnan(tmp))
        error('How do you have NaN in tmp?')
    end
    tmp_xv = reshape(tmp,[Nbad Nlu]) ;
    
    % Report initial areas in debug cell
    if do_debug
        tmp_db = tmp_xv(I_db,:) ;
        fprintf(['   Before: ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], tmp_db) ;
        tmp_db_before = tmp_db ;
    end
    
    % Set negatives to zero
    tmp_xv(tmp_xv<0) = 0 ;
    tmp_xSum = sum(tmp_xv,2) ;
    
    % Report updated areas in debug cell
    if do_debug
        tmp_db = tmp_xv(I_db,:) ;
        fprintf(['        0: ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], tmp_db) ;
    end
    
    j = 0 ;
    while max(abs(tmp_xSum-1)) > 3*eps
        j = j+1 ;
        if j>50
            error('Possible infinite loop in fixing tiny negative areas!')
        end
        tmp_xv = tmp_xv ./ repmat(tmp_xSum,[1 Nlu]) ;
        tmp_xSum = sum(tmp_xv,2) ;
        
        if do_debug
            tmp_db = tmp_xv(I_db,:) ;
            fprintf(['         %d: ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], j, tmp_db) ;
        end
    end
    
    % Report updated areas in debug cell
    if do_debug
        tmp_db = tmp_xv(I_db,:) ;
        fprintf(['    After: ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], tmp_db) ;
        fprintf(['     Diff: ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], tmp_db - tmp_db_before) ;
    end
    
    out_YXv(isBad_YXv) = tmp_xv(:) .* tmp_landArea ;
    
    % Where differences are very small, zero them out
    smallDiff_YXv = abs(out_YXv - in_YXv) < 1e-6 ;
    out_YXv(smallDiff_YXv) = in_YXv(smallDiff_YXv) ;

    if do_debug
        fprintf([' in_YXv(%d,%d,:) = ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], ...
            dbI, dbJ, shiftdim(in_YXv(dbI,dbJ,:)))
        fprintf(['out_YXv(%d,%d,:) = ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], ...
            dbI, dbJ, shiftdim(out_YXv(dbI,dbJ,:)))
        fprintf(['             diff = ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], ...
            shiftdim(out_YXv(dbI,dbJ,:)) - shiftdim(in_YXv(dbI,dbJ,:)))
        disp(' ')
    end
end

% Small negatives are fine; just zero them out.
smallNeg_YXv = out_YXv<0 & out_YXv>-outPrec/2;
if any(any(any(smallNeg_YXv)))
    total_negs = sum(out_YXv(smallNeg_YXv)) ;
	 if abs(total_negs) >= 1
        warning('Zeroing out %0.1g total small negatives', total_negs)
    end
	 out_YXv(smallNeg_YXv) = 0 ;
end


end
