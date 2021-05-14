function e2p_check_correct_zeros(data_gvt, which_file, varNames, ...
    bl_or_fu, getbasenamei, varargin)
% If which_file is yield, returns error if no positive cells found for one
% or more crops.
% If which_file is gsirrigation, returns error if no positive cells found
% for one or more irrigated crops, or if any positive cells are found for
% one or more rainfed crops.

warn_in_test_gt0_is_good = false ;
missing_N = {} ;
if ~isempty(varargin)
    if ~isempty(varargin(1))
        warn_in_test_gt0_is_good = varargin{1} & strcmp(which_file, 'gsirrigation') ;
    end
    if length(varargin) > 1
        missing_N = varargin{2} ;
        if length(varargin) > 2
            error('At most two optional arguments allowed (warn_in_test_gt0_is_good, missing_N)')
        end
    end
end

basenamesi = getbasenamei(varNames) ;

if strcmp(which_file, 'yield')
    test_gt0_is_good(data_gvt, warn_in_test_gt0_is_good, varNames, bl_or_fu)
elseif strcmp(which_file, 'gsirrigation')
    test_isirrig = @(x) strcmp(x(end),'i') ;
    isirrig = cellfun(test_isirrig, basenamesi) ;
    data_gvt_rf = data_gvt(:,~isirrig,:) ;
    data_gvt_ir = data_gvt(:,isirrig,:) ;
    test_gt0_is_bad(data_gvt_rf)
    test_gt0_is_good(data_gvt_ir, warn_in_test_gt0_is_good, varNames, bl_or_fu)
else
    error('which_file (%s) not recognized', which_file)
end


end


function test_gt0_is_good(data_gvt, do_warn, varNames, bl_or_fu)

for t = 1:size(data_gvt, 3)
    if any(~any(data_gvt(:,:,t)>0,1))
        Ifound = find(~any(data_gvt(:,:,t)>0,1)) ;
        Nfound = length(Ifound) ;
        msg = sprintf('%s: No IWD found for %d/%d variables at timestep %d', ...
            bl_or_fu, Nfound, size(data_gvt, 2), t) ;
        for ii = 1:Nfound
            thisI = Ifound(ii) ;
            msg = sprintf('%s\n    %d: %s', msg, thisI, varNames{thisI}) ;
        end
        if do_warn
            warning('%s\n  Make sure rainfed and irrigated yields are the same!!', msg)
            error(['You need to finish dealing with this situation. The problem is that ', ...
                'you might have irrigation in spring wheat but not winter wheat (or ', ...
                'vice versa). This means that you shouldn''t have included the missing ', ...
                'one in the max_wheat yield choice, but at the time you didn''t know ', ...
                'that, because you hadn''t loaded irrigation data yet'])
        else
            error('%s: No positive-yield cells found for %d/%d variables (first %d, %s) at timestep %d', ...
                bl_or_fu, Nfound, size(data_gvt, 2), Ifound(1), varNames{Ifound(1)}, t)
        end
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