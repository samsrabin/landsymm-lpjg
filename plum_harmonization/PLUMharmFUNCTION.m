function PLUMharmFUNCTION(varargin)

PLUM_in_toptop = {...
%                   'SSP1.v10.s1' ;
%                   'SSP3.v10.s1' ;
%                   'SSP4.v10.s1' ;
%                   'SSP5.v10.s1' ;
%                   'SSP1.v11.s1' ;
%                   'SSP3.v11.s1' ;
%                   'SSP4.v11.s1' ;
%                   'SSP5.v11.s1' ;
%                   'SSP1.v11.s1test' ;
%                   'ssp11/SSP5/s3' ;
                  'SSP1.v12.s1' ;
                  'SSP3.v12.s1' ;
                  'SSP4.v12.s1' ;
                  'SSP5.v12.s1' ;
                  } ;

% Replace PLUM_in_toptop with input, if provided
if ~isempty(varargin)
    PLUM_in_toptop = varargin{1} ;
    if length(varargin) > 1
        error('At most one optional argument (PLUM_in_toptop) can be provided.')
    end
end
if ~iscellstr(PLUM_in_toptop)
    error('PLUM_in_toptop must be a cell array of strings!')
end

% Determine which system you're on
tmp = pwd ;
if strcmp(tmp(1:5),'/User')
    onMac = true ;
elseif strcmp(tmp(1:5),'/pfs/') || strcmp(tmp(1:5),'/home')
    onMac = false ;
else
    error('What system are you on?')
end
clear tmp

if onMac
    PLUM_in_toptop = strcat('/Users/Shared/PLUM/PLUM_outputs_for_LPJG/',PLUM_in_toptop) ;
    addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;
    PLUMharm_top = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/plum_harmonization/' ;
    inDir_protectedAreas = '/Users/Shared/PLUM/input/protected_areas/' ;
else
    PLUM_in_toptop = strcat('/home/fh1-project-lpjgpi/lr8247/PLUM/input/PLUMouts_2011-2100/',PLUM_in_toptop) ;
    addpath(genpath( '/home/fh1-project-lpjgpi/lr8247/paper02-matlab-work')) ;
    PLUMharm_top = '/home/fh1-project-lpjgpi/lr8247/plum_harmonization/' ;
    inDir_protectedAreas = '/home/fh1-project-lpjgpi/lr8247/PLUM/input/protected_areas/' ;
    addpath(genpath('/home/fh1-project-lpjgpi/lr8247/matlab-general/'))
    addpath(genpath('/home/fh1-project-lpjgpi/lr8247/matlab-general-fromshared/'))
    addpath(genpath('/home/fh1-project-lpjgpi/lr8247/lpj-guess-crop-calibration/'))
end

% Do harmonization
PLUMharm


end