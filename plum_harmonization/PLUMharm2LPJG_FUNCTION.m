function PLUMharm2LPJG_FUNCTION(varargin)

% Directories for harmonized PLUM outputs
dirList = {...
%               'SSP1.v10.s1.harm' ;
%               'SSP3.v10.s1.harm' ;
%               'SSP4.v10.s1.harm' ;
%               'SSP5.v10.s1.harm' ;
              'SSP1.v12.s1.harm' ;
              'SSP3.v12.s1.harm' ;
              'SSP4.v12.s1.harm' ;
              'SSP5.v12.s1.harm' ;
              } ;
          
% Replace dirList with input, if provided
if ~isempty(varargin)
    dirList = varargin{1} ;
    if length(varargin) > 1
        error('At most one optional argument (dirList) can be provided.')
    end
end
if ~iscellstr(dirList)
    error('dirList must be a cell array of strings!')
end

% Determine which system you're on
tmp = pwd ;
if strcmp(tmp(1:5),'/User')
    onMac = true ;
elseif strcmp(tmp(1:5),'/pfs/')
    onMac = false ;
else
    error('What system are you on?')
end
clear tmp
if onMac
    addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work/')
    landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
    lpjg_in_file = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851/yield.out.gz' ;
else
    addpath(genpath('/home/fh1-project-lpjgpi/lr8247/paper02-matlab-work')) ;
    addpath(genpath('/home/fh1-project-lpjgpi/lr8247/matlab-general/'))
    addpath(genpath('/home/fh1-project-lpjgpi/lr8247/matlab-general-fromshared/'))
    addpath(genpath('/home/fh1-project-lpjgpi/lr8247/lpj-guess-crop-calibration/'))
    landarea_file = '/home/fh1-project-lpjgpi/lr8247/PLUM/input/LUH2/supporting/staticData_quarterdeg.nc' ;
    dirList = strcat('/home/fh1-project-lpjgpi/lr8247/PLUM/input/',dirList) ;
    lpjg_in_file = '/home/fh1-project-lpjgpi/lr8247/PLUM/trunk_runs_out/LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851/yield.out.gz' ;
end

% Do processing
PLUMharm2LPJG


end
