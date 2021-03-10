function PLUMharm2LPJG_FUNCTION(varargin)

% Directories for harmonized PLUM outputs
dirList = {...
%               'SSP1.v10.s1.harm' ;
%               'SSP3.v10.s1.harm' ;
%               'SSP4.v10.s1.harm' ;
%               'SSP5.v10.s1.harm' ;
%               'SSP1.v12.s1.harm' ;
%               'SSP3.v12.s1.harm' ;
%               'SSP4.v12.s1.harm' ;
%               'SSP5.v12.s1.harm' ;
%     'SSP1/s1.harm' ;
    'SSP2/s1.harm' ;
    'SSP3/s1.harm' ;
    'SSP4/s1.harm' ;
    'SSP5/s1.harm' ;
              } ;
base_year = 2010 ;
y1 = 2011 ;
yN = 2100 ;
yStep = 1 ;

% dirList = {...
%           'halfearth/HEoct/baseline/s1.harm';
%           'halfearth/HEoct/halfearth/s1.harm';
%           } ;
% base_year = 2010 ; %#ok<*NASGU>
% y1 = 2011 ;
% yN = 2060 ;
% yStep = 1 ;

do_gzip = true ;
          
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

% Determine which system you're on and set up.
thisSystem = get_system_name() ;
if strcmp(thisSystem, 'ssr_mac')
    addpath(genpath('/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/MATLAB_work')) ;
    plumharm_repo_path = '/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/plum_harmonization' ;
elseif strcmp(thisSystem, 'ssr_keal')
    addpath(genpath('/pd/data/lpj/sam/paper02-matlab-work')) ;
    plumharm_repo_path = '/pd/data/lpj/sam/PLUM/plum_harmonization' ;
else
    error('thisSystem not recognized: %s', thisSystem)
end
addpath(genpath(plumharm_repo_path))

% Do processing
PLUMharm2LPJG


end
