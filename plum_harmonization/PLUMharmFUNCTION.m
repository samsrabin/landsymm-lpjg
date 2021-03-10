function PLUMharmFUNCTION(varargin)

% dirList = {...
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
%                   'SSP1.v12.s1' ;
%                   'SSP3.v12.s1' ;
%                   'SSP4.v12.s1' ;
%                   'SSP5.v12.s1' ;
%                   } ;
% base_year = 2010 ;
% year1 = 2011 ;
% yearN = 2100 ;
% fruitveg_sugar_2oil = false ;
% % baseline_ver = 1 ;
% % baseline_ver = 2 ;   % Based on remap_v6
% baseline_ver = 3 ;   % Based on remap_v6p7

% dirList = {...
%                   'halfearth/HEoct/baseline/s1';
%                   'halfearth/HEoct/halfearth/s1';
%                   } ;
% base_year = 2010 ;
% year1 = 2011 ;
% yearN = 2060 ;
% fruitveg_sugar_2oil = false ;
% baseline_ver = 4 ;   % Based on remap_v8b(2oil, if fruitveg_sugar_2oil true)

dirList = {...
%     'SSP1/s1' ;
    'SSP2/s1' ;
    'SSP3/s1' ;
    'SSP4/s1' ;
    'SSP5/s1' ;
    } ;
base_year = 2010 ;
year1 = 2011 ;
yearN = 2100 ;
fruitveg_sugar_2oil = false ;
baseline_ver = 4 ;   % Based on remap_v8b(2oil, if fruitveg_sugar_2oil true)

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

% Determine which system you're on and set up
thisSystem = get_system_name() ;
if strcmp(thisSystem, 'ssr_mac')
    addpath(genpath('/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/MATLAB_work')) ;
    plumharm_repo_path = '/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/plum_harmonization/' ;
elseif strcmp(thisSystem, 'ssr_keal')
    addpath(genpath('/pd/data/lpj/sam/paper02-matlab-work')) ;
    plumharm_repo_path = '/pd/data/lpj/sam/PLUM/plum_harmonization' ;
else
    error('thisSystem not recognized: %s', thisSystem)
end
addpath(genpath(plumharm_repo_path))

% Do harmonization
PLUMharm


end