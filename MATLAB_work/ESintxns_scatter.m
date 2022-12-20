function hs = ESintxns_scatter( ...
    xdata, ydata, ...
    varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parsing arguments %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

is1x1num = @(x) numel(x)==1 & isnumeric(x) ;
is1x1int = @(x) numel(x)==1 & isint(x) ;
% is1x1logical = @(x) numel(x)==1 & (x==0 | x==1) ;
% iscellorempty = @(x) iscell(x) | isempty(x) ;
iscolor = @(x) (isnumeric(x) & isequal(size(shiftdim(x)),[1 3])) ...
    | (ischar(x) & length(x)==1) ;
isAlpha = @(x) is1x1num(x) & x>=0 & x<=1 ;

p = inputParser ;
addRequired(p,'xdata',@isnumeric) ;
addRequired(p,'ydata',@isnumeric) ;

% Points
addParameter(p,'marker','.',@ischar) ;
addParameter(p,'markerSize',30,is1x1int) ; % Size in one-d points
addParameter(p,'markerColor','b',iscolor) ;
addParameter(p,'markerAlpha',1,isAlpha) ;

% Parse inputs
parse(p,...
    xdata, ydata, ...
    varargin{:});
pFields = fieldnames(p.Results) ;
Nfields = length(pFields) ;
for f = 1:Nfields
    thisField = pFields{f} ;
    if ~exist(thisField,'var')
        eval([thisField ' = p.Results.' thisField ' ;']) ;
    end
    clear thisField
end ; clear f
clear p


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The actual function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

hs = scatter(xdata, ydata, ...
    markerSize^2, markerColor, marker, ...
    'MarkerEdgeColor', markerColor, ...
    'MarkerFaceColor', markerColor, ...
    'MarkerEdgeAlpha', markerAlpha, ...
    'MarkerFaceAlpha', markerAlpha) ;


end