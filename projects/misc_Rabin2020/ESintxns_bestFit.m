function hl = ESintxns_bestFit( ...
    xdata, ydata, bestFitType, ...
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
% isAlpha = @(x) is1x1num(x) & x>=0 & x<=1 ;

p = inputParser ;
addRequired(p,'xdata',@isnumeric) ;
addRequired(p,'ydata',@isnumeric) ;
addRequired(p,'bestFitType',@ischar) ;
addParameter(p,'lineStyle','-',@ischar) ;
addParameter(p,'lineWidth',2,is1x1int) ; % Size in one-d points
addParameter(p,'lineColor','r',iscolor) ;

% Parse inputs
parse(p,...
    xdata, ydata, bestFitType, ...
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

if any(isnan(xdata) | isnan(ydata))
    error('xdata and ydata must have no NaN values in order for best-fit calculations to work!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The actual function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on

if strcmp(bestFitType,'linear')
    
    xlim = [min(xdata) max(xdata)] ;
    coef_fit = polyfit(xdata, ydata, 1) ;
    y_fit = polyval(coef_fit, xlim) ;
    hl = plot(xlim, y_fit, ...
        'LineStyle', lineStyle, ...
        'LineWidth', lineWidth, ...
        'Color', lineColor) ;
    
    mdl = fitlm(xdata,ydata) ;
    bestfit_string = sprintf( ...
        'y = %0.3fx + %0.3f, r^2 = %0.3f, p = %0.3g', ...
        mdl.Coefficients.Estimate('x1'), ...
        mdl.Coefficients.Estimate('(Intercept)'), ...
        mdl.Rsquared.Ordinary, ...
        mdl.Coefficients.pValue('x1')) ;
    legend(hl,bestfit_string,'Location','SouthOutside')
    
elseif strcmp(bestFitType,'splines')
    
    ylims = get(gca,'YLim') ;
    
    [f,gof,out] = fit(xdata, ydata,'smoothingspline','SmoothingParam',1e-6/length(xdata));
    hl = plot(f) ;
    set(hl, ...
        'LineStyle', lineStyle, ...
        'LineWidth', lineWidth, ...
        'Color', lineColor) ;
    set(gca,'YLim',ylims)
    
else
    hold off
    error('bestFitType "%s" not recognized!', bestFitType) ;
end

hold off



end