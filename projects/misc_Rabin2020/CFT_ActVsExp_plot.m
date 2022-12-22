function CFT_ActVsExp_plot(...
    cell_yr_act, cell_yr_exp, ...
    yearList_future, ...
    CFTnames, units, title_prefix, stdLegend, do_caps, ...
    varargin)

is1x1num = @(x) numel(x)==1 & isnumeric(x) ;
is1x1int = @(x) numel(x)==1 & isint(x) ;
is1x1logical = @(x) numel(x)==1 & (x==0 | x==1) ;

p = inputParser ;
addRequired(p,'cell_yr_act',@iscell) ;
addRequired(p,'cell_yr_exp',@iscell) ;
addRequired(p,'yearList_future',@isnumeric) ;
addRequired(p,'CFTnames',@iscellstr) ;
addRequired(p,'units',@ischar) ;
addRequired(p,'title_prefix',@ischar) ;
addRequired(p,'stdLegend',@iscellstr) ;
addParameter(p,'conv_fact',1,is1x1num) ;
addParameter(p,'fontSize',14,is1x1num) ;
addParameter(p,'spacing',0.1,@isnumeric) ;   % [vert, horz]
addParameter(p,'Nsmth',1,is1x1int) ;
addParameter(p,'lineWidth',2,is1x1num) ;
addParameter(p,'figure_pos',figurePos,@isnumeric) ;

% Parse inputs
parse(p,...
    cell_yr_act, cell_yr_exp, yearList_future, ...
    CFTnames, units, title_prefix, stdLegend, ...
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


% Setup
if length(cell_yr_act) ~= length(cell_yr_exp)
    error('length(cell_yr_act) ~= length(cell_yr_exp)')
else
    Nplots = length(cell_yr_act) ;
end
if Nplots <= 6
    ny = 2 ;
    nx = 3 ;
elseif Nplots <= 8
    ny = 2 ;
    nx = 4 ;
    elseif Nplots == 9
    ny = 3 ;
    nx = 3 ;
else
    error(['Set ny and nx for Nplots=' num2str(Nplots)]) ; 
end

% Avoid subscripting where you intend an underscore
CFTnames = strrep(CFTnames,'_','\_') ;


% Make figure
figure('Position',figure_pos,'Color','w') ;
for p = 1:Nplots
    subplot_tight(ny,nx,p,spacing)
    plot(yearList_future,movmean(100*(cell_yr_act{p}-cell_yr_exp{p})./cell_yr_exp{p},Nsmth),'LineWidth',lineWidth)
    hold on
    plot(get(gca,'XLim'),[0 0],'--k')
    hold off
    legend(stdLegend, ...
        'Location','Best') ;
    set(gca,'FontSize',fontSize)
    xlabel('Year')
    ylabel(units)
    ht = title([title_prefix ': ' CFTnames{p}]) ;
    letterlabel_align0(char(p + 64),ht,do_caps) ;
end




end
