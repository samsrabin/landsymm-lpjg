function file_suffix = ...
    CFT_timeseries_plot(...
    cell_bl, cell_yr, cell_fao, ...
    yearList_baseline, yearList_future, yearList_fao, rebase, ...
    CFTnames, units, title_prefix, legend_in, do_caps, skip3rdColor, ...
    varargin)

is1x1num = @(x) numel(x)==1 & isnumeric(x) ;
is1x1int = @(x) numel(x)==1 & isint(x) ;
is1x1logical = @(x) numel(x)==1 & (x==0 | x==1) ;
iscellorempty = @(x) iscell(x) | isempty(x) ;

p = inputParser ;
addRequired(p,'cell_bl',iscellorempty) ;
addRequired(p,'cell_yr',@iscell) ;
addRequired(p,'cell_fao',iscellorempty) ;
addRequired(p,'yearList_baseline',@isnumeric) ;
addRequired(p,'yearList_future',@isnumeric) ;
addRequired(p,'yearList_fao',@isnumeric) ;
addRequired(p,'rebase',is1x1logical) ;
addRequired(p,'CFTnames',@iscellstr) ;
addRequired(p,'units',@ischar) ;
addRequired(p,'title_prefix',@ischar) ;
addRequired(p,'legend_in',@iscellstr) ;
addRequired(p,'do_caps',is1x1int) ;
addRequired(p,'skip3rdColor',is1x1logical) ;
addParameter(p,'conv_fact',1,is1x1num) ;
addParameter(p,'fontSize',14,is1x1num) ;
addParameter(p,'spacing',0.1,@isnumeric) ;   % [vert, horz]
addParameter(p,'ignYrs',0,is1x1int) ;
addParameter(p,'Nsmth',1,is1x1int) ;
addParameter(p,'lineWidth',2,is1x1num) ;
addParameter(p,'figure_pos',figurePos,@isnumeric) ;
addParameter(p,'legend_loc','NorthWest',@ischar) ;
addParameter(p,'dashed_zero_line',false,is1x1logical) ;
addParameter(p,'plum_area_adjustment',1,@isnumeric) ;
addParameter(p,'lpjg_area_adjustment',1,@isnumeric) ;
addParameter(p,'fao_linestyle','--',@isstr) ;

% Parse inputs
parse(p,...
    cell_bl, cell_yr, cell_fao, yearList_baseline, yearList_future, yearList_fao, rebase, ...
    CFTnames, units, title_prefix, legend_in, do_caps, skip3rdColor, ...
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
if ~isempty(cell_bl) && length(cell_bl) ~= length(cell_yr)
    error('length(cell_bl) ~= length(cell_yr)')
else
    Nplots = length(cell_yr) ;
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

if numel(lpjg_area_adjustment)==1
    lpjg_area_adjustment = repmat(lpjg_area_adjustment,size(yearList_baseline)) ;
end


% Make figure
figure('Position',figure_pos,'Color','w') ;
for p = 1:Nplots
    if ~isempty(cell_bl)
        this_cell_bl = cell_bl{p} ;
    else
        this_cell_bl = [] ;
    end
    units_toUse = units ;
    [ts_future, title_suffix, file_suffix] = ...
        rebase_future2baseline(rebase, Nsmth, this_cell_bl, cell_yr{p}, ignYrs, yearList_future) ;
    subplot_tight(ny,nx,p,spacing)
    h_list = [] ;
    if ~isempty(cell_bl)
        h1 = plot(yearList_baseline,...
            movmean((cell_bl{p}*conv_fact).*lpjg_area_adjustment,...
            Nsmth),'-k','LineWidth',lineWidth) ;
        h_list = [h_list ; h1(1)] ;
        hold on
    end
    if ~isempty(cell_fao)
        [~,IA] = intersect(yearList_baseline,yearList_fao) ;
        h2 = plot(yearList_fao,movmean((cell_fao{p}*conv_fact).*lpjg_area_adjustment(IA),Nsmth),[fao_linestyle 'k'],'LineWidth',lineWidth) ;
        h_list = [h_list ; h2(1)] ;
        if isempty(cell_bl)
            hold on
        end
    end
    if skip3rdColor
        set(gca,'ColorOrderIndex',1) ;
        h3 = plot(yearList_future,conv_fact*ts_future(:,1),'LineWidth',lineWidth) ;
        h4 = plot(yearList_future,conv_fact*ts_future(:,2),'LineWidth',lineWidth) ;
        set(h4,'Color',[255 159 56]/255) ;
        set(gca,'ColorOrderIndex',4) ;
        h5 = plot(yearList_future,conv_fact*ts_future(:,3),'LineWidth',lineWidth) ;
        h_list = [h_list ; h3 ; h4 ; h5] ;
    else
        h3 = plot(yearList_future,ts_future*conv_fact*plum_area_adjustment,'LineWidth',lineWidth) ;
        h_list = [h_list ; h3] ;
    end
    % Specify handles for legend, as described at http://www.mathworks.com/matlabcentral/answers/127613#comment_211148
    legend(h_list,legend_in,'Location',legend_loc,'AutoUpdate','off') ;
    if ~isempty(cell_bl) || ~isempty(cell_fao)
        hold off
    end
    if dashed_zero_line && min(get(gca,'YLim'))<0 && max(get(gca,'YLim'))>0
        hold on
        xlims = get(gca,'XLim') ;
        plot(xlims,[0 0],'--k') ;
        hold off
    end
    set(gca,'FontSize',fontSize)
    xlabel('Year')
    if plum_area_adjustment ~= 1
%         title_suffix = [title_suffix ' (PLUM\times' num2str(plum_area_adjustment) ')'] ;
        units_toUse = [units_toUse ' (PLUM\times' num2str(plum_area_adjustment) ')'] ;
        file_suffix = [file_suffix '_plumAdj' num2str(plum_area_adjustment)] ;
    end
    if lpjg_area_adjustment ~= 1
%         title_suffix = [title_suffix ' (blAdj)'] ;
        units_toUse = [units_toUse ' (blAdj)'] ;
        file_suffix = [file_suffix '_blAdj'] ;
    end
    ylabel(units_toUse)
    ht = title([title_prefix ': ' CFTnames{p} title_suffix]) ;
    letterlabel_align0(char(p + 64),ht,do_caps) ;
end




end
