%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make time series of livestock demand and crop usage %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisDir = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP5.v12.s1' ;
do_save = true ;


%% Setup

% Get (and create) output directory
out_dir = sprintf('%s/figures', thisDir) ;
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

% Get run identifier, stripping periods to avoid issues with La
[~, name, ext] = fileparts(thisDir) ;
run_id = [name ext] ;
run_id = strrep(run_id,'.','') ;


%% Import

demand = readtable(sprintf('%s/demand.txt', thisDir)) ;
fbs = readtable(sprintf('%s/fbs.txt', thisDir)) ;


%% Make figure: Livestock demand, crop usae

% Options %%%%%%%%%%%%%%%%%%%%
lineWidth = 3 ;
fontSize = 14 ;
thisPos = [1    33   720   772] ;
spacing = [0.075 0.08] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',thisPos,'Color','w') ;

subplot_tight(2,1,1,spacing)
plot(unique(demand.Year), ...
    [demand.Amount_Mt_(strcmp(demand.Commodity,'ruminants')) ...
    demand.Amount_Mt_(strcmp(demand.Commodity,'monogastrics'))], ...
    'LineWidth', lineWidth) ;
legend({'Ruminants', 'Monogastrics'},'Location','Best')
% xlabel('Year')
ylabel('Demand (Mt)')
title('Livestock demand')
set(gca,'FontSize',fontSize)

subplot_tight(2,1,2,spacing)
tmp_rfeed = nansum(reshape(fbs.RuminantsFeed,[9 91]),1)' ;
tmp_mfeed = nansum(reshape(fbs.MonogastricsFeed,[9 91]),1)' ;
tmp_foodEn = nansum(reshape(fbs.FoodAnd1stGen,[9 91]),1)' ;
plot(unique(demand.Year), ...
    [tmp_rfeed ...
    tmp_mfeed ...
    tmp_foodEn], ...
    'LineWidth', lineWidth) ;
legend({'Ruminants feed', 'Monogastrics feed', 'Food & Biofuels'},'Location','Best')
% xlabel('Year')
ylabel('Usage (Mt)')
title('Crop usage')
set(gca,'FontSize',fontSize)

if do_save
    disp('Saving...')
    out_file = sprintf('%s/livestock_and_feed_%s.pdf', out_dir, run_id) ;
    export_fig(out_file, '-r150')
    disp('Done.')
end


%% Make figure: Crop usage as fraction of livestock demand

% Options %%%%%%%%%%%%%%%%%%%%
lineWidth = 3 ;
fontSize = 14 ;
thisPos = [337   293   720   296] ;
spacing = [0.075 0.08] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',thisPos,'Color','w') ;

plot(unique(demand.Year), ...
    100*tmp_rfeed./demand.Amount_Mt_(strcmp(demand.Commodity,'ruminants')), ...
    'LineWidth', lineWidth)
% xlabel('Year')
ylabel('Percentage (%)')
title('Share of ruminant food provided by feed crops (SSP5-85)')
set(gca,'FontSize',fontSize)

if do_save
    disp('Saving...')
    out_file = sprintf('%s/ruminant_feed_share_%s.pdf', out_dir, run_id) ;
    export_fig(out_file, '-r150')
    disp('Done.')
end



%% Save figure





