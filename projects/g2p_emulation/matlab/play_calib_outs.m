%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% How do the calibration outputs look? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisModel = 'LPJmL' ;


%% Import

% Import
in_dir = '/Users/Shared/GGCMI2PLUM/emulator/Sam/outputs_calib' ;
filename = sprintf('%s/%s.out', in_dir, thisModel) ;
yield = lpjgu_matlab_readTable_then2map(filename, 'force_mat_save', false, 'force_mat_nosave', true) ;

% Get crop names
crop_List = {} ;
for c = 1:length(yield.varNames)
    
end


%% Plot

spacing = [0.025 0.025] ;

figure('Color','w','Position',figurePos) ;

for c = 1:12
    subplot_tight(3,4,c,spacing) ;
    pcolor(yield.maps_YXv(:,:,c)) ; shading flat
    axis equal tight off
    colorbar
    title(yield.varNames{c})
end