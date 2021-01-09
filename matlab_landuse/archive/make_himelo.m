

in_file = '/Volumes/Reacher/G2P/inputs/LU/remaps_v7b/cropfracs.remapv7b.txt' ;
% in_file = '/Volumes/Reacher/G2P/inputs/LU/remaps_v7b/nfert.remapv7b.txt' ;



%% Do it

cd '/Users/sam/Documents/git_repos/g2p_emulation/matlab_landuse'
addpath(genpath(pwd))

del_crops = {'Miscanthus','Miscanthusi'} ;
        
Nlevels = [10 60 200 1000] ;
        

tic
removeCols_addZeros(in_file, del_crops, Nlevels) ;
disp(toc_hms(toc))
