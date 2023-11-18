function harmDirs = PLUMharm_get_harmDirs(plumDirs, fruitveg_sugar_2oil, combineCrops)
% SSR 2023-11-18: UNTESTED!

harmDirs = cell(length(plumDirs), 1) ;
for d = 1:length(plumDirs)
    harmDir = [plumDir '.harm'] ;
    harmDir = get_harm_dir(harmDir, fruitveg_sugar_2oil, combineCrops) ;
    harmDirs{d} = harmDir ;
end

end