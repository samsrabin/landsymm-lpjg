function harmDirs = PLUMharm_get_harmDirs(plumDirs, fruitveg_sugar_2oil, combineCrops, allow_unveg)

harmDirs = cell(length(plumDirs), 1) ;
for d = 1:length(plumDirs)
    plumDir = plumDirs{d} ;
    harmDir = [plumDir '.harm'] ;
    harmDir = get_harm_dir(harmDir, fruitveg_sugar_2oil, combineCrops, allow_unveg) ;
    harmDirs{d} = harmDir ;
end

end