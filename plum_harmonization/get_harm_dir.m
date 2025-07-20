function harmDir = get_harm_dir(inDir, fruitveg_sugar_2oil, combineCrops, allow_unveg)
% Append option strings to directory name as needed.

harmDir = inDir ;

if fruitveg_sugar_2oil
    harmDir = [harmDir '.fvs2oil'] ;
end
if combineCrops
    harmDir = [harmDir '.combineCrops'] ;
end
if allow_unveg
    harmDir = [harmDir '.allow_unveg'] ;
end

end