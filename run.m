%run this m file only once
Path1= cd;  Path1=[Path1,'\data1'];path(Path1,path);
Path1 = cd; Path1 = [Path1, '\MainRoutine']; path(Path1, path);
Path1 = cd; Path1 = [Path1, '\HoughTransform']; path(Path1, path);
Path1 = cd; Path1 = [Path1, '\CorrectShape']; path(Path1, path);
Path1 = cd; Path1 = [Path1, '\LevelSet']; path(Path1, path);
%Path1 = cd; Path1 = [Path1, '\createVideo']; path(Path1, path);


mex LevelSet\TwoDlevelSet.cpp

Main.m