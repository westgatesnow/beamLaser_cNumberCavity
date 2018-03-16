%Import simulation parameters
dataTable = readtable('input.txt');
dataVariables = dataTable.Variables;
%Define all the simulation parameters one by one;
%   Should change the parameters according to different parameter setups.
dt = str2double(dataVariables{1,2});
tmax = str2double(dataVariables{2,2});
nTrajectory = str2double(dataVariables{3,2});
nstore = str2double(dataVariables{4,2});
yWall = str2double(dataVariables{5,2});
sigmaXX = str2double(dataVariables{6,2});
sigmaXZ = str2double(dataVariables{7,2});
transitTime = str2double(dataVariables{8,2});
sigmaPX = str2double(dataVariables{9,2});
sigmaPY = str2double(dataVariables{10,2});
sigmaPZ = str2double(dataVariables{11,2});
density = str2double(dataVariables{12,2});
rabi = str2double(dataVariables{13,2});
kappa = str2double(dataVariables{14,2});
name = dataVariables{15,2};



