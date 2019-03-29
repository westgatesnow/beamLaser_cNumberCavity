%cNumberCavity
%Import simulation parameters
dataTable = readtable('input.txt');
dataVariables = dataTable.Variables;
%Define all the simulation parameters one by one;
%   Should change the parameters according to different parameter setups.
%   CAN rewrite this function with genvarname()????????
dt = str2double(dataVariables{1,2});
tmax = str2double(dataVariables{2,2});
nStore = str2double(dataVariables{3,2});
nTrajectory = str2double(dataVariables{4,2});
nBin = str2double(dataVariables{5,2});
yWall = str2double(dataVariables{6,2});
lambda = str2double(dataVariables{7,2});
deltaZ = str2double(dataVariables{8,2});
deltaPz = str2double(dataVariables{9,2});
transitTime = str2double(dataVariables{10,2});
density = str2double(dataVariables{11,2});
rabi = str2double(dataVariables{12,2});
kappa = str2double(dataVariables{13,2});
invT2 = str2double(dataVariables{14,2});
controlType = dataVariables{15,2};
name = dataVariables{16,2};
% pois = dataVariables{17,2};



