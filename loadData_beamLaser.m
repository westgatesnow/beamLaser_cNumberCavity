%Initialization
clear; close all; clc;

%Get parameters
getParam_beamLaser;
projectdir = name;
cd(projectdir);

%Load data.
fprintf('Loading data... This may take several minutes.\n')
intensity = load('intensity.dat');
inversionAve = load('inversionAve.dat');
sxFinal = load('sxFinal.dat');
syFinal = load('syFinal.dat');
szFinal = load('szFinal.dat');
qMatrix = load('qMatrix.dat');
pMatrix = load('pMatrix.dat');
nTimeStep = tmax/dt;

%Add the previous directory into path
addpath ~/Desktop/codes/beamLaser_Proj/cNumberCavity/;