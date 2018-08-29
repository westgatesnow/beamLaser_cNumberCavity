%cNumberCavity
%Initialization
clear; close all; 
%clc;

%Get parameters
getParam_beamLaser;
projectdir1 = controlType;
cd(projectdir1);
projectdir2 = name;
cd(projectdir2);

%Load data.
fprintf('Loading data... This may take several minutes.\n')
nAtom = load('nAtom.dat');
intensity = load('intensity.dat');
inversionAve = load('inversionAve.dat');
qMatrix = load('qMatrix.dat');
pMatrix = load('pMatrix.dat');
spinSpinCorAve_re = load('spinSpinCorAve_re.dat');
spinSpinCorAve_im = load('spinSpinCorAve_im.dat');
spinSpinCor_re = load('spinSpinCor_re.dat');
spinSpinCor_im = load('spinSpinCor_im.dat');
sxFinal = load('sxFinal.dat');
syFinal = load('syFinal.dat');
szFinal = load('szFinal.dat');
sxMatrix = load('sxMatrix.dat');
syMatrix = load('syMatrix.dat');
szMatrix = load('szMatrix.dat');

nTimeStep = tmax/dt;

%Add the previous directory into path
addpath ~/Desktop/codes/beamLaser_Proj/cNumberCavity/;