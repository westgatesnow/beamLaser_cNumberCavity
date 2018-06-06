%cNumberCavity
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
% spinSpinCorAve_re = load('spinSpinCorAve_re.dat');
% spinSpinCorAve_im = load('spinSpinCorAve_im.dat');
sxFinal = load('sxFinal.dat');
syFinal = load('syFinal.dat');
szFinal = load('szFinal.dat');
qMatrix = load('qMatrix.dat');
pMatrix = load('pMatrix.dat');
% spinSpinCor_re = load('spinSpinCor_re.dat');
% spinSpinCor_im = load('spinSpinCor_im.dat');
nTimeStep = tmax/dt;
% NBIN = size(spinSpinCor_re,1);

%Add the previous directory into path
addpath ~/Desktop/codes/beamLaser_Proj/cNumberCavity/;