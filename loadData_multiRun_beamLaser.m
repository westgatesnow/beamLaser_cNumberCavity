%cNumberCavity
%This file is used to load the data from running "multirun.sh"

%Initialization
clear; close all; clc;
addpath ~/Desktop/codes/beamLaser_Proj/cNumberCavity/;

%bash variables; 
%tau
nMaxTau = 100;
initTau = 0.1;
intervalTau = 0.1;
%nAtom
nMaxNAtom = 1;
initNAtom = 100;
intervalNAtom = 100;

%get parameters from the input.txt in the parent directory
%all the parameters are universal except for transitTime, density, and name
getParam_beamLaser;
nTimeStep = tmax/dt;

%define empty data structure for variables
intensity = zeros(nMaxTau,nMaxNAtom,nstore);
inversionAve = zeros(nMaxTau,nMaxNAtom,nstore);
spinSpinCorAve = zeros(nMaxTau,nMaxNAtom,nstore);
szFinal = zeros(nMaxTau,nMaxNAtom,nTimeStep);
qMatrix = zeros(nMaxTau,nMaxNAtom,nTrajectory,nstore);
pMatrix = zeros(nMaxTau,nMaxNAtom,nTrajectory,nstore);

%load all the data needed
timeRequired = 3.2*nMaxTau*nMaxNAtom;
printWords1 = ['Begin loading. Time required: ' , num2str(timeRequired), ' seconds.'];
disp(printWords1);
tic;
for i = 1:nMaxTau
    for j = 1:nMaxNAtom
        %define the name of the directory
        tau = initTau+intervalTau*(i-1);
        nAtom = initNAtom+intervalNAtom*(j-1);
        filename = ['pois_tau0', num2str(tau,'%.1f'), '_g2_k40_N', num2str(nAtom)];
        cd(filename);
        %Load data.
        intensity(i,j,:) = load('intensity.dat');
        inversionAve(i,j,:) = load('inversionAve.dat');
        %spinSpinCorAve(i,j,:) = load('spinSpinCorAve_re.dat');
        %spinSpinCorAve_im = load('spinSpinCorAve_im.dat');
        %sxFinal = load('sxFinal.dat');
        %syFinal = load('syFinal.dat');
        szFinal(i,j,:) = load('szFinal.dat');
        qMatrix(i,j,:,:) = load('qMatrix.dat');
        pMatrix(i,j,:,:) = load('pMatrix.dat');
        %spinSpinCor_re = load('spinSpinCor_re.dat');
        %spinSpinCor_im = load('spinSpinCor_im.dat');  
        cd ..;
    end
    %Print out info
    percentage = i/nMaxTau*100;
    if mod(percentage,10) == 0
        printWords2 = ['Finish loading ', num2str(percentage),'%.'];
        disp(printWords2);
    end
end
toc;