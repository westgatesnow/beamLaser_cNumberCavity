%cNumberCavity
%This file is used to load the data from running "multirun.sh"

%Initialization
clear; close all; clc;
addpath ~/Desktop/codes/beamLaser_Proj/cNumberCavity/;

%get parameters from the input.txt in the parent directory
%all the parameters are universal except for transitTime, density, and name
getParam_beamLaser;
nTimeStep = tmax/dt;
gc = rabi^2/kappa;

%bash variables; 

%tau
%tauList1 = 0.05 ~ 0.95;
nMaxTau1 = 19;%19
initTau1 = 0.05;%0.05/0.60
intervalTau1 = 0.05;
tauList1 = (initTau1+(0:nMaxTau1-1)*intervalTau1)*gc;
%tauList = 1:nMaxTau;%in the unit of tau, if initTau = intervalTau;
%tauList2 = 1.0 ~ 10.5
nMaxTau2 = 20;%20/0
initTau2 = 1.0;
intervalTau2 = 0.5;
tauList2 = (initTau2+(0:nMaxTau2-1)*intervalTau2)*gc;
%for convenience
nMaxTau = nMaxTau1+nMaxTau2;
tauList = [tauList1,tauList2];

%nAtomAve
nMaxNAtom = 1;%10
initNAtom = 100;
intervalNAtom = 100;
nAtomAveList =initNAtom+(0:nMaxNAtom-1)*intervalNAtom;

%define empty data structure for variables
nAtom = zeros(nMaxTau,nMaxNAtom,nstore);
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
for j = 1:nMaxNAtom
    for i = 1:nMaxTau1
        %define the name of the directory
        tau1 = initTau1+intervalTau1*(i-1);
        nAtomAve = initNAtom+intervalNAtom*(j-1);
        filename = ['pois_tau', num2str(tau1,'%.2f'), '_g2_k40_N', num2str(nAtomAve)];
        cd(filename);
        %Load data.
        nAtom(i,j,:) = load('nAtom.dat');
        intensity(i,j,:) = load('intensity.dat');
        inversionAve(i,j,:) = load('inversionAve.dat');
        spinSpinCorAve(i,j,:) = load('spinSpinCorAve_re.dat');
        %spinSpinCorAve_im = load('spinSpinCorAve_im.dat');
        %sxFinal = load('sxFinal.dat');
        %syFinal = load('syFinal.dat');
        szFinal(i,j,:) = load('szFinal.dat');
        qMatrix(i,j,:,:) = load('qMatrix.dat');
        pMatrix(i,j,:,:) = load('pMatrix.dat');
        %spinSpinCor_re = load('spinSpinCor_re.dat');
        %spinSpinCor_im = load('spinSpinCor_im.dat');  
        cd ..;
            %Print out info
        percentage = (i+(j-1)*nMaxTau)/nMaxTau/nMaxNAtom*100;
        printWords2 = ['Finish loading ', num2str(percentage,'%.1f'),'%...'];
        disp(printWords2);
    end
    for i = 1:nMaxTau2
        %define the name of the directory
        tau2 = initTau2+intervalTau2*(i-1);
        nAtomAve = initNAtom+intervalNAtom*(j-1);
        filename = ['pois_tau', num2str(tau2,'%.1f'), '_g2_k40_N', num2str(nAtomAve)];
        cd(filename);
        %Load data.
        nAtom(i+nMaxTau1,j,:) = load('nAtom.dat');
        intensity(i+nMaxTau1,j,:) = load('intensity.dat');
        inversionAve(i+nMaxTau1,j,:) = load('inversionAve.dat');
        spinSpinCorAve(i+nMaxTau1,j,:) = load('spinSpinCorAve_re.dat');
        %spinSpinCorAve_im = load('spinSpinCorAve_im.dat');
        %sxFinal = load('sxFinal.dat');
        %syFinal = load('syFinal.dat');
        szFinal(i+nMaxTau1,j,:) = load('szFinal.dat');
        qMatrix(i+nMaxTau1,j,:,:) = load('qMatrix.dat');
        pMatrix(i+nMaxTau1,j,:,:) = load('pMatrix.dat');
        %spinSpinCor_re = load('spinSpinCor_re.dat');
        %spinSpinCor_im = load('spinSpinCor_im.dat');  
        cd ..;
            %Print out info
        percentage = (i+nMaxTau1+(j-1)*nMaxTau)/nMaxTau/nMaxNAtom*100;
        printWords2 = ['Finish loading ', num2str(percentage,'%.1f'),'%...'];
        disp(printWords2);
    end
end
toc;