%cNumberCavity
%This file is used to load the data from running "multirun_beamLaser.sh"

%Initialization
clear; close all; clc;
addpath ~/Desktop/codes/beamLaser_Proj/cNumberCavity/;

%get parameters from the input.txt in the controlType directory
prompt = "What is the controlType?\n"; 
controlType = input(prompt,'s');
cd(controlType);
getParam_beamLaser;
nTimeStep = tmax/dt;
gc = rabi^2/kappa;

%bash variables; 

%tau
% %tauList1 = 0.5 ~ 1.4;
% nMaxTau1 = 10;%
% initTau1 = 0.1;%
% intervalTau1 = 0.1;
% tauList1 = (initTau1+(0:nMaxTau1-1)*intervalTau1)*gc;
% %tauList = 1:nMaxTau;%in the unit of tau, if initTau = intervalTau;
% %tauList2 = 1.0 ~ 10.5
% nMaxTau2 = 0;%20/0
% initTau2 = 1.0;
% intervalTau2 = 0.5;
% tauList2 = (initTau2+(0:nMaxTau2-1)*intervalTau2)*gc;
% %for convenience
% nMaxTau = nMaxTau1+nMaxTau2;
% tauList = [tauList1,tauList2];
%Can use uniform distri in log scale??????

nMaxTau = 10;%
initTau = 0.1;%
intervalTau = 0.1;
tauList = (initTau+(0:nMaxTau-1)*intervalTau);%dim = 1*nMaxTau

%dens
nMaxDens = 10;%
initDens = 50;
intervalDens = 50;
densList =(initDens+(0:nMaxDens-1)*intervalDens);%dim = 1*densList

%define empty data structure for variables
nAtom = zeros(nMaxTau,nMaxDens,nStore);
intensity = zeros(nMaxTau,nMaxDens,nStore);
inversionAve = zeros(nMaxTau,nMaxDens,nStore);
spinSpinCorAve = zeros(nMaxTau,nMaxDens,nStore);
szFinal = zeros(nMaxTau,nMaxDens,nTimeStep);
qMatrix = zeros(nMaxTau,nMaxDens,nTrajectory,nStore);
pMatrix = zeros(nMaxTau,nMaxDens,nTrajectory,nStore);

%load all the data needed
timeRequired = 3.2*nMaxTau*nMaxDens;
printWords1 = ['Begin loading. Time required: ' , num2str(timeRequired), ' seconds.'];
disp(printWords1);
tic;
for j = 1:nMaxDens
    for i = 1:nMaxTau
        %define the name of the directory
        tau = initTau+intervalTau*(i-1);
        dens = initDens+intervalDens*(j-1);
        filename = ['pois_dt0.005_dZ0.5_dPz1.0_tau', ...
            num2str(tau,'%.1f'), '_nBin20_nAtom',num2str(dens),'_g3_k90_yWall1.0', ];
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
        percentage = (i+(j-1)*nMaxTau)/nMaxTau/nMaxDens*100;
        printWords2 = ['Finish loading ', num2str(percentage,'%.1f'),'%...'];
        disp(printWords2);
    end
%     for i = 1:nMaxTau2
%         %define the name of the directory
%         tau2 = initTau2+intervalTau2*(i-1);
%         dens = initDens+intervalDens*(j-1);
%         filename = ['pois_tau', num2str(tau2,'%.1f'), '_g2_k40_N', num2str(dens)];
%         cd(filename);
%         %Load data.
%         nAtom(i+nMaxTau1,j,:) = load('nAtom.dat');
%         intensity(i+nMaxTau1,j,:) = load('intensity.dat');
%         inversionAve(i+nMaxTau1,j,:) = load('inversionAve.dat');
%         spinSpinCorAve(i+nMaxTau1,j,:) = load('spinSpinCorAve_re.dat');
%         %spinSpinCorAve_im = load('spinSpinCorAve_im.dat');
%         %sxFinal = load('sxFinal.dat');
%         %syFinal = load('syFinal.dat');
%         szFinal(i+nMaxTau1,j,:) = load('szFinal.dat');
%         qMatrix(i+nMaxTau1,j,:,:) = load('qMatrix.dat');
%         pMatrix(i+nMaxTau1,j,:,:) = load('pMatrix.dat');
%         %spinSpinCor_re = load('spinSpinCor_re.dat');
%         %spinSpinCor_im = load('spinSpinCor_im.dat');  
%         cd ..;
%             %Print out info
%         percentage = (i+nMaxTau1+(j-1)*nMaxTau)/nMaxTau/nMaxDens*100;
%         printWords2 = ['Finish loading ', num2str(percentage,'%.1f'),'%...'];
%         disp(printWords2);
%     end
end
cd ..;
toc;