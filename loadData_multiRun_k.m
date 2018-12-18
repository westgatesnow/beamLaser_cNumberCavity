%cNumberCavity
%This file is used to load the data from running "multiRun_k.sh"

%Initialization
clear; close all; clc;
addpath ~/Desktop/codes/beamLaser_Proj/cNumberCavity/;

%get parameters from the input.txt in the parent directory
%all the parameters are universal except for kappa, dens and name
getParam_beamLaser;
nTimeStep = tmax/dt;

%bash variables; 

%kappa
%kappaList1 = 0.0 ~ 100.0;
nMaxKappa = 30;%
initKappa = -1.0;%
intervalKappa = 0.1;
kappaList = (10.^(initKappa+(0:nMaxKappa-1)*intervalKappa))';
%tauList = 1:nMaxKappa;%in the unit of tau, if initKappa = intervalKappa;
% nMaxKappa2 = 0;%20/0
% initKappa2 = 1.0;
% intervalKappa2 = 1.0;
% kappaList2 = (initKappa2+(0:nMaxKappa2-1)*intervalKappa2);
% for convenience
% nMaxKappa = nMaxKappa+nMaxKappa2;
% kappaList = [kappaList1,kappaList2];

%define empty data structure for variables
intensity = zeros(nMaxKappa,nStore);
inversionAve = zeros(nMaxKappa,nStore);
spinSpinCorAve = zeros(nMaxKappa,nStore);
szFinal = zeros(nMaxKappa,nTimeStep);
qMatrix = zeros(nMaxKappa,nTrajectory,nStore);
pMatrix = zeros(nMaxKappa,nTrajectory,nStore);

%load all the data needed
timeRequired = 3.2*nMaxKappa;
printWords1 = ['Begin loading. Time required: ' , num2str(timeRequired), ' seconds.'];
disp(printWords1);
tic;

cd('kappa');
    for i = 1:nMaxKappa
        kappa = initKappa+intervalKappa*(i-1);
        filename = ['dt0.001_dZ0_dPz0_tau1.0_nBin20_dens100_g3_kexp',num2str(kappa,'%.1f'),];
        cd(filename);
        %Load data.
        intensity(i,:) = load('intensity.dat');
        inversionAve(i,:) = load('inversionAve.dat');
        spinSpinCorAve(i,:) = load('spinSpinCorAve_re.dat');
        szFinal(i,:) = load('szFinal.dat');
        qMatrix(i,:,:) = load('qMatrix.dat');
        pMatrix(i,:,:) = load('pMatrix.dat');
        cd ..;
        %Print out info
        percentage = (i)/nMaxKappa*100;
        printWords2 = ['Finish loading ', num2str(percentage,'%.1f'),'%...'];
        disp(printWords2);  
    end
cd ..;
toc;