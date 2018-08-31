%cNumberCavity
%This file is used to load the data from running "multirun_pz.sh"

%Initialization
clear; close all; clc;
addpath ~/Desktop/codes/beamLaser_Proj/cNumberCavity/;

%get parameters from the input.txt in the parent directory
%all the parameters are universal except for pz, dens and name
getParam_beamLaser;
nTimeStep = tmax/dt;
gc = rabi^2/kappa;

%bash variables; 

%pz
%pzList1 = 1.0 ~ 11.0;
nMaxPz1 = 10;%
initPz1 = 1.0;%
intervalPz1 = 1.0;
pzList1 = (initPz1+(0:nMaxPz1-1)*intervalPz1);
%tauList = 1:nMaxPz;%in the unit of tau, if initPz = intervalPz;
% nMaxPz2 = 20;%20/0
% initPz2 = 1.0;
% intervalPz2 = 0.5;
% tauList2 = (initPz2+(0:nMaxPz2-1)*intervalPz2)*gc;
%for convenience
% nMaxPz = nMaxPz1+nMaxPz2;
% tauList = [tauList1,tauList2];

% %dens
% nMaxDens = 1;%10
% initDens = 100;
% intervalDens = 100;
% densAveList =initDens+(0:nMaxDens-1)*intervalDens;

%define empty data structure for variables
dens = zeros(nMaxPz,nStore);
intensity = zeros(nMaxPz,nStore);
inversionAve = zeros(nMaxPz,nStore);
spinSpinCorAve = zeros(nMaxPz,nStore);
szFinal = zeros(nMaxPz,nTimeStep);
qMatrix = zeros(nMaxPz,nTrajectory,nStore);
pMatrix = zeros(nMaxPz,nTrajectory,nStore);

%load all the data needed
timeRequired = 3.2*nMaxPz*nMax;
printWords1 = ['Begin loading. Time required: ' , num2str(timeRequired), ' seconds.'];
disp(printWords1);
tic;

    for i = 1:nMaxPz1
        %define the name of the directory
        tau1 = initPz1+intervalPz1*(i-1);
        filename = ['dt0.005_dZ0.5_dPz',num2str(pz1,'%.2f'), ...
            '_tau1.0_nBin30_dens100_g3_k90_yWall1.0'];
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
        percentage = (i+(j-1)*nMaxPz)/nMaxPz*100;
        printWords2 = ['Finish loading ', num2str(percentage,'%.1f'),'%...'];
        disp(printWords2);
    end
toc;