%Should input the name of the data folder as a variable. Change later.
cd test;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
clear; close all; clc;

%load data.
intensity = load('intensity.dat');
inversion = load('inversion.dat');
qMatrix = load('qMatrix.dat');
pMatrix = load('pMatrix.dat');

%Plot I and Sz. Can add labels later.
figure(1);
subplot(2,1,1);
plot(intensity);
subplot(2,1,2);
plot(inversion);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g(1)function.

%Take 10*transitTime as the steady state time for now. DO LATER.
%Should load input.dat here. DO LATER
dt = 0.01;
tmax = 20;
transitTime = 1;
nTrajectory = 1000;
nTimeStep = tmax/dt;

steadyMultiplier = 10;
t0 = steadyMultiplier*transitTime;
n0 = ceil(t0/dt);

q = qMatrix(:,n0:nTimeStep);
p = pMatrix(:,n0:nTimeStep);

%dim of g1 function vector
m = size(q,2);

realG1Pos = (q(:,1)'*q+p(:,1)'*p)/nTrajectory/4;
imagG1Pos = (p(:,1)'*q-q(:,1)'*p)/nTrajectory/4;

figure(2);
subplot(2,1,1);
plot(realG1Pos);
subplot(2,1,2);
plot(imagG1Pos);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spectrum

%time reversal and complex
realG1 = [fliplr(realG1Pos),realG1Pos(2:end)];
imagG1 = [-fliplr(imagG1Pos),imagG1Pos(2:end)];
G1 = realG1+1i*imagG1;
figure(3);
subplot(2,1,1);
plot(realG1);
subplot(2,1,2);
plot(imagG1);

S = fft(ifftshift(G1));
figure(4);
subplot(2,1,1);
plot(real(S));
subplot(2,1,2);
plot(imag(S));

s = real(S);
fSamp = 1/dt;
L = size(s,2);
f = fSamp*(-L/2+1:L/2)/L;
figure(5);
plot(f,fftshift(s));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit
axis = linspace(0,20,length(s));

cd ..;