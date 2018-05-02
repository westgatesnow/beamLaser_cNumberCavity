%cNumberCavity
%This file is used to load the data from running "multirun_transitTime.sh"

%Initialization
clear; close all; clc;

%bash variables; 
%We are varying transitTime (tau) here, with tau*dens = nAtomAve. With
%   nAtomAve fixed, dens is passively changed when actively changing tau.
nMax = 100;
tauInit = 0.1;
tauInterval = 0.1;
nAtomAve = 100;



