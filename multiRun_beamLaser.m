%cNumberCavity
%This code needs to be run after running "loadData_multiRun_beamLaser.m".

%set the steadyMultiplier
steadyMultiplier = 5;%This value can be varied if needed. 5 is empirical.
t0 = steadyMultiplier*transitTime;
n0_nstore = ceil(t0/tmax*nstore);
n0_nTimeStep = ceil(t0/tmax*nTimeStep);

%nAtomSS
plot_nAtom_multiRun_beamLaser;
fprintf('Intracavity atom number shown. Press enter to see next.\n');
pause;

%intensitySS
plot_intensity_multiRun_beamLaser;
fprintf('Intensity shown. Press enter to see population next.\n');
pause;

%inversionAveSS and szFinalSS
plot_sz_multiRun_beamLaser;
fprintf('Population inversion shown. Press enter to see next.\n');
pause;

% %intensitySS and szFinalSS
% plot_IvsSz_multiRun_beamLaser;
% fprintf('Comparison between I and szFinal shown. Press enter to see spin-spin correlations.\n');
% pause;

%ssCorSS
plot_ssCor_multiRun_beamLaser;
fprintf('Spin-spin correlations shown. Press enter to see next.\n');
pause;

%linewidth
plot_linewidth_multiRun_beamLaser;
fprintf('Linewidth shown. Press enter to see next.\n');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g(2)
%g(2)function

%g(2) at 0 (take the last data point as our sample data)
% temp1 = sum((qMatrix(:,nstore).^2+pMatrix(:,nstore).^2-2).^2)/nTrajectory;
% temp2 = sum((qMatrix(:,nstore).^2+pMatrix(:,nstore).^2-2))/nTrajectory;
% g2=(temp1-4*temp2)/temp2^2;
% formatSpec = 'The g2(0) value is %f \n';
% fprintf(formatSpec, g2);
% fprintf('Program paused. Press enter to continue.\n');