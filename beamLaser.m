%cNumberCavity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot nAtom.
figure(1);
set(gca,'FontSize',20);
plot(linspace(0,tmax,nStore)/transitTime, nAtom);
xlabel('t/T','FontSize', 20);
ylabel('N(t)');

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot I and szFinal. Can add labels later.
figure(2);
set(gca,'FontSize',20);
subplot(2,1,1);
plot(linspace(0,tmax,nStore)/transitTime, intensity);
xlabel('t/T','FontSize', 20);
ylabel('\kappa \langle a^+(t) a(t) \rangle');

subplot(2,1,2);
hold on;
scatter(linspace(0,tmax,nTimeStep)/transitTime, szFinal, 5, 'filled');
plot(linspace(0,tmax,nStore)/transitTime, 1-2*intensity/density);
hold off;
xlabel('t/T','FontSize', 20);
ylabel('\langle j^z(t) \rangle_o');

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%comparison to adiabatic elimination
figure(12);
%q and qAE

set(gca,'FontSize',20);
subplot(2,1,1);
hold on;
plot(linspace(0,tmax,nStore)/transitTime, mean(qMatrix,1));
plot(linspace(0,tmax,nStore)/transitTime, mean(qAEMatrix,1));
hold off;
xlabel('t/T','FontSize', 20);
ylabel('q');

Deltat=0.1;
%Deltat defined above is impirical and used to acount for the right Delta t
%that should used for the adiabatic elimination
subplot(2,1,2);
hold on;
plot(linspace(0,tmax,nStore)/transitTime, mean(qMatrix.^2,1));
plot(linspace(0,tmax,nStore)/transitTime, mean(qAEMatrix.^2,1)+4/kappa/Deltat);
hold off;
xlabel('t/T','FontSize', 20);
ylabel('q^2');

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot sx and sy final.
% figure(3);
% 
% subplot(2,1,1);
% scatter(linspace(0,tmax,nTimeStep)/transitTime, sxFinal, 5, 'filled');
% xlabel('t/T','FontSize', 20);
% ylabel('\langle j^x(t)\rangle');
% 
% subplot(2,1,2);
% scatter(linspace(0,tmax,nTimeStep)/transitTime, syFinal, 5, 'filled');
% xlabel('t/T','FontSize', 20);
% ylabel('\langle j^y(t) \rangle');
% 
% fprintf('Program paused. Press enter to continue.\n');
% pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot szAve and szMatrix.
figure(4);

subplot(2,1,1);
scatter(linspace(0,tmax,nStore)/transitTime, inversionAve, 5, 'filled');
xlabel('t/T','FontSize', 20);
ylabel('\langle j^z(t)\rangle');

subplot(2,1,2);
plot(1:nBin, szMatrix(:,end));%, 5, 'filled');
xlabel('nBin','FontSize', 20);
ylabel('\langle j^z(t) \rangle');

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot sxMatrix and syMatrix.
% figure(5);
% 
% subplot(2,1,1);
% scatter(1:nBin, sxMatrix(:,end), 20, 'filled');
% xlabel('nBin','FontSize', 20);
% ylabel('\langle j^x(t) \rangle');
% 
% subplot(2,1,2);
% scatter(1:nBin, syMatrix(:,end), 20, 'filled');
% xlabel('nBin','FontSize', 20);
% ylabel('\langle j^y(t) \rangle');
% 
% fprintf('Program paused. Press enter to continue.\n');
% pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot spinSpinCorAve.
figure(6);
set(gca,'FontSize',20);
subplot(2,1,1);
plot(linspace(0,tmax,nStore)/transitTime,spinSpinCorAve_re);
axis([0 inf 0 0.15]);
xlabel('t/T','FontSize', 20);
ylabel('Re[\langle \sigma^+_i(t) \sigma^-_j(t) \rangle]');

subplot(2,1,2);
plot(linspace(0,tmax,nStore)/transitTime, spinSpinCorAve_im);
xlabel('t/T','FontSize', 20);
ylabel('Im[\langle \sigma^+_i(t) \sigma^-_j(t) \rangle]');

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot spinSpinCor.
figure(7);
set(gca,'FontSize',20);
subplot(2,1,1);
s_re = reshape(spinSpinCor_re,nBin,nBin,[]);
[X, Y] = meshgrid(1:nBin, 1:nBin);
surf(X, Y, s_re(:,:,end),s_re(:,:,end));
colorbar;
rotate3d on;

subplot(2,1,2);
s_im = reshape(spinSpinCor_im,nBin,nBin,[]);
surf(X, Y, s_im(:,:,end),s_im(:,:,end));
colorbar;
rotate3d on;

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g(1)function.

%Take steadyMultiplier*transitTime as the steady state time for now. DO LATER.
steadyMultiplier = 15;
t0 = steadyMultiplier*transitTime;
n0 = ceil(t0/tmax*nStore);

q = qMatrix(:,n0:nStore);
p = pMatrix(:,n0:nStore);

%dim of g1 function vector
m = size(q,2);

realG1Pos = (q(:,1)'*q+p(:,1)'*p)/nTrajectory/4;
imagG1Pos = (p(:,1)'*q-q(:,1)'*p)/nTrajectory/4;

realG1 = [linspace(0,(tmax-t0)/transitTime,size(realG1Pos,2));realG1Pos]';
save realG1.dat realG1 -ascii;

figure(8);
subplot(2,1,1);
plot(linspace(0,(tmax-t0)/transitTime,size(realG1Pos,2)),realG1Pos);
xlabel('t/T','FontSize', 20);
ylabel('Re[g^{(1)}(t)]');
subplot(2,1,2);
plot(linspace(0,(tmax-t0)/transitTime,size(realG1Pos,2)),imagG1Pos);
xlabel('t/T','FontSize', 20);
ylabel('Im[g^{(1)}(t)]');

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spectrum

%Method 1: fit the realG1Pos to exponential decay
x = linspace(0,(tmax-t0)/transitTime,size(realG1Pos,2))';
f = coeffvalues(fit(x,realG1Pos','exp1','startpoint',[1,-1])); 
linewidth = -f(2)/pi;
gammac = rabi^2/kappa*transitTime;
ratio_gc = linewidth/gammac;
printWords = ['The linewidth is ', num2str(linewidth,'%.2f'),...
    ' transitTime inverse or ', num2str(ratio_gc,'%.2f'),' gammac.'];
disp(printWords);

%Method 2: plot the spectra
%time reversal and complex
realG1Plot = [fliplr(realG1Pos),realG1Pos(2:end)];
imagG1Plot = [-fliplr(imagG1Pos),imagG1Pos(2:end)];
G1Plot = realG1Plot+1i*imagG1Plot;
%figure(9);
%subplot(2,1,1);
%plot(realG1);
%subplot(2,1,2);
%plot(imagG1);

S = fft(ifftshift(G1Plot));
% figure(10);
% subplot(2,1,1);
% plot(real(S));
% subplot(2,1,2);
% plot(imag(S));

timeInterval = tmax/nStore;
df = 1/timeInterval;
%smoothdata() equivalent to decreasing total time
s = real(S);
%Normalize s such that the peak is at 1;
s = fftshift(s)/max(s);
L = size(s,2);
f = df*(-L/2+1:L/2)/L;
spectra = [f*transitTime; s]';

%Save spectrum up to k transitTime^-1
k = 100;
effectivePosition = find(abs(spectra(:,1)) < k);
effectiveSpectra = spectra(effectivePosition,:);
save spectra.dat effectiveSpectra -ascii;

%plot the spectrum
figure(11);
plot(effectiveSpectra(:,1),effectiveSpectra(:,2));
set(gca,'FontSize',20);
title('Spectrum of the Beam Laser','FontSize',20);
xlabel('Frequency/T^{-1}','FontSize', 16);

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g(2)function

%g(2) at 0 (take the last data point as our sample data)
q4 = sum(qMatrix(:,nStore).^4)/nTrajectory;
p4 = sum(pMatrix(:,nStore).^4)/nTrajectory;
q2p2 = sum(qMatrix(:,nStore).^2.*pMatrix(:,nStore).^2)/nTrajectory;
q2 = sum(qMatrix(:,nStore).^2)/nTrajectory;
p2 = sum(pMatrix(:,nStore).^2)/nTrajectory;
g2=(q4+p4+2*q2p2-8*(q2+p2)+8)/(q2+p2-2)^2;
formatSpec = 'The g2(0) value is %f \n';
fprintf(formatSpec, g2);
fprintf('Program paused. Press enter to continue.\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
