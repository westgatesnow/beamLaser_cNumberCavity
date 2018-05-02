%cNumberCavity
%Plot I and Sz. Can add labels later.
figure(1);
set(gca,'FontSize',20);
subplot(2,1,1);
plot(linspace(0,tmax,nstore)/transitTime,intensity);
xlabel('t/T','FontSize', 20);
ylabel('\kappa \langle a^+(t) a(t) \rangle');

subplot(2,1,2);
scatter(linspace(0,tmax,nTimeStep)/transitTime, szFinal, 5, 'filled');
xlabel('t/T','FontSize', 20);
ylabel('\langle j^z(t) \rangle');

%Plot sx and sy final.
% figure(2);
% 
% subplot(2,1,1);
% scatter(linspace(0,tmax,nTimeStep)/transitTime, sxFinal, 5, 'filled');
% xlabel('t/T','FontSize', 20);
% ylabel('\langle j^x(t) \rangle');
% 
% subplot(2,1,2);
% scatter(linspace(0,tmax,nTimeStep)/transitTime, syFinal, 5, 'filled');
% xlabel('t/T','FontSize', 20);
% ylabel('\langle j^y(t) \rangle');

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spinSpinCor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g(1)function.

%Take steadyMultiplier*transitTime as the steady state time for now. DO LATER.
steadyMultiplier = 10;
t0 = steadyMultiplier*transitTime;
n0 = ceil(t0/tmax*nstore);

q = qMatrix(:,n0:nstore);
p = pMatrix(:,n0:nstore);

%dim of g1 function vector
m = size(q,2);

realG1Pos = (q(:,1)'*q+p(:,1)'*p)/nTrajectory/4;
imagG1Pos = (p(:,1)'*q-q(:,1)'*p)/nTrajectory/4;

realG1 = [linspace(0,(tmax-t0)/transitTime,size(realG1Pos,2));realG1Pos]';
save realG1.dat realG1 -ascii;

figure(3);
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

%time reversal and complex
realG1Plot = [fliplr(realG1Pos),realG1Pos(2:end)];
imagG1Plot = [-fliplr(imagG1Pos),imagG1Pos(2:end)];
G1Plot = realG1Plot+1i*imagG1Plot;
%figure(4);
%subplot(2,1,1);
%plot(realG1);
%subplot(2,1,2);
%plot(imagG1);

S = fft(ifftshift(G1Plot));
% figure(5);
% subplot(2,1,1);
% plot(real(S));
% subplot(2,1,2);
% plot(imag(S));

timeInterval = tmax/nstore;
df = 1/timeInterval;
%smoothdata() equivalent to decreasing total time
s = real(S);
%Normalize s such that the peak is at 1;
s = fftshift(s)/max(s);
L = size(s,2);
f = df*(-L/2+1:L/2)/L;
spectra = [f*transitTime; s]';

%Save spectrum up to k transitTime^-1
k = 3;
effectivePosition = find(abs(spectra(:,1)) < k);
effectiveSpectra = spectra(effectivePosition,:);
save spectra.dat effectiveSpectra -ascii;

%plot the spectrum
figure(6);
plot(effectiveSpectra(:,1),effectiveSpectra(:,2));
set(gca,'FontSize',20);
title('Spectrum of the Beam Laser','FontSize',20);
xlabel('Frequency/T^{-1}','FontSize', 16);

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g(2)function

%g(2) at 0 (take the last data point as our sample data)
temp1 = sum((qMatrix(:,nstore).^2+pMatrix(:,nstore).^2-2).^2)/nTrajectory;
temp2 = sum((qMatrix(:,nstore).^2+pMatrix(:,nstore).^2-2))/nTrajectory;
g2=(temp1-4*temp2)/temp2^2;
formatSpec = 'The g2(0) value is %f \n';
fprintf(formatSpec, g2);
fprintf('Program paused. Press enter to continue.\n');