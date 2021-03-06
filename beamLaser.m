%cNumberCavity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREPARATIONS
nTimeStep = tmax/dt;
tList = linspace(0,tmax,nStore)/transitTime;%in units of transitTime
tScatterList = linspace(0,tmax,nTimeStep)/transitTime;%in units of transitTime

%Take steadyMultiplier*transitTime as the steady state time. 
%This is empirical for now. Improve later?
steadyMultiplier = 10;

t0 = steadyMultiplier*transitTime;
n0 = ceil(t0/tmax*nStore);

%Solve the analytical solution
aValue = getAValue(transitTime,kappa,rabi,density);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot nAtom.
figure(1);
set(gca,'FontSize',20);
hold on;
plot(tList, nAtom);
plot(tList, density*transitTime*ones(1,nStore));
hold off;
xlabel('t/T','FontSize', 20);
ylabel('N(t)');

nAtomPrint = mean(nAtom(n0:end));
formatSpec = 'The steady-state nAtom is %d. \n';
fprintf(formatSpec, nAtomPrint);
fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot I and szFinal. Can add labels later.
figure(2);
set(gca,'FontSize',20);
subplot(2,1,1);
hold on;
plot(tList, intensity/kappa);
plot(tList, aValue^2*ones(1,nStore));
hold off;
xlabel('t/T','FontSize', 20);
ylabel('\langle a^+(t) a(t) \rangle');

subplot(2,1,2);
hold on;
scatter(tScatterList, szFinal, 5, 'filled');
plot(tList, 1-2*intensity/density);
plot(tList, 1-2*kappa*aValue^2*ones(1,nStore)/density);
hold off;
xlabel('t/T','FontSize', 20);
ylabel('\langle j^z(t) \rangle_o');

intensityPrint = mean(intensity(n0:end));
szFinalPrint = mean(szFinal(n0:end), 'omitnan');
formatSpec = 'The steady-state intensity is %4.4f. \nThe steady-state szFinal is %4.4f.\n';
fprintf(formatSpec, intensityPrint, szFinalPrint);
fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%comparison to adiabatic elimination
% figure(22);
% %q and qAE
% 
set(gca,'FontSize',20);
subplot(2,1,1);
hold on;
plot(tList, mean(qMatrix,1));
% plot(tList, mean(-JyMatrix/kappa,1));
hold off;
xlabel('t/T','FontSize', 20);
ylabel('q');

Deltat=0;
%Deltat defined above is impirical and used to acount for the right Delta t
%that should used for the adiabatic elimination
subplot(2,1,2);
hold on;
plot(tList, mean(qMatrix.^2,1));
% plot(tList, mean((-JyMatrix/kappa).^2,1));%+4/kappa/Deltat);
hold off;
xlabel('t/T','FontSize', 20);
ylabel('q^2');

subplot(3,1,3);
JySq = mean((JyMatrix).^2,1);
plot(linspace(0,tmax,nStore)/transitTime, JySq);
xlabel('t/T','FontSize', 20);
ylabel('q^2');
JySqPrint = mean(JySq(n0:end));
formatSpec = 'The steady-state JySq is %4.4f. \n';
fprintf(formatSpec, JySqPrint);

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot sx and sy final.
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
hold on;
plot(tList, inversionAve);%, 5, 'filled');
plot(tList, sinc(rabi*aValue*transitTime/pi)*ones(1,nStore));%sinc has pi
hold off;
xlabel('t/T','FontSize', 20);
ylabel('\langle j^z(t)\rangle');

subplot(2,1,2);
hold on;
plot(1:nBin, szMatrix(:,end));%, 5, 'filled');
plot(1:nBin, cos(rabi*aValue*transitTime*(1:nBin)/nBin));
hold off;
xlabel('nBin','FontSize', 20);
ylabel('\langle j^z(t) \rangle');

inversionAvePrint = mean(inversionAve(n0:end));
formatSpec = 'The steady-state szAve is %4.4f.\n';
fprintf(formatSpec, inversionAvePrint);
fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot sxMatrix and syMatrix.
figure(5);

subplot(2,1,1);
scatter(1:nBin, mean(sxMatrix(:,n0:nStore),2), 20, 'filled');
xlabel('nBin','FontSize', 20);
ylabel('\langle j^x(t) \rangle');

subplot(2,1,2);
scatter(1:nBin, mean(syMatrix(:,n0:nStore),2), 20, 'filled');
xlabel('nBin','FontSize', 20);
ylabel('\langle j^y(t) \rangle');

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot spinSpinCorAve.
figure(6);
set(gca,'FontSize',20);
% subplot(2,1,1);
hold on;
plot(tList,spinSpinCorAve_re);
plot(tList,1/4*sinc(rabi*transitTime*aValue/pi/2)^2*ones(1,nStore));
hold off;
% axis([0 inf 0 0.20]);
xlabel('t/T','FontSize', 20);
ylabel('Re[\langle \sigma^+_i(t) \sigma^-_j(t) \rangle]');

% subplot(2,1,2);
% plot(linspace(0,tmax,nStore)/transitTime, spinSpinCorAve_im);
% xlabel('t/T','FontSize', 20);
% ylabel('Im[\langle \sigma^+_i(t) \sigma^-_j(t) \rangle]');

ssCorRePrint = mean(spinSpinCorAve_re(n0:end));
formatSpec = 'The steady-state ssCorRe is %4.4f.\n';
fprintf(formatSpec, ssCorRePrint);
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
[X, Y] = meshgrid(1:nBin, 1:nBin);
s_re_model = 1/4*sin(rabi*aValue*transitTime*X/nBin)...
    .*sin(rabi*aValue*transitTime*Y/nBin);
surf(X, Y, s_re_model, s_re_model);
colorbar;
rotate3d on;

% subplot(2,1,2);
% s_im = reshape(spinSpinCor_im,nBin,nBin,[]);
% surf(X, Y, s_im(:,:,end),s_im(:,:,end));
% colorbar;
% rotate3d on;

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Movie of space-space correlation over time
%This program uses the definition of s_re and X, Y, Z.
% movie_ssCor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single-point g(1) and decorrelation

q = qMatrix(:,n0:nStore);
p = pMatrix(:,n0:nStore);

%dim of g1 function vector
m = size(q,2);

realG1Pos = (q(:,1)'*q+p(:,1)'*p)/nTrajectory/4;
imagG1Pos = (p(:,1)'*q-q(:,1)'*p)/nTrajectory/4;

realG1 = [linspace(0,(tmax-t0)/transitTime,m);realG1Pos]';
save realG1.dat realG1 -ascii;

figure(8);
subplot(2,1,1);
plot(linspace(0,(tmax-t0)/transitTime,m),realG1Pos)
xlabel('t/T','FontSize', 20);
ylabel('Re[g^{(1)}(t)]');
subplot(2,1,2);
plot(linspace(0,(tmax-t0)/transitTime,m),imagG1Pos)
xlabel('t/T','FontSize', 20);
ylabel('Im[g^{(1)}(t)]');

fprintf('Program paused. Press enter to continue.\n');
pause;

% check <jz q q(0)> = <jz><q q(0)> 
% jz = JzMatrix(:,n0:nStore);
% left = mean(jz.*q.*q(:,1),1);
% right = mean(jz,1).*mean(q(:,1).*q,1)+mean(q,1).*mean(q(:,1).*jz,1)...
%     +mean(q(:,1),1).*mean(jz.*q)-2*mean(jz,1).*mean(q,1).*mean(q(:,1),1);
% 
% ratioLeft = left./(mean(q(:,1).*q,1));
% ratioRight = (nAtomPrint-1)/nAtomPrint*mean(jz,1)-1;
% 
% leftSave = [linspace(0,(tmax-t0)/transitTime,size(left,2));left]';
% save left.dat leftSave -ascii;
% 
% figure(81);
% subplot(2,1,1);
% hold on;
% plot(linspace(0,(tmax-t0)/transitTime,size(left,2)),left);
% plot(linspace(0,(tmax-t0)/transitTime,size(left,2)),right);
% hold off;
% xlabel('t/T','FontSize', 20);
% subplot(2,1,2);
% hold on;
% plot(linspace(0,(tmax-t0)/transitTime,size(left,2)),ratioLeft);
% plot(linspace(0,(tmax-t0)/transitTime,size(left,2)),ratioRight);
% hold off;
% xlabel('t/T','FontSize', 20);

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single-point spectrum

%Method 1: fit the realG1Pos to exponential decay
x = linspace(0,tmax-t0,m)';

%Normalize realG1Pos
realG1PosNorm = realG1Pos./max(realG1Pos);

f = coeffvalues(fit(x,realG1PosNorm','exp1','startpoint',[1,-0.1])); 
linewidth1 = -f(2)/pi;
gammac = rabi^2/kappa;
ratio_tau1 = linewidth1*transitTime;
ratio_gc1 = linewidth1/gammac;

figure(9);
hold on;
yCalc1 = f(1)*exp(f(2).*x);
plot(x,realG1PosNorm);
plot(x,yCalc1,'--')
hold off;

printWords1 = ['The linewidth is ', num2str(ratio_tau1,'%.4f'),...
    ' transitTime inverse or ', num2str(ratio_gc1,'%.4f'),' gammac.'];
disp(printWords1);
fprintf('Program paused. Press enter to continue.\n');
pause;
%Method 2: fit the log realG1Pos to linear
x2 = [ones(m,1) x];
y2 = real(log(realG1Pos'));
b = x2\y2;

figure(91);
hold on;
yCalc2 = x2*b;
plot(x, log(realG1Pos))
plot(x,yCalc2)
% legend('Data','Slope','Slope & Intercept','Location','best');
hold off;

linewidth2 = -b(2)/pi;
ratio_tau2 = linewidth2*transitTime;
ratio_gc2 = linewidth2/gammac;
printWords = ['The linewidth is ', num2str(ratio_tau2,'%.4f'),...
    ' transitTime inverse or ', num2str(ratio_gc2,'%.4f'),' gammac.'];
disp(printWords);
fprintf('Program paused. Press enter to continue.\n');
pause;
%Method 3: plot the spectra
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
effectivePosition = abs(spectra(:,1)) < k;
effectiveSpectra = spectra(effectivePosition,:);
save spectra.dat effectiveSpectra -ascii;

%plot the spectrum
figure(92);
plot(effectiveSpectra(:,1),effectiveSpectra(:,2));
set(gca,'FontSize',20);
title('Spectrum of the Beam Laser','FontSize',20);
xlabel('Frequency/T^{-1}','FontSize', 16);

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%multi-point g(1) function.

step = 100;
nStep = (nStore-n0)/step;
mulRealG1Pos = zeros(1,step+1);
mulImagG1Pos = zeros(1,step+1);
mulLeft = zeros(1,step+1);
mulRight = zeros(1,step+1);
mulRatioLeft = zeros(1,step+1);
mulRatioRight = zeros(1,step+1);
for i=1:nStep
    mulq = qMatrix(:,n0+step*(i-1):n0+step*i);
    mulp = pMatrix(:,n0+step*(i-1):n0+step*i);
    %dim of g1 function vector
    mulRealG1Pos = mulRealG1Pos + (mulq(:,1)'*mulq+p(:,1)'*mulp)/nTrajectory/4;
    mulImagG1Pos = mulImagG1Pos + (mulp(:,1)'*mulq-mulq(:,1)'*mulp)/nTrajectory/4;
    
    jz = JzMatrix(:,n0+step*(i-1):n0+step*i);
    mulLeft = mulLeft + mean(jz.*mulq.*mulq(:,1),1);
    mulRight = mulRight + mean(jz,1).*mean(mulq(:,1).*mulq,1);
%     +mean(q,1).*mean(q(:,1).*jz,1)...
%         +mean(q(:,1),1).*mean(jz.*q)-2*mean(jz,1).*mean(q,1).*mean(q(:,1),1);

    mulRatioLeft = mulRatioLeft + mean(jz.*mulq.*mulq(:,1),1)./(mean(mulq(:,1).*mulq,1));
    mulRatioRight = mulRatioRight + (nAtomPrint-1)/nAtomPrint*mean(jz,1)-1;
end
mulRealG1Pos = mulRealG1Pos/nStep;
mulImagG1Pos = mulImagG1Pos/nStep;
mulLeft = mulLeft/nStep;
mulRight = mulRight/nStep;
mulRatioLeft = mulRatioLeft/nStep;
mulRatioRight = mulRatioRight/nStep;

figure(10);
subplot(2,1,1);
plot(linspace(0,tmax/nStore*step/transitTime,step+1),mulRealG1Pos)
xlabel('t/T','FontSize', 20);
ylabel('Re[g^{(1)}(t)]');
subplot(2,1,2);
plot(linspace(0,tmax/nStore*step/transitTime,step+1),mulImagG1Pos)
xlabel('t/T','FontSize', 20);
ylabel('Im[g^{(1)}(t)]');

fprintf('Program paused. Press enter to continue.\n');
pause;

figure(101);
subplot(2,1,1);
hold on;
plot(linspace(0,tmax/nStore*step/transitTime,step+1),mulLeft);
plot(linspace(0,tmax/nStore*step/transitTime,step+1),mulRight);
hold off;
xlabel('t/T','FontSize', 20);
subplot(2,1,2);
hold on;
plot(linspace(0,tmax/nStore*step/transitTime,step+1),mulRatioLeft);
plot(linspace(0,tmax/nStore*step/transitTime,step+1),mulRatioRight);
hold off;
xlabel('t/T','FontSize', 20);

fprintf('Program paused. Press enter to continue.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%multiple-point spectrum

%Method 1: fit the mulRealG1Pos to exponential decay
mulx = linspace(0,tmax/nStore*step,step+1)';
mulRealG1PosNorm = mulRealG1Pos/max(mulRealG1Pos);
mulf = coeffvalues(fit(mulx,mulRealG1PosNorm','exp1','startpoint',[1,-1])); 
mulLinewidth1 = -mulf(2)/pi;
mulRatio_tau1 = mulLinewidth1*transitTime;
mulRatio_gc1 = mulLinewidth1/gammac;

figure(11);
hold on;
mulyCalc1 = mulf(1)*exp(mulf(2).*mulx);%??????
plot(mulx,mulRealG1PosNorm);
plot(mulx,mulyCalc1,'--')
hold off;

printWords1 = ['The linewidth is ', num2str(mulRatio_tau1,'%.4f'),...
    ' transitTime inverse or ', num2str(mulRatio_gc1,'%.4f'),' gammac.'];
disp(printWords1);
fprintf('Program paused. Press enter to continue.\n');
pause;
%Method 2: fit the log realG1Pos to linear
mulx2 = [ones(step+1,1) mulx];
muly2 = real(log(mulRealG1Pos'));
mulb = mulx2\muly2;

figure(110);
hold on;
mulyCalc2 = mulx2*mulb;
plot(mulx, log(mulRealG1Pos))
plot(mulx,mulyCalc2)
% legend('Data','Slope','Slope & Intercept','Location','best');
hold off;

mulLinewidth2 = -mulb(2)/pi;
mulRatio_tau2 = mulLinewidth2*transitTime;
mulRatio_gc2 = mulLinewidth2/gammac;
printWords = ['The linewidth is ', num2str(mulRatio_tau2,'%.4f'),...
    ' transitTime inverse or ', num2str(mulRatio_gc2,'%.4f'),' gammac.'];
disp(printWords);
fprintf('Program paused. Press enter to continue.\n');
pause;
%Method 3: plot the spectra
%time reversal and complex
realG1Plot = [fliplr(mulRealG1Pos),mulRealG1Pos(2:end)];
imagG1Plot = [-fliplr(mulImagG1Pos),mulImagG1Pos(2:end)];
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
figure(111);
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
formatSpec = 'The g2(0) value is %.2f \n';
fprintf(formatSpec, g2);
fprintf('Program paused. Press enter to continue.\n');