%Plot I and Sz. Can add labels later.
figure(1);
set(gca,'FontSize',20);
subplot(3,1,1);
plot(linspace(0,tmax,nstore)/transitTime,intensity);
xlabel('t/T','FontSize', 20);
ylabel('\kappa \langle a^+(t) a(t) \rangle');

subplot(3,1,2);
plot(linspace(0,tmax,nstore)/transitTime,inversionAve);
xlabel('t/T','FontSize', 20);
ylabel('\langle j^z(t) \rangle');

subplot(3,1,3);
plot(linspace(0,tmax,nstore)/transitTime,inversionFinal);
xlabel('t/T','FontSize', 20);
ylabel('\langle j^z(t) \rangle');

fprintf('Program paused. Press enter to continue.\n');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g(1)function.

%Take steadyMultiplier*transitTime as the steady state time for now. DO LATER.
steadyMultiplier = 6;
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
realG1 = [fliplr(realG1Pos),realG1Pos(2:end)];
imagG1 = [-fliplr(imagG1Pos),imagG1Pos(2:end)];
G1 = realG1+1i*imagG1;
%figure(3);
%subplot(2,1,1);
%plot(realG1);
%subplot(2,1,2);
%plot(imagG1);

S = fft(ifftshift(G1));
% figure(4);
% subplot(2,1,1);
% plot(real(S));
% subplot(2,1,2);
% plot(imag(S));

df = 1/dt;
s = real(S);
L = size(s,2);
f = df*(-L/2+1:L/2)/L;
figure(5);
plot(f*transitTime,fftshift(s));
set(gca,'FontSize',20);
title('Spectrum of the Beam Laser','FontSize',20);
xlabel('Frequency/T^{-1}','FontSize', 16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit
%axis = linspace(0,20,length(s));

cd ..;