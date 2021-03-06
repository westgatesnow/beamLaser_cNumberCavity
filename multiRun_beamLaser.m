%cNumberCavity
%This code needs to be run after running "loadData_multiRun_beamLaser.m".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREPARATIONS

%analytical solutions
aValueList = zeros(nMaxTau,nMaxDens);
for j = 1:nMaxDens
    for i = 1:nMaxTau
        aValueList(i,j) = getAValue(tauList(i),kappa,rabi,densList(j));
    end
end

%set the steadyMultiplier
steadyMultiplier = 10;%This value can be varied if needed. 5 is empirical.
t0 = steadyMultiplier*transitTime;
n0_nStore = ceil(t0/tmax*nStore);
n0_nTimeStep = ceil(t0/tmax*nTimeStep);

%get the steady-state observables (mean and std):
%nAtomSS
nAtomSS = nAtom(:,:,n0_nStore:nStore);
nAtomSS_mean = mean(nAtomSS, 3);
nAtomSS_std = std(nAtomSS, 0, 3);

%intensitySS
intensitySS = intensity(:,:,n0_nStore:nStore);%in units of gc
intensitySS_mean = mean(intensitySS, 3);
intensitySS_std = std(intensitySS, 0, 3);

%inversionAveSS 
inversionAveSS = inversionAve(:,:,n0_nStore:nStore);
inversionAveSS_mean = mean(inversionAveSS, 3);
inversionAveSS_std = std(inversionAveSS, 0, 3);

%szFinalSS
szFinalSS = szFinal(:,:,n0_nTimeStep:nTimeStep);
szFinalSS_mean = mean(szFinalSS, 3, 'omitnan');
szFinalSS_std = std(szFinalSS, 0, 3, 'omitnan');

%spinSpinCorAve
spinSpinCorAveSS = spinSpinCorAve(:,:,n0_nStore:nStore);
spinSpinCorAveSS_mean = mean(spinSpinCorAveSS, 3);
spinSpinCorAveSS_std = std(spinSpinCorAveSS, 0, 3);

%linewidth
fMatrix = zeros(nMaxTau,nMaxDens);
%loop over all tau's and nAtom's;
%get corresponding parameters from loadData_multiRun_beamLaser.m;
for j = 1:nMaxDens
    for i = 1:nMaxTau
        %get current transitTime
        tau = tauList(i);
        %get q and p for a single simulation
        q(:,:) = qMatrix(i,j,:,n0_nStore:nStore);
        p(:,:) = pMatrix(i,j,:,n0_nStore:nStore);
        %dim of g1 function vector
        realG1Pos = (q(:,1)'*q+p(:,1)'*p)/nTrajectory/4;
        %imagG1Pos = (p(:,1)'*q-q(:,1)'*p)/nTrajectory/4;
        x = linspace(0,(tmax-t0),size(realG1Pos,2))';
        %x = linspace(0,(tmax-t0)*gc/tau,size(realG1Pos,2))';for tau^-1unit
        %x = linspace(0,(tmax-t0),size(realG1Pos,2))';for gc unit
        f = coeffvalues(fit(x,realG1Pos','exp1','startpoint',[1,-0.5])); 
        fMatrix(i,j) = -f(2)/pi;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3-D PLOTS

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %nAtom
% figure(21);
% set(gca,'FontSize',20);
% [X, Y] = meshgrid(tauList, nAtomAveList);
% surf(X, Y, nAtomSS_mean', nAtomSS_mean');%for some reason, need to use transpose
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('\langle N \rangle','FontSize', 20);
% title('Steady-state Intracavity Atom Number');
% colorbar;
% rotate3d on;
% 
% fprintf('Intracavity atom number shown. Press enter to see next.\n');
% pause;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %intensity
% figure(22);
% set(gca,'FontSize',20);
% [X, Y] = meshgrid(tauList, nAtomAveList);
% surf(X, Y, intensitySS_mean', intensitySS_mean');%for some reason, need to use transpose
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('\langle N \rangle','FontSize', 20);
% title('Steady-state Mean Photon Emission Rate');
% colorbar;
% rotate3d on;
% 
% fprintf('Intensity shown. Press enter to see population next.\n');
% pause;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %szFinal
% figure(23);
% set(gca,'FontSize',20);
% [X, Y] = meshgrid(tauList, nAtomAveList);
% surf(X, Y, szFinalSS_mean', szFinalSS_mean');%for some reason, need to use transpose
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('\langle N \rangle','FontSize', 20);
% title('Steady-state Final Population Inversion');
% colorbar;
% rotate3d on;
% 
% fprintf('Final population inversion shown. Press enter to see next.\n');
% pause;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %inversionAve
% figure(24);
% set(gca,'FontSize',20);
% [X, Y] = meshgrid(tauList, nAtomAveList);
% surf(X, Y, inversionAveSS_mean', inversionAveSS_mean');
% %for some reason, need to use transpose
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('\langle N \rangle','FontSize', 20);
% title('Steady-state Mean Population Inversion');
% colorbar;
% rotate3d on;
% 
% fprintf('Mean population inversion shown. Press enter to see next.\n');
% pause;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %ssCorSS
% %THIS IS NOT A GOOG OBSERVABLE FOR BEAM LASER
% % plot_ssCor_multiRun_beamLaser;
% % fprintf('Spin-spin correlations shown. Press enter to see next.\n');
% % pause;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %linewidth
% figure(25);
% [X, Y] = meshgrid(tauList, nAtomAveList);
% surf(X, Y, fMatrix',fMatrix');
% set(gca,'FontSize',20);
% set(gca,'zscale','log')
% % axis(f23,'square')
% % axis(f23,'tight', 'manual')
% % set(f23,'zlim',axi_lim)
% %for some reason, need to use transpose
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('\langle N \rangle','FontSize', 20);
% zlabel('\Delta \nu/\Gamma_c','FontSize', 20);
% title('Steady-state Linewidth');
% %make a logscale colorbar
% % colormap(jet);
% cb = colorbar();
% cb.Ruler.Scale = 'log';
% cb.Ruler.MinorTick = 'on';
% rotate3d on;
% 
% fprintf('Linewidth shown. Press enter to see 2-D plots.\n');
% pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2-D PLOTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %2-D plots (v.s. tau)
% %set the wanted nAtom
% nAtomAve_want = 100;
% 
% %get the wanted observables
% %nAtomSS_want
% nAtomSS_mean_want = nAtomSS_mean(:,(nAtomAveList == nAtomAve_want));
% nAtomSS_std_want = nAtomSS_std(:,(nAtomAveList == nAtomAve_want));
% 
% %intensitySS_want
% intensitySS_mean_want = intensitySS_mean(:,(nAtomAveList == nAtomAve_want));
% intensitySS_std_want = intensitySS_std(:,(nAtomAveList == nAtomAve_want));
% 
% %inversionAveSS_want
% inversionAveSS_mean_want = inversionAveSS_mean(:,(nAtomAveList == nAtomAve_want));
% inversionAveSS_std_want = inversionAveSS_std(:,(nAtomAveList == nAtomAve_want));
% 
% %szFinalSS_want
% szFinalSS_mean_want = szFinalSS_mean(:,(nAtomAveList == nAtomAve_want));
% szFinalSS_std_want = szFinalSS_std(:,(nAtomAveList == nAtomAve_want));
% 
% %linewidth_want
% fMatrix_want = fMatrix(:,(nAtomAveList == nAtomAve_want));
% 
% %PLOTS
% %nAtom
% figure(1);
% set(gca,'FontSize',20);
% subplot(2,1,1);
% plot(tauList,nAtomSS_mean_want);
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('N');
% subplot(2,1,2);
% plot(tauList,nAtomSS_std_want.^2);
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('\Delta N^2');
% 
% fprintf('Press enter to see next.\n');
% pause;
% 
% %intensity
% figure(2);
% h2 = errorbar(tauList, intensitySS_mean_want, intensitySS_std_want, 'o-');
% %used to logplot with xbar
% % set(get(h2,'Parent'), 'XScale', 'log');
% hold on;
% plot(tauList, nAtomSS_mean_want'/2./tauList.*(1.-szFinalSS_mean_want'));
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('I/\Gamma_c');
% hold off;
% 
% fprintf('Press enter to see next.\n');
% pause;
% 
% %inversionAve and szFinal
% figure(3);
% set(gca,'FontSize',20);
% subplot(2,1,1);
% h31 = errorbar(tauList, inversionAveSS_mean_want, inversionAveSS_std_want);
% hold on;
% %used to logplot with xbar
% % set(get(h31,'Parent'), 'XScale', 'log');
% ida = sqrt(intensitySS_mean_want');
% plot(tauList, sinc(ida.*tauList));
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('\langle j^z \rangle_e');
% hold off;
% subplot(2,1,2);
% h32 = errorbar(tauList, szFinalSS_mean_want, szFinalSS_std_want);
% hold on;
% % set(get(h32,'Parent'), 'XScale', 'log');
% plot(tauList,cos(ida.*tauList).*exp(-fMatrix_want'/2.*tauList));
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('\langle j^z_{out} \rangle_e');
% hold off;
% 
% fprintf('Press enter to see next.\n');
% pause;
% 
% %linewidth
% figure(4);
% set(gca,'FontSize',20);
% plot(tauList,fMatrix_want);
% hold on;
% % semilogx(tauList,fMatrix_want);
% % plot(tauList,nAtomAve_want*inversionAveSS_mean_want/2);
% % plot(tauList,ida.*cot(ida.*tauList/2));
% xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
% ylabel('\Delta \nu/\Gamma_c','FontSize', 20);
% hold off;
% 
% fprintf('2-D plots v.s. tau shown. Press enter to see next.\n');
% pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2-D plots (v.s. dens)
%set the wanted tau
tau_want = 0.50;
% tau_want = tau_want*gc;%in unit of gc^-1

%get the wanted observables
%nAtomSS_want
nAtomSS_mean_want = nAtomSS_mean((tauList == tau_want),:);
nAtomSS_std_want = nAtomSS_std((tauList == tau_want),:);

%intensitySS_want
intensitySS_mean_want = intensitySS_mean((tauList == tau_want),:);
intensitySS_std_want = intensitySS_std((tauList == tau_want),:);

%inversionAveSS_want
inversionAveSS_mean_want = inversionAveSS_mean((tauList == tau_want),:);
inversionAveSS_std_want = inversionAveSS_std((tauList == tau_want),:);

%szFinalSS_want
szFinalSS_mean_want = szFinalSS_mean((tauList == tau_want),:);
szFinalSS_std_want = szFinalSS_std((tauList == tau_want),:);

%spinSpinCorAve_want
spinSpinCorAveSS_mean_want = spinSpinCorAveSS_mean((tauList == tau_want),:);
spinSpinCorAveSS_std_want = spinSpinCorAveSS_std((tauList == tau_want),:);

%linewidth_want
fMatrix_want = fMatrix((tauList == tau_want),:);

%PLOTS
%nAtom
figure(11);
set(gca,'FontSize',20);
subplot(2,1,1);
hold on;
plot(densList,nAtomSS_mean_want);
plot(densList,densList.*tau_want);
xlabel('Density','FontSize', 20);
ylabel('N');
hold off;
subplot(2,1,2);
hold on;
plot(densList,nAtomSS_std_want.^2);
plot(densList,densList.*tau_want);
xlabel('Density','FontSize', 20);
ylabel('\Delta N^2');
pause;

%intensity
figure(12);
hold on;
errorbar(densList/gc, intensitySS_mean_want/gc, intensitySS_std_want/gc, 'o-');
% %used to logplot with xbar
% h12 = errorbar(densList/gc, intensitySS_mean_want/gc, intensitySS_std_want/gc, 'o-');
% set(get(h12,'Parent'), 'XScale', 'log');
plot(densList/gc, kappa/gc*aValueList.^2);
% semilogx(densList, iValueList/gc);
% loglog(densList, iValueList/gc);
plot(densList/gc, densList/gc);
xlabel('Density/\Gamma_c','FontSize', 20);
ylabel('I/\Gamma_c');%in units of gc
pause;

%inversionAve and szFinal
figure(13);
set(gca,'FontSize',20);
subplot(2,1,1);
hold on;
h131 = errorbar(densList/gc, inversionAveSS_mean_want, inversionAveSS_std_want);
% %used to logplot with xbar
% h131 = errorbar(densList/gc, inversionAveSS_mean_want, inversionAveSS_std_want);
% set(get(h131,'Parent'), 'XScale', 'log');
plot(densList/gc, sinc(aValueList*rabi*tau_want/pi));
% semilogx(densList/gc, sinc(sqrt(iValueList*gc)*tau_want/pi));
xlabel('Density/\Gamma_c','FontSize', 20);
ylabel('\langle j^z \rangle_e');
hold off;
subplot(2,1,2);
hold on;
errorbar(densList/gc, szFinalSS_mean_want, szFinalSS_std_want);
% %used to logplot with xbar
% h132 = errorbar(densList/gc, szFinalSS_mean_want, szFinalSS_std_want);
% set(get(h132,'Parent'), 'XScale', 'log');
plot(densList/gc, cos(aValueList*rabi*tau_want));
% semilogx(densList/gc, cos(sqrt(iValueList*gc)*tau_want));
xlabel('Density/\Gamma_c','FontSize', 20);
ylabel('\langle j^z_{out} \rangle_e');
hold off;
pause;

%spinSpinCorAve
figure(14);
hold on;
h12 = errorbar(densList/gc, spinSpinCorAveSS_mean_want, ...
     spinSpinCorAveSS_std_want, 'o-');
%used to logplot with xbar
% h12 = errorbar(densList, spinSpinCorAveSS_mean_want/gc, ...
%     spinSpinCorAveSS_std_want/gc, 'o-');
% set(get(h12,'Parent'), 'XScale', 'log');
plot(densList/gc, sinc(aValueList*rabi*tau_want/2/pi).^2/4);
% semilogx(densList/gc, sinc(sqrt(iValueList*gc)*tau_want/2/pi)^2/4);
% loglog(densList/gc, iValueList/gc);
xlabel('Density/\Gamma_c','FontSize', 20);
ylabel('ssCor/\Gamma_c');%in units of gc
pause;

%linewidth
figure(15);
set(gca,'FontSize',20);
subplot(2,1,1);
% semilogy(densList, fMatrix_want/gc,'-+');%in units of gc
loglog(densList/gc, fMatrix_want/gc,'-+');%in units of tau^-1
xlabel('Density/\Gamma_c','FontSize', 20);
ylabel('\Delta \nu/\Gamma_c','FontSize', 20);
subplot(2,1,2);
% semilogy(densList, fMatrix_want*tau_want,'-+');%in units of tau^-1
% semilogx(densList,fMatrix_want);
loglog(densList/gc, fMatrix_want*tau_want,'-+');%in units of tau^-1
xlabel('Density/\Gamma_c','FontSize', 20);
ylabel('\Delta \nu/\tau^{-1}','FontSize', 20);

fprintf('2-D plots v.s. nAtomAve shown. Press enter to see next.\n');
pause;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %g(2)
% %g(2)function
% 
% %g(2) at 0 (take the last data point as our sample data)
% % temp1 = sum((qMatrix(:,nStore).^2+pMatrix(:,nStore).^2-2).^2)/nTrajectory;
% % temp2 = sum((qMatrix(:,nStore).^2+pMatrix(:,nStore).^2-2))/nTrajectory;
% % g2=(temp1-4*temp2)/temp2^2;
% % formatSpec = 'The g2(0) value is %f \n';
% % fprintf(formatSpec, g2);
% % fprintf('Program paused. Press enter to continue.\n');