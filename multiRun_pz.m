%cNumberCavity
%This code needs to be run after running "loadData_multiRun_pz.m".

%set the steadyMultiplier
steadyMultiplier = 15;%This value can be varied if needed. 5 is empirical.
t0 = steadyMultiplier*transitTime;
n0_nStore = ceil(t0/tmax*nStore);
n0_nTimeStep = ceil(t0/tmax*nTimeStep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nAtomSS
% plot_nAtom_multiRun_beamLaser;
% fprintf('Intracavity atom number shown. Press enter to see next.\n');
% pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%intensitySS
%we first look at the two extreme cases
% figure(3);
% set(gca,'FontSize',20);
% subplot(2,1,1);
% intensity_i(:) = intensity(1,:);
% plot(linspace(0,tmax,nStore)/transitTime,intensity_i);
% xlabel('t/T','FontSize', 20);
% ylabel('\kappa \langle a^+(t) a(t) \rangle');
% subplot(2,1,2);
% intensity_f(:) = intensity(end,:);
% plot(linspace(0,tmax,nStore)/transitTime,intensity_f);
% xlabel('t/T','FontSize', 20);
% ylabel('\kappa \langle a^+(t) a(t) \rangle');

%2-D plots
%get E[I_s] and V[I_s]
intensitySS = intensity(:,n0_nStore:nStore)/gc;%in units of gc
intensitySS_mean = mean(intensitySS, 2);
intensitySS_std = std(intensitySS, 0, 2);
%for 2-D plots with error bars (v.s. pz)
figure(4);
h4 = errorbar(pzList/gc, intensitySS_mean, intensitySS_std, 'o-');%used to logplot with xbar
set(get(h4,'Parent'), 'XScale', 'log');
xlabel('kv/\Gamma_c','FontSize', 20);
ylabel('I/\Gamma_c');
%for 3-D plots
%figure(5)

fprintf('Intensity shown. Press enter to see population next.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inversionAveSS and szFinalSS

%get E[inversionAve_s] and V[inversionAve_s]
inversionAveSS = inversionAve(:,n0_nStore:nStore);
inversionAveSS_mean = mean(inversionAveSS, 2);
inversionAveSS_std = std(inversionAveSS, 0, 2);
%get E[szFinal_s] and V[szFinal_s]
szFinalSS = szFinal(:,n0_nTimeStep:nTimeStep);
szFinalSS_mean = mean(szFinalSS, 2, 'omitnan');
szFinalSS_std = std(szFinalSS, 0, 2, 'omitnan');

%2-D plots (v.s. pz)
figure(6);
set(gca,'FontSize',20);
subplot(2,1,1);
h61 = errorbar(pzList/gc, inversionAveSS_mean, inversionAveSS_std);%used to logplot with xbar
set(get(h61,'Parent'), 'XScale', 'log');
xlabel('kv/\Gamma_c','FontSize', 20);
ylabel('\langle j^z \rangle_e');
subplot(2,1,2);
h62 = errorbar(pzList/gc, szFinalSS_mean, szFinalSS_std);
set(get(h62,'Parent'), 'XScale', 'log');
xlabel('kv/\Gamma_c','FontSize', 20);
ylabel('\langle j^z_{out} \rangle_e');

%for 3-D plots
%figure(7)

fprintf('Population inversion shown. Press enter to see next.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %intensitySS and szFinalSS
% figure(8);
% yyaxis left
% h4 = errorbar(pzList/gc, intensitySS_mean, intensitySS_std, 'o-');%used to logplot with xbar
% set(get(h4,'Parent'), 'XScale', 'log');
% xlabel('kv/\Gamma_c','FontSize', 20);
% ylabel('I/\Gamma_c');
% yyaxis right
% h62 = errorbar(pzList/gc, szFinalSS_mean, szFinalSS_std);
% set(get(h62,'Parent'), 'XScale', 'log');
% xlabel('kv/\Gamma_c','FontSize', 20);
% ylabel('\langle j^z_{out} \rangle_e');
% 
% fprintf('Comparison between I and szFinal shown. Press enter to see spin-spin correlations.\n');
% pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ssCorSS

% %2-D plots
% %get E[ssCor_s] and V[ssCor_s]
% ssCorSS = spinSpinCorAve(:,n0_nStore:nStore);
% ssCorSS_mean = mean(ssCorSS, 2);
% ssCorSS_std = std(ssCorSS, 0, 2);
% %for 2-D plots with error bars (v.s. pz)
% figure(10);
% h8 = errorbar(pzList/gc, ssCorSS_mean, ssCorSS_std, 'o-');%used to logplot with xbar
% set(get(h8,'Parent'), 'XScale', 'log');
% axis([0 inf 0 0.005]);
% xlabel('kv/\Gamma_c','FontSize', 20);
% ylabel('\langle \sigma_i^+ \sigma_j^-\rangle');

%for 3-D plots
%figure(11)
% fprintf('Spin-spin correlations shown. Press enter to see next.\n');
% pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%linewidth

linewidth = pzList*0;
%loop over all pz's;
%get corresponding parameters from loadData_multiRun_beamLaser.m;
for i = 1:nMaxPz
    %get q and p for a single simulation
    clear q;
    clear p;
    q(:,:) = qMatrix(i,:,n0_nStore:nStore);
    p(:,:) = pMatrix(i,:,n0_nStore:nStore);
    %dim of g1 function vector
    realG1Pos = (q(:,1)'*q+p(:,1)'*p)/nTrajectory/4;
    %imagG1Pos = (p(:,1)'*q-q(:,1)'*p)/nTrajectory/4;
    x = linspace(0,(tmax-t0)/transitTime,size(realG1Pos,2))';
    f = coeffvalues(fit(x,realG1Pos','exp1','startpoint',[1,-1]));
    linewidth(i) = -f(2)/pi/gc/transitTime;
end

%for 2-D plots (v.s. pz)
figure(12);
set(gca,'FontSize',20);
% fMatrixMannual = [1.07, 1.22, 1.22, 1.35, 1.97, 2.22, ...
%     3.16, 3.18, 3.65, 6.84, 17.36, ];
% hold on;
semilogx(pzList/gc,linewidth);
% plot(pzList/gc, pzList/gc);
% hold off;
% plot(pzList/gc,fMatrix(:,1));
xlabel('kv/\Gamma_c','FontSize', 20);
ylabel('\Delta \nu/\Gamma_c','FontSize', 20);

%for 3-D plots
%figure(13)

fprintf('Linewidth shown. Press enter to see next.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g(2)
%g(2)function
%g2Mannual = [0.84, 1.54, 1.47, 1.47, 0.91, 1.78. 0.58, ...
% 2.67]
%g(2) at 0 (take the last data point as our sample data)
% temp1 = sum((qMatrix(:,nStore).^2+pMatrix(:,nStore).^2-2).^2)/nTrajectory;
% temp2 = sum((qMatrix(:,nStore).^2+pMatrix(:,nStore).^2-2))/nTrajectory;
% g2=(temp1-4*temp2)/temp2^2;
% formatSpec = 'The g2(0) value is %f \n';
% fprintf(formatSpec, g2);
% fprintf('Program paused. Press enter to continue.\n');