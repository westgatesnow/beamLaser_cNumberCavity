%cNumberCavity
%plot_IvsSz_multiRun_beamLaser

%for 2-D plots (v.s. tau)
figure(8);
set(gca,'FontSize',20);
%plot intensitySS
subplot(3,1,1);
h81 = errorbar(tauList, intensitySS_mean, intensitySS_std, 'o-');%used to logplot with xbar
set(get(h81,'Parent'), 'XScale', 'log');
xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
ylabel('I/\Gamma_c');
%plot nAtomAve./tauList/2.*(1-szFinalSS)
intensitySS_predicted = nAtomAve./tauList'/2.*(1-szFinalSS_mean);
subplot(3,1,2);
h82 = errorbar(tauList, intensitySS_predicted, nAtomAve./tauList/2.*szFinalSS_std','o-');
set(get(h82,'Parent'), 'XScale', 'log');
xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
ylabel('\I`/\Gamma_c');
subplot(3,1,3);
h83 = semilogx(tauList, (intensitySS_predicted-intensitySS_mean)./intensitySS_mean,'o-');
xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
ylabel('percentage');
%for 3-D plots
%figure(9)
