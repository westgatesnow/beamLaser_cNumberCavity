%cNumberCavity
%plot_nAtom_multiRun_beamLaser

%2-D plots:
%get E[N_s] and V[N_s]
nAtomSS = nAtom(:,:,n0_nstore:nstore);
nAtomSS_mean = mean(nAtomSS, 3);
nAtomSS_std = std(nAtomSS, 0, 3);
%for 2-D plots (v.s. tau)
if nMaxNAtom == 1
    figure(1);
    set(gca,'FontSize',20);
    subplot(2,1,1);
    plot(tauList,nAtomSS_mean);
    xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
    ylabel('N');
    subplot(2,1,2);
    plot(tauList,nAtomSS_std);
    xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
    ylabel('\Delta N');
end
%for 2-D plots (v.s. nAtomAve)
if nMaxTau == 1
    figure(1);
    set(gca,'FontSize',20);
    subplot(2,1,1);
    plot(nAtomAveList,nAtomSS_mean);
    xlabel('Mean Atom Number','FontSize', 20);
    ylabel('N');
    subplot(2,1,2);
    plot(nAtomAveList,nAtomSS_std.^2);
    xlabel('Mean Atom Number','FontSize', 20);
    ylabel('\Delta N^2');
end

%3-D plots
%figure(2)
