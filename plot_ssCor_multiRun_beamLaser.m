%cNumberCavity
%plot_ssCor_multiRun_beamLaser

%2-D plots
%get E[ssCor_s] and V[ssCor_s]
ssCorSS = spinSpinCorAve(:,:,n0_nstore:nstore);
ssCorSS_mean = mean(ssCorSS, 3);
ssCorSS_std = std(ssCorSS, 0, 3);
%for 2-D plots with error bars (v.s. tau)
if nMaxNAtom == 1
    figure(10);
    h8 = errorbar(tauList, ssCorSS_mean, ssCorSS_std, 'o-');%used to logplot with xbar
    set(get(h8,'Parent'), 'XScale', 'log');
    axis([0 inf 0 0.15]);
    xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
    ylabel('\langle \sigma_i^+ \sigma_j^-\rangle');
end
%for 2-D plots with error bars (v.s. nAtomAve)
if nMaxTau == 1
    figure(10);
    h8 = errorbar(nAtomAveList, ssCorSS_mean, ssCorSS_std, 'o-');%used to logplot with xbar
    set(get(h8,'Parent'), 'XScale', 'log');
    axis([0 inf 0 0.15]);
    xlabel('Mean Atom Number','FontSize', 20);
    ylabel('\langle \sigma_i^+ \sigma_j^-\rangle');
end

%for 3-D plots
%figure(11)