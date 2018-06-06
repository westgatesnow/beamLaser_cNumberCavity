%cNumberCavity
%plot_sz_multiRun_beamLaser

%get E[inversionAve_s] and V[inversionAve_s]
inversionAveSS = inversionAve(:,:,n0_nstore:nstore);
inversionAveSS_mean = mean(inversionAveSS, 3);
inversionAveSS_std = std(inversionAveSS, 0, 3);
%get E[szFinal_s] and V[szFinal_s]
szFinalSS = szFinal(:,:,n0_nTimeStep:nTimeStep);
szFinalSS_mean = mean(szFinalSS, 3, 'omitnan');
szFinalSS_std = std(szFinalSS, 0, 3, 'omitnan');

%2-D plots (v.s. tau)
if nMaxNAtom == 1
    figure(6);
    set(gca,'FontSize',20);
    subplot(2,1,1);
    h61 = errorbar(tauList, inversionAveSS_mean, inversionAveSS_std);%used to logplot with xbar
    set(get(h61,'Parent'), 'XScale', 'log');
    xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
    ylabel('\langle j^z \rangle_e');
    subplot(2,1,2);
    h62 = errorbar(tauList, szFinalSS_mean, szFinalSS_std);
    set(get(h62,'Parent'), 'XScale', 'log');
    xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
    ylabel('\langle j^z_{out} \rangle_e');
end
%2-D plots (v.s. nAtomAve)
if nMaxTau == 1
    figure(6);
    set(gca,'FontSize',20);
    subplot(2,1,1);
    h61 = errorbar(nAtomAveList, inversionAveSS_mean, inversionAveSS_std);%used to logplot with xbar
    set(get(h61,'Parent'), 'XScale', 'log');
    xlabel('Mean Atom Number','FontSize', 20);
    ylabel('\langle j^z \rangle_e');
    subplot(2,1,2);
    h62 = errorbar(nAtomAveList, szFinalSS_mean, szFinalSS_std);
    set(get(h62,'Parent'), 'XScale', 'log');
    xlabel('Mean Atom Number','FontSize', 20);
    ylabel('\langle j^z_{out} \rangle_e');
end

%for 3-D plots
%figure(7)
