%cNumberCavity
%plot_intensity_multiRun_beamLaser

%we first look at the two extreme cases
% figure(3);
% set(gca,'FontSize',20);
% subplot(2,1,1);
% intensity_i(:) = intensity(1,end,:);
% plot(linspace(0,tmax,nstore)/transitTime,intensity_i);
% xlabel('t/T','FontSize', 20);
% ylabel('\kappa \langle a^+(t) a(t) \rangle');
% subplot(2,1,2);
% intensity_f(:) = intensity(end,1,:);
% plot(linspace(0,tmax,nstore)/transitTime,intensity_f);
% xlabel('t/T','FontSize', 20);
% ylabel('\kappa \langle a^+(t) a(t) \rangle');

%2-D plots
%get E[I_s] and V[I_s]
intensitySS = intensity(:,:,n0_nstore:nstore)/gc;%in units of gc
intensitySS_mean = mean(intensitySS, 3);
intensitySS_std = std(intensitySS, 0, 3);
%for 2-D plots with error bars (v.s. tau)
if nMaxNAtom == 1
    figure(4);
    h4 = errorbar(tauList, intensitySS_mean, intensitySS_std, 'o-');%used to logplot with xbar
    set(get(h4,'Parent'), 'XScale', 'log');
    xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
    ylabel('I/\Gamma_c');
end
%for 2-D plots with error bars (v.s. nAtomAve)
if nMaxTau == 1
    figure(4);
    h4 = errorbar(nAtomAveList, intensitySS_mean, intensitySS_std, 'o-');%used to logplot with xbar
    set(get(h4,'Parent'), 'XScale', 'log');
    xlabel('Mean Atom Number','FontSize', 20);
    ylabel('I/\Gamma_c');
end
%for 3-D plots
%figure(5)
