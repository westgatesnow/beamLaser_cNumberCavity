%cNumberCavity
%plot_linewidth_multiRun_beamLaser

fMatrix = zeros(nMaxTau,nMaxNAtom);
%loop over all tau's and nAtom's;
%get corresponding parameters from loadData_multiRun_beamLaser.m;
for j = 1:nMaxNAtom
    for i = 1:nMaxTau
        %get current transitTime
        if i <= nMaxTau1
            tau = tauList1(i);
        else
            tau = tauList2(i-nMaxTau1);
        end
        %get q and p for a single simulation
        q(:,:) = qMatrix(i,j,:,n0_nstore:nstore);
        p(:,:) = pMatrix(i,j,:,n0_nstore:nstore);
        %dim of g1 function vector
        realG1Pos = (q(:,1)'*q+p(:,1)'*p)/nTrajectory/4;
        %imagG1Pos = (p(:,1)'*q-q(:,1)'*p)/nTrajectory/4;
        x = linspace(0,(tmax-t0)*gc/tau,size(realG1Pos,2))';
        %x = linspace(0,(tmax-t0)*gc/tau,size(realG1Pos,2))';for tau^-1unit
        %x = linspace(0,(tmax-t0),size(realG1Pos,2))';for gc unit
        f = coeffvalues(fit(x,realG1Pos','exp1','startpoint',[0.1,-0.1])); 
        fMatrix(i,j) = -f(2)/pi;
    end
end

%for 2-D plots (v.s. tau)
if nMaxNAtom == 1
    figure(12);
    set(gca,'FontSize',20);
    semilogx(tauList,fMatrix(:,1));
    xlabel('\tau/\Gamma_c^{-1}','FontSize', 20);
    ylabel('\Delta \nu/\tau^{-1}','FontSize', 20);
end
%for 2-D plots (v.s. tau)
if nMaxTau == 1
    figure(12);
    set(gca,'FontSize',20);
    semilogx(tauList,fMatrix(1,:));
    xlabel('Mean Atom Number','FontSize', 20);
    ylabel('\Delta \nu/\tau^{-1}','FontSize', 20);
end


%for 3-D plots
%figure(13)
