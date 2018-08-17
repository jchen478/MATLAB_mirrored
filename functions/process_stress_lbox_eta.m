function [r, eta, etaE] = process_stress_lbox_eta( filePrefix, ithStat,sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0)
%PROCESS_INTENSITY Find stats for intensity
r = 0;
eta = 0;
etaE = 0;
if (exist([filePrefix,'Stress_tensor.txt'], 'file') ~= 0)
    [stress_strain, sigxz_out] = read_stress([filePrefix,'Stress_tensor.txt']);
    % define lboxArr to calculate stress
    [box_strain, sidex] = read_box([filePrefix,'Lbox.txt']);
    lboxArr = find_lboxArr(stress_strain, box_strain, sidex);
    % find relevant parameters for stress calculations
    nPt = length(stress_strain);
    [Seff,rp,nL3,L,gamma] = pCalc(nfib*ones(nPt,1),...
        nseg*ones(nPt,1), rps, kb, lboxArr, a, EY, Imom, eta0);
    
    % Dimensionalize stresses
    [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
        stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
        zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
    
    etarel = sigma./gamma;
    % Non-dimensionalize
    sigmap_nondim = sigmap.*L.^4/EY/Imom;
    
    %% calculate interval averages
    r = round(intervals(box_strain,sidex),0);
    %     r(end-1) = 1156;
    %     r(end) = 1256
    % if simulation has not finished
    if round(stress_strain(end),1) ~= round(box_strain(end),1)
        stress_strain(end)
        return
    end
    
    %     figure()
    %     plot(stress_strain, smooth(etarel,50))
    %     xlim([0 stress_strain(end)])
    
    %     sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r');
    eta_stat = interval_average(stress_strain,etarel,r');
    
    eta = eta_stat(ithStat,1);
    etaE = eta_stat(ithStat,2);
    
%     figure()
%     hold on
%     title(filePrefix)
%     plot(stress_strain, smooth(etarel),'color',rgb('MediumBlue'));
%     xlabel('$\gamma$')
%     ylabel('$\eta$')
    
end
% if (exist([filePrefix,'extra_Stress_tensor.txt'], 'file') ~= 0)
%     [stress_strain_extra, sigxz_out_extra] = read_stress([filePrefix,'extra_Stress_tensor.txt']);
%     % define lboxArr to calculate stress
% %     [box_strain, sidex] = read_box([filePrefix,'Lbox.txt']);
% %     lboxArr = find_lboxArr(stress_strain_extra, box_strain, sidex);
%     % find relevant parameters for stress calculations
%     nPt = length(stress_strain_extra);
%     [Seff,rp,nL3,L,gamma] = pCalc(nfib*ones(nPt,1),...
%         nseg*ones(nPt,1), rps, kb, lboxArr(end)*ones(nPt,1), a, EY, Imom, eta0);
%     % Dimensionalize stresses
%     [sigmap_extra,sigma_extra, std_sigmap_extra, std_sigma_extra, N1_extra, N2_extra] = ...
%         stressDim(sigxz_out_extra, zeros(nPt,1), zeros(nPt,1), ...
%         zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
%     
%     etarel_extra = sigma_extra./gamma;
%     
%     figure()
%     hold on
%     title(filePrefix)
%     plot(stress_strain_extra+stress_strain(end), smooth(etarel_extra),'color',rgb('Crimson'));
%     xlabel('$\gamma$')
%     ylabel('$\eta$')
    
%     if (stress_strain(end) < 150)
%         return;
%     end
%     eta_stat = interval_average(stress_strain,etarel,[100 stress_strain(end)]);
%     eta = eta_stat(1,1);
%     etaE = eta_stat(1,2);
    
% end
end

