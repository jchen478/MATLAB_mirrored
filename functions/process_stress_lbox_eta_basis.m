function [r, eta, etaE] = process_stress_lbox_eta_basis( filePrefix, basisStrain, ithStat, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0)
%PROCESS_INTENSITY Find stats for intensity
r = 0; 

eta = 0; 
etaE = 0; 
if (exist([filePrefix,'Stress_tensor.txt'], 'file') ~= 0)
    [stress_strain, sigxz_out] = read_stress([filePrefix,'Stress_tensor.txt']);
    
    nPt = length(stress_strain);
    [Seff,rp,nL3,L,gamma] = pCalc(nfib*ones(nPt,1),...
        nseg*ones(nPt,1), rps, kb, sidex*ones(nPt,1), a, EY, Imom, eta0);
    
    % Dimensionalize stresses
    [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
        stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
        zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
    
    etarel = sigma./gamma;
    % Non-dimensionalize
    sigmap_nondim = sigmap.*L.^4/EY/Imom;
    if (stress_strain(end) <= basisStrain)
        return
    end
    r = [basisStrain stress_strain(end)];
%     sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r');
    eta_stat = interval_average(stress_strain,etarel,r');
    
    eta = eta_stat(ithStat,1);
    etaE = eta_stat(ithStat,2);
end

end

