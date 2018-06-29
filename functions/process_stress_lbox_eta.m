function [r, eta] = process_stress_lbox_eta( filePrefix, ithStat,sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0)
%PROCESS_INTENSITY Find stats for intensity
r = 0; 
eta = 0; 
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
    % if simulation has not finished
    if round(stress_strain(end),1) ~= round(box_strain(end),1)
        stress_strain(end)
        return
    end

    sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r');
    eta_stat = interval_average(stress_strain,etarel,r');
    
    eta = eta_stat(ithStat,1);
end

end

