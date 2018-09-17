function [block_strain, block_etarel] = block_stress_lbox_eta_basis( filePrefix, blockSize, ithStat, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0)
%PROCESS_INTENSITY Find stats for intensity
block_strain = 0;
block_etarel = 0;
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
  
    block_strain = block_average(stress_strain, blockSize);
    block_etarel = block_average(etarel, blockSize);
    
end

end

