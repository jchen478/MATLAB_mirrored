function [r, block_strain, block_etarel] = block_stress_lbox_eta( filePrefix, blockSize, ithStat,sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0)
% return block average of data with specified block size
r = 0;
block_strain = [];
block_etarel = [];
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
    % if simulation has not finished
    if round(stress_strain(end),1) ~= round(box_strain(end),1)
        stress_strain(end)
        return
    end
    
    r = round(intervals(box_strain,sidex),0);
    
end

%% discard regions not needed
stress_strain (stress_strain <= r(4)) = 0;
ind = 0;
for i=length(stress_strain):-1:0
    if stress_strain(i) == 0
        ind = i;
        break;
    end
end
stress_strain(1:ind) = [];
etarel(1:ind) = [];
stress_strain = stress_strain - stress_strain(1);

%% read in extra files
if (exist([filePrefix,'extra_Stress_tensor.txt'], 'file') ~= 0)
    [stress_strain_extra, sigxz_out_extra] = read_stress([filePrefix,'extra_Stress_tensor.txt']);
    % find relevant parameters for stress calculations
    nPt = length(stress_strain_extra);
    [Seff,rp,nL3,L,gamma] = pCalc(nfib*ones(nPt,1),...
        nseg*ones(nPt,1), rps, kb, lboxArr(end)*ones(nPt,1), a, EY, Imom, eta0);
    % Dimensionalize stresses
    [sigmap_extra,sigma_extra, std_sigmap_extra, std_sigma_extra, N1_extra, N2_extra] = ...
        stressDim(sigxz_out_extra, zeros(nPt,1), zeros(nPt,1), ...
        zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
    
    etarel_extra = sigma_extra./gamma;
    
    if (length(stress_strain) == 0)
        return;
    end
end

%% concatenate results with extra files
% if (exist([filePrefix,'extra_Stress_tensor.txt'], 'file') ~= 0)
%     stress_strain_extra = stress_strain_extra + stress_strain(end);
%     strain_total = [stress_strain ; stress_strain_extra];
%     etarel_total = [etarel ; etarel_extra];
% else
    strain_total = stress_strain;
    etarel_total = etarel;
% end
% 
block_strain = block_average(strain_total, blockSize);
block_etarel = block_average(etarel_total, blockSize);

end


