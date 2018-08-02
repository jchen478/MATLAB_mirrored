[stress_strain, sigxz_out] = read_stress(['Stress_tensor.txt']);
[NC_strain, NC] = read_NC(['Number_of_Contacts.txt']);

close all
figure(1)
plot(stress_strain, sigxz_out)
figure(2)
plot(NC_strain, NC)


% % find relevant parameters for stress calculations
% nPt = length(stress_strain);
% [Seff,rpdum,nL3,L,gamma] = pCalc(nfib*ones(nPt,1),...
%     nseg*ones(nPt,1), rps, kb,...
%     sidex*ones(nPt,1), a, EY, Imom, eta0);
% 
% % Dimensionalize stresses
% [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
%     stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
%     zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
% 
% % Calculate viscosity
% etarel = sigma./gamma;