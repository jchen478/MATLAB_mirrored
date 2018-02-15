function [sigmap,sigma, se_sigmap, se_sigma, N1, N2] = stressDim(sigmapstar, se_sigmapstar, N1star, N2star, nseg, rps, eta0, gamma, nL3)
% Obtain dimensional stress from simulation output
%   Detailed explanation goes here

sigmap = sigmapstar*pi.*nL3./(6*nseg.^3.*log(2*rps))*eta0.*gamma;
se_sigmap = se_sigmapstar*pi.*nL3./(6*nseg.^3.*log(2*rps))*eta0.*gamma;
N1 =  N1star*pi.*nL3./(6*nseg.^3.*log(2*rps))*eta0.*gamma;
N2 =  N2star*pi.*nL3./(6*nseg.^3.*log(2*rps))*eta0.*gamma;
sigma = sigmap + eta0*gamma; 
se_sigma = se_sigmap; 

end

