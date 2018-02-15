function [Seff,rp,nL3,L,gamma] = pCalc(nfib, nseg, rps, kb, side, a, EY, Imom, eta0)
% Calculate relevant simulation parameters

Seff = pi.*kb./nseg.^4;
rp = rps.*nseg;
nL3 = nfib .* (2.*rp./side).^3;

L = 2*rps.*nseg*a;
gamma = EY*Imom ./ (Seff.*L.^4*eta0);

end

