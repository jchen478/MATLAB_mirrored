function [ I ] = intensityCalc( rx, ry, rz, side, nbin, nsegTot, rps )
%INTENSITYCALC Summary of this function goes here
%   Detailed explanation goes here

dbin = side/nbin;
bin = zeros(nbin,nbin,nbin);
Vseg = 2*pi*rps; 
Vcell = dbin*dbin*dbin; 
phiavg = 2*pi*rps*nsegTot / (side*side*side); 

%%{
% place fiber segments into bin
xloc = floor(rx/dbin) + 1;
yloc = floor(ry/dbin) + 1;
zloc = floor(rz/dbin) + 1;
for mi=1:nsegTot
    bin(xloc(mi),yloc(mi),zloc(mi)) = bin(xloc(mi),yloc(mi),zloc(mi)) + 1; 
end
% calculate local volume fraction
bin = bin * Vseg/Vcell;

if (max(max(max(bin))) > 1)
    display('Local volume fraction exceeded 1')
end

% calculate intensity
bin = (bin - phiavg).^2; 
I = 1/nbin^3 * sum(sum(sum(bin))) / (phiavg*(1-phiavg)); 
%}

end

