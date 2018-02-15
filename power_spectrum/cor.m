%% correlation function
function [cAB] = cor(A,B)
% Calculates the 2D cross correlation for two 2D fields A and B
% A and B must be the same size
% Correlation normalized by RMS values and mesh size

N=size(A);
cAB=fftshift(ifft2(conj(fft2(A)).*fft2(B)))./(sqrt(mean2(A.^2))*sqrt(mean2(B.^2))*N(1)*N(2));

end