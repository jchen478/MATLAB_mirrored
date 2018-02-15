clc; 
close all;

mu_stat = 1.0;
mu_kin = 0.0; 
alpha = sqrt(mu_stat*(mu_stat-mu_kin));
nfric = 100; 
eps = 1e-4; 

npt = 300; 
dU = linspace(-5,5,npt)';

mu = (-mu_kin*dU.*sqrt(dU.^2+eps/nfric^2)-2*alpha/nfric*dU)./(dU.^2+1/nfric^2); 

plot(dU,mu)