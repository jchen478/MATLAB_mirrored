clc;
clear;
close all;

% Define common parameters
a = 1.60E-05;       % radius (m)
Imom = pi*a^4/4;    % area moment (m^4)
EY = 1e9;           % Young's modulus (Pa m^2)
eta0 = 1;           % fluid viscosity (Pa s)

name = ['mu0.txt'];
File = fopen(name,'r');
data = fscanf(File,'%f',[14 Inf])';
fclose(File);

% Assign variable names to data
nfib = data(:,1); 
nseg = data(:,2); 
rps = data(:,3); 
kb = data(:,4); 
mu = data(:,5); 
side = data(:,6); 
sigmapstar = data(:,7); 
se_sigmapstar = data(:,8); 
N1star = data(:,9); 
N2star = data(:,10); 
nc = data(:,11); 
se_nc = data(:,12); 
elas = data(:,13); 
se_elas = data(:,14); 

% Dimensionalization calculations
[Seff,rp,nL3,L,gamma] = pCalc(nfib, nseg, rps, kb, side, a, EY, Imom, eta0);
[sigmap,sigma, std_sigmap, se_sigma, N1, N2] =...
    stressDim(sigmapstar, se_sigmapstar, N1star, N2star, ...
    nseg, rps, eta0, gamma, nL3);

% Calculate viscosity
etarel = sigma./gamma;
se_etarel = se_sigma./gamma;

% Plotting
figure('Units','Inches','Position',[0.5 4 4.5 3.5])
hold on
box on
title('\boldmath{$\mu_{stat} = 0$}')
plot(nL3,etarel-1,'-x')
xlabel('$nL^3$')
ylabel('$\eta_{sp}$')