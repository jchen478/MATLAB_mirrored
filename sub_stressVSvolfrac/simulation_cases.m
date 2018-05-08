% specify simulation cases for stressVSvolfrac
clc; 
clear;
close all;

% path
dataPath = 'data_stressVSvolfrac/';

% case
shape = 'theta6';
nfib = 1280;
nseg = 5; 
a = 1.60E-05;       % radius (m)
Imom = pi*a^4/4;    % area moment (m^4)
EY = 1e9;           % Young's modulus (Pa m^2)
eta0 = 1;           % fluid viscosity (Pa s)
rps = 15;
kb = 10;

% parameter space
nL3Arr = [20, 40, 60, 80, 100];
muArr = [0, 5, 10, 20];
sidex = [600 476.2 416 378 351];

% volume fraction
volfracArr = nfib*nseg*2*rps*pi./sidex.^3*100;

nMu = length(muArr);
nNL3 = length(nL3Arr); 

% legends
muLegend = cell(nMu,1);

for i=1:nMu
    muLegend{i} = ['$\mu = $ ',num2str(muArr(i))];
end