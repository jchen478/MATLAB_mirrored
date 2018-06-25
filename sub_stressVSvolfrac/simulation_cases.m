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
muArr = [0, 5, 10, 15];
attArr = [0, 9, 20, 50];
% attArr = [0, 9, 20, 30, 50];
sidex = [600 476.2 416 378 351];

% volume fraction
volfracArr = nfib*nseg*2*rps*pi./sidex.^3*100;

% weight percent
weightfracArr = 1.6*volfracArr;

nMu = length(muArr);
nNL3 = length(nL3Arr); 
nAtt = length(attArr); 

muC = caseArr('$\mu$',muArr,1);
attC = caseArr('$A_N$',attArr,2);
volfracC = caseArr('$\phi (\%)$',round(volfracArr,2),3);
weightfracArr = caseArr('$C (w/w\%)$',round(weightfracArr,2),3);