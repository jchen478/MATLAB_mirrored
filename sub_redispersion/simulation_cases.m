% specify simulation cases for redispersion
clc; 
clear;
close all;

% path
dataPath = '../data_all_redispersion/';
dataPathBasis = '../data_basis/';

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
sidex = 600; 
volfrac = 2*pi*rps*nseg*nfib/sidex^3; 

% parameter space
% aArr = [2, 3, 4, 5];
aArr = [2 3 4 5];
muArr = [0, 5, 10, 15];
attArr = [0, 9, 20, 30, 35, 50];

nMu = length(muArr);
nA = length(aArr); 
nAtt = length(attArr);

muC = caseArr('$\mu$',muArr,1);
attC = caseArr('$A_N$',attArr,2);
aC = caseArr('$r_{\phi}$',aArr,3);


