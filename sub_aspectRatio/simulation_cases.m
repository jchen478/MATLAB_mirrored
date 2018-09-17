% specify simulation cases for redispersion
clc; 
clear;
close all;

% path
% dataPath = 'data_aspectRatio/';
dataPath = '../data_all_redispersion/';
dataPathBasis = '../data_basis/';

% case
shape = 'theta6';
nseg = 5; 
a = 1.60E-05;       % radius (m)
Imom = pi*a^4/4;    % area moment (m^4)
EY = 1e9;           % Young's modulus (Pa m^2)
eta0 = 1;           % fluid viscosity (Pa s)
sidex = 600; 

% parameter space
% kb = caseArr('$k_b$',[50.625 10 3.164],4);
% rps = caseArr('$r_{ps}$',[10 15 20],4);
% rp = caseArr('$r_p$',[50 75 100],1);
% mu = caseArr('$\mu$',[5 10 15],2);
% att = caseArr('$A_N$',[0 9 20 50],3);

nfib = caseArr('$N_{fib}$',[1920 1280],4);
kb = caseArr('$k_b$',[50.625 10 ],4);
rps = caseArr('$r_{ps}$',[10 15 ],4);
rp = caseArr('$r_p$',[50 75 ],1);
mu = caseArr('$\mu$',[5 10 15],2);
att = caseArr('$A_N$',[0 9 20 30 50],3);

% mu = caseArr('$\mu$',[0 5 10 15],2);
% att = caseArr('$A_N$',[0 9 20 30 50],3);


