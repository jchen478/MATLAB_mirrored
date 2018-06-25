% specify simulation cases for redispersion
clc; 
clear;
close all;

% path
dataPath = 'data_kinetic/';

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

% parameter space
muArr = [15];
attArr = [0];
mukinArr = [0 0.2 0.4 0.6 0.8 1.0 5.0];

muC = caseArr('$\mu$',muArr,1);
attC = caseArr('$A_N$',attArr,2);
mukinC = caseArr('$\mu_{kin}$',mukinArr,3);


