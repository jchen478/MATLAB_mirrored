% specify simulation cases for redispersion
clc; 
clear;
close all;

% path
dataPath = 'data_redispersion/';

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
aArr = [2, 3, 4, 5, 6, 7];
muArr = [5, 10];

nMu = length(muArr);
nA = length(aArr); 