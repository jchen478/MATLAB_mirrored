% specify simulation cases for redispersion
clc; 
clear;
close all;

% path
dataPath = 'data_redispersion/';
% dataPath = 'data_aniso/aniso_';

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
aArr = [2, 3, 4, 5];
muArr = [0, 5, 10, 15];
% attArr = [0, 9, 20, 30, 50];
attArr = [0, 9, 20, 50];
% aniso set
% aArr = [2, 3, 4, 5];
% muArr = [0, 5];
% attArr = [0, 9, 20];

nMu = length(muArr);
nA = length(aArr); 
nAtt = length(attArr);

muC = caseArr('$\mu$',muArr,1);
attC = caseArr('$A_N$',attArr,2);
aC = caseArr('$a$',aArr,3);


