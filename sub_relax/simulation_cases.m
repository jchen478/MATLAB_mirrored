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
sidex = 600; 
rps = 15;
kb = 10;

% parameter space
shear = caseArr('case',{'rand_shear_shear','shear_shear',...
    'shear_noshear','noshear_shear','noshear_noshear'},1);
shear.legend = {'0: Rand SS','1: SS','2: SN',...
    '3: NS','4: NN'};
shearplot = caseArr('case',[0 1 2 3 4],1);
shearplot.legend = shear.legend;
mu = caseArr('$\mu$',[5 10 15],2);
att = caseArr('$A_N$',[0 9 20 50],3);


