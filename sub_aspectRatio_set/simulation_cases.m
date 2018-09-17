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
% rp = caseArr('$r_p$',50:5:100,1);
% cond = caseArr('set ',{'mu5_att0','mu5_att20','mu5_att50','mu15_att20'},2);
% cond.legend = {'$(\mu,A_N) = (5,0)$','$(\mu,A_N) = (5,20)$','$(\mu,A_N) = (5,50)$','$(\mu,A_N) = (15,20)$'};
% 
% condMu = caseArr('mu',{'5','5','5','15'},3);
% condAtt = caseArr('att',{'0','20','50','20'},3);

rp = caseArr('$r_p$',50:5:100,1);
cond = caseArr('set ',{'mu5_att0','mu5_att20','mu15_att20'},2);
cond.legend = {'$(\mu,A_N) = (5,0)$','$(\mu,A_N) = (5,20)$','$(\mu,A_N) = (15,20)$'};

condMu = caseArr('mu',{'5','5','15'},3);
condAtt = caseArr('att',{'0','20','20'},3);

nfib = caseArr('$N_{fib}$',[1920 1760 ...
    1600 1472 ...
    1376 1280 ...
    1216 1120 ...
    1056 1024 ...
    960],3);
kb = caseArr('$k_b$',[50.63 34.58 ...
    24.41 17.73 ...
    13.18 10.0 ...
    7.72 6.06 ...
    4.82 3.88 ...
    3.16],3);
rps = caseArr('$r_{ps}$',10:1:20,3);

