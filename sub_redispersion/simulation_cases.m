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
aArr = [2, 3, 4, 5];
muArr = [0, 5, 10, 15];
attArr = [0, 9, 20, 30, 50];

nMu = length(muArr);
nA = length(aArr); 
nAtt = length(attArr);

muColorArr = {rgb('Crimson'), rgb('DarkOrange'), ...
    rgb('DarkGreen'), rgb('MediumBlue')};
aColorArr = {rgb('Crimson'), rgb('DarkOrange'), ...
    rgb('DarkGreen'), rgb('MediumBlue')};

muLegend = cell(nMu,1);
aLegend = cell(nA,1);
for i=1:nMu
    muLegend{i} = ['$\mu = $ ',num2str(muArr(i))];
end
for i=1:nA
    aLegend{i} = ['$a = $ ',num2str(aArr(i))];
end