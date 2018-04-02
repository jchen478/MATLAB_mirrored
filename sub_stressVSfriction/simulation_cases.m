%{
%%%%%%%%%%%%%%%%%%%%% U-shaped fibers %%%%%%%%%%%%%%%%%%%%%
% fileNameArr = {'theta0'};thetaArr = 0;
% fileNameArr = {'theta1'}; thetaArr = 1;
% fileNameArr = {'theta3'}; thetaArr = 3;
% fileNameArr = {'theta6'}; thetaArr = 6;

nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200 1293];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('Olive') rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue') rgb('Plum') rgb('Purple') };
%}
%%{
%%%%%%%%%%%%%%%%%%%%% Helical fibers %%%%%%%%%%%%%%%%%%%%% 
nfibArr = [160 240  320 640 1280 3200 6400];
lboxArr = [300 343.4 378 476.2 600 814.3  1026];
muArr = [0 1 2 3 4 5 10 15 20];
fileNameArr = {'helical'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('DarkGreen') rgb('LightSkyBlue') rgb('Plum')};
%}

% fiber dimensions
a = 1.60E-05;       % radius (m)
Imom = pi*a^4/4;    % area moment (m^4)
EY = 1e9;           % Young's modulus (Pa m^2)
eta0 = 1;           % fluid viscosity (Pa s)

% nTheta=1;
nMu = length(muArr);
nLbox = length(lboxArr);
nNfib = length(nfibArr);
muLegendArr = cell(nMu,1);
thetaNfibLegendArr = cell(nNfib,1);


for i=1:nMu
    muLegendArr{i} = ['$\mu =$ ',num2str(muArr(i))];
end
if strcmpi(fileNameArr,'helical')
        for j=1:nNfib
            thetaNfibLegendArr{j} = ['$(\theta_{eq},\phi_{eq},N_{fib}) =$ (0.8, 0.7, ',num2str(nfibArr(j)),')'];
        end
else
        for j=1:nNfib
            thetaNfibLegendArr{j} = ['$(\theta_{eq},N_{fib}) =$ (0.',num2str(thetaArr(i)),', ',num2str(nfibArr(j)),')'];
        end
end