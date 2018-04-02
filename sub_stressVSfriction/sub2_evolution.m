%%%
%%% Plot evolution of properties
%%% -- properties include nc, sigP
%%%

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%%%%%%%%%%%%%%%%%%%% U-shaped fibers %%%%%%%%%%%%%%%%%%%%%
% fileNameArr = {'theta0'};
% fileNameArr = {'theta1'};
% fileNameArr = {'theta3'}; thetaArr = 3;
fileNameArr = {'theta6'}; thetaArr = 6;

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


% nfibArr = [160 240  320 640 1280 ];
% lboxArr = [300 343.4 378 476.2 600];
% muArr = [0 1 2 3 4 5 10 15 20];

thetaArr = [3];
fileNameArr = {'helical'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('DarkGreen') rgb('LightSkyBlue') rgb('Plum')};
%}

xLowLim = 1500;
xUpLim = 2000;

% fiber dimensions
a = 1.60E-05;       % radius (m)
Imom = pi*a^4/4;    % area moment (m^4)
EY = 1e9;           % Young's modulus (Pa m^2)
eta0 = 1;           % fluid viscosity (Pa s)
nseg = 5;           % number of segments
rps = 15;
kb = 10;

nTheta = length(thetaArr);
nMu = length(muArr);
nLbox = length(lboxArr);
nNfib = length(nfibArr);
muLegendArr = cell(nMu,1);
thetaNfibLegendArr = cell(nTheta*nNfib,1);
markersize = 50;

for i=1:nMu
    muLegendArr{i} = ['$\mu =$ ',num2str(muArr(i))];
end
if strcmpi(fileNameArr,'helical')
    for i=1:nTheta
        for j=1:nNfib
            thetaNfibLegendArr{(i-1)*nNfib+j} = ['$(\theta_{eq},\phi_{eq},N_{fib}) =$ (0.8, 0.7, ',num2str(nfibArr(j)),')'];
        end
    end
else
    for i=1:nTheta
        for j=1:nNfib
            thetaNfibLegendArr{(i-1)*nNfib+j} = ['$(\theta_{eq},N_{fib}) =$ (0.',num2str(thetaArr(i)),', ',num2str(nfibArr(j)),')'];
        end
    end
end

for j=1:nNfib
    figure('Units','Inches','Position',[3 3 3.5 6.5])
    set(gca,'XMinorTick','on')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NC

dataPath = '../data_stressVSfriction/NC/';
% dataPath = '../data_stressVSfriction/NC_replicate1/';
% dataPath = '../data_stressVSfriction/NC_replicate2/';
%%{
for j=1:nNfib
    figure(j)
    subplot(2,1,1)
    set(gca,'XMinorTick','on')
    hold on;
    ylabel('$N_c$')
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_NC_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File = fopen(name,'r');
        data = fscanf(File,'%f',[5 Inf])';
        fclose(File);
        diff = data(2,1) - data(1,1);
        % make sure that strain is increasing
        for ii=2:length(data)
            if (data(ii,1) < data(ii-1,1))
                data(ii,1) = data(ii-1,1) + diff;
            end
        end
        % block average
        block = 25;
        npt = floor(length(data)/block);
        B = zeros(npt,1);
        Bt = zeros(npt,1);
        ind = 1;
        for ii=1:block:length(data)-block+1
            B(ind) = mean(data(ii:ii+block-1,4));
            Bt(ind) = mean(data(ii:ii+block-1,1));
            ind = ind + 1;
        end
        plot(Bt,B,'-o','color',colorArr{i},...
            'MarkerFaceColor',colorArr{i},...
            'MarkerEdgeColor',colorArr{i},...
            'MarkerSize',2)
    end
    box on
    xlabel('$\gamma$')
    xlim([xLowLim xUpLim])
    legend(muLegendArr,'location','best')
    title([fileNameArr{1},' $N_{fib} =$ ',num2str(nfibArr(j))])
end
%}

% % Stress
dataPath = '../data_stressVSfriction/Stress/';
% dataPath = '../data_stressVSfriction/Stress_replicate1/';
% dataPath = '../data_stressVSfriction/Stress_replicate2/';

for j=1:nNfib
    figure(j)
    subplot(2,1,2)
    hold on;
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_stress_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File = fopen(name,'r');
        data = fscanf(File,'%f',[7 Inf])';
        fclose(File);
        diff = data(2,1) - data(1,1);
        nStep = length(data);
        % make sure that strain is increasing
        for ii=2:nStep
            if (data(ii,1) < data(ii-1,1))
                data(ii,1) = data(ii-1,1) + diff;
            end
        end
        
        nfib = nfibArr(j)*ones(nStep,1);
        nsegArr = nseg*ones(nStep,1);
        side = lboxArr(j)*ones(nStep,1);
        % calculate relevant parameters
        [Seff,rp,nL3,L,gamma] = pCalc(nfib, nsegArr, rps, kb, side, a, EY, Imom, eta0);
        
        % Dimensionalize stresses
        [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
            stressDim(data(:,4), zeros(nStep,1), zeros(nStep,1), ...
            zeros(nStep,1), nseg, rps, eta0, gamma, nL3);
        
        
        % Non-dimensionalize
        sigmap_nondim = sigmap.*L.^4/EY/Imom;
        strain = data(:,1);
        
        % block average
        block = 25;
        npt = floor(length(sigmap_nondim)/block);
        B = zeros(npt,1);
        Bt = zeros(npt,1);
        ind = 1;
        for ii=1:block:length(sigmap_nondim)-block+1
            B(ind) = mean(sigmap_nondim(ii:ii+block-1));
            Bt(ind) = mean(strain(ii:ii+block-1));
            ind = ind + 1;
        end
        plot(Bt,B,'-o','color',colorArr{i},...
            'MarkerFaceColor',colorArr{i},...
            'MarkerEdgeColor',colorArr{i},...
            'MarkerSize',2)
    end
    box on
    xlabel('$\gamma$')
    ylabel('$\sigma_{p,xz} L^4/ E_Y I$')
    xlim([xLowLim xUpLim])
    legend(muLegendArr,'location','best')
    title([fileNameArr{1},' $N_{fib} =$ ',num2str(nfibArr(j))])
end

