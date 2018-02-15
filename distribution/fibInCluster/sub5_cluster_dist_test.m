%%%
%%% Plot averaged values
%%% -- properties include intensity, nc, sigP, eta, N1, N2, Eelas
%%%

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%%%%%%%%%%%%%%%%%%%% Theta1 fibers %%%%%%%%%%%%%%%%%%%%%
nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200 1293];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
thetaArr = 1;
fileNameArr = {'theta1'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('Olive') rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue') rgb('Plum') rgb('Purple') };
%}

%{
%%%%%%%%%%%%%%%%%%%%% Theta3 fibers %%%%%%%%%%%%%%%%%%%%%
nfibArr = [160 240 320 640 1280 3200 6400 10240];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
thetaArr = 3;
fileNameArr = {'theta3'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('Olive') rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue') rgb('Plum') rgb('Purple') };
%}

%{
%%%%%%%%%%%%%%%%%%%%% Theta6 fibers %%%%%%%%%%%%%%%%%%%%%
nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200 1293];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
thetaArr = 6;
fileNameArr = {'theta6'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('Olive') rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue') rgb('Plum') rgb('Purple') };
%}

%{
%%%%%%%%%%%%%%%%%%%%% Helical fibers %%%%%%%%%%%%%%%%%%%%%
nfibArr = [160 240  320 640 1280 3200 6400];
lboxArr = [300 343.4 378 476.2 600 814.3  1026];
muArr = [0 1 2 3 4 5 10 15 20];
thetaArr = [3];
fileNameArr = {'helical'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('DarkGreen') rgb('LightSkyBlue') rgb('Plum')};
%}

%%{
%%%%%%%%%%%%%%%%%%%%% Straight fibers %%%%%%%%%%%%%%%%%%%%%
nfibArr = [160 240 320 640 1280 3200 6400];
lboxArr = [300 343.4 378 476.2 600 814.3 1026];
muArr = [0 1 2 3 4 5 10 15 20];
thetaArr = 0;
fileNameArr = {'theta0'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('DarkGreen') rgb('LightSkyBlue') rgb('Plum')};
%}

figStart = 1;
nMu = length(muArr);
nNfib = length(nfibArr);
muLegendArr = cell(nMu,1);
nfibLegendArr = cell(nNfib,1);
for j=1:nMu
    muLegendArr{j} = ['$\mu = $ ',num2str(muArr(j))];
end
for j=1:nNfib
    nfibLegendArr{j} = ['$N_{fib} = $ ',num2str(nfibArr(j))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read files and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
% Plot for each Nfib sum of bins
figure('units','normalized','outerposition',[0.05 0.1 0.95 0.75])
for i=1:nNfib
    subplot(2,4,i)
    title(['$N_{fib} =$ ', num2str(nfibArr(i))])
    hold on
    for j=1:nMu
        name = ['data_nfibInCluster/',fileNameArr{1},'_Sdist_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'.txt'];
        data = csvread(name);
        data(:,1:4) = [];
        Sbin = size(data,2);
        SbinArr = 2:Sbin+1;
        Shist = mean(data,1);
        plot(SbinArr,Shist,'-.o','MarkerSize',10, 'linewidth',2.5, 'color',colorArr{j})
    end
end

subplot(2,4,1)
hold on
legend(muLegendArr,'location','best')
for i=1:nNfib
    subplot(2,4,i)
    hold on
    set(gca,'fontsize',16)
    xlabel('$S$')
    ylabel('$n_S$')
    xlim([2 inf])
end
figStart = figStart + 1;
%}

%%{
% Plot for histogram and compare nfib
for j=1:nMu
   figure('units','normalized','outerposition',[0.2 0.2 0.5 0.5])
    hold on
    for i=1:nNfib        
        name = ['data_nfibInCluster/',fileNameArr{1},'_Sdist_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'.txt'];
        data = csvread(name);
        data(:,1:4) = [];
        Sbin = size(data,2);
        SbinArr = 2:Sbin+1;
        Shist = mean(data,1);
        plot(SbinArr,Shist/nfibArr(i),'-.o','MarkerSize',10, 'linewidth',2.5)      
        title(['$\mu = $ ',num2str(muArr(j))])

    end
end
%
for i=figStart:figStart+nMu-1
    figure(i)
    xlabel('$S$')
    ylabel('$n_S / N_{fib}$')
    set(gca,'yscale','log')
    legend(nfibLegendArr,'location','best')
end
figStart = figStart + nMu;
%}