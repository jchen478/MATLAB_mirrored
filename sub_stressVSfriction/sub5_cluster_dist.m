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

%%{
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

%{
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
figStart = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read files and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
% Plot for each Nfib sum of bins
figure('units','normalized','outerposition',[0.05 0.1 0.95 0.75])
for i=1:nNfib
    subplot(2,5,i)
    title(['$N_{fib} =$ ', num2str(nfibArr(i))])
    hold on
    for j=1:nMu
        name = ['../data_stressVSfriction/SDist/',fileNameArr{1},'_Sdist_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'.txt'];
        data = csvread(name);
        data(:,1:4) = [];
        Sbin = size(data,2);
        SbinArr = 2:Sbin+1;
        Shist = mean(data,1);
        plot(SbinArr,Shist,'-.o','MarkerSize',10, 'linewidth',2.5, 'color',colorArr{j})
    end
end

dataFile = {['fig5_cluster_dist/',fileNameArr{1},'_nfib_subplot']};

subplot(2,5,1)
hold on
legend(muLegendArr,'location','best')
for i=1:nNfib
    subplot(2,5,i)
    hold on
    set(gca,'fontsize',16)
    xlabel('$S$')
    ylabel('$n_S$')
    xlim([2 inf])
end
figStart = figStart + 1;
print(dataFile{1},'-dpng')
%}

dataFile = cell(nMu,1); 
%%{
% Plot for histogram and compare nfib
for j=1:nMu
    figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6])
    dataFile{j} = ['fig5_cluster_dist/',fileNameArr{1},'_nS_norm_mu',num2str(muArr(j))]; 
    hold on
    for i=1:nNfib
        name = ['../data_stressVSfriction/SDist/',fileNameArr{1},'_Sdist_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'.txt'];
        data = csvread(name);
        data(:,1:4) = [];
        Sbin = size(data,2);
        SbinArr = 2:Sbin+1;
        Shist = mean(data,1);
        scatter(SbinArr,Shist/nfibArr(i),'filled',...
            'MarkerFaceColor',colorArr{i},...
            'MarkerEdgeColor',colorArr{i})
        title(['$\mu = $ ',num2str(muArr(j))])
        
    end
end
%
ind = 1;
for i=figStart:figStart+nMu-1
    figure(i)
    hold on
    xlabel('$S$')
    ylabel('$<n_S>/N_{fib}$')
    set(gca,'yscale','log','fontsize',18')
    legend(thetaNfibLegendArr,'location','best')
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end
figStart = figStart + nMu;
%}

name = ['sub5_',fileNameArr{1}]; 
save(name)