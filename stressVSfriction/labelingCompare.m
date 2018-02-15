clc;
clear;
close all;

nfibArr = [160 1280 12800];
nNfib = length(nfibArr);

data_160 = load('160_meanClusterSize.mat');
data_1280 = load('1280_meanClusterSize.mat');
data_12800 = load('12800_meanClusterSize.mat');

mu160 = data_160.CountAverage(:,1);
mu1280 = data_1280.CountAverage(:,1);
mu12800 = data_12800.CountAverage(:,1);
S160 = data_160.CountAverage(:,2);
S1280 = data_1280.CountAverage(:,2);
S12800 = data_12800.CountAverage(:,2);
se_S160 = data_160.CountAverage(:,3);
se_S1280 = data_1280.CountAverage(:,3);
se_S12800 = data_12800.CountAverage(:,3);

nMu = length(mu160);

color = {rgb('DarkGreen'),rgb('DarkOrange'),rgb('MediumBlue')};

%% S vs mu
figure(1)
hold on;
box on;
errorbar(mu160,S160,se_S160/2,'-.o',...
    'MarkerFaceColor',color{1},'Linewidth',2.5,...
    'color',color{1});
errorbar(mu1280,S1280,se_S1280/2,'-.o',...
    'MarkerFaceColor',color{2},'Linewidth',2.5,...
    'color',color{2});
errorbar(mu12800,S12800,se_S1280/2,'-.o',...
    'MarkerFaceColor',color{3},...
    'Linewidth',2.5,'color',color{3});

xlabel('$\mu$')
ylabel('$\bar{S}$')

legendArr = cell(nNfib,1);
for i=1:nNfib
    legendArr{i} = (['$N_{fib} =$ ',num2str(nfibArr(i))]);
end
legend(legendArr{:},'location','best')

xlim([0 25])
ylim([0 inf])

%% S vs nfib
plotMarker=['-o' '-s' '-*' '-o' '-s' '-*' '-o' '-s' '-*' '-o' '-s' '-*' ];
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed') rgb('DarkOrange') rgb('Gold') rgb('Lime') rgb('Olive') rgb('DarkGreen') rgb('SkyBlue') rgb('MediumBlue') rgb('Orchid') rgb('Purple')};

figure(2)
hold on
box on
Snfib = zeros(nNfib,1);
se_Snfib = zeros(nNfib,1);
for i=1:nMu    
    Snfib(1) = S160(i);
    Snfib(2) = S1280(i);
    Snfib(3) = S12800(i);
    se_Snfib(1) = se_S160(i);
    se_Snfib(2) = se_S1280(i);
    se_Snfib(3) = se_S12800(i);
    plot(nfibArr,Snfib,'-o','color',colorArr{i},'Linewidth',2.5);
%     errorbar(nfibArr,Snfib,se_Snfib/2,scatterMarker(i),'Linewidth',2.5);
end
legendArr = cell(nMu,1);
for i=1:nMu
    legendArr{i} = (['$\mu =$ ',num2str(mu160(i))]);
end
legend(legendArr{:},'location','best')
% set(gca,'xscale','log','yscale','log')
xlabel('$N_{fib}$')
ylabel('$\bar{S}$')

figure(3)
hold on
box on
Snfib = zeros(nNfib,1);
se_Snfib = zeros(nNfib,1);
LboxArr = [300 600 1293]; 
for i=1:nMu    
    Snfib(1) = S160(i);
    Snfib(2) = S1280(i);
    Snfib(3) = S12800(i);
    se_Snfib(1) = se_S160(i);
    se_Snfib(2) = se_S1280(i);
    se_Snfib(3) = se_S12800(i);
    plot(LboxArr,Snfib,'-o','color',colorArr{i},'Linewidth',2.5);
%     errorbar(nfibArr,Snfib,se_Snfib/2,scatterMarker(i),'Linewidth',2.5);
end
legendArr = cell(nMu,1);
for i=1:nMu
    legendArr{i} = (['$\mu =$ ',num2str(mu160(i))]);
end
legend(legendArr{:},'location','best')
% set(gca,'xscale','log','yscale','log')
xlabel('$L_{box}$')
ylabel('$\bar{S}$')

figure(4)
hold on
box on
Snfib = zeros(nNfib,1);
se_Snfib = zeros(nNfib,1);
LboxArr = [300 600 1293]; 
for i=1:nMu    
    Snfib(1) = S160(i);
    Snfib(2) = S1280(i);
    Snfib(3) = S12800(i);
    se_Snfib(1) = se_S160(i);
    se_Snfib(2) = se_S1280(i);
    se_Snfib(3) = se_S12800(i);
    plot(150./LboxArr,Snfib,'-o','color',colorArr{i},'Linewidth',2.5);
%     errorbar(nfibArr,Snfib,se_Snfib/2,scatterMarker(i),'Linewidth',2.5);
end
legendArr = cell(nMu,1);
for i=1:nMu
    legendArr{i} = (['$\mu =$ ',num2str(mu160(i))]);
end
legend(legendArr{:},'location','bestoutside')
% set(gca,'xscale','log','yscale','log')
xlabel('$L/L_{box}$')
ylabel('$\bar{S}$')
title('\bf{Mean cluster count}')
ylim([0 inf])

figure(5)
hold on
box on
Snfib = zeros(nNfib,1);
se_Snfib = zeros(nNfib,1);
for i=1:nMu    
    Snfib(1) = S160(i);
    Snfib(2) = S1280(i);
    Snfib(3) = S12800(i);
    se_Snfib(1) = se_S160(i);
    se_Snfib(2) = se_S1280(i);
    se_Snfib(3) = se_S12800(i);
    plot(1./nfibArr,Snfib,'-o','color',colorArr{i},'Linewidth',2.5);
%     errorbar(nfibArr,Snfib,se_Snfib/2,scatterMarker(i),'Linewidth',2.5);
end
legendArr = cell(nMu,1);
for i=1:nMu
    legendArr{i} = (['$\mu =$ ',num2str(mu160(i))]);
end
legend(legendArr{:},'location','best')
% set(gca,'xscale','log','yscale','log')
xlabel('$1/N_{fib}$')
ylabel('$\bar{S}$')