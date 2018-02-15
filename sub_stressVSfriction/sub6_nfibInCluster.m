%%%
%%% Plot total number of fibers with at least one contact
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

figStart = 1;
nMu = length(muArr);
nNfib = length(nfibArr);
nfibLegendArr = cell(nNfib,1);

thetaNfibLegendArr = cell(nNfib,1);
if strcmpi(fileNameArr,'helical')
    for j=1:nNfib
        thetaNfibLegendArr{j} = ['$(\theta_{eq},\phi_{eq},N_{fib}) =$ (0.8, 0.7, ',num2str(nfibArr(j)),')'];
    end
else
    for j=1:nNfib
        thetaNfibLegendArr{j} = ['$(\theta_{eq},N_{fib}) =$ (0.',num2str(thetaArr(1)),', ',num2str(nfibArr(j)),')'];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read files and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot of the number of fibers in a contact
NfibContact = zeros(nMu,nNfib);

for i=1:nNfib
    for j=1:nMu
        name = ['../data_stressVSfriction/nfibInCluster/',fileNameArr{1},'_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'_nfibInCluster.txt'];
        File = fopen(name,'r');
        data = fscanf(File,'%d',[2 Inf])';
        fclose(File);
        NfibContact(j,i) = mean(data(:,2));
    end
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
ylabel('$N_{fib}$ in contact')
hold on
for i=1:nNfib
    plot(muArr,NfibContact(:,i),'-.o','MarkerSize',10,'Linewidth',2.5)
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
ylabel('Normalized $N_{fib}$ in contact')
hold on
for i=1:nNfib
    plot(muArr,NfibContact(:,i)/nfibArr(i),'-.o','MarkerSize',10,'Linewidth',2.5)
end

dataFile = cell(2,1);
for ii=1:2
    dataFile{ii} = ['fig6_nfibInCluster/',fileNameArr{1},'_nfib'];
    for j=1:nNfib
        dataFile{ii} = [ dataFile{ii},num2str(nfibArr(j)),'_'];
    end
end
dataFile{1} = [ dataFile{1},'nfibInCluster'];
dataFile{2} = [ dataFile{2},'norm_nfibInCluster'];

for i=figStart:figStart + 1
    figure(i)
    legend(thetaNfibLegendArr,'location','bestoutside')
    xlabel('$\mu$')
    xlim([0 inf])
    print(dataFile{i},'-dpng')
end

name = ['sub6_',fileNameArr{1}]; 
save(name)