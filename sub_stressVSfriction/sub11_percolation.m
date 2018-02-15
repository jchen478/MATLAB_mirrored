%%%
%%% Plot transient percolation results
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

%%{
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
nfibArr = [160 240 320 640 1280 3200 ];
lboxArr = [300 343.4 378 476.2 600 814.3];
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
% nfibArr = [160 240 320 640 1280 3200];
% lboxArr = [300 343.4 378 476.2 600 814.3];
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

dataFile = cell(2,1);
for ii=1:2
    dataFile{ii} = ['fig11_percolation/',fileNameArr{1},'_nfib'];
    for j=1:nNfib
        dataFile{ii} = [ dataFile{ii},num2str(nfibArr(j)),'_'];
    end
end
dataFile{1} = [ dataFile{1},'nP_strain'];
dataFile{2} = [ dataFile{2},'nP_averaged'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot Percolation vs. strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perc_averaged = zeros(nMu, nNfib); 
dataPath = '../data_stressVSfriction/Percolation/';

for j=1:nNfib
    figure('units','normalized','outerposition',[0.2 0.2 0.5 0.8])
    hold on
    title(fileNameArr{1});
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_perc_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File=fopen(name,'r');
        data = fscanf(File,'%f',[2 Inf])';
        fclose(File);
        t = data(:,1);
        perc = data(:,2);
        perc_averaged(i,j) = mean(perc);
        scatter(t,perc,'filled')
    end
end

box on
xlabel('$\gamma$')
ylabel('Number of percolating structures')
legend(muLegendArr,'location','bestoutside')
set(gca,'fontsize',16)
print(dataFile{1},'-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot Percolation averaged over strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
ylabel('$<N_{p}>$ (over frame)'); 
hold on
for i=1:nNfib
    plot(muArr,perc_averaged(:,i),'-.o','MarkerSize',10,'Linewidth',2.5)
end

legend(thetaNfibLegendArr,'location','bestoutside')
xlabel('$\mu$')
xlim([0 inf])
print(dataFile{2},'-dpng')
