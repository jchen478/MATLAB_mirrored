%%%
%%% Plot transient diffusivities and msd
%%%

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%%%%%%%%%%%%%%%%%%%%% U-shaped fibers %%%%%%%%%%%%%%%%%%%%%
% fileNameArr = {'theta0'}; thetaArr = 0;
% fileNameArr = {'theta1'}; thetaArr = 1;
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
thetaArr = [3];
fileNameArr = {'helical'};
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate diffusivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dyy = zeros(nMu, nNfib); 
Dzz = zeros(nMu, nNfib); 
conf_Dyy = zeros(nMu, nNfib); 
conf_Dzz = zeros(nMu, nNfib); 

strain_eliminat = 50; 

dataPath = '../data_stressVSfriction/MSD/';

for j=1:nNfib
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_MSD_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File=fopen(name,'r');
        data = fscanf(File,'%f',[3 Inf])';
        fclose(File);
        t = data(:,1);
        MSDy = data(:,2);
        MSDz = data(:,3);
        N = length(t); 
        if ( N < 5 ) 
            continue; 
        end
%         t(1:strain_eliminat*2) = [];
%         MSDy(1:strain_eliminat*2) = [];
%         MSDz(1:strain_eliminat*2) = [];         
        [Dyy(i,j), Dyyconf] = regress(MSDy,t);
        [Dzz(i,j), Dzzconf] = regress(MSDz,t);
        conf_Dyy(i,j) = Dyyconf(2)-Dyyconf(1);
        conf_Dzz(i,j) = Dzzconf(2)-Dzzconf(1);       
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  File storage name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataFile = cell(2,1);
for ii=1:2
    dataFile{ii} = ['fig7_diffusivity/',fileNameArr{1},'_nfib'];
    for j=1:nNfib
        dataFile{ii} = [ dataFile{ii},num2str(nfibArr(j)),'_'];
    end
end
dataFile{1} = [ dataFile{1},'Dyy_vs_mu'];
dataFile{2} = [ dataFile{2},'Dzz_vs_mu'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot vorticity diffusivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','Inches','Position',[3 3 3.5 6.5])
subplot(2,1,1)
hold on
for j=1:nNfib
    errorbar(muArr,Dyy(:,j),conf_Dyy(:,j)/2,...
        '-.o','MarkerSize',5, 'linewidth',2.5,...
        'MarkerFaceColor',colorArr{j}, ...
        'MarkerEdgeColor',colorArr{j},...
        'color',colorArr{j});
end

box on
xlim([0 inf])
set(gca,'XMinorTick','on')
xlabel('$\mu$')
ylabel('Vorticity $D_{yy}$')
legend(thetaNfibLegendArr,'location','bestoutside')
print(dataFile{1},'-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot gradient diffusivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
hold on
for j=1:nNfib   
    errorbar(muArr,Dzz(:,j),conf_Dzz(:,j)/2,...
        '-.o','MarkerSize',5, 'linewidth',2.5,...
        'MarkerFaceColor',colorArr{j}, ...
        'MarkerEdgeColor',colorArr{j},...
        'color',colorArr{j});
end

box on
xlim([0 inf])
set(gca,'XMinorTick','on')
xlabel('$\mu$')
ylabel('Gradient $D_{zz}$')
legend(thetaNfibLegendArr,'location','bestoutside')
print(dataFile{2},'-dpng')

name = ['sub7_',fileNameArr{1}]; 
save(name)
