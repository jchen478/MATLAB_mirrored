%%%
%%% Plot transient diffusivities and msd
%%%

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
%%%%%%%%%%%%%%%%%%%%% U-shaped fibers %%%%%%%%%%%%%%%%%%%%%
% fileNameArr = {'theta0'}; thetaArr = 0;
% fileNameArr = {'theta1'}; thetaArr = 1;
fileNameArr = {'theta3'}; thetaArr = 3;
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

%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot MSDs
%%%%%%%%%%%%%%%%%%%%%%%%

dataPath = '../data_stressVSfriction/MSD/';

for j=1:nNfib
    figure('units','normalized','outerposition',[0.2 0.2 0.5 0.8])
    hold on
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_MSD_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File=fopen(name,'r');
        data = fscanf(File,'%f',[3 Inf])';
        fclose(File);
        t = data(:,1);
        MSDy = data(:,2);
        MSDz = data(:,3);
        
        subplot(2,1,1)
        title(['$N_{fib} = $ ',num2str(nfibArr(j))])
        hold on
        plot(t,MSDy)
        ylabel('vorticity MSD')
        
        subplot(2,1,2)
        hold on
        plot(t,MSDz)
        ylabel('gradient MSD')
    end
end

for i=figStart:figStart+nNfib-1
    figure(i)
    for j=1:2
        subplot(2,1,j)
        hold on
        box on
        xlim([10 inf])
        xlabel('$\gamma$')
        legend(muLegendArr,'location','bestoutside')
        set(gca,'fontsize',16)
    end
end

figStart = figStart + nNfib; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate moving regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nStrain = 3000;
window = 2; 
Dyy = zeros(nStrain,nMu,nNfib);
Dzz = zeros(nStrain,nMu,nNfib);

for j=1:nNfib   
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_MSD_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File=fopen(name,'r');
        data = fscanf(File,'%f',[3 Inf])';
        fclose(File);
        t = data(:,1);
        MSDy = data(:,2);
        MSDz = data(:,3);
        N = length(MSDy); 
        for k=1:N-window+1
            s1 = k;
            s2 = k+window-1;
            Dyy(k,i,j) = regress(MSDy(s1:s2),t(s1:s2));
            Dzz(k,i,j) = regress(MSDz(s1:s2),t(s1:s2));
        end
       
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot transient diffusivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strain=0:nStrain-1; 
strain = strain*0.5; 
for j=1:nNfib
    figure('units','normalized','outerposition',[0.2 0.2 0.5 0.8])
    hold on
    for i=1:nMu
        subplot(2,1,1)
        title(['$N_{fib} = $ ',num2str(nfibArr(j))])
        hold on
        plot(strain, Dyy(:,i,j))
        ylabel('Transiet D vorticity')
        
        subplot(2,1,2)
        hold on
        plot(strain, Dzz(:,i,j))
        ylabel('Transiet D gradient')
    end
end

for i=figStart:figStart+nNfib-1
    figure(i)
    for j=1:2
        subplot(2,1,j)
        hold on
        box on
        xlim([10 inf])
        xlabel('$\gamma$')
        legend(muLegendArr,'location','bestoutside')
        set(gca,'fontsize',16)
    end
end

figStart = figStart + nNfib; 

