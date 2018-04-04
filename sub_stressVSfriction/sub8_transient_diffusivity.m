%%%
%%% Plot transient diffusivities and msd
%%%

clc;
clear;
close all;

% Define common parameters
simulation_cases;

% data path
dataPath = '../data_stressVSfriction/MSD/';

% array allocation
start_msd = 1600;
nStrain = 4000;
window = 2; 
strain = zeros(nStrain,nMu,nNfib);
Dyy = zeros(nStrain,nMu,nNfib);
Dzz = zeros(nStrain,nMu,nNfib);

% figure info
figStart = 1;

%%  Plot MSDs
for j=1:nNfib
    figure('Units','Inches','Position',[3 3 6.5 6.5])
    set(gca,'XMinorTick','on')
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
        xlim([start_msd inf])
        xlabel('$\gamma$')
        legend(muLegendArr,'location','bestoutside')
    end
end

figStart = figStart + nNfib; 

%%  Calculate moving regression
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
        strain(1:N,i,j) = t;
    end
end

%%  Plot transient diffusivity
for j=1:nNfib
    figure('Units','Inches','Position',[3 3 6.5 6.5])
    set(gca,'XMinorTick','on')
    hold on
    for i=1:nMu
        subplot(2,1,1)
        title(['$N_{fib} = $ ',num2str(nfibArr(j))])
        hold on
        plot(strain(:,i,j), Dyy(:,i,j))
        ylabel('Transiet D vorticity')
        
        subplot(2,1,2)
        hold on
        plot(strain(:,i,j), Dzz(:,i,j))
        ylabel('Transiet D gradient')
    end
end

for i=figStart:figStart+nNfib-1
    figure(i)
    for j=1:2
        subplot(2,1,j)
        hold on
        box on
        xlim([start_msd inf])
        xlabel('$\gamma$')
        legend(muLegendArr,'location','bestoutside')
    end
end

figStart = figStart + nNfib; 
