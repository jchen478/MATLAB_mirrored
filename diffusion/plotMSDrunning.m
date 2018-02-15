close all;
clc;
clear;

colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed') rgb('Orange') rgb('Gold')...
    rgb('Lime')...
    rgb('Olive') rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue')...
    rgb('Plum') rgb('Purple') };

%% common parameters
%{
nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200 1293];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
thetaArr = [1 6];
%}

%%{
nfibArr = [160 240 320 640  3200 6400 ];
lboxArr = [300 343.4 378 476.2  814.3 1026 ];
muArr = [0 1 2 3 4 5 10 15 20];
theta = [0];
rpFiber = 75;
naverage = 20;
%}


% Test case
%{
nfibArr = [640];
lboxArr = [344.7];
muArr = [0];
thetaArr = [1 2 3 5 6];
%}


nFig = 1;
nMu = length(muArr);
nNfib = length(nfibArr);

muLegendArr = cell(nMu,1);

for i=1:nMu
    muLegendArr{i} = ['$\mu = $',num2str(muArr(i))]; 
end

%% calculate Dyy
nStrain = 7200;
Dyy = zeros(nStrain,nMu,nNfib);

dataPath = '../data_stressVSfriction/MSD/';

for j=1:nNfib
    for i=1:nMu
        
        name = [dataPath,'theta',num2str(theta),'_MSD_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File=fopen(name,'r');
        data = fscanf(File,'%f',[3 Inf])';
        fclose(File);
        t = data(:,1);
        MSDy = data(:,2);
        MSDz = data(:,3);
        N = length(t);
        if(N > nStrain)
            disp(['Need to allocate more space for strain, N = ',num2str(N)])
        end
        
        % compute running slope
        for k=1:N-naverage+1
            s1 = k; 
            s2 = k+naverage-1; 
            Dyy(k,i,j) = regress(MSDy(s1:s2),t(s1:s2));
        end
        Dyyplot = Dyy(1:N-naverage+1,i,j); 
        tplot = t(1:N-naverage+1); 
        %%{
        figure(j)
        hold on
        plot(tplot,Dyyplot,'color',colorArr{i},'linewidth',2);
        title(['$\theta_{eq}$',num2str(theta),' ',num2str(nfibArr(j)),' ',num2str(muArr(i))]);
        %}
    end
end

for j=1:nNfib
    figure(j)
    hold on
    xlabel('$\gamma$')
    ylabel('$D_{yy}(\gamma)$')
    legend(muLegendArr,'location','bestoutside')
    title(['$N_{fib} = $',num2str(nfibArr(j)), ' Vorticity'])
    set(gca,'fontsize',18)
end

%% calculate Dzz
Dzz = zeros(nStrain,nMu,nNfib);

dataPath = '../data_stressVSfriction/MSD/';

for j=1:nNfib
    for i=1:nMu
        
        name = [dataPath,'theta',num2str(theta),'_MSD_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File=fopen(name,'r');
        data = fscanf(File,'%f',[3 Inf])';
        fclose(File);
        t = data(:,1);
        MSDy = data(:,2);
        MSDz = data(:,3);
        N = length(t);
        if(N > nStrain)
            disp(['Need to allocate more space for strain, N = ',num2str(N)])
        end
        
        % compute running slope
        for k=1:N-naverage+1
            s1 = k; 
            s2 = k+naverage-1; 
            Dzz(k,i,j) = regress(MSDz(s1:s2),t(s1:s2));
        end
        Dzzplot = Dzz(1:N-naverage+1,i,j); 
        tplot = t(1:N-naverage+1); 
        %%{
        figure(j+nNfib)
        hold on
        plot(tplot,Dzzplot,'color',colorArr{i},'linewidth',2);
        title(['$\theta_{eq}$',num2str(theta),' ',num2str(nfibArr(j)),' ',num2str(muArr(i))]);
        %}
    end
end

for j=1:nNfib
    figure(j+nNfib)
    hold on
    xlabel('$\gamma$')
    ylabel('$D_{zz}(\gamma)$')
    legend(muLegendArr,'location','bestoutside')
    title(['$N_{fib} = $',num2str(nfibArr(j)), ' Gradient'])
    set(gca,'fontsize',18)
end