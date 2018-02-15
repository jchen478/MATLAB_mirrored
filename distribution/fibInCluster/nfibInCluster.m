%%{
clc
clear
close all
%}

%% Common parameters
% nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
% muArr = [0 1 2 3 4 5 7 10 15 17 20 23];

nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];

nfibArr = [160 240 320 640 1280 3200 6400];
muArr = [0 1 2 3 4 5 10 15 20];

% nfibArr = [160 240 320 640 1280 3200 6400  ];
% muArr = [0 1 2 3 4 5 10 15 20];
strainArr = 200:200:1500;
figStart = 1;

nMu = length(muArr);
nNfib = length(nfibArr);
nStrain = length(strainArr);

strainLegendArr = cell(nStrain,1);
muLegendArr = cell(nMu,1);
nfibLegendArr = cell(nNfib,1);
for j=1:nMu
    muLegendArr{j} = ['$\mu = $ ',num2str(muArr(j))];
end
for j=1:nStrain
    strainLegendArr{j} = ['$\gamma = $ ',num2str(strainArr(j))];
end
for j=1:nNfib
    nfibLegendArr{j} = ['$N_{fib} = $ ',num2str(nfibArr(j))];
end

%{
% Find number of clusters exceeding threshold
nCluster2 = zeros(nMu,nNfib);
nCluster5 = zeros(nMu,nNfib);
nCluster10 = zeros(nMu,nNfib);
nCluster15 = zeros(nMu,nNfib);

for i=1:nNfib
    for j=1:nMu
        name = ['data_nfibInCluster/theta1_Sdist_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'.txt'];
        data = csvread(name);
        
        % delete data before steady configs
        data(1:20,:) = [];
        ndata = size(data,1);
        
        % 2
        nCluster2(j,i) = sum(data(:,4))/ndata;
        
        % delete non histogram info
        data(:,1:4) = [];
        
        % 5
        if (size(data,2) <= 3)
            continue;
        end
        data(:,1:3) = [];
        nCluster5(j,i) = sum(sum(data))/ndata;
        
        % 10
        if (size(data,2) <= 5)
            continue;
        end
        data(:,1:5) = [];
        nCluster10(j,i) = sum(sum(data))/ndata;
        
        % 15
        if (size(data,2) <= 5)
            continue;
        end
        data(:,1:5) = [];
        nCluster15(j,i) = sum(sum(data))/ndata;
    end
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
title('Number of cluster $(S \ge 2)$')
ylabel('$n_S$')
% set(gca,'yscale','log')
hold on
for i=1:nNfib
    plot(muArr,nCluster2(:,i),'-.o')
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
title('Number of cluster $(S \ge 5)$')
ylabel('$n_S$')
% set(gca,'yscale','log')
hold on
for i=1:nNfib
    plot(muArr,nCluster5(:,i),'-.o')
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
title('Number of cluster $(S \ge 10)$')
ylabel('$n_S$')
% set(gca,'yscale','log')
hold on
for i=1:nNfib
    plot(muArr,nCluster10(:,i),'-.o')
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
title('Number of cluster $(S \ge 15)$')
ylabel('$n_S$')
% set(gca,'yscale','log')
hold on
for i=1:nNfib
    plot(muArr,nCluster15(:,i),'-.o')
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
title('Number of cluster $(S \ge 2)$')
ylabel('$n_S/N_{fib}$')
hold on
for i=1:nNfib
    plot(muArr,nCluster2(:,i)/nfibArr(i),'-.o')
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
title('Number of cluster $(S \ge 5)$')
ylabel('$n_S/N_{fib}$')
hold on
for i=1:nNfib
    plot(muArr,nCluster5(:,i)/nfibArr(i),'-.o')
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
title('Number of cluster $(S \ge 10)$')
ylabel('$n_S/N_{fib}$')
hold on
for i=1:nNfib
    plot(muArr,nCluster10(:,i)/nfibArr(i),'-.o')
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
title('Number of cluster $(S \ge 15)$')
ylabel('$n_S/N_{fib}$')
hold on
for i=1:nNfib
    plot(muArr,nCluster15(:,i)/nfibArr(i),'-.o')
end

for i=figStart:figStart + 7
    figure(i)
    legend(nfibLegendArr,'location','best')
    xlabel('$\mu$')
    
end
figStart = figStart + 8;
%}

%{
% Plot for each Nfib sum of bins
for i=1:nNfib
    figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
    title(['$N_{fib} =$ ', num2str(nfibArr(i))])
    hold on
    for j=1:nMu
        name = ['data_nfibInCluster/theta1_Sdist_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'.txt'];
        data = csvread(name);
        data(:,1:4) = [];
        Sbin = size(data,2);
        SbinArr = 2:Sbin+1;
        Shist = mean(data,1);
        plot(SbinArr,Shist,'-.o')
    end
end

for i=figStart:figStart + nNfib - 1
    figure(i)
    legend(muLegendArr,'location','best')
    xlabel('$S$')
    ylabel('$n_S$')
     xlim([5 inf])
end
figStart = figStart + nNfib;
%}

%%{
% Plot for histogram and compare nfib
figure('units','normalized','outerposition',[0.05 0.05 0.85 0.85])
hold on
for i=1:nNfib
    for j=1:nMu
        subplot(4,3,j)
        hold on
        name = ['data_nfibInCluster/theta1_Sdist_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'.txt'];
        data = csvread(name);
        data(:,1:4) = [];
        Sbin = size(data,2);
        SbinArr = 2:Sbin+1;
        Shist = mean(data,1);
        plot(SbinArr,Shist/nfibArr(i),'-.o')
        set(gca,'yscale','log')
        title(['$\mu = $ ',num2str(muArr(j))])
%         xlim([2 15])
%         ylim([1e-6 1e-1])
        xlabel('$S$')
        ylabel('$n_S / N_{fib}$')
    end
end
%
for i=figStart:figStart
    figure(i)
    subplot(4,3,1)
%     legend(nfibLegendArr,'location','best')
    for j=1:nMu
        subplot(4,3,j)
        xlabel('$S$')
        ylabel('$n_S$')
        xlim([2 80])
    end
end
figStart = figStart + 1;
%}

%{
% Plot for each Nfib sum of bins multiply by bin value
for i=1:nNfib
    figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
    title(['$N_{fib} =$ ', num2str(nfibArr(i))])
    hold on
    for j=1:nMu
        name = ['data_nfibInCluster/theta1_Sdist_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'.txt'];
        data = csvread(name);
        data(:,1:4) = [];
        Sbin = size(data,2);
        SbinArr = 2:Sbin+1;
        ShistNorm = sum(data(26:end,:),1).*SbinArr;
        plot(SbinArr,ShistNorm,'-.o')
    end
end

for i=figStart:figStart + nNfib - 1
    figure(i)
    legend(muLegendArr,'location','best')
    xlabel('$S$')
    ylabel('$n_S \times S$')
    xlim([5 inf])
end
figStart = figStart + nNfib;
%}

%%{
% Plot of the number of fibers in a contact
NfibContact = zeros(nMu,nNfib);

for i=1:nNfib
    for j=1:nMu
        name = ['data_nfibInCluster/helical_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'_nfibInCluster.txt'];
        File = fopen(name,'r');
        data = fscanf(File,'%d',[2 Inf])';
        fclose(File);
        NfibContact(j,i) = mean(data(:,2));
    end
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
ylabel('$N_{fib}$ in contact')
title('Number of fibers in a contact')
hold on
for i=1:nNfib
    plot(muArr,NfibContact(:,i),'-.o')
end

figure('units','normalized','outerposition',[0.2 0.2 0.5 0.6])
ylabel('Normalized $N_{fib}$ in contact')
title('Normalized number of fibers in a contact')
hold on
for i=1:nNfib
    plot(muArr,NfibContact(:,i)/nfibArr(i),'-.o')
end

for i=figStart:figStart + 1
    figure(i)
    legend(nfibLegendArr,'location','bestoutside')
    xlabel('$\mu$')
    xlim([0 inf])
end
figStart = figStart + 2;
%}

%{
% Plot of distribution evolution
j = 11;
figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9])
for i=1:nNfib
    subplot(3,3,i)
    hold on
    xlim([5 inf])
    xlabel('$S$')
    ylabel('$n_S/N_{fib}$')
    set(gca,'fontsize',18)
    title(['$N_{fib} =$ ', num2str(nfibArr(i)),' $\mu =$ ',num2str(muArr(j))],'fontsize',18)
    name = ['data_nfibInCluster/theta1_Sdist_nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'.txt'];
    data = csvread(name);
    n = 1;
    k = 1;
    while (n <= nStrain)
        if (data(k,3) == strainArr(n))
            n = n + 1;
            k = k + 1;
        else
            data(k,:) = [];
        end
    end
    data(nStrain+1:end,:) = [];
    data(:,1:4) = [];
    Sbin = size(data,2);
    SbinArr = 2:Sbin+1;
    Shist = data;
    for k=1:nStrain
        plot(SbinArr,Shist(k,:)/nfibArr(i),'-.o')
        xlim([2 20])
        ylim([0 0.06])
    end
    
end


for i=figStart:figStart
    figure(i)
    subplot(3,3,1)
    legend(strainLegendArr,'location','best','fontsize',14)
end
figStart = figStart + 1;
%}