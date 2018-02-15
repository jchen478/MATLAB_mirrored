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
nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200 1293];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
thetaArr = [1];
rpFiber = 75;
naverage = 1000;
%}

%{
% Straight fiber parameters
nfibArr = [160 240 320 640 1280 3200 6400  ];
lboxArr = [300 343.4 378 476.2 600 814.3 1026  ];
muArr = [0 1 2 3 4 5 10 15 20];
thetaArr = [0];
rpFiber = 75;
naverage = 1000;
%}

% Test case
%{
nfibArr = [640];
lboxArr = [344.7];
muArr = [0];
thetaArr = [1 2 3 5 6];
%}


nFig = 1;
nTheta = length(thetaArr);
nMu = length(muArr);
nNfib = length(nfibArr);

thetaNfibLegendArr = cell(nTheta*nNfib,1);
muLegendArr = cell(nMu,1);

for i=1:nMu
    muLegendArr{i} = ['$\mu =$ ',num2str(muArr(i))];
end
for i=1:nTheta
    for j=1:nNfib
        thetaNfibLegendArr{(i-1)*nNfib+j} = ['$(\theta_{eq},N_{fib}) =$ (0.',num2str(thetaArr(i)),', ',num2str(nfibArr(j)),')'];
    end
end

%% calculate Diffusivitys
Dyy = zeros(nMu,nNfib,nTheta);
Dzz = zeros(nMu,nNfib,nTheta);
se_Dyy = zeros(nMu,nNfib,nTheta);
se_Dzz = zeros(nMu,nNfib,nTheta);

dataPath = '../data_stressVSfriction/MSD/';
% dataPath = 'rp20/';

for k=1:nTheta
    for j=1:nNfib
        for i=1:nMu
            name = [dataPath,'theta',num2str(thetaArr(k)),'_MSD_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File=fopen(name,'r');
            data = fscanf(File,'%f',[3 Inf])';
            fclose(File);
            t = data(:,1);
            MSDy = data(:,2);
            MSDz = data(:,3);
            %{
            figure()
            plot(t,MSDy,t,MSDz);
            legend('y','z','location','best')
            title(['$\theta_{eq}$',num2str(thetaArr(k)),' ',num2str(nfibArr(j)),' ',num2str(muArr(i))]);
            %}
            N = length(t);
            n = N-naverage;
            t(1:n) = [];
            MSDy(1:n) = [];
            MSDz(1:n) = [];
            [Dyy(i,j,k), Dyyconf] = regress(MSDy,t);
            [Dzz(i,j,k), Dzzconf] = regress(MSDz,t);
            se_Dyy(i,j,k) = Dyyconf(2)-Dyyconf(1);
            se_Dzz(i,j,k) = Dzzconf(2)-Dzzconf(1);
        end
    end
end



%% plotting diffusion with mu as x-axis

dataFile = cell(1,1);
for i=1:1
    dataFile{i} = '../stressVSfriction/dynamic_figures/theta';
    for k=1:nTheta
        dataFile{i} = [dataFile{i},num2str(thetaArr(k)),'_'];
    end
    dataFile{i} = [dataFile{i},'nfib'];
    for j=1:nNfib
        dataFile{i} = [dataFile{i},num2str(nfibArr(j)),'_'];
    end
end

dataFile{1} = [dataFile{1},'D_vs_mu'];
figure('units','normalized','outerposition',[0.1 0.1 0.7 0.8])
subplot(2,1,1)
hold on
title('\textbf{Vorticity}','fontsize',18)
ylabel('$D_{yy}/b^2\dot{\gamma}$')
for k=1:nTheta
    for j=1:nNfib
        errorbar(muArr,Dyy(:,j,k),se_Dyy(:,j,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

subplot(2,1,2)
hold on
title('\textbf{Gradient}','fontsize',18)
ylabel('$D_{zz}/b^2\dot{\gamma}$')
for k=1:nTheta
    for j=1:nNfib
        errorbar(muArr,Dzz(:,j,k),se_Dzz(:,j,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure(nFig)
for j=1:2
    subplot(2,1,j)
    hold on
    box on
    xlabel('$\mu$')
    xlim([0 inf])
    legend(thetaNfibLegendArr,'location','bestoutside','fontsize',18)
    set(gca,'fontsize',18)
end
%     savefig([dataFile{nFig},'.fig']);
print(dataFile{nFig},'-dpng')
nFig = nFig + 1;

%% plotting diffusion with Nfib as x-axis
%{
dataFile = cell(nTheta,1);
for i=1:nTheta
    dataFile{i} = ['../stressVSfriction/dynamic_figures/theta',num2str(thetaArr(i)),'_nfib'];
    for j=1:nNfib
        dataFile{i} = [dataFile{i},num2str(nfibArr(j)),'_'];
    end
    dataFile{i} = [dataFile{i},'D_vs_nfib'];
end

for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.02 0.7 0.98])
    subplot(2,1,1)
    hold on
    title(['Vorticity, $\theta_{eq} = 0.$',num2str(thetaArr(k))],'fontsize',18)
    ylabel('$D_{yy}/b^2\dot{\gamma}$')
    for j=1:nMu
        errorbar(nfibArr,Dyy(j,:,k),se_Dyy(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'color',colorArr{j},...
            'MarkerEdgeColor',colorArr{j})
    end
end

for k=1:nTheta
    figure(nFig+k-1)
    subplot(2,1,2)
    hold on
    title(['Gradient, $\theta_{eq} = 0.$',num2str(thetaArr(k))],'fontsize',18)
    ylabel('$D_{zz}/b^2\dot{\gamma}$')
    for j=1:nMu
        errorbar(nfibArr,Dzz(j,:,k),se_Dzz(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'color',colorArr{j},...
            'MarkerEdgeColor',colorArr{j})
    end
end

ind = 1;
for i=nFig:nFig+nTheta-1
    figure(i)
    for j=1:2
        subplot(2,1,j)
        hold on
        box on
        xlabel('$N_{fib}$')
        xlim([0 inf])
        legend(muLegendArr,'location','bestoutside','fontsize',18)
    end
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end
nFig = nFig + nTheta;
%}

%% plotting diffusion with L/Lbox as x-axis
%{
dataFile = cell(nTheta,1);
for i=1:nTheta
    dataFile{i} = ['../stressVSfriction/dynamic_figures/theta',num2str(thetaArr(i)),'_nfib'];
    for j=1:nNfib
        dataFile{i} = [dataFile{i},num2str(nfibArr(j)),'_'];
    end
    dataFile{i} = [dataFile{i},'D_vs_overLbox'];
end

for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.02 0.7 0.98])
    subplot(2,1,1)
    hold on
    title(['Vorticity, $\theta_{eq} = 0.$',num2str(thetaArr(k))],'fontsize',18)
    ylabel('$D_{yy}/b^2\dot{\gamma}$')
    for j=1:nMu
        errorbar(2*rpFiber./lboxArr,Dyy(j,:,k),se_Dyy(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'color',colorArr{j},...
            'MarkerEdgeColor',colorArr{j})
    end
end

for k=1:nTheta
    figure(nFig+k-1)
    subplot(2,1,2)
    hold on
    title(['Gradient, $\theta_{eq} = 0.$',num2str(thetaArr(k))],'fontsize',18)
    ylabel('$D_{zz}/b^2\dot{\gamma}$')
    for j=1:nMu
        errorbar(2*rpFiber./lboxArr,Dzz(j,:,k),se_Dzz(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'color',colorArr{j},...
            'MarkerEdgeColor',colorArr{j})
    end
end

ind = 1;
for i=nFig:nFig+nTheta-1
    figure(i)
    for j=1:2
        subplot(2,1,j)
        hold on
        box on
        xlabel('$L/L_{box}$')
        xlim([0 inf])
        legend(muLegendArr,'location','bestoutside','fontsize',18)
    end
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end
nFig = nFig + nTheta;
%}

%% plotting diffusivity dependence on overLbox with mu as x-axis
% Dyy_slope = zeros(nMu,1,nTheta);
% Dzz_slope = zeros(nMu,1,nTheta);
% se_Dyy_slope = zeros(nMu,1,nTheta);
% se_Dzz_slope = zeros(nMu,1,nTheta);
% % find the slope
% for k=1:nTheta
%     for i=1:nMu
%         [Dyy_slope(i,1,k),Dyyconf] = regress(Dyy(i,:,k)',2*rpFiber./lboxArr');
%         [Dzz_slope(i,1,k),Dzzconf] = regress(Dzz(i,:,k)',2*rpFiber./lboxArr');
%         se_Dyy_slope(i,1,k) = Dyyconf(2)-Dyyconf(1);
%         se_Dzz_slope(i,1,k) = Dzzconf(2)-Dzzconf(1);
%     end
%     figure('units','normalized','outerposition',[0.1 0.1 0.6 0.6])
%     hold on
%     title(['$\theta_{eq} = 0.$',num2str(thetaArr(k))]);
%     errorbar(muArr,Dyy_slope(:,1,k),se_Dyy_slope(:,1,k)/2,...
%         '-.o','MarkerSize',10, 'linewidth',2.5);
%     errorbar(muArr,Dzz_slope(:,1,k),se_Dzz_slope(:,1,k)/2,...
%         '-.o','MarkerSize',10, 'linewidth',2.5);
% end
% 
% for i=nFig:nFig+nTheta-1
%     figure(i)
%     box on
%     legend('$\dot{D_{yy}}/b^2\dot{\gamma}$','$\dot{D_{zz}}/b^2\dot{\gamma}$','location','best')
%     xlabel('$\mu$')
%     ylabel('$\dot{D}/b^2\dot{\gamma}$')
% end
% nFig = nFig + nTheta;


%% plotting diffusivity comparison between two directions with mu as x-axis
%{
for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
        hold on
        plot(muArr,Dyy(:,j,k),'-.o');
        plot(muArr,Dzz(:,j,k),'-.o');
        title(['$(\theta_{eq},N_{fib}) = (0.$',num2str(thetaArr(k)),', ',num2str(nfibArr(j)),'$)$'])
    end
end
for i=nFig:nFig+nTheta*nNfib-1
    figure(i)
    box on
    legend('$D_{yy}$','$D_{zz}$','location','best')
    xlabel('$\mu$')
    ylabel('$D$')
end
nFig = nFig + nTheta*nNfib;
%}

%% plotting diffusivity comparison between two directions with theta as x-axis
%{
for j=1:nNfib
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    hold on
    Dy = reshape(Dyy(:,j,:),[nTheta,1]);
    Dz = reshape(Dzz(:,j,:),[nTheta,1]);
    scatter(thetaArr/10,Dy,40,'filled','markerFaceColor','black','markerEdgeColor','black','linewidth',2.0);
    scatter(thetaArr/10,Dz,40,'markerEdgeColor','black','linewidth',2.0);
    legend('$D_{yy}/b^2\dot{\gamma}$','$D_{zz}/b^2\dot{\gamma}$','location','best')
    xlabel('$\theta_{eq}$')
    ylabel('$D/b^2\dot{\gamma}$')
%     title(['$N_{fib} =$ ',num2str(nfibArr(j))])
    ylim([0 0.05])
    xlim([0 0.7])
end
%}