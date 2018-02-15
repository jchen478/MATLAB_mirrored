%%%
%%% Plot cluster size related properties
%%% -- properties include maxDim
%%%

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%%%%%%%%%%%%%%%%%%%% Theta1 fibers %%%%%%%%%%%%%%%%%%%%% 
nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200 1293];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
thetaArr = [1];
fileNameArr = {'theta1'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('Olive') rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue') rgb('Plum') rgb('Purple') };
scatterMarker=['o' 'v' '>' '<' '^' 's' 'd' 'v' '>' '<' '^' 's' 'd' 'o' ];
%}

%{
%%%%%%%%%%%%%%%%%%%%% Theta3 fibers %%%%%%%%%%%%%%%%%%%%% 
nfibArr = [160 240 320 640 1280 3200 6400 10240];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
thetaArr = [3];
fileNameArr = {'theta3'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('Olive') rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue') rgb('Plum') rgb('Purple') };
scatterMarker=['o' 'v' '>' '<' '^' 's' 'd' 'v' '>' '<' '^' 's' 'd' 'o' ];
%}

%%{
%%%%%%%%%%%%%%%%%%%%% Theta6 fibers %%%%%%%%%%%%%%%%%%%%% 
nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200 1293];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
thetaArr = [6];
fileNameArr = {'theta6'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('Olive') rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue') rgb('Plum') rgb('Purple') };
scatterMarker=['o' 'v' '>' '<' '^' 's' 'd' 'v' '>' '<' '^' 's' 'd' 'o' ];
%}

%{
%%%%%%%%%%%%%%%%%%%%% Helical fibers %%%%%%%%%%%%%%%%%%%%% 
nfibArr = [160 240  320 640 1280 3200 6400];
lboxArr = [300 343.4 378 476.2 600 814.3  1026];
muArr = [0 1 2 3 4 5 10 15 20];
thetaArr = [0];
fileNameArr = {'helical'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('DarkGreen') rgb('LightSkyBlue') rgb('Plum')};
scatterMarker=['o' 'v' '>' '<' '^' 's' 'd' 'v' '>' '<' '^' 's' 'd' 'o' ];
%}

%{
%%%%%%%%%%%%%%%%%%%%% Straight fibers %%%%%%%%%%%%%%%%%%%%% 
nfibArr = [160 240 320 640 1280 3200 6400  ];
lboxArr = [300 343.4 378 476.2 600 814.3 1026  ];
muArr = [0 1 2 3 4 5 10 15 20];
thetaArr = [0];
fileNameArr = {'theta0'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('DarkGreen') rgb('LightSkyBlue') rgb('Plum')};
scatterMarker=['o' 'v' '>' '<' '^' 's' 'd' 'v' '>' '<' '^' 's' 'd' 'o' ];
%}

%{
% Straight fiber parameters
nfibArr = [160 240 320 640 1280 3200 6400  ];
lboxArr = [300 343.4 378 476.2 600 814.3 1026  ];
muArr = [0 1 2 3 4 5 10 15 20];
thetaArr = [0];
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('DarkGreen') rgb('LightSkyBlue') rgb('Plum')};
scatterMarker=['o' 'v' '>' '<' '^' 's' 'd' 'v' '>' '<' '^' 's' 'd' 'o' ];
%}

% fiber dimensions
a = 1.60E-05;       % radius (m)
Imom = pi*a^4/4;    % area moment (m^4)
EY = 1e9;           % Young's modulus (Pa m^2)
eta0 = 1;           % fluid viscosity (Pa s)

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
%%% maxDim evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
dataPath = '../data_stressVSfriction/maxDim/';
dataFile = cell(nTheta*nNfib,1);

for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.2 0.2 0.6 0.9])
        hold on;
        dataFile{(k-1)*nNfib+j} = ['fig4_cluster_dim/',fileNameArr{1},'_nfib',num2str(nfibArr(j)),'_maxDim_vs_strain'];
        maxStrain = 0;
        subplot(3,1,1)
        hold on;
        for i=1:nMu
            % read files
            name = [dataPath,fileNameArr{1},'_maxDim_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File = fopen(name,'r');
            data = fscanf(File,'%f',[4 Inf])';
            fclose(File);            
            
            % plot maxDim evolution with strain
            subplot(3,1,1)
            hold on
            ylabel('Max $x$ size')
            scatter(data(:,1),data(:,2),markersize,scatterMarker(i),...
                'MarkerFaceColor',colorArr{i},...
                'MarkerEdgeColor',colorArr{i})
            subplot(3,1,2)
            hold on
            ylabel('Max $y$ size')
            scatter(data(:,1),data(:,3),markersize,scatterMarker(i),...
                'MarkerFaceColor',colorArr{i},...
                'MarkerEdgeColor',colorArr{i})
            subplot(3,1,3)
            hold on
            ylabel('Max $z$ size')
            scatter(data(:,1),data(:,4),markersize,scatterMarker(i),...
                'MarkerFaceColor',colorArr{i},...
                'MarkerEdgeColor',colorArr{i})
            if(data(end,1) > maxStrain)
                maxStrain = data(end,1);
            end
        end
        for ii=1:3
            subplot(3,1,ii)
            hold on
            line([0 maxStrain],[lboxArr(j) lboxArr(j)],'linewidth',2,'color',rgb('OrangeRed'))
        end
    end
end

ind = 1;
for i=figStart:figStart+nNfib*nTheta-1
    figure(i)
    subplot(3,1,1)
    set(gca,'XTick',[]);
    subplot(3,1,2)
    set(gca,'XTick',[]);
    subplot(3,1,3)
    xlabel('$\gamma$')
    for j=1:3
        hold on
        subplot(3,1,j)
        ylim([0 inf])
        xlim([1000 inf])
    end
%     print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + nNfib*nTheta;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% maxDim evolution for each mu
%%% Calculate maxDim average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
maxDimx = zeros(nMu,nNfib,nTheta);
maxDimy = zeros(nMu,nNfib,nTheta);
maxDimz = zeros(nMu,nNfib,nTheta);

dataPath = '../data_stressVSfriction/maxDim/';
dataFile = cell(nTheta*nNfib,1);

for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.1 0.1 0.9 1])
        dataFile{(k-1)*nNfib+j} = ['fig4_cluster_dim/',fileNameArr{1},'_nfib',num2str(nfibArr(j)),'_maxDim_vs_strain_vs_mu'];
        for i=1:nMu
            
            % read files
            name = [dataPath,fileNameArr{1},'_maxDim_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File = fopen(name,'r');
            data = fscanf(File,'%f',[4 Inf])';
            fclose(File);
            % plot maxDim evolution with strain
            subplot(4,3,i)
            hold on
            scatter(data(:,1),data(:,2),markersize,'filled','MarkerFaceColor',rgb('Crimson'),'MarkerEdgeColor',rgb('Crimson'))
            scatter(data(:,1),data(:,3),markersize,'filled','MarkerFaceColor',rgb('MediumBlue'),'MarkerEdgeColor',rgb('MediumBlue'))
            scatter(data(:,1),data(:,4),markersize,'filled','MarkerFaceColor',rgb('DarkOrange'),'MarkerEdgeColor',rgb('DarkOrange'))
            maxDimx(i,j,k) = mean(data(:,2));
            maxDimy(i,j,k) = mean(data(:,3));
            maxDimz(i,j,k) = mean(data(:,4));
        end
    end
end

ind = 1;
for i=figStart:figStart+nNfib*nTheta-1
    figure(i)
    for j=1:nMu
        subplot(4,3,j)
        hold on
        box on
        ylim([0 inf])
        xlim([1000 inf])
        ylabel('$L_{max}$')
        xlabel('$\gamma$')
        title(['$\mu =$ ',num2str(muArr(j))])
        legend('$x$','$y$','$z$','location','bestoutside')
        set(gca,'fontsize',20)
    end
%     print(dataFile{ind},'-dpng')
    ind = ind + 1;
end
figStart = figStart + nNfib*nTheta;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% maxDim vs. mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
dataFile = cell(3*nTheta,1);
for i=1:3
    for k=1:nTheta
        dataFile{(i-1)*nTheta+k} = ['fig4_cluster_dim/',fileNameArr{1},'_nfib'];
        for j=1:nNfib
            dataFile{(i-1)*nTheta+k} = [ dataFile{(i-1)*nTheta+k},num2str(nfibArr(j)),'_'];
        end
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{1} = [dataFile{1},'maxDimx_vs_mu'];
ylabel('MaxDim $x$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,maxDimx(:,j,k),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{2} = [dataFile{2},'maxDimy_vs_mu'];
ylabel('MaxDim $y$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,maxDimy(:,j,k),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{3} = [dataFile{3},'maxDimz_vs_mu'];
ylabel('MaxDim $z$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,maxDimz(:,j,k),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

ind = 1;
for i=figStart:figStart+3-1
    figure(i)
    box on
    xlabel('$\mu$')
    xlim([0 inf])
    legend(thetaNfibLegendArr,'location','bestoutside')
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + 3;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% maxDim / Lbox vs. mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
dataFile = cell(3*nTheta,1);
for i=1:3
    for k=1:nTheta
        dataFile{(i-1)*nTheta+k} = ['fig4_cluster_dim/',fileNameArr{1},'_nfib'];
        for j=1:nNfib
            dataFile{(i-1)*nTheta+k} = [ dataFile{(i-1)*nTheta+k},num2str(nfibArr(j)),'_'];
        end
    end
end


figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('MaxDim')
dataFile{1} = [dataFile{1},'maxDimxFrac_vs_mu'];
ylabel('MaxDim $x$ / $L_{box}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,maxDimx(:,j,k)/lboxArr(j),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('MaxDim')
dataFile{2} = [dataFile{2},'maxDimyFrac_vs_mu'];
ylabel('MaxDim $y$ / $L_{box}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,maxDimy(:,j,k)/lboxArr(j),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('MaxDim')
dataFile{3} = [dataFile{3},'maxDimzFrac_vs_mu'];
ylabel('MaxDim $z$ / $L_{box}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,maxDimz(:,j,k)/lboxArr(j),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

ind = 1;
for i=figStart:figStart+3-1
    figure(i)
    box on
    xlabel('$\mu$')
    xlim([0 inf])
    legend(thetaNfibLegendArr,'location','bestoutside')
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + 3;
%}

name = ['sub4_',fileNameArr{1}]; 
save(name)