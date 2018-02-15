%%%
%%% Plot cluster count related properties
%%% -- properties include S, nCluster 
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

%{
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

nTheta = length(thetaArr);
nMu = length(muArr);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S vs. strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = '../data_stressVSfriction/meanS/';
S = zeros(nMu,nNfib,nTheta);
se_S = zeros(nMu,nNfib,nTheta);

% storing file names
dataFile = cell(nTheta*nNfib,1);

for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.1 0.1 0.5 0.6])
        hold on;
        dataFile{(k-1)*nNfib+j} = ['fig3_cluster_count/',fileNameArr{1},'_nfib',num2str(nfibArr(j)),'_S_vs_strain'];
        for i=1:nMu
            % read files
            name = [dataPath,fileNameArr{1},'_meanS_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File = fopen(name,'r');
            data = fscanf(File,'%f',[4 Inf])';
            fclose(File);
            % remove friction coeffient column
            data(:,1) = [];
            % number of sample points
            nstrain = length(data(:,1));           
            S_mu = data(:,3)./data(:,2);
            % plot S evolution with strain
            scatter(data(:,1),S_mu,markersize,scatterMarker(i),...
                'MarkerFaceColor',colorArr{i},...
                'MarkerEdgeColor',colorArr{i})
            % calculate average
            Savg = mean(S_mu);
            se_Savg = std(S_mu)/sqrt(length(S_mu));
            S(i,j,k) = Savg;
            se_S(i,j,k) = se_Savg;
        end
    end
end

ind = 1;
for i=figStart:figStart+nNfib*nTheta-1
    figure(i)
    box on
    xlabel('$\gamma$')
    ylabel('$S$')
    ylim([2 inf])
    xlim([1000 inf])
    set(gca,'yscale','log');
    legend(muLegendArr,'location','bestoutside')
%     print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

% spreadfigures();
figStart = figStart + nNfib*nTheta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nCluster vs. strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = '../data_stressVSfriction/nCluster/';
dataFile = cell(nTheta*nNfib,1);

% array used to store average values
nCluster = zeros(nMu,nNfib,nTheta);
se_nCluster = zeros(nMu,nNfib,nTheta);

for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.1 0.1 0.5 0.6])
        hold on;
        dataFile{(k-1)*nNfib+j} = ['fig3_cluster_count/',fileNameArr{1},'_nfib',num2str(nfibArr(j)),'_nCluster_vs_strain'];
        for i=1:nMu
            % read files
            name = [dataPath,fileNameArr{1},'_nCluster_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File = fopen(name,'r');
            data = fscanf(File,'%f',[2 Inf])';
            fclose(File);
            % plot nCluster evolution with strain
            scatter(data(:,1),data(:,2),markersize,scatterMarker(i),...
                'MarkerFaceColor',colorArr{i},...
                'MarkerEdgeColor',colorArr{i})           
            % calculate average
            nClusteravg = mean(data(:,2));
            se_nClusteravg = std(data(:,2))/sqrt(length(data(:,2)));
            nCluster(i,j,k) = nClusteravg;
            se_nCluster(i,j,k) = se_nClusteravg;
        end
    end
end

ind = 1;
for i=figStart:figStart+nNfib*nTheta-1
    figure(i)
    box on
    xlabel('$\gamma$')
    ylabel('$N_{cluster}$')
    ylim([0 inf])
    xlim([1000 inf])
%     set(gca,'yscale','log');
    legend(muLegendArr,'location','bestoutside')
%     print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

% spreadfigures();
figStart = figStart + nNfib*nTheta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% average values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plotting averaged values with mu as x-axis

dataFile = cell(2,1);
for i=1:2
    dataFile{i} = ['fig3_cluster_count/',fileNameArr{1},'_nfib'];
    for j=1:nNfib
        dataFile{i} = [dataFile{i},num2str(nfibArr(j)),'_'];
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{1} = [dataFile{1},'S_vs_mu'];
ylabel('$<S>$')
hold on
for k=1:nTheta
    for j=1:nNfib
        errorbar(muArr,S(:,j,k),se_S(:,j,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    end
end
set(gca,'yscale','log')

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{2} = [dataFile{2},'nClusterFrac_vs_mu'];
ylabel('$<N_{cluster}>/N_{fib}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        errorbar(muArr,nCluster(:,j,k)/nfibArr(j),se_nCluster(:,j,k)/2./nfibArr(j),...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

ind = 1;
for i=figStart:figStart+1
    figure(i)
    box on
    xlabel('$\mu$')
    xlim([0 inf])
    legend(thetaNfibLegendArr,'location','bestoutside')
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + 2;

name = ['sub3_',fileNameArr{1}]; 
save(name)