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
Dyy = zeros(nMu, nNfib); 
Dzz = zeros(nMu, nNfib); 
conf_Dyy = zeros(nMu, nNfib); 
conf_Dzz = zeros(nMu, nNfib); 

% figure info
figStart = 1;

%% calculate diffusivity
for j=1:nNfib
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_MSD_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File=fopen(name,'r');
        data = fscanf(File,'%f',[3 Inf])';
        fclose(File);
        t = data(:,1);
        MSDy = data(:,2);
        MSDz = data(:,3);      
        if length(t) < 100
            continue;
        end
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
figure('Units','Inches','Position',[3 3 6.5 6.5])
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
