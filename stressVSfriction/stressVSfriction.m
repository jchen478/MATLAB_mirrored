%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot all data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%{
% common parameters for theta6 and theta1
nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200 1293];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
thetaArr = [1];
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('Olive') rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue') rgb('Plum') rgb('Purple') };
%}


%%{
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
nAvg = 4; % only used for cluster data
markersize = 50;

for i=1:nMu
    muLegendArr{i} = ['$\mu =$ ',num2str(muArr(i))];
end
for i=1:nTheta
    for j=1:nNfib
        thetaNfibLegendArr{(i-1)*nNfib+j} = ['$(\theta_{eq},N_{fib}) =$ (0.',num2str(thetaArr(i)),', ',num2str(nfibArr(j)),')'];
    end
end
figStart = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S vs. strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = '../data_stressVSfriction/meanCount_min2/';
S = zeros(nMu,nNfib,nTheta);
se_S = zeros(nMu,nNfib,nTheta);

% storing file names
dataFile = cell(nTheta*nNfib,1);

for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.1 0.1 0.5 0.6])
        hold on;
        title(['$(\theta_{eq},N_{fib}) = (0.$',num2str(thetaArr(k)),', ',num2str(nfibArr(j)),'$)$'])
        dataFile{(k-1)*nNfib+j} = ['dynamic_figures/theta',num2str(thetaArr(k)),'_nfib',num2str(nfibArr(j)),'_S_vs_strain'];
        for i=1:nMu
            % read files
            name = [dataPath,'theta',num2str(thetaArr(k)),'_meanCount_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File = fopen(name,'r');
            dataUnsorted = fscanf(File,'%f',[4 Inf])';
            fclose(File);
            % remove friction coeffient column
            dataUnsorted(:,1) = [];
            % number of sample points
            nstrain = length(dataUnsorted(:,1));
            % sort data
            data = zeros(size(dataUnsorted));
            [data(:,1), sortI] = sort(dataUnsorted(:,1));
            for ii=1:nstrain
                data(ii,2) = dataUnsorted(sortI(ii),2);
                data(ii,3) = dataUnsorted(sortI(ii),3);
            end
            S_mu = zeros(nstrain,1);
            for ii=1:nstrain
                if(data(ii,2) == 0)
                    continue;
                end
                S_mu(ii) = data(ii,3)/data(ii,2);
            end
            % plot S evolution with strain
            scatter(data(:,1),S_mu,markersize,scatterMarker(i),...
                'MarkerFaceColor',colorArr{i},...
                'MarkerEdgeColor',colorArr{i})
            % calculate average
            Savg = mean(S_mu((nstrain-nAvg+1):end));
            se_Savg = std(S_mu((nstrain-nAvg+1):end))/sqrt(nAvg);
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
    ylim([0 inf])
    xlim([0 inf])
    set(gca,'yscale','log');
    legend(muLegendArr,'location','bestoutside')
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

% spreadfigures();
figStart = figStart + nNfib*nTheta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max dimension vs. strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataPath = '../data_stressVSfriction/clusterSize_min2/';
dataFile = cell(nTheta*nNfib,1);

for k=1:nTheta
    for j=1:nNfib
        
        figure('units','normalized','outerposition',[0.2 0.2 0.6 0.9])
        hold on;
        dataFile{(k-1)*nNfib+j} = ['dynamic_figures/theta',num2str(thetaArr(k)),'_nfib',num2str(nfibArr(j)),'_maxDim_vs_strain'];
        maxStrain = 0;
        subplot(3,1,1)
        hold on;
        title(['$(\theta_{eq},N_{fib}) = (0.$',num2str(thetaArr(k)),', ',num2str(nfibArr(j)),'$)$'])
        for i=1:nMu
            % read files
            name = [dataPath,'theta',num2str(thetaArr(k)),'_clusterSize_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File = fopen(name,'r');
            dataUnsorted = fscanf(File,'%f',[4 Inf])';
            fclose(File);
            % number of sample points
            nstrain = length(dataUnsorted(:,1));
            % sort data
            data = zeros(size(dataUnsorted));
            [data(:,1), sortI] = sort(dataUnsorted(:,1));
            for ii=1:nstrain
                data(ii,2:end) = dataUnsorted(sortI(ii),2:end);
            end
            % plot S evolution with strain
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
        %         box on
        ylim([0 inf])
        xlim([0 inf])
    end
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

% spreadfigures();
figStart = figStart + nNfib*nTheta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max dimension vs. strain / one plot for each mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxDimx = zeros(nMu,nNfib,nTheta);
maxDimy = zeros(nMu,nNfib,nTheta);
maxDimz = zeros(nMu,nNfib,nTheta);

dataPath = '../data_stressVSfriction/clusterSize_min2/';
dataFile = cell(nTheta*nNfib,1);

for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.1 0.1 0.8 1])
        dataFile{(k-1)*nNfib+j} = ['dynamic_figures/theta',num2str(thetaArr(k)),'_nfib',num2str(nfibArr(j)),'_maxDim_vs_strain_vs_mu'];
        %         title(['$(\theta_{eq},N_{fib}) = (0.$',num2str(thetaArr(k)),', ',num2str(nfibArr(j)),'$)$'])
        for i=1:nMu
            
            % read files
            name = [dataPath,'theta',num2str(thetaArr(k)),'_clusterSize_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File = fopen(name,'r');
            dataUnsorted = fscanf(File,'%f',[4 Inf])';
            fclose(File);
            % number of sample points
            nstrain = length(dataUnsorted(:,1));
            % sort data
            data = zeros(size(dataUnsorted));
            [data(:,1), sortI] = sort(dataUnsorted(:,1));
            for ii=1:nstrain
                data(ii,2:end) = dataUnsorted(sortI(ii),2:end);
            end
            % plot S evolution with strain
            subplot(4,3,i)
            hold on
            scatter(data(:,1),data(:,2),markersize,'filled')
            scatter(data(:,1),data(:,3),markersize,'filled')
            scatter(data(:,1),data(:,4),markersize,'filled')
            begin=length(data(:,2))-nAvg+1;
            maxDimx(i,j,k) = mean(data(begin:end,2));
            maxDimy(i,j,k) = mean(data(begin:end,3));
            maxDimz(i,j,k) = mean(data(begin:end,4));
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
        xlim([0 inf])
        ylabel('$L_{max}$')
        xlabel('$\gamma$')
        title(['$\mu =$ ',num2str(muArr(j))])
        legend('$x$','$y$','$z$','location','bestoutside')
        set(gca,'fontsize',20)
    end
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

% spreadfigures();
figStart = figStart + nNfib*nTheta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nCluster vs. strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = '../data_stressVSfriction/clusterStrain_min2/';
dataFile = cell(nTheta*nNfib,1);

% array used to store average values
nCluster = zeros(nMu,nNfib,nTheta);
se_nCluster = zeros(nMu,nNfib,nTheta);

for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.1 0.1 0.5 0.6])
        hold on;
        title(['$(\theta_{eq},N_{fib}) = (0.$',num2str(thetaArr(k)),', ',num2str(nfibArr(j)),'$)$'])
        dataFile{(k-1)*nNfib+j} = ['dynamic_figures/theta',num2str(thetaArr(k)),'_nfib',num2str(nfibArr(j)),'_nCluster_vs_strain'];
        for i=1:nMu
            % read files
            name = [dataPath,'theta',num2str(thetaArr(k)),'_clusterStrain_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File = fopen(name,'r');
            dataUnsorted = fscanf(File,'%f',[2 Inf])';
            fclose(File);
            % number of sample points
            nstrain = length(dataUnsorted(:,1));
            % sort data
            data = zeros(size(dataUnsorted));
            [data(:,1), sortI] = sort(dataUnsorted(:,1));
            for ii=1:nstrain
                data(ii,2) = dataUnsorted(sortI(ii),2);
            end
            % plot nCluster evolution with strain
            scatter(data(:,1),data(:,2),markersize,scatterMarker(i),...
                'MarkerFaceColor',colorArr{i},...
                'MarkerEdgeColor',colorArr{i})
            % calculate average
            nClusteravg = mean(data((nstrain-nAvg+1):end,2));
            se_nClusteravg = std(data((nstrain-nAvg+1):end,2))/sqrt(nAvg);
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
    xlim([0 inf])
    set(gca,'yscale','log');
    legend(muLegendArr,'location','bestoutside')
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

% spreadfigures();
figStart = figStart + nNfib*nTheta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% average values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = '../data_stressVSfriction/averagedValues/';

% pre-averaged values
particle_stress = zeros(nMu,nNfib,nTheta);
se_particle_stress = zeros(nMu,nNfib,nTheta);
relative_viscosity = zeros(nMu,nNfib,nTheta);
se_relative_viscosity = zeros(nMu,nNfib,nTheta);
norm1 = zeros(nMu,nNfib,nTheta);
norm2 = zeros(nMu,nNfib,nTheta);
number_of_contacts = zeros(nMu,nNfib,nTheta);
Eelas = zeros(nMu,nNfib,nTheta);
se_Eelas = zeros(nMu,nNfib,nTheta);


% find average
for k=1:nTheta
    for j=1:nNfib
        % read files
        name = [dataPath,'theta',num2str(thetaArr(k)),'_nfib',num2str(nfibArr(j)),'.txt'];
        File = fopen(name,'r');
        data = fscanf(File,'%f',[28 Inf])';
        fclose(File);
        
        % define variables
        nfib = data(:,1);
        nseg = data(:,2);
        rps = data(:,3);
        kb = data(:,4);
        side = data(:,6);
        sigmapstar = data(:,7);
        N1star = data(:,8);
        N2star = data(:,9);
        nc = data(:,10);
        se_sigmapstar = data(:,18);
        se_nc = data(:,19);
        elas = data(:,27);
        se_elas = data(:,28);
        
        % calculate relevant parameters
        [Seff,rp,nL3,L,gamma] = pCalc(nfib, nseg, rps, kb, side, a, EY, Imom, eta0);
        
        % Dimensionalize stresses
        [sigmap,sigma, se_sigmap, se_sigma, N1, N2] = stressDim(sigmapstar, se_sigmapstar, N1star, N2star, nseg, rps, eta0, gamma, nL3);
        
        % Calculate viscosity
        etarel = sigma./gamma;
        se_etarel = se_sigma./gamma;
        
        particle_stress(:,j,k) = sigmap;
        relative_viscosity(:,j,k) = etarel;
        norm1(:,j,k) = N1;
        norm2(:,j,k) = N2;
        se_particle_stress(:,j,k) = se_sigmap;
        se_relative_viscosity(:,j,k) = se_etarel;
        number_of_contacts(:,j,k) = nc;
        Eelas(:,j,k) = elas;
        se_Eelas(:,j,k) =  se_elas;
        
    end
end

%% Intensity fit
intensity = zeros(nMu,nNfib,nTheta);
dbin = 20;
for k=1:nTheta
    for j=1:nNfib
        intensity(:,j,k) = intensityFitFunc(thetaArr(k),nfibArr(j), dbin, lboxArr(j), nMu, j+nNfib*(k-1));
    end
end



%%{
dataFile = cell(1,1);
for i=1:1
    dataFile{i} = 'dynamic_figures/theta';
    for k=1:nTheta
        dataFile{i} = [dataFile{i},num2str(thetaArr(k)),'_'];
    end
    dataFile{i} = [dataFile{i},'nfib'];
    for j=1:nNfib
        dataFile{i} = [dataFile{i},num2str(nfibArr(j)),'_'];
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{1} = [dataFile{1},'I_vs_mu'];
title(['Intensity $d_{bin} = $ ',num2str(dbin)])
ylabel('$I$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,intensity(:,j,k),'-.o','MarkerSize',10,'Linewidth',2.5);
    end
end

ind = 1;
for i=figStart:figStart+1-1
    figure(i)
    box on
    xlabel('$\mu$')
    xlim([0 inf])
    legend(thetaNfibLegendArr,'location','bestoutside')
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + 1;

%}

%% plotting averaged values with mu as x-axis

dataFile = cell(9,1);
for i=1:9
    dataFile{i} = 'dynamic_figures/theta';
    for k=1:nTheta
        dataFile{i} = [dataFile{i},num2str(thetaArr(k)),'_'];
    end
    dataFile{i} = [dataFile{i},'nfib'];
    for j=1:nNfib
        dataFile{i} = [dataFile{i},num2str(nfibArr(j)),'_'];
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{1} = [dataFile{1},'I_vs_mu'];
title(['Intensity $d_{bin} = $',num2str(dbin)])
ylabel('$I$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,intensity(:,j,k),'-.o','MarkerSize',10,'Linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('Particle stress contribution')
dataFile{2} = [dataFile{2},'sigP_vs_mu'];
ylabel('$\sigma_{p,xz} L^4/ E_Y I$')
hold on
for k=1:nTheta
    for j=1:nNfib
        errorbar(muArr,particle_stress(:,j,k).*L.^4/EY/Imom,...
            se_particle_stress(:,j,k).*L.^4/EY/Imom/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('First normal stress difference')
dataFile{3} = [dataFile{3},'norm1_vs_mu'];
ylabel('$N_1 L^4/ E_Y I$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,norm1(:,j,k).*L.^4/EY/Imom,'-.o','MarkerSize',10,'Linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('Second normal stress difference')
dataFile{4} = [dataFile{4},'norm2_vs_mu'];
ylabel('$N_2 L^4/ E_Y I$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,norm2(:,j,k).*L.^4/EY/Imom,'-.o','MarkerSize',10,'Linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('Relative viscosity')
dataFile{5} = [dataFile{5},'etarel_vs_mu'];
ylabel('$\eta_{rel}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        errorbar(muArr,relative_viscosity(:,j,k),se_relative_viscosity(:,j,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('Number of contacts per fiber')
dataFile{6} = [dataFile{6},'nc_vs_mu'];
ylabel('$N_c$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,number_of_contacts(:,j,k),...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('Cluster size')
dataFile{7} = [dataFile{7},'S_vs_mu'];
ylabel('$\bar{S}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        errorbar(muArr,S(:,j,k),se_S(:,j,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    end
end
set(gca,'yscale','log')

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('Cluster count')
dataFile{8} = [dataFile{8},'nCluster_vs_mu'];
ylabel('$\bar{N}_{cluster}/N_{fib}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        errorbar(muArr,nCluster(:,j,k)./nfibArr(j),se_nCluster(:,j,k)/2./nfibArr(j),...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    end
end
% set(gca,'yscale','log')

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{9} = [dataFile{9},'Eelas_vs_mu'];
title('Elastic energy stored per fiber')
ylabel('$<E_{elastic}>/ 8\pi\eta_0\dot{\gamma}l^3$')
hold on
for k=1:nTheta
    for j=1:nNfib
        errorbar(muArr,Eelas(:,j,k),se_Eelas(:,j,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    end
end


ind = 1;
for i=figStart:figStart+9-1
    figure(i)
    box on
    xlabel('$\mu$')
    xlim([0 inf])
    legend(thetaNfibLegendArr,'location','bestoutside')
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + 9;

%% plotting averaged values with nfib as x-axis
dataFile = cell(9*nTheta,1);
for i=1:9
    for k=1:nTheta
        dataFile{(i-1)*nTheta+k} = ['dynamic_figures/theta',num2str(thetaArr(k)),'_nfib'];
        for j=1:nNfib
            dataFile{(i-1)*nTheta+k} = [ dataFile{(i-1)*nTheta+k},num2str(nfibArr(j)),'_'];
        end
    end
end

ind = 1;
for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'I_vs_nfib'];
    ind = ind+1;
    title(['Intensity with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$I$')
    hold on
    for j=1:nMu
        plot(nfibArr,intensity(j,:,k),'-.o',...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j},...
            'MarkerSize',10,'Linewidth',2.5);
    end
end


for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'sigP_vs_nfib'];
    ind = ind+1;
    title(['Particle stress contribution with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$\sigma_{p,xz} L^4/ E_Y I$')
    hold on
    for j=1:nMu
        errorbar(nfibArr,particle_stress(j,:,k)*L(1).^4/EY/Imom,...
            se_particle_stress(j,:,k)*L(1).^4/EY/Imom/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end



for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'norm1_vs_nfib'];
    ind = ind+1;
    title(['First normal stress difference with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$N_1 L^4/ E_Y I$')
    hold on
    for j=1:nMu
        plot(nfibArr,norm1(j,:,k)*L(1).^4/EY/Imom,'-.o',...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j},...
            'MarkerSize',10,'Linewidth',2.5);
    end
end

for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'norm2_vs_nfib'];
    ind = ind+1;
    title(['Second normal stress difference with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$N_2 L^4/ E_Y I$')
    hold on
    for j=1:nMu
        plot(nfibArr,norm2(j,:,k)*L(1).^4/EY/Imom,'-.o',...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j},...
            'MarkerSize',10,'Linewidth',2.5);
    end
end


for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'etarel_vs_nfib'];
    ind = ind+1;
    title(['Relative viscosity with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$\eta_{rel}$')
    hold on
    for j=1:nMu
        errorbar(nfibArr,relative_viscosity(j,:,k),se_relative_viscosity(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end

for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'nc_vs_nfib'];
    ind = ind+1;
    title(['Number of contacts per fiber with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$N_c$')
    hold on
    for j=1:nMu
        plot(nfibArr,number_of_contacts(j,:,k),...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end


for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'S_vs_nfib'];
    ind = ind+1;
    title(['Cluster size with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$\bar{S}$')
    set(gca,'yscale','log')
    hold on
    for j=1:nMu
        errorbar(nfibArr,S(j,:,k),se_S(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end


for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'nCluster_vs_nfib'];
    ind = ind+1;
    title(['Cluster count with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$\bar{N}_{cluster}$')
    set(gca,'yscale','log')
    hold on
    for j=1:nMu
        errorbar(nfibArr,nCluster(j,:,k),se_nCluster(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end

for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'Eelas_vs_nfib'];
    ind = ind+1;
    title(['Elastic energy with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$<E_{elastic}>/ 8\pi\eta_0\dot{\gamma}l^3$')
    hold on
    for j=1:nMu
        errorbar(nfibArr,Eelas(j,:,k),se_Eelas(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end

ind = 1;
for i=figStart:figStart+9*nTheta-1
    figure(i)
    box on
    xlabel('$N_{fib}$')
    xlim([0 inf])
    legend(muLegendArr,'location','bestoutside')
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + 9*nTheta;

%% plotting averaged values with L/Lbox as x-axis
dataFile = cell(9*nTheta,1);
for i=1:9
    for k=1:nTheta
        dataFile{(i-1)*nTheta+k} = ['dynamic_figures/theta',num2str(thetaArr(k)),'_nfib'];
        for j=1:nNfib
            dataFile{(i-1)*nTheta+k} = [ dataFile{(i-1)*nTheta+k},num2str(nfibArr(j)),'_'];
        end
    end
end

ind = 1;
for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'I_vs_overLbox'];
    ind = ind+1;
    title(['Intensity with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$I$')
    hold on
    for j=1:nMu
        plot(2*rpFiber./lboxArr,intensity(j,:,k),'-.o',...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j},...
            'MarkerSize',10,'Linewidth',2.5);
    end
end


for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'sigP_vs_overLbox'];
    ind = ind+1;
    title(['Particle stress contribution with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$\sigma_{p,xz} L^4/ E_Y I$')
    hold on
    for j=1:nMu
        errorbar(2*rpFiber./lboxArr,particle_stress(j,:,k)*L(1).^4/EY/Imom,...
            se_particle_stress(j,:,k)*L(1).^4/EY/Imom/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end



for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'norm1_vs_overLbox'];
    ind = ind+1;
    title(['First normal stress difference with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$N_1 L^4/ E_Y I$')
    hold on
    for j=1:nMu
        plot(2*rpFiber./lboxArr,norm1(j,:,k)*L(1).^4/EY/Imom,'-.o',...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j},...
            'MarkerSize',10,'Linewidth',2.5);
    end
end

for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'norm2_vs_overLbox'];
    ind = ind+1;
    title(['Second normal stress difference with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$N_2 L^4/ E_Y I$')
    hold on
    for j=1:nMu
        plot(2*rpFiber./lboxArr,norm2(j,:,k)*L(1).^4/EY/Imom,'-.o',...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j},...
            'MarkerSize',10,'Linewidth',2.5);
    end
end


for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'etarel_vs_overLbox'];
    ind = ind+1;
    title(['Relative viscosity with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$\eta_{rel}$')
    hold on
    for j=1:nMu
        errorbar(2*rpFiber./lboxArr,relative_viscosity(j,:,k),se_relative_viscosity(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end


for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'nc_vs_overLbox'];
    ind = ind+1;
    title(['Number of contacts per fiber with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$N_c$')
    hold on
    for j=1:nMu
        plot(2*rpFiber./lboxArr,number_of_contacts(j,:,k),...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end



for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'S_vs_overLbox'];
    ind = ind+1;
    title(['Cluster size with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$\bar{S}$')
    set(gca,'yscale','log')
    hold on
    for j=1:nMu
        errorbar(2*rpFiber./lboxArr,S(j,:,k),se_S(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end


for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'nCluster_vs_overLbox'];
    ind = ind+1;
    title(['Cluster count with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$\bar{N}_{cluster}$')
    set(gca,'yscale','log')
    hold on
    for j=1:nMu
        errorbar(2*rpFiber./lboxArr,nCluster(j,:,k),se_nCluster(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end

for k=1:nTheta
    figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
    dataFile{ind} = [dataFile{ind},'Eelas_vs_overLbox'];
    ind = ind+1;
    title(['Elastic energy per fiber with $\theta_{eq} = 0.$',num2str(thetaArr(k))])
    ylabel('$<E_{elastic}>/ 8\pi\eta_0\dot{\gamma}l^3$')
    hold on
    for j=1:nMu
        errorbar(2*rpFiber./lboxArr,Eelas(j,:,k),se_Eelas(j,:,k)/2,...
            '-.o','MarkerSize',10, 'linewidth',2.5,...
            'MarkerFaceColor',colorArr{j},...
            'MarkerEdgeColor',colorArr{j},...
            'color',colorArr{j});
    end
end

ind = 1;
for i=figStart:figStart+9*nTheta-1
    figure(i)
    box on
    xlabel('$L/L_{box}$')
    xlim([0 inf])
    legend(muLegendArr,'location','bestoutside')
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + 9*nTheta;

%% plotting maxDim with mu as x-axis
dataFile = cell(3*nTheta,1);
for i=1:3
    for k=1:nTheta
        dataFile{(i-1)*nTheta+k} = ['dynamic_figures/theta',num2str(thetaArr(k)),'_nfib'];
        for j=1:nNfib
            dataFile{(i-1)*nTheta+k} = [ dataFile{(i-1)*nTheta+k},num2str(nfibArr(j)),'_'];
        end
    end
end


figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{1} = [dataFile{1},'maxDimx_vs_mu'];
ylabel('Max $L_{floc,x}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,maxDimx(:,j,k),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{2} = [dataFile{2},'maxDimy_vs_mu'];
ylabel('Max $L_{floc,y}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,maxDimy(:,j,k),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
dataFile{3} = [dataFile{3},'maxDimz_vs_mu'];
ylabel('Max $L_{floc,z}$')
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

%% plotting maxDim fraction with mu as x-axis
dataFile = cell(3*nTheta,1);
for i=1:3
    for k=1:nTheta
        dataFile{(i-1)*nTheta+k} = ['dynamic_figures/theta',num2str(thetaArr(k)),'_nfib'];
        for j=1:nNfib
            dataFile{(i-1)*nTheta+k} = [ dataFile{(i-1)*nTheta+k},num2str(nfibArr(j)),'_'];
        end
    end
end


figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('MaxDim')
dataFile{1} = [dataFile{1},'maxDimxFrac_vs_mu'];
ylabel('$MaxDim_{x} / L_{box}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,maxDimx(:,j,k)/lboxArr(j),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('MaxDim')
dataFile{2} = [dataFile{2},'maxDimyFrac_vs_mu'];
ylabel('$MaxDim_{y} / L_{box}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,maxDimy(:,j,k)/lboxArr(j),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('MaxDim')
dataFile{3} = [dataFile{3},'maxDimzFrac_vs_mu'];
ylabel('$MaxDim_{z} / L_{box}$')
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
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + 3;

%% S / Nfib vs. mu
dataFile = cell(1,1);
for i=1:1
    dataFile{i} = 'dynamic_figures/theta';
    for k=1:nTheta
        dataFile{i} = [dataFile{i},num2str(thetaArr(k)),'_'];
    end
    dataFile{i} = [dataFile{i},'nfib'];
    for j=1:nNfib
        dataFile{i} = [dataFile{i},num2str(nfibArr(j)),'_'];
    end
end

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('Mean cluster size normalized by system size')
dataFile{1} = [dataFile{1},'Sfrac_vs_mu'];
ylabel('$S/N_{fib}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        plot(muArr,S(:,j,k)/nfibArr(j),'-.o','MarkerSize',10, 'linewidth',2.5);
    end
end

ind = 1;
for i=figStart:figStart+1-1
    figure(i)
    box on
    xlabel('$\mu$')
    xlim([0 inf])
    legend(thetaNfibLegendArr,'location','bestoutside')
    %     savefig([dataFile{ind},'.fig']);
    print(dataFile{ind},'-dpng')
    ind = ind + 1;
end
figStart = figStart + 1;