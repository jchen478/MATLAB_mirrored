%%%
%%% Plot averaged values
%%% -- properties include intensity, nc, sigP, eta, N1, N2, Eelas
%%%

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%%%%%%%%%%%%%%%%%%%% U-shaped fibers %%%%%%%%%%%%%%%%%%%%%
% fileNameArr = {'theta0'};thetaArr = 0;
% fileNameArr = {'theta1'}; thetaArr = 1;
% fileNameArr = {'theta3'}; thetaArr = 3;
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

% replicate info
nReplicate = 2; 
replicate_flag = zeros(nNfib,1); 

figStart = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataPath = '../data_stressVSfriction/averagedValues/';

% allocate space
particle_stress = zeros(nMu,nNfib,nTheta);
std_particle_stress = zeros(nMu,nNfib,nTheta);
relative_viscosity = zeros(nMu,nNfib,nTheta);
std_relative_viscosity = zeros(nMu,nNfib,nTheta);
norm1 = zeros(nMu,nNfib,nTheta);
std_norm1 = zeros(nMu,nNfib,nTheta);
norm2 = zeros(nMu,nNfib,nTheta);
std_norm2 = zeros(nMu,nNfib,nTheta);
number_of_contacts = zeros(nMu,nNfib,nTheta);
std_number_of_contacts = zeros(nMu,nNfib,nTheta);
Eelas = zeros(nMu,nNfib,nTheta);
std_Eelas = zeros(nMu,nNfib,nTheta);

sigmapstar = zeros(nMu,nReplicate+1); 
N1star = zeros(nMu,nReplicate+1); 
N2star = zeros(nMu,nReplicate+1); 
nc = zeros(nMu,nReplicate+1); 
elas = zeros(nMu,nReplicate+1); 

% read in averaged
for k=1:nTheta
    for j=1:nNfib
        % read files
        name = [dataPath,fileNameArr{k},'_nfib',num2str(nfibArr(j)),'.txt'];
        File = fopen(name,'r');
        data = fscanf(File,'%f',[28 Inf])';
        fclose(File);
        
        % define variables
        nfib = data(:,1);
        nseg = data(:,2);
        rps = data(:,3);
        kb = data(:,4);
        side = data(:,6);
        sigmapstar(:,1) = data(:,7);
        N1star(:,1) = data(:,8);
        N2star(:,1) = data(:,9);
        nc(:,1) = data(:,10);
        elas(:,1) = data(:,27);
        
        for r=1:nReplicate
            replicate_name = ['../data_stressVSfriction/averagedValues_replicate',num2str(r),'/',fileNameArr{k},'_nfib',num2str(nfibArr(j)),'.txt'];
            if exist(replicate_name, 'file') == 2
                % if replicates exist, then average the values and find std
                replicate_flag(j) = r;
                File = fopen(replicate_name,'r');
                data = fscanf(File,'%f',[28 Inf])';
                fclose(File);
                sigmapstar(:,r+1) = data(:,7);
                N1star(:,r+1) = data(:,8);
                N2star(:,r+1) = data(:,9);
                nc(:,r+1) = data(:,10);
                elas(:,r+1) = data(:,27);
            end
        end
        
        % average
        average_range = replicate_flag(j) + 1;
        
        std_sigmapstar = std(sigmapstar(:,1:average_range),0,2);
        sigmapstar(:,1) = mean(sigmapstar(:,1:average_range),2);
        
        std_N1star = std(N1star(:,1:average_range),0,2);
        N1star(:,1) = mean(N1star(:,1:average_range),2);
        
        std_N2star = std(N2star(:,1:average_range),0,2);
        N2star(:,1) = mean(N2star(:,1:average_range),2);
        
        std_nc = std(nc(:,1:average_range),0,2);
        nc(:,1) = mean(nc(:,1:average_range),2);
        
        std_elas = std(elas(:,1:average_range),0,2);
        elas(:,1) = mean(elas(:,1:average_range),2);
        
        % calculate relevant parameters
        [Seff,rp,nL3,L,gamma] = pCalc(nfib, nseg, rps, kb, side, a, EY, Imom, eta0);
        
        % Dimensionalize stresses
        [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = stressDim(sigmapstar(:,1), std_sigmapstar, N1star(:,1), N2star(:,1), nseg, rps, eta0, gamma, nL3);
        [dum, dum1, dum2, dum3, std_N1, std_N2] = stressDim(sigmapstar(:,1), std_sigmapstar, std_N1star, std_N2star, nseg, rps, eta0, gamma, nL3);
        
        % Calculate viscosity
        etarel = sigma./gamma;
        std_etarel = std_sigma./gamma; 
        
        particle_stress(:,j,k) = sigmap;
        std_particle_stress(:,j,k) = std_sigmap;
        relative_viscosity(:,j,k) = etarel;
        std_relative_viscosity(:,j,k) = std_etarel;
        norm1(:,j,k) = N1;
        std_norm1(:,j,k) = std_N1;
        norm2(:,j,k) = N2;
        std_norm2(:,j,k) = std_N2;
        number_of_contacts(:,j,k) = nc(:,1);
        std_number_of_contacts(:,j,k) = std_nc;
        Eelas(:,j,k) = elas(:,1);
        std_Eelas(:,j,k) = std_elas;
      
        
    end
end

% Find intensity through fit function
intensity = zeros(nMu,nNfib,nTheta);
std_intensity = zeros(nMu,nNfib,nTheta);
dbin = 20;
I = zeros(nMu,nReplicate); 
%{
for k=1:nTheta
    for j=1:nNfib
        intensityDataPath = '../data_stressVSfriction/intensityBin/';
        I(:,1) = intensityFitFunc(fileNameArr{k},...
            intensityDataPath, nfibArr(j), dbin, lboxArr(j), nMu);
        for r=1:replicate_flag(j)
            intensityDataPath = ['../data_stressVSfriction/intensityBin_replicate',num2str(r),'/'];            
            I(:,1+r) = intensityFitFunc(fileNameArr{k},...
                intensityDataPath, nfibArr(j), dbin, lboxArr(j), nMu);            
        end
        
        % average
        average_range = replicate_flag(j) + 1;
        
        std_I = std(I(:,1:average_range),0,2);
        I(:,1) = mean(I(:,1:average_range),2);
        
        intensity(:,j,k) = I(:,1);
        std_intensity(:,j,k) = std_I;
    end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot P vs. mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{

dataFile = cell(7,1);
for i=1:7
    dataFile{i} = 'fig15_averaged_values_replicates/';
    for k=1:nTheta
        dataFile{i} = [dataFile{i},fileNameArr{k},'_'];
    end
    dataFile{i} = [dataFile{i},'nfib'];
    for j=1:nNfib
        dataFile{i} = [dataFile{i},num2str(nfibArr(j)),'_'];
    end
end

figure('Units','Inches','Position',[0.5 0.5 6.5 3.5])
dataFile{1} = [dataFile{1},'I_vs_mu'];
title(['Intensity $d_{bin} = $',num2str(dbin)])
ylabel('$I$')
hold on
for k=1:nTheta
    for j=1:nNfib        
         if replicate_flag(j) == 0
             plot(muArr,intensity(:,j,k),...
                 '-.x','MarkerSize',10,'Linewidth',2.5')
        else
            errorbar(muArr,intensity(:,j,k),...
            std_intensity(:,j,k)/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
        end    
    end
end

figure('Units','Inches','Position',[0.5 4 6.5 3.5])
title('Particle stress contribution')
dataFile{2} = [dataFile{2},'sigP_vs_mu'];
ylabel('$\sigma_{p,xz} L^4/ E_Y I$')
hold on
for k=1:nTheta
    for j=1:nNfib
        if replicate_flag(j) == 0
             plot(muArr,particle_stress(:,j,k).*L.^4/EY/Imom,...
                '-.x','MarkerSize',10,'Linewidth',2.5')
        else
            errorbar(muArr,particle_stress(:,j,k).*L.^4/EY/Imom,...
            std_particle_stress(:,j,k).*L.^4/EY/Imom/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
        end        
    end
end

figure('Units','Inches','Position',[7 0.5 6.5 3.5])
title('First normal stress difference')
dataFile{3} = [dataFile{3},'norm1_vs_mu'];
ylabel('$N_1 L^4/ E_Y I$')
hold on
for k=1:nTheta
    for j=1:nNfib       
         if replicate_flag(j) == 0
            plot(muArr,norm1(:,j,k).*L.^4/EY/Imom,...
                '-.x','MarkerSize',10,'Linewidth',2.5')
        else
            errorbar(muArr,norm1(:,j,k).*L.^4/EY/Imom,...
            std_norm1(:,j,k).*L.^4/EY/Imom/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
        end        
    end
end

figure('Units','Inches','Position',[7 4 6.5 3.5])
title('Second normal stress difference')
dataFile{4} = [dataFile{4},'norm2_vs_mu'];
ylabel('$N_2 L^4/ E_Y I$')
hold on
for k=1:nTheta
    for j=1:nNfib
        if replicate_flag(j) == 0
            plot(muArr,norm2(:,j,k).*L.^4/EY/Imom,...
             '-.x','MarkerSize',10,'Linewidth',2.5')
        else
            errorbar(muArr,norm2(:,j,k).*L.^4/EY/Imom,...
                std_norm2(:,j,k).*L.^4/EY/Imom/2,...
               '-.o','MarkerSize',8, 'linewidth',2.5);
        end
    end
end

figure('Units','Inches','Position',[13.5 0.5 6.5 3.5])
title('Relative viscosity')
dataFile{5} = [dataFile{5},'etarel_vs_mu'];
ylabel('$\eta_{rel}$')
hold on
for k=1:nTheta
    for j=1:nNfib
        if replicate_flag(j) == 0
            plot(muArr,relative_viscosity(:,j,k),...
               '-.x','MarkerSize',10,'Linewidth',2.5')
        else
            errorbar(muArr,relative_viscosity(:,j,k),...
                std_relative_viscosity(:,j,k)/2,...
                 '-.o','MarkerSize',8, 'linewidth',2.5);
        end
    end
end

figure('Units','Inches','Position',[13.5 4 6.5 3.5])
title('Number of contacts per fiber')
dataFile{6} = [dataFile{6},'nc_vs_mu'];
ylabel('$N_c$')
hold on
for k=1:nTheta
    for j=1:nNfib
        if replicate_flag(j) == 0
            plot(muArr,number_of_contacts(:,j,k),...
              '-.x','MarkerSize',10,'Linewidth',2.5')
        else
            errorbar(muArr,number_of_contacts(:,j,k),...
                std_number_of_contacts(:,j,k)/2,...
               '-.o','MarkerSize',8, 'linewidth',2.5);
        end
    end
end

figure('Units','Inches','Position',[7 7 6.5 3.5])
dataFile{7} = [dataFile{7},'Eelas_vs_mu'];
title('Elastic energy stored per fiber')
ylabel('$<E_{elastic}>/ 8\pi\eta_0\dot{\gamma}l^3$')
hold on
for k=1:nTheta
    for j=1:nNfib        
        if replicate_flag(j) == 0
             plot(muArr,Eelas(:,j,k),...
             '-.x','MarkerSize',10,'Linewidth',2.5')
        else
             errorbar(muArr,Eelas(:,j,k),std_Eelas(:,j,k)/2,...
             '-.o','MarkerSize',8, 'linewidth',2.5);
        end
    end
end

ind = 1;
for i=figStart:figStart+6
    figure(i)
    box on
    xlabel('$\mu$')
    xlim([0 inf])
    legend(thetaNfibLegendArr,'location','bestoutside')
    h = figure(i); 
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,dataFile{ind},'-dpdf','-r0')
    ind = ind + 1;
end

figStart = figStart + 7;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot P vs. nfib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

dataFile = cell(7*nTheta,1);
for i=1:7
    for k=1:nTheta
        dataFile{(i-1)*nTheta+k} = ['fig1_averaged_values/',fileNameArr{k},'_nfib'];
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
    dataFile{ind} = [dataFile{ind},'Eelas_vs_nfib'];
    ind = ind+1;
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
for i=figStart:figStart+7*nTheta-1
    figure(i)
    box on
    xlabel('$N_{fib}$')
    xlim([0 inf])
    legend(muLegendArr,'location','bestoutside')
    %     savefig([dataFile{ind},'.fig']);
    %     print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + 7*nTheta;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot P vs. L/Lbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
dataFile = cell(7*nTheta,1);
for i=1:7
    for k=1:nTheta
        dataFile{(i-1)*nTheta+k} = ['fig1_averaged_values/',fileNameArr{k},'_nfib'];
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
    dataFile{ind} = [dataFile{ind},'Eelas_vs_overLbox'];
    ind = ind+1;
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
for i=figStart:figStart+7*nTheta-1
    figure(i)
    box on
    xlabel('$L/L_{box}$')
    xlim([0 inf])
    legend(muLegendArr,'location','bestoutside')
    %     savefig([dataFile{ind},'.fig']);
    %     print(dataFile{ind},'-dpng')
    ind = ind + 1;
end

figStart = figStart + 7*nTheta;
%}

save(fileNameArr{1})


name = ['sub15_',fileNameArr{1}]; 
save(name)

