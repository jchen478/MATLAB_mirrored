%%%
%%% Plot averaged values
%%% -- properties include intensity, nc, sigP, eta, N1, N2, Eelas
%%%

clc;
clear;
close all;

% Define common parameters
simulation_cases;

% data path
dataPath = '../data_stressVSfriction/averagedValues/';

% replicate info
nReplicate = 2;
replicate_flag = zeros(nNfib,1);

% figure info
figStart = 1;

% array allocation
particle_stress = zeros(nMu,nNfib);
std_particle_stress = zeros(nMu,nNfib);
relative_viscosity = zeros(nMu,nNfib);
std_relative_viscosity = zeros(nMu,nNfib);
norm1 = zeros(nMu,nNfib);
std_norm1 = zeros(nMu,nNfib);
norm2 = zeros(nMu,nNfib);
std_norm2 = zeros(nMu,nNfib);
number_of_contacts = zeros(nMu,nNfib);
std_number_of_contacts = zeros(nMu,nNfib);
Eelas = zeros(nMu,nNfib);
std_Eelas = zeros(nMu,nNfib);
sigmapstar = zeros(nMu,nReplicate+1);
N1star = zeros(nMu,nReplicate+1);
N2star = zeros(nMu,nReplicate+1);
nc = zeros(nMu,nReplicate+1);
elas = zeros(nMu,nReplicate+1);

%% read in averaged values
for j=1:nNfib
    % read files
    name = [dataPath fileNameArr{1},'_nfib',num2str(nfibArr(j)),'.txt'];
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
        replicate_name = ['../data_stressVSfriction/averagedValues_replicate',num2str(r),'/',fileNameArr{1},'_nfib',num2str(nfibArr(j)),'.txt'];
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
    
    % assign to overall arrays
    particle_stress(:,j) = sigmap;
    std_particle_stress(:,j) = std_sigmap;
    relative_viscosity(:,j) = etarel;
    std_relative_viscosity(:,j) = std_etarel;
    norm1(:,j) = N1;
    std_norm1(:,j) = std_N1;
    norm2(:,j) = N2;
    std_norm2(:,j) = std_N2;
    number_of_contacts(:,j) = nc(:,1);
    std_number_of_contacts(:,j) = std_nc;
    Eelas(:,j) = elas(:,1);
    std_Eelas(:,j) = std_elas;
end

% Find intensity through fit function
intensity = zeros(nMu,nNfib);
std_intensity = zeros(nMu,nNfib);
dbin = 20;
I = zeros(nMu,nReplicate);
%%{
for j=1:nNfib
    intensityDataPath = '../data_stressVSfriction/intensityBin/';
    I(:,1) = intensityFitFunc(fileNameArr{1},...
        intensityDataPath, nfibArr(j), dbin, lboxArr(j), nMu);
    for r=1:replicate_flag(j)
        intensityDataPath = ['../data_stressVSfriction/intensityBin_replicate',num2str(r),'/'];
        I(:,1+r) = intensityFitFunc(fileNameArr{1},...
            intensityDataPath, nfibArr(j), dbin, lboxArr(j), nMu);
    end
    
    % average
    average_range = replicate_flag(j) + 1;
    
    std_I = std(I(:,1:average_range),0,2);
    I(:,1) = mean(I(:,1:average_range),2);
    
    intensity(:,j) = I(:,1);
    std_intensity(:,j) = std_I;
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot P vs. mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{

dataFile = cell(7,1);
for i=1:7
    dataFile{i} = 'fig1_averaged_values/';
    
    dataFile{i} = [dataFile{i},fileNameArr{1},'_'];
    
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
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,intensity(:,j),...
            '-.x','MarkerSize',10,'Linewidth',2.5')
    else
        errorbar(muArr,intensity(:,j),...
            std_intensity(:,j)/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
    end
end

figure('Units','Inches','Position',[0.5 4 6.5 3.5])
title('Particle stress contribution')
dataFile{2} = [dataFile{2},'sigP_vs_mu'];
ylabel('$\sigma_{p,xz} L^4/ E_Y I$')
hold on
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,particle_stress(:,j).*L.^4/EY/Imom,...
            '-.x','MarkerSize',10,'Linewidth',2.5')
    else
        errorbar(muArr,particle_stress(:,j).*L.^4/EY/Imom,...
            std_particle_stress(:,j).*L.^4/EY/Imom/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
    end
end

figure('Units','Inches','Position',[7 0.5 6.5 3.5])
title('First normal stress difference')
dataFile{3} = [dataFile{3},'norm1_vs_mu'];
ylabel('$N_1 L^4/ E_Y I$')
hold on
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,norm1(:,j).*L.^4/EY/Imom,...
            '-.x','MarkerSize',10,'Linewidth',2.5')
    else
        errorbar(muArr,norm1(:,j).*L.^4/EY/Imom,...
            std_norm1(:,j).*L.^4/EY/Imom/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
    end
end

figure('Units','Inches','Position',[7 4 6.5 3.5])
title('Second normal stress difference')
dataFile{4} = [dataFile{4},'norm2_vs_mu'];
ylabel('$N_2 L^4/ E_Y I$')
hold on
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,norm2(:,j).*L.^4/EY/Imom,...
            '-.x','MarkerSize',10,'Linewidth',2.5')
    else
        errorbar(muArr,norm2(:,j).*L.^4/EY/Imom,...
            std_norm2(:,j).*L.^4/EY/Imom/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
    end
end

figure('Units','Inches','Position',[13.5 0.5 6.5 3.5])
title('Relative viscosity')
dataFile{5} = [dataFile{5},'etarel_vs_mu'];
ylabel('$\eta_{rel}$')
hold on
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,relative_viscosity(:,j),...
            '-.x','MarkerSize',10,'Linewidth',2.5')
    else
        errorbar(muArr,relative_viscosity(:,j),...
            std_relative_viscosity(:,j)/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
    end
end

figure('Units','Inches','Position',[13.5 4 6.5 3.5])
title('Number of contacts per fiber')
dataFile{6} = [dataFile{6},'nc_vs_mu'];
ylabel('$N_c$')
hold on
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,number_of_contacts(:,j),...
            '-.x','MarkerSize',10,'Linewidth',2.5')
    else
        errorbar(muArr,number_of_contacts(:,j),...
            std_number_of_contacts(:,j)/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
    end
end

figure('Units','Inches','Position',[7 7 6.5 3.5])
dataFile{7} = [dataFile{7},'Eelas_vs_mu'];
title('Elastic energy stored per fiber')
ylabel('$<E_{elastic}>/ 8\pi\eta_0\dot{\gamma}l^3$')
hold on
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,Eelas(:,j),...
            '-.x','MarkerSize',10,'Linewidth',2.5')
    else
        errorbar(muArr,Eelas(:,j),std_Eelas(:,j)/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
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

% save(fileNameArr{1})
% name = ['sub1_',fileNameArr{1}];
% save(name)

