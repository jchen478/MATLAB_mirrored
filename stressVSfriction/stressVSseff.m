%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot all data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

%% common parameters
%%{
nfibArr = [160 1280 3200 6400 10240];
lboxArr = [300 600 814.3 1026 1200];
muArr = [0 1 2 3 4 5 10 15 20];
seffArr = [0.001 0.01 0.05];
seffNameArr = {'02','2','10'}; 
rpFiber = 75;
%}

%%{
nfibArr = [160 1280 3200];
lboxArr = [300 600 814.3];
seffArr = [0.001 0.01 0.05];
seffNameArr = {'02' '2' '10'}; 
%}

%{
nfibArr = [160 1280 3200 6400];
lboxArr = [300 600 814.3 1026];
seffArr = [0.001 0.05 ];
seffNameArr = {'02' '10'}; 

%}



%% plotting specs
scatterMarker=['o' 'v' '>' '<' '^' 's' 'd' 'v' '>' '<' '^' 's' 'd' 'o' ];
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed') rgb('Orange') rgb('Gold')...
    rgb('Lime') rgb('DarkGreen') rgb('LightSkyBlue') rgb('Plum')};

%% fiber dimensions
a = 1.60E-05;       % radius (m)
Imom = pi*a^4/4;    % area moment (m^4)
EY = 1e9;           % Young's modulus (Pa m^2)
eta0 = 1;           % fluid viscosity (Pa s)

nMu = length(muArr);
nLbox = length(lboxArr);
nNfib = length(nfibArr);
nSeff = length(seffArr); 
muLegendArr = cell(nMu,1);
nfibLegendArr = cell(nNfib,1);
nAvg = 5; % only used for cluster data
markersize = 50;
for i=1:nMu
    muLegendArr{i} = ['$\mu =$ ',num2str(muArr(i))];
end
for j=1:nNfib
    nfibLegendArr{j} = ['$N_{fib} =$ ',num2str(nfibArr(j))];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find averaged values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
% allocate arrays for pre-averaged values
particle_stress = zeros(nMu,nNfib,nSeff);
se_particle_stress = zeros(nMu,nNfib,nSeff);
relative_viscosity = zeros(nMu,nNfib,nSeff);
se_relative_viscosity = zeros(nMu,nNfib,nSeff);
norm1 = zeros(nMu,nNfib,nSeff);
norm2 = zeros(nMu,nNfib,nSeff);
number_of_contacts = zeros(nMu,nNfib,nSeff);
Eelas = zeros(nMu,nNfib,nSeff);
se_Eelas = zeros(nMu,nNfib,nSeff);

% read into array
dataPath = '../data_stressVSseff/';

for j=1:nNfib
    for k=1:nSeff
        % read files
        name = [dataPath,'theta6_nfib',num2str(nfibArr(j)),'_kb',seffNameArr{k},'.txt']; 
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
%}

%% Normalize intensity

%{
dataPath = '../data_stressVSfriction/averagedValues/';
intensityBase = zeros(1,nNfib,1);
for j=1:nNfib
    % read files
    name = [dataPath,'theta6_nfib',num2str(nfibArr(j)),'_baseI.txt'];
    File = fopen(name,'r');
    data = fscanf(File,'%f',[1 Inf])';
    fclose(File);
    intensityBase(1,j,1) = mean(data);
    
end
%}

%{
intensityNorm = zeros(nMu,nNfib,nSeff);
se_intensityNorm = zeros(nMu,nNfib,nSeff);
for j=1:nNfib 
    base = intensityBase(1,j,1);
    intensityNorm(:,j,:) = (intensity(:,j,:)-base)/base;
    se_intensityNorm(:,j,:) = se_intensity(:,j,:)/base;
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot averaged values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{ Stress
figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9])
for i=1:nMu
    subplot(3,3,i)
    hold on
    for j=1:nNfib
        plot(1./seffArr,squeeze(particle_stress(i,j,:)*L(1).^4/EY/Imom),'-o')
    end
end
for i=1:nMu
    subplot(3,3,i)
    hold on
    title(['$\mu = $',num2str(muArr(i))],'fontsize',16);
    legend(nfibLegendArr,'location','best','fontsize',16);
    xlabel('$1/S_{eff}$')
    ylabel('$\sigma_{p,xz} L^4/ E_Y I$')
    set(gca,'fontsize',16)
end
%}
%%{ Number of contacts
figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9])
for i=1:nMu
    subplot(3,3,i)
    hold on
    for j=1:nNfib
        plot(1./seffArr,squeeze(number_of_contacts(i,j,:)),'-o')
    end
end
for i=1:nMu
    subplot(3,3,i)
    hold on
    title(['$\mu = $',num2str(muArr(i))],'fontsize',16);
    legend(nfibLegendArr,'location','best','fontsize',16);
    xlabel('$1/S_{eff}$')
    ylabel('$N_{C}$')
    set(gca,'fontsize',16)
end
%}

%{
% Intensity
figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9])
for i=1:nMu
    subplot(3,3,i)
    hold on
    for j=1:nNfib
        plot(1./seffArr,squeeze(intensity(i,j,:)),'-o')
    end
end
for i=1:nMu
    subplot(3,3,i)
    hold on
    title(['$\mu = $',num2str(muArr(i))],'fontsize',16);
    legend(nfibLegendArr,'location','best','fontsize',16);
    xlabel('$1/S_{eff}$')
    ylabel('$I$')
    set(gca,'fontsize',16)
end
%}

%{ 
% NormalizedIntensity
figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9])
for i=1:nMu
    subplot(3,3,i)
    hold on
    for j=1:nNfib
        plot(1./seffArr,squeeze(intensityNorm(i,j,:)),'-o')
    end
end
for i=1:nMu
    subplot(3,3,i)
    hold on
    title(['$\mu = $',num2str(muArr(i))],'fontsize',16);
    legend(nfibLegendArr,'location','best','fontsize',16);
    xlabel('$1/S_{eff}$')
    ylabel('$I_{norm}$')
    set(gca,'fontsize',16)
end
%}
%%{ Elastic energy
figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9])
for i=1:nMu
    subplot(3,3,i)
    hold on
    for j=1:nNfib
        plot(1./seffArr,squeeze(Eelas(i,j,:)),'-o')
    end
end
for i=1:nMu
    subplot(3,3,i)
    hold on
    title(['$\mu = $',num2str(muArr(i))],'fontsize',16);
    legend(nfibLegendArr,'location','best','fontsize',16);
    xlabel('$1/S_{eff}$')
    ylabel('$<E_{elastic}>/ 8\pi\eta_0\dot{\gamma}l^3$')
    set(gca,'fontsize',16)
end
%}