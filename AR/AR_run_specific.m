close all;
clc;
clear;

%% parameters
a = 1.60E-05;
Imom = pi*a^4/2;
EY = 1e9; 
eta0 = 1; 

%% read files
filename={'run007.txt'};
sizeA = [24 Inf];

File = fopen(filename{1},'r');
data = fscanf(File,'%f',sizeA)';
fclose(File);

%% define variables
nfib = data(:,1);
nseg = data(:,2); 
rps = data(:,3);
kb = data(:,4);
mu = data(:,5); 
side = data(:,6); 
tauxzstar = data(:,7);
se_tauxzstar = data(:,8);
I = data(:,9);
se_I = data(:,10);
Sx = data(:,11);
se_Sx = data(:,12);
Sy = data(:,13);
se_Sy = data(:,14);
Sz = data(:,15);
se_Sz = data(:,16);
Vx = data(:,17);
se_Vx = data(:,18);
Vy = data(:,19);
se_Vy = data(:,20);
Vz = data(:,21);
se_Vz = data(:,22);
N1star = data(:,23);
N2star = data(:,24);

%% simulation parameters
Seff = pi.*kb./nseg.^4; 
rp = rps.*nseg; 
nL3 = nfib .* (2.*rp./side).^3; 

L = 2*rps.*nseg*a;
gamma = EY*Imom ./ (Seff.*L.^4*eta0); 

% dimensionalize 
tauxz = tauxzstar*pi.*nL3./(6*nseg.^3.*log(2*rps))*eta0.*gamma;
se_tauxz = se_tauxzstar*pi.*nL3./(6*nseg.^3.*log(2*rps))*eta0.*gamma;
N1 =  N1star*pi.*nL3./(6*nseg.^3.*log(2*rps))*eta0.*gamma;
N2 =  N2star*pi.*nL3./(6*nseg.^3.*log(2*rps))*eta0.*gamma;

sigxz = tauxz + eta0*gamma;
se_sigxz = se_tauxz;

etasp = tauxz./gamma;
se_etasp = se_tauxz./gamma;

etarel = etasp+1;
se_etarel = se_tauxz./gamma;

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 24;
linewidth = 1.5;

%% plotting
color = {rgb('Crimson'),rgb('MediumBlue'),rgb('DarkOliveGreen')};

figure()
hold on;
errorbar(1./Seff,etasp,se_etasp/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color{2},'Linewidth',linewidth,'color',color{2});
set(gca,'xscale','log','yscale','log')
ylabel('\eta_{sp}')

figure()
hold on;
errorbar(1./Seff,etarel,se_etarel/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color{2},'Linewidth',linewidth,'color',color{2});
set(gca,'xscale','log','yscale','log')
ylabel('\eta_{rel}')

figure()
hold on;
errorbar(1./Seff,sigxz.*L.^4/EY/Imom,se_sigxz/2.*L.^4/EY/Imom,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color{2},'Linewidth',linewidth,'color',color{2});
set(gca,'xscale','log','yscale','log')
ylabel('\it{\sigma_{xz} L^4/ E_Y I}')

figure()
hold on;
errorbar(1./Seff,I,se_I/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color{2},'Linewidth',linewidth,'color',color{2});
set(gca,'xscale','log')
ylabel('\it{}I')

figure()
hold on;
errorbar(1./Seff,Sx,se_Sx/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color{2},'Linewidth',linewidth,'color',color{2});
set(gca,'xscale','log')
ylabel('\it{}S_x')

figure()
hold on;
errorbar(1./Seff,Sy,se_Sy/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color{2},'Linewidth',linewidth,'color',color{2});
set(gca,'xscale','log')
ylabel('\it{}S_y')

figure()
hold on;
errorbar(1./Seff,Sz,se_Sz/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color{2},'Linewidth',linewidth,'color',color{2});
set(gca,'xscale','log')
ylabel('\it{}S_z')

figure()
hold on;
plot(1./Seff,N1/EY,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color{2},'Linewidth',linewidth,'color',color{2});
set(gca,'xscale','log')
ylabel('\it{}N_1/E_Y')

figure()
hold on;
plot(1./Seff,N2/EY,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color{2},'Linewidth',linewidth,'color',color{2});
set(gca,'xscale','log')
ylabel('\it{}N_2/E_Y') 

%% Figure formatting
for i=1:9
    figure(i)
    hold on
    box on
    xlabel('1/S_{eff}')
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    xlim([0 inf])
end
