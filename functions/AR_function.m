function AR_function(filename,color)

%% parameters
a = 1.60E-05;
Imom = pi*a^4/4;
EY = 1e9;
eta0 = 1;

%% read files
sizeA = [26 Inf];

File = fopen(filename,'r');
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
N1star = data(:,8);
N2star = data(:,9);
nc = data(:,10);
I = data(:,11);
Sx = data(:,12);
Sy = data(:,13);
Sz = data(:,14);
Vx = data(:,15);
Vy = data(:,16);
Vz = data(:,17);
se_tauxzstar = data(:,18);
se_nc = data(:,19);
se_I = data(:,20);
se_Sx = data(:,21);
se_Sy = data(:,22);
se_Sz = data(:,23);
se_Vx = data(:,24);
se_Vy = data(:,25);
se_Vz = data(:,26);

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

%% second stiffness scale
scale2 = eta0*gamma/EY;

%% Figure formatting parameters
markersize = 12;
tickx = 1.5;
fontsize = 24;
linewidth = 2.5;

%% plotting

figure(1)
hold on;
errorbar(1./Seff,etasp,se_etasp/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
% set(gca,'xscale','log','yscale','log')
set(gca,'xscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}\eta_{sp}')
title('Specific viscosity')


figure(2)
hold on;
errorbar(1./Seff,etarel,se_etarel/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}\eta_{rel}')
title('Relative viscosity')

figure(3)
hold on;
errorbar(1./Seff,sigxz.*L.^4/EY/Imom,se_sigxz/2.*L.^4/EY/Imom,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{\sigma_{xz} L^4/ E_Y I}')
title('Bulk stress')

figure(4)
hold on;
errorbar(1./Seff,I,se_I/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}I')

figure(5)
hold on;
errorbar(1./Seff,Sx,se_Sx/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}S_x')

figure(6)
hold on;
errorbar(1./Seff,Sy,se_Sy/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}S_y')

figure(7)
hold on;
errorbar(1./Seff,Sz,se_Sz/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}S_z')

figure(8)
hold on;
plot(1./Seff,N1/EY,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}N_1 / E_Y')

figure(9)
hold on;
plot(1./Seff,N2/EY,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}N_2 / E_Y')

figure(10)
hold on;
errorbar(scale2,etarel,se_etarel/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{}\eta_{rel}')
title('Relative viscosity')

figure(11)
hold on;
errorbar(scale2,sigxz.*L.^4/EY/Imom,se_sigxz/2.*L.^4/EY/Imom,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{\sigma_{xz} L^4/ E_Y I}')
title('Bulk stress')

figure(12)
hold on;
errorbar(scale2,I,se_I/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{}I')

figure(13)
hold on;
errorbar(scale2,Sy,se_Sy/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{}S_y')

figure(14)
hold on;
errorbar(scale2,etasp,se_etasp/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
% set(gca,'xscale','log','yscale','log')
set(gca,'xscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{}\eta_{sp}')
title('Specific viscosity')

figure(15)
hold on;
plot(scale2,N1/EY,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{}N_1 / E_Y')

figure(16)
hold on;
plot(scale2,N2/EY,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{}N_2 / E_Y')

figure(17)
hold on;
errorbar(scale2,etarel,se_etarel/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{}\eta_{rel}')
title('Relative viscosity')

figure(18)
hold on;
errorbar(scale2,tauxz.*L.^4/EY/Imom,se_tauxz/2.*L.^4/EY/Imom,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
% set(gca,'xscale','log','yscale','log')
set(gca,'xscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{\sigma_{xz,p} L^4/ E_Y I}')
title('Particle stress contribution')

figure(19)
hold on;
errorbar(1./Seff,tauxz.*L.^4/EY/Imom,se_tauxz/2.*L.^4/EY/Imom,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
% set(gca,'xscale','log','yscale','log')
set(gca,'xscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{\sigma_{xz,p} L^4/ E_Y I}')
title('Particle stress contribution')

figure(20)
hold on;
errorbar(1./Seff,tauxz./(eta0*gamma),se_tauxz/2./(eta0*gamma),'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
% set(gca,'xscale','log','yscale','log')
set(gca,'xscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{\sigma_{xz,p} / \eta_0 d\gamma/dt}')
title('Particle stress contribution')

figure(21)
hold on;
errorbar(1./Seff,nc,se_nc/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}N_C')
title('Number of contacts per fiber')

end