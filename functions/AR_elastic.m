function AR_elastic(filename,color)

%% parameters
a = 1.60E-05;
Imom = pi*a^4/4;
EY = 1e9;
eta0 = 1;

%% read files
sizeA = [12 Inf];

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
Ebstar = data(:,7);
Etstar = data(:,8);
Etotstar = data(:,9);
se_Ebstar = data(:,10);
se_Etstar = data(:,11);
se_Etotstar = data(:,12);

%% simulation parameters
Seff = pi.*kb./nseg.^4;
rp = rps.*nseg;
nL3 = nfib .* (2.*rp./side).^3;

l = rps*a;
L = 2*rps.*nseg*a;
gamma = EY*Imom ./ (Seff.*L.^4*eta0);

%% second stiffness scale
scale2 = eta0*gamma/EY;

%% Dimensionalize
Etot = Etotstar/2*8*pi*eta0.*l.^3.*gamma;
se_Etot = se_Etotstar/2*8*pi*eta0.*l.^3.*gamma;

%% Figure formatting parameters
markersize = 12;
tickx = 1.5;
fontsize = 24;
linewidth = 2.5;

%% plotting

figure(1)
subplot(2,1,1)
hold on;
errorbar(1./Seff,Ebstar,se_Ebstar/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}E_{b,elastic}*')
title('Bending elastic energy')

subplot(2,1,2)
hold on;
errorbar(1./Seff,Etstar,se_Etstar/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}E_{t,elastic}*')
title('Twisting elastic energy')

figure(2)
subplot(2,1,1)
hold on;
errorbar(scale2,Ebstar,se_Ebstar/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{}E_{b,elastic}*')
title('Bending elastic energy')

subplot(2,1,2)
hold on;
errorbar(scale2,Etstar,se_Etstar/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{}E_{t,elastic}*')
title('Twisting elastic energy')

figure(3)
hold on;
errorbar(1./Seff,Etotstar,se_Etotstar/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}E_{tot,elastic}*')
title('Total elastic energy')

figure(4)
hold on;
errorbar(scale2,Etotstar,se_Etotstar/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}\mu d\gamma/dt / E_Y')
ylabel('\it{}E_{tot,elastic}*')
title('Total elastic energy')

figure(5)
hold on;
errorbar(1./Seff,Etot,se_Etot/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
set(gca,'xscale','log','yscale','log')
xlabel('\it{}1/S_{eff}')
ylabel('\it{}E_{tot,elastic} (J)')
title('Total elastic energy')

end