function [] = relax_yield( filename, color)

%%%%%% Defining Stress Symbols
%
%      sigmapstar 
%           output from simulation (nondimensionalized by eta0 * gamma)
%      sigmap
%           dimensionalized particle stress
%      sigma 
%           dimensionalized bulk stress
%      sigmastar
%           dimensionless bulk stress (by eta0 * gamma)
%

%% read files
File = fopen(filename,'r');
data = fscanf(File,'%f',[26 Inf])';
fclose(File);

%% define variables
nfib = data(:,1);
nseg = data(:,2);
rps = data(:,3);
kb = data(:,4);
mu = data(:,5);
side = data(:,6);
sigmapstar = data(:,7);
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
se_sigmapstar = data(:,18);
se_nc = data(:,19);
se_I = data(:,20);
se_Sx = data(:,21);
se_Sy = data(:,22);
se_Sz = data(:,23);
se_Vx = data(:,24);
se_Vy = data(:,25);
se_Vz = data(:,26);

%% fiber dimensions
a = 1.60E-05;       % radius (m)
Imom = pi*a^4/4;    % area moment (m^4)
EY = 1e9;           % Young's modulus (Pa m^2)
eta0 = 1;           % fluid viscosity (Pa s)

%% calculate relevant parameters
[Seff,rp,nL3,L,gamma] = pCalc(nfib, nseg, rps, kb, side, a, EY, Imom, eta0);

%% Dimensionalize stresses
[sigmap,sigma, se_sigmap, se_sigma, N1, N2] = stressDim(sigmapstar, se_sigmapstar, N1star, N2star, nseg, rps, eta0, gamma, nL3); 

%% Normalize intensity
I_normalized = I/I(1);
se_I_normalized = se_I/I(1);

%% Plotting parameters
markersize = 10; 
linewidth = 2.5; 

%% Plotting 

%% intensity I
figure(1)
hold on;
ylabel('\it{}I')
errorbar(1./Seff,I,se_I/2,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
title('Intensity')

%% normalized I
figure(2)
hold on;
ylabel('\it{}I_{normalized}')
errorbar(1./Seff,I_normalized,se_I_normalized/2,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
title('Normalized intensity')

%% first normal difference
figure(3)
hold on;
ylabel('\it{N_1 L^4/ E_Y I}')
plot(1./Seff,N1.*L.^4/EY/Imom,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
title('First normal stress')

%% second normal difference
figure(4)
hold on;
ylabel('\it{N_2 L^4/ E_Y I}')
plot(1./Seff,N2.*L.^4/EY/Imom,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
title('Second normal stress')

%% Particle stress contribution
figure(5)
hold on;
ylabel('\it{\sigma_{xz,p} L^4/ E_Y I}')
errorbar(1./Seff,sigmap.*L.^4/EY/Imom,se_sigmap.*L.^4/EY/Imom/2,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
title('Particle stress contribution')

%% Fitting parameters
x = 1./Seff; 
y = sigmap.*L.^4/EY/Imom; 

%% Linear fit
linFit = polyfit(x,y,1);
y_linFit = linFit(1)*x+linFit(2);

m = linFit(1)
b = linFit(2)

figure(6)
hold on
scatter(x,y,60,'filled','MarkerEdgeColor',color,'MarkerFaceColor',color)
plot(x,y_linFit,'-x','linewidth',linewidth,'markersize',10,'color',color)
ylabel('\it{\sigma_{xz,p} L^4/ E_Y I}')

RSSlin = (y-y_linFit).^2; 

figure(7)
title('Normalized squared residual')
hold on

plot(x,RSSlin/max(RSSlin),'-o','linewidth',linewidth,'markersize',10,'color',color)
ylabel('RS')

sum(RSSlin)

%% Casson model 
%  sqrt(tau) = sqrt(tau_0) + sqrt(eta0*gamma)
% x = sqrt(1./Seff); 
% y = sqrt(sigmap.*L.^4/EY/Imom); 
% 
% linFit = polyfit(x,y,1);
% y_linFit = linFit(1)*x+linFit(2);
% y_linFit = y_linFit.^2;
% 
% linFit(1)
% linFit(2)^2
% 
% x = 1./Seff; 
% y = sigmap.*L.^4/EY/Imom; 
% 
% figure(6)
% hold on
% scatter(x,y,60,'filled','MarkerEdgeColor',color,'MarkerFaceColor',color)
% plot(x,y_linFit,'-x','linewidth',linewidth,'markersize',10,'color',color)
% ylabel('\it{\sigma_{xz,p} L^4/ E_Y I}')
% 
% RSSlin = (y-y_linFit).^2; 
% 
% figure(7)
% title('Normalized squared residual')
% hold on
% 
% plot(x,RSSlin/max(RSSlin),'-o','linewidth',linewidth,'markersize',10,'color',color)
% ylabel('RS')
% 
% sum(RSSlin)

%% Power-law

% x = log( 1./Seff);
% y = log(sigmap.*L.^4/EY/Imom);
% 
% linFit = polyfit(x,y,1);
% y_linFit = linFit(1)*x+linFit(2);
% y_linFit = exp(y_linFit);
% 
% n = linFit(1)
% k = exp(linFit(2))
% 
% x = 1./Seff; 
% y = sigmap.*L.^4/EY/Imom; 
% 
% figure(6)
% hold on
% scatter(x,y,60,'filled','MarkerEdgeColor',color,'MarkerFaceColor',color)
% plot(x,y_linFit,'-x','linewidth',linewidth,'markersize',10,'color',color)
% ylabel('\it{\sigma_{xz,p} L^4/ E_Y I}')
% 
% RSSlin = (y-y_linFit).^2; 
% 
% figure(7)
% title('Normalized squared residual')
% hold on
% 
% plot(x,RSSlin/max(RSSlin),'-o','linewidth',linewidth,'markersize',10,'color',color)
% ylabel('RS')
% 
% sum(RSSlin)

end

