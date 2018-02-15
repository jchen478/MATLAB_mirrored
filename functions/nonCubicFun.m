function [] = nonCubicFun( filename, color, sidex, sidey, sidez)

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
Seff = pi.*kb./nseg.^4;
rp = rps.*nseg;
nL3 = nfib .* (2.*rp).^3/(sidex*sidey*sidez);

L = 2*rps.*nseg*a;
gamma = EY*Imom ./ (Seff.*L.^4*eta0);


%% Dimensionalize stresses
[sigmap,sigma, se_sigmap, se_sigma, N1, N2] = stressDim(sigmapstar, se_sigmapstar, N1star, N2star, nseg, rps, eta0, gamma, nL3); 

%% Calculate viscosity
etasp = sigmap./gamma;
se_etasp = se_sigmap./gamma;

etarel = sigma./gamma;
se_etarel = se_sigma./gamma;

%% Normalize intensity
I_normalized = I/I(1);
se_I_normalized = se_I/I(1);

%% Plotting parameters
markersize = 10; 
linewidth = 2.5; 

%% Plotting 
figure(1)
hold on;
errorbar(mu,sigma.*L.^4/EY/Imom,se_sigma.*L.^4/EY/Imom/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',color,'Linewidth',linewidth,'color',color);
ylabel('\it{\sigma_{xz} L^4/ E_Y I}')
title('Bulk stress')

%% intensity I
figure(2)
hold on;
ylabel('\it{}I')
errorbar(mu,I,se_I/2,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
title('Intensity')

%% normalized I
figure(3)
hold on;
ylabel('\it{}I_{normalized}')
errorbar(mu,I_normalized,se_I_normalized/2,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
title('Normalized intensity')

%% first normal difference
figure(4)
hold on;
ylabel('\it{N_1 L^4/ E_Y I}')
plot(mu,N1.*L.^4/EY/Imom,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
title('First normal stress')

%% second normal difference
figure(5)
hold on;
ylabel('\it{N_2 L^4/ E_Y I}')
plot(mu,N2.*L.^4/EY/Imom,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
title('Second normal stress')

%% Particle stress contribution
figure(6)
hold on;
ylabel('\it{\sigma_{xz,p} L^4/ E_Y I}')
errorbar(mu,sigmap.*L.^4/EY/Imom,se_sigmap.*L.^4/EY/Imom/2,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
title('Particle stress contribution')

%% Specific viscosity
figure(7)
hold on;
errorbar(mu,etasp,se_etasp/2,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
ylabel('\it{}\eta_{sp}')
title('Specific viscosity')

%% Relative viscosity
figure(8)
hold on;
errorbar(mu,etarel,se_etarel/2,'-.o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
ylabel('\it{}\eta_{rel}')
title('Relative viscosity')

end

