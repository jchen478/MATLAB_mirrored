clc; 
close all; 
clear;

%% parameters

% simulation
nfib = 1280;
kb = 1; 
nseg = 5;
rps = 15;
strain = 2000;
dt = 0.0001; 
Seff = pi*kb/nseg^4; 
over_Seff = 1/Seff;
Lx_init = 600;
fconc = 0.9987;
fexpd = 1.0006;
nconc = 500;
nexpd = 1030;
eqconc = 1;
eqexpd = 0.8;
startconc = 50;
startexpd = 850;

% data write
config_write = 5000;
stress_write = 2000;
time_stress = stress_write*dt:stress_write*dt:strain;
time_config = 0:config_write*dt:strain;

% time smooth
window = 50;

%% open file
% Stress
File = fopen('Stress_tensor.txt','r');
data = fscanf(File,'%f',[7 Inf])';
fclose(File);
sigxz = data(:,4); 
n_stress = length(sigxz); 

% Entropy
File = fopen('Shannon_parameters.txt','r');
textscan(File, '%[^\n]',1);
dx = fscanf(File,'%f',1);
fclose(File);


File = fopen('Shannon_Entropy.txt','r');
data = fscanf(File,'%f',[2 Inf])';
fclose(File);
entropy = data(2:end,2); 
n_entropy = length(entropy); 

% Box size
% File = fopen('Box_Size.txt','r');
% Box = fscanf(File,'%f',[4 Inf])';
% fclose(File);
% Lx = Box(:,2);
% Ly = Box(:,3);
% Lz = Box(:,4);

Lx = zeros(n_stress,1); 

for i=1:n_stress
    if time_stress(i) < startconc
        Lx(i) = Lx_init;
    elseif time_stress(i) >= startconc && time_stress(i) < startconc+(nconc-1)*eqconc
        power = ceil((time_stress(i) - startconc) / eqconc);
        Lx(i) = Lx_init*fconc^power;
    elseif time_stress(i) >= startconc+(nconc-1)*eqconc && time_stress(i) < startexpd
        Lx(i) = Lx_init*fconc^nconc;
    elseif time_stress(i) >= startexpd && time_stress(i) < startexpd+nexpd*eqexpd
        power = ceil((time_stress(i) - startexpd) / eqexpd);
        Lx(i) = Lx_init*fconc^nconc*fexpd^power;
    elseif time_stress(i) >= startexpd+nexpd*eqexpd
        Lx(i) = Lx_init*fconc^nconc*fexpd^nexpd;
    end
end
Ly = Lx;
Lz = Lx; 

%% concentration calculation
nL3 = nfib./(Lx.*Ly.*Lz)*(2*rps*nseg)^3; 

%% dimensionalize
sigxz_total = sigxz*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff + over_Seff;
sigxz_total_ts = zeros(n_stress,1); 

%% time-smooth
for j=1:n_stress-window    
% for j=1:5
    sum = 0;
    for k=1:window        
        sum = sum + sigxz_total(j+k);     
    end
    sum = sum/window;
    sigxz_total_ts(j) = sum;
end
for j=n_stress-window+1:n_stress
    sigxz_total_ts(j) = mean(sigxz_total(j:end));
end

%% plotting
figure('color','white')
% stress
% subplot(3,1,1)
% hold on;
% box on
% h = plot(time_stress,sigxz_total,'color',rgb('MediumBlue'),'linewidth',3);
% hfriction.LineStyle = '-.';
% h.LineStyle = '-';
% ylabel('\it{\sigma_{xz} L^4/ E_Y I}')
% xlabel('\gamma')
% title('Stress')
% set(gca,'fontsize',18,'linewidth',2,'fontname','Times New Roman')
% set(gca,'YMinorTick','on','XMinorTick','on')

% time-smoothed stress
subplot(3,1,1)
hold on
h = plot(time_stress,sigxz_total_ts,'color',rgb('MediumBlue'),'linewidth',3);
hfriction.LineStyle = '-.';
h.LineStyle = '-';
ylabel('$\bar{\it{\sigma}}_{xz}\it{L^4/ E_Y I}$','interpreter','latex')
xlabel('\gamma')
title(['Time-smoothed Stress (Window = ',num2str(window),')'])
set(gca,'fontsize',18,'linewidth',2,'fontname','Times New Roman')
set(gca,'YMinorTick','on','XMinorTick','on')

% entropy
subplot(3,1,2)
hold on
h = plot(time_config,entropy,'color',rgb('MediumBlue'),'linewidth',3);
hfriction.LineStyle = '-.';
h.LineStyle = '-';
ylabel('$\it{H}$','interpreter','latex')
xlabel('\gamma')
title(['Shannon Entropy (dx = ',num2str(dx,'%.2f'),')'])
set(gca,'fontsize',18,'linewidth',2,'fontname','Times New Roman')
set(gca,'YMinorTick','on','XMinorTick','on')

% concentration
subplot(3,1,3)
hold on;
box on
h = plot(time_stress,nL3,'color',rgb('MediumBlue'),'linewidth',3);
hfriction.LineStyle = '-.';
h.LineStyle = '-';
ylabel('\it{nL^3}')
xlabel('\gamma')
% title('Concentration')
set(gca,'fontsize',18,'linewidth',2,'fontname','Times New Roman')
set(gca,'YMinorTick','on','XMinorTick','on')