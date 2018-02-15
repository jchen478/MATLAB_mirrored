clc;
clear;
close all;

casename={'theta04mu15'};

%% Simulation parameters
nfib = 1280;
nseg = 5; 
kb = 5; 
rps = 15; 
rp = 2*nseg*rps; 
Seff = pi*kb/nseg^4;
dt = 0.0001;
strain = 3000;
stress_write = 2000;
n_stress = strain/dt/stress_write;
time_stress = stress_write*dt:stress_write*dt:strain;
config_write = 5000;
time_config = 0:config_write*dt:strain;

%% Expansion parameters
Lx_init = 600;
fconc = 0.9986;
fexpd = 1.0007;
startconc = 106;
startexpd = 1100;
nconc = 494;
nexpd = 985;
eqconc = 1;
eqexpd = 1;

%% Dimensional parameters
a = 1.60E-05;
Imom = pi*a^4/2;
EY = 1e9; 
eta0 = 1; 

%% Dimensional variables
L = 2*rps*nseg*a;
gamma = EY*Imom / (Seff.*L^4*eta0); 

%% Stress
File = fopen('Stress_tensor_low_nL3.txt','r');
data = fscanf(File,'%f',[7 Inf])';
fclose(File);
tauxzstar(:,1) = data(:,4); 

File = fopen('Stress_tensor.txt','r');
data = fscanf(File,'%f',[7 Inf])';
fclose(File);
tauxzstar(:,2) = data(:,4); 

%% intensity
File = fopen('ISV_low_nL3.txt','r');
data = fscanf(File,'%f',[8 Inf])';
fclose(File);

tINS = data(:,1);
I(:,1) = data(:,2);
Sx(:,1) = data(:,3);
Sy(:,1) = data(:,4);
Sz(:,1) = data(:,5);
Vx(:,1) = data(:,6);
Vy(:,1) = data(:,7);
Vz(:,1) = data(:,8);

File = fopen('ISV.txt','r');
data = fscanf(File,'%f',[8 Inf])';
fclose(File);

I(:,2) = data(:,2);
Sx(:,2) = data(:,3);
Sy(:,2) = data(:,4);
Sz(:,2) = data(:,5);
Vx(:,2) = data(:,6);
Vy(:,2) = data(:,7);
Vz(:,2) = data(:,8);

%% Box size
Lx(:,2) = zeros(n_stress,1); 
for i=1:n_stress
    if time_stress(i) < startconc
        Lx(i,2) = Lx_init;
    elseif time_stress(i) >= startconc && time_stress(i) < startconc+(nconc-1)*eqconc
        power = ceil((time_stress(i) - startconc) / eqconc);
        Lx(i,2) = Lx_init*fconc^power;
    elseif time_stress(i) >= startconc+(nconc-1)*eqconc && time_stress(i) < startexpd
        Lx(i,2) = Lx_init*fconc^nconc;
    elseif time_stress(i) >= startexpd && time_stress(i) < startexpd+nexpd*eqexpd
        power = ceil((time_stress(i) - startexpd) / eqexpd);
        Lx(i,2) = Lx_init*fconc^nconc*fexpd^power;
    elseif time_stress(i) >= startexpd+nexpd*eqexpd
        Lx(i,2) = Lx_init*fconc^nconc*fexpd^nexpd;
    end
end
Lx(:,1) = Lx_init*fconc^nconc*fexpd^nexpd*ones(n_stress,1);
Ly = Lx; 
Lz = Lx; 

%% concentration calculation
% low_nL3
nL3(:,1) = nfib./(Lx(:,1).*Ly(:,1).*Lz(:,1))*(2*rps*nseg)^3; 

% expansionWithBox
nL3(:,2) = nfib./(Lx(:,2).*Ly(:,2).*Lz(:,2))*(2*rps*nseg)^3; 

%% dimensionalize
tauxz = tauxzstar*pi.*nL3./(6*nseg.^3.*log(2*rps))*eta0.*gamma;
sigxz = tauxz + eta0*gamma;
eta = tauxz/gamma;
etasp = eta/eta0; 

%% plotting
title={'nL3','eta','stress','intensity'};
colorname={rgb('MediumBlue'),rgb('Crimson')};

for i=1:2
    
    figure(1)
    hold on
    plot(time_stress,nL3(:,i),'linewidth',2,'color',colorname{i})
    ylabel('\it{}nL3')
    xlim([0 3000])
    
    figure(2)
    hold on
    plot(time_stress,etasp(:,i),'linewidth',2,'color',colorname{i})
    ylabel('\it{}\eta_{sp}')
    xlim([0 3000])
    
    figure(3)
    hold on
    plot(time_stress,sigxz(:,i)*L^4/EY/Imom,'linewidth',2,'color',colorname{i})
    ylabel('\it{}\sigma_{xz}L^4 / E_Y I')
    xlim([0 3000])
    
    figure(4)
    hold on
    plot(tINS,I(:,i),'linewidth',2,'color',colorname{i});
    ylabel('\it{}I')
    xlim([0 915])
    
    figure(5)
    subplot(2,1,1)
    hold on
    plot(time_stress,sigxz(:,i)*L^4/EY/Imom,'linewidth',2,'color',colorname{i})
    set(gca,'yscale','log')
    ylabel('\it{}\sigma_{xz}L^4 / E_Y I')
    xlim([0 3000])
    
    subplot(2,1,2)
    hold on
    plot(tINS,I(:,i),'linewidth',2,'color',colorname{i});
    ylabel('\it{}I')
    xlim([0 915])
    
end


%% formatting
tickx = 1.5;
viewPoint = [-12 16];
name={'Never concentrated','Redispersion'};
% for i=1:length(title)   
%     figure(i)
%     hold on;
%     box on;
%     xlabel('\gamma')
%     set(gca,'fontsize',40,'linewidth',2,'fontname','Times New Roman')
%     set(gca,'YMinorTick','on','XMinorTick','on','linewidth',2)
%     set(gca,'ticklength',tickx*get(gca,'ticklength'))
%     set(gcf, 'color','white')
%     legend(name{:},'location','best');
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
%     print(title{i},'-dpng');
%     savefig(title{i});
% end

for i=1:2
    figure(5)
    subplot(2,1,i)
    hold on;
    box on;
    xlabel('\gamma')
    set(gca,'fontsize',24,'linewidth',2,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on','linewidth',2)
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    legend(name{:},'location','best');
end
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])

print(casename{1},'-dpng');
savefig(casename{1});