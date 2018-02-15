clc
clear
close all

Theta = [0 0.1 0.2 0.3 0.4 0.5 0.6];
Mu = [0 1 3 5 10 15 20];

%% Simulation parameters
nfib = 1280;
nseg = 5;
rps = 15; 
kb = 5; 
rp = nseg*rps; 
Seff = pi*kb/nseg^4; 
side = 598.31; 

%% Dimensional parameters
a = 1.60E-05;
Imom = pi*a^4/2;
EY = 1e9; 
eta0 = 1; 

L = 2*rps.*nseg*a;
gamma = EY*Imom ./ (Seff.*L.^4*eta0); 
nL3 = nfib*(2*rp/side)^3; 

%% file spec
filename={'low_nL3.txt'}; 
name={'Never concentrated'}; 

ncases = length(filename); 
ndata = length(Theta)*length(Mu); 
sizeA = [18 ndata];

%% allocate arrays
theta = zeros(ndata,ncases); 
mu = zeros(ndata,ncases); 
nfib = zeros(ndata,ncases); 
% nseg = zeros(ndata,ncases); 
% rps = zeros(ndata,ncases); 
kb = zeros(ndata,ncases); 
side = zeros(ndata,ncases); 
I = zeros(ndata,ncases); 
se_I = zeros(ndata,ncases);
I_normalized = zeros(ndata,ncases); 
se_I_normalized = zeros(ndata,ncases);
Sx = zeros(ndata,ncases); 
se_Sx = zeros(ndata,ncases);
Sy = zeros(ndata,ncases); 
se_Sy = zeros(ndata,ncases);
Sz = zeros(ndata,ncases); 
se_Sz = zeros(ndata,ncases);
Vx = zeros(ndata,ncases); 
se_Vx = zeros(ndata,ncases);
Vy = zeros(ndata,ncases); 
se_Vy = zeros(ndata,ncases);
Vz = zeros(ndata,ncases); 
se_Vz = zeros(ndata,ncases);
tauxzstar = zeros(ndata,ncases); 
se_tauxzstar = zeros(ndata,ncases);
N1star = zeros(ndata,ncases); 
N2star = zeros(ndata,ncases);

%% read files
for i=1:ncases
    File = fopen(filename{i},'r');
    data = fscanf(File,'%f',sizeA)';
    fclose(File);
    theta(:,i) = data(:,1);
    mu(:,i) = data(:,2);
    nfib(:,i) = data(:,3);
%     nseg(:,i) = data(:,4);
%     rps(:,i) = data(:,5);
    kb(:,i) = data(:,6);
    side(:,i) = data(:,7);
    I(:,i) = data(:,8);
    Sx(:,i) = data(:,9);
    Sy(:,i) = data(:,10);
    Sz(:,i) = data(:,11);
    Vx(:,i) = data(:,12);
    Vy(:,i) = data(:,13);
    Vz(:,i) = data(:,14);
    tauxzstar(:,i) = data(:,15); 
    se_tauxzstar(:,i) = data(:,16); 
    N1star(:,i) = data(:,17); 
    N2star(:,i) = data(:,18); 
end

%% dimensionalize
tauxz = tauxzstar*pi*nL3/(6*nseg^3*log(2*rps))*eta0*gamma;
se_tauxz = se_tauxzstar*pi*nL3/(6*nseg^3*log(2*rps))*eta0*gamma;
N1 =  N1star*pi*nL3/(6*nseg^3*log(2*rps))*eta0*gamma;
N2 =  N2star*pi*nL3/(6*nseg^3*log(2*rps))*eta0*gamma;

sigxz = tauxz + eta0*gamma;
se_sigxz = se_tauxz;

eta = tauxz/gamma;
se_eta = se_tauxz/gamma;

%% plot
title={'intensity','eta','stress'};
subname={'Sx','Sy','Sz','Vx','Vy','Vz'};

[meshX,meshY]=meshgrid(Theta,Mu);

figure(1)
surf(meshX,meshY,reshape(I,length(Mu),length(Theta)),'linestyle','none')
zlabel('\it{}I')

figure(2)
surf(meshX,meshY,reshape(eta,length(Mu),length(Theta)),'linestyle','none')
zlabel('\eta_{sp}')

figure(3)
surf(meshX,meshY,reshape(sigxz*L^4/EY/Imom,length(Mu),length(Theta)),'linestyle','none')
zlabel('\sigma_{xz} L^4 / E_Y I')

figure(4)
subplot(2,3,1)
hold on
surf(meshX,meshY,reshape(Sx,length(Mu),length(Theta)),'linestyle','none')
zlabel('\it{}Sx')

subplot(2,3,2)
hold on
surf(meshX,meshY,reshape(Sy,length(Mu),length(Theta)),'linestyle','none')
zlabel('\it{}Sy')

subplot(2,3,3)
hold on
surf(meshX,meshY,reshape(Sz,length(Mu),length(Theta)),'linestyle','none')
zlabel('\it{}Sz')

subplot(2,3,4)
hold on
surf(meshX,meshY,reshape(Vx,length(Mu),length(Theta)),'linestyle','none')
zlabel('\it{}Vx')

subplot(2,3,5)
hold on
surf(meshX,meshY,reshape(Vy,length(Mu),length(Theta)),'linestyle','none')
zlabel('\it{}Vy')

subplot(2,3,6)
hold on
surf(meshX,meshY,reshape(Vz,length(Mu),length(Theta)),'linestyle','none')
zlabel('\it{}Vz')

%% figure format
tickx = 1.5;
viewPoint = [-12 16];
for i=1:length(title)   
    figure(i)
    hold on;
    box on;
    xlabel('\theta_{eq}')
    ylabel('\mu_{stat}')
    set(gca,'fontsize',40,'linewidth',2,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    legend(name{:},'location','best');
    xlim([0 0.6])
    ylim([0 20])
    colorbar('eastoutside')
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
    print(title{i},'-dpng');
    savefig(title{i});
end
figure(4)
for i=1:length(subname)   
    subplot(2,3,i);
    hold on;
    box on;
    xlabel('\theta_{eq}')
    ylabel('\mu_{stat}')
    set(gca,'fontsize',22,'linewidth',2,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    legend(name{:},'location','best');
    xlim([0 0.6])
    ylim([0 20])
    colorbar('eastoutside')
    view(viewPoint);
end
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
print('lengthscale','-dpng');
savefig('lengthscale');