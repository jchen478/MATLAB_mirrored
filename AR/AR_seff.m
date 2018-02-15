close all;
clc;
clear;

%% parameters
a = 1.60E-05;
Imom = pi*a^4/2;
EY = 1e9; 
eta0 = 1; 

%% read files
filename={'run002.txt'};
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

eta = tauxz./gamma
se_eta = se_tauxz./gamma;

%% plotting
color = {rgb('Crimson'),rgb('MediumBlue'),rgb('DarkOliveGreen')};

nstiff = 2; 

rp = reshape(rp,length(rp)/nstiff,nstiff);
eta = reshape(eta,length(eta)/nstiff,nstiff);
se_eta = reshape(se_eta,length(se_eta)/nstiff,nstiff);
N1 = reshape(N1,length(N1)/nstiff,nstiff);

for i=1:nstiff

figure(1)
subplot(2,1,1)
hold on;
xlabel('\it{}r_p')
ylabel('\it{}\eta_{sp}')
errorbar(rp(:,i),eta(:,i),se_eta(:,i)/2,'-.o','MarkerSize',10,'MarkerEdgeColor',color{i},'Linewidth',2.5,'color',color{i});

subplot(2,1,2)
hold on;
xlabel('\it{}r_p')
ylabel('\it{}N_1 / E_Y')
plot(rp(:,i),N1(:,i)/EY,'-.o','MarkerSize',10,'MarkerEdgeColor',color{i},'Linewidth',2.5,'color',color{i});

end


%% figure formatting
tickx = 2;
name = {['\it{}\eta d\gamma/dt / E_Y = '],['\it{}rp_s \rm{}= 10.0']};
for i=1:2 
    
%     figure(i)
    subplot(2,1,i)
    hold on;
    box on;    
    set(gca,'fontsize',24,'linewidth',2,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
    legend(name{:}, 'location','northwest');
%     print(title{i},'-dpng');
%     savefig(title{i});
    
end
    print('case1','-dpng');
    savefig('case1');
    