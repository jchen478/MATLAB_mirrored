clc; 
close all;

% figure(4)
% rgb chart

%% Parameters
nseg = 5;
rps = 15;
EYI = 1.128e-10; % N m2
L = 3.45e-3; % m
color = {rgb('Maroon'),rgb('DeepPink'),rgb('DarkOrange'),rgb('DarkGreen'),rgb('MediumBlue'),rgb('DarkMagenta')};

%% 320
File = fopen('320.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f ';
sizeA = [10 Inf];
stress = flipud(fscanf(File,formatSpec,sizeA)');
fclose(File);

% nfib nL3 kb volfrac sigxz stdev sem contact 1 2

sigxz_p = stress(:,5);
kb = stress(:,3); 
volfrac = stress(:,4);
nL3 = stress(:,2); 
SEM = stress(:,7);

Seff = pi*kb/nseg^4; 
over_Seff = 1./Seff;

sigxz_total = sigxz_p*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff + over_Seff;
SEM = SEM*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff/2; 

j = 0;
S = cell(2*length(kb)/3,1);
sig0 = zeros(length(kb)/3,1); 
sig0_phi = zeros(length(kb)/3,1); 

for i=1:3:length(kb)
    
    S{2*j+1} = ['\Phi = ', num2str(volfrac(i),'%.3f')];
    j = j+1; 
    
end

figure(1)
box on;
hold on;
set(gcf,'Color','white');
title('nfib = 320')

regressX = linspace(0,max(over_Seff),400);
j = 1; 

for i=1:3:length(kb)
    
  
    s(j) = scatter(over_Seff(i:i+2),sigxz_total(i:i+2),100,'MarkerEdgeColor',color{j},'Linewidth',2.0); 
    over_Seff(i:i+2);
    sigxz_total(i:i+2);
    fit = polyfit(over_Seff(i:i+2),sigxz_total(i:i+2),1);
    sigxz_regress = fit(1)*regressX+fit(2); 
    h(j) = plot(regressX,sigxz_regress,'color',color{j},'linewidth',2);
    sig0(j) = fit(2); 
    sig0_phi(j) = volfrac(i); 
    S{2*j} = ['\sigma_{xz} = ', num2str(fit(1),3), ' / S_{eff} + ', num2str(fit(2),'%.2f')];
    j = j+1; 
    
end

s(1).Marker = 'd';
s(2).Marker = 'o';
s(3).Marker = '^';
s(4).Marker = 'v';


h(1).LineStyle = '-';
h(2).LineStyle = '--';
h(3).LineStyle = ':';
h(4).LineStyle = '-.';

legend({'\it{}\Phi\rm{} = 0.013',...
' \it{}\sigma_{xz}^*\rm{} = (5.09\pm0.30) / \it{}S_{eff}\rm{} + (156\pm8)',...
'\it{}\Phi\rm{} = 0.007',...
' \it{}\sigma_{xz}^*\rm{} = (1.24\pm0.06) / \it{}S_{eff}\rm{} + (12\pm2)',...
'\it{}\Phi\rm{} = 0.004',...
' \it{}\sigma_{xz}^*\rm{} = (1.13\pm0.01) / \it{}S_{eff}\rm{} + (0.8\pm0.2)',...
'\it{}\Phi\rm{} = 0.003',...
' \it{}\sigma_{xz}^*\rm{} = (1.070\pm0.003) / \it{}S_{eff}\rm{} - (0.026\pm0.07)'},...
'Location','northwest');

xlabel('1/\it{}S_{eff}')
ylabel('\it{}\sigma_{xz} L^4 / E_Y I')
set(gca,'fontsize',24,'fontName','Times New Roman','linewidth',2)
ylim([0 800])

sig0err = [5.125 0.313 0.0636 0.0097];

figure(2)
hold on;    
E1 = errorbar(sig0_phi,sig0,sig0err/2,'v','MarkerEdgeColor',color{1},'Linewidth',2.5,'color',color{1});
sig0_phi_log = log(sig0_phi);
sig0_log = log(sig0);
corr = polyfit(sig0_phi_log,sig0_log,1);
sigxz_regressX = linspace(sig0_phi(1), sig0_phi(4),200);
alpha = exp(corr(2));
beta = corr(1);
sigxz_regress = alpha*sigxz_regressX.^beta; 
h2(1) = plot(sigxz_regressX,sigxz_regress,'linewidth',3,'color',color{1});

%% 640
File = fopen('640.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f ';
sizeA = [10 Inf];
stress = flipud(fscanf(File,formatSpec,sizeA)');
fclose(File);

% nfib nL3 kb volfrac sigxz stdev sem contact 1 2

sigxz_p = stress(:,5);
kb = stress(:,3); 
volfrac = stress(:,4);
nL3 = stress(:,2); 
SEM = stress(:,7);

Seff = pi*kb/nseg^4; 
over_Seff = 1./Seff;

sigxz_total = sigxz_p*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff + over_Seff;
SEM = SEM*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff/2; 

j = 0;
S = cell(2*length(kb)/3,1);
sig0 = zeros(length(kb)/3,1); 
sig0_phi = zeros(length(kb)/3,1); 

for i=1:3:length(kb)
    
    S{2*j+1} = ['\Phi = ', num2str(volfrac(i),'%.3f')];
    j = j+1; 
    
end

figure(3)
box on;
hold on;
set(gcf,'Color','white');
title('nfib = 640')

regressX = linspace(0,max(over_Seff),400);
j = 1; 

for i=1:3:length(kb)
    
  
    s(j) = scatter(over_Seff(i:i+2),sigxz_total(i:i+2),100,'MarkerEdgeColor',color{j},'Linewidth',2.0); 
    over_Seff(i:i+2);
    sigxz_total(i:i+2);
    fit = polyfit(over_Seff(i:i+2),sigxz_total(i:i+2),1);
    sigxz_regress = fit(1)*regressX+fit(2); 
    h(j) = plot(regressX,sigxz_regress,'color',color{j},'linewidth',2);
    sig0(j) = fit(2); 
    sig0_phi(j) = volfrac(i); 
    S{2*j} = ['\sigma_{xz} = ', num2str(fit(1),3), ' / S_{eff} + ', num2str(fit(2),'%.2f')];
    j = j+1; 
    
end

s(1).Marker = 'd';
s(2).Marker = 'o';
s(3).Marker = '^';
s(4).Marker = 'v';


h(1).LineStyle = '-';
h(2).LineStyle = '--';
h(3).LineStyle = ':';
h(4).LineStyle = '-.';

legend({'\it{}\Phi\rm{} = 0.013',...
' \it{}\sigma_{xz}^*\rm{} = (5.09\pm0.30) / \it{}S_{eff}\rm{} + (156\pm8)',...
'\it{}\Phi\rm{} = 0.007',...
' \it{}\sigma_{xz}^*\rm{} = (1.24\pm0.06) / \it{}S_{eff}\rm{} + (12\pm2)',...
'\it{}\Phi\rm{} = 0.004',...
' \it{}\sigma_{xz}^*\rm{} = (1.13\pm0.01) / \it{}S_{eff}\rm{} + (0.8\pm0.2)',...
'\it{}\Phi\rm{} = 0.003',...
' \it{}\sigma_{xz}^*\rm{} = (1.070\pm0.003) / \it{}S_{eff}\rm{} - (0.026\pm0.07)'},...
'Location','northwest');

xlabel('1/\it{}S_{eff}')
ylabel('\it{}\sigma_{xz} L^4 / E_Y I')
set(gca,'fontsize',24,'fontName','Times New Roman','linewidth',2)
ylim([0 800])

sig0err = [17.09 1.436 0.199 0.073];

figure(2) 
hold on;    
E2 = errorbar(sig0_phi,sig0,sig0err/2,'v','MarkerEdgeColor',color{2},'Linewidth',2.5,'color',color{2});
sig0_phi_log = log(sig0_phi);
sig0_log = log(sig0);
corr = polyfit(sig0_phi_log,sig0_log,1);
sigxz_regressX = linspace(sig0_phi(1), sig0_phi(4),200);
alpha = exp(corr(2));
beta = corr(1);
sigxz_regress = alpha*sigxz_regressX.^beta; 
h2(2) = plot(sigxz_regressX,sigxz_regress,'linewidth',3,'color',color{2});

%% 1280
File = fopen('1280.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f ';
sizeA = [10 Inf];
stress = flipud(fscanf(File,formatSpec,sizeA)');
fclose(File);

% nfib nL3 kb volfrac sigxz stdev sem contact 1 2

sigxz_p = stress(:,5);
kb = stress(:,3); 
volfrac = stress(:,4);
nL3 = stress(:,2); 
SEM = stress(:,7);

Seff = pi*kb/nseg^4; 
over_Seff = 1./Seff;

sigxz_total = sigxz_p*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff + over_Seff;
SEM = SEM*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff/2; 

j = 0;
S = cell(2*length(kb)/3,1);
sig0 = zeros(length(kb)/3,1); 
sig0_phi = zeros(length(kb)/3,1); 

for i=1:3:length(kb)
    
    S{2*j+1} = ['\Phi = ', num2str(volfrac(i),'%.3f')];
    j = j+1; 
    
end

regressX = linspace(0,max(over_Seff),400);
j = 1; 

for i=1:3:length(kb)
    
    fit = polyfit(over_Seff(i:i+2),sigxz_total(i:i+2),1);
    sigxz_regress = fit(1)*regressX+fit(2); 
   
    sig0(j) = fit(2); 
    sig0_phi(j) = volfrac(i); 
    sigxz_total(i:i+2);
    S{2*j} = ['\sigma_{xz} = ', num2str(fit(1),3), ' / S_{eff} + ', num2str(fit(2),3)];
    j = j+1; 
    
end

sig0err = [21.73 1.05 0.174 0.034];

figure(2)
hold on;    
E3 = errorbar(sig0_phi,sig0,sig0err/2,'v','MarkerEdgeColor',color{3},'Linewidth',2.5,'color',color{3});
sig0_phi_log = log(sig0_phi);
sig0_log = log(sig0);
corr = polyfit(sig0_phi_log,sig0_log,1);
sigxz_regressX = linspace(sig0_phi(1), sig0_phi(4),200);
alpha = exp(corr(2));
beta = corr(1);
sigxz_regress = alpha*sigxz_regressX.^beta; 
h2(3) = plot(sigxz_regressX,sigxz_regress,'linewidth',3,'color',color{3});

%% 2560
File = fopen('2560.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f ';
sizeA = [10 Inf];
stress = flipud(fscanf(File,formatSpec,sizeA)');
fclose(File);

% nfib nL3 kb volfrac sigxz stdev sem contact 1 2

sigxz_p = stress(:,5);
kb = stress(:,3); 
volfrac = stress(:,4);
nL3 = stress(:,2); 
SEM = stress(:,7);

Seff = pi*kb/nseg^4; 
over_Seff = 1./Seff;

sigxz_total = sigxz_p*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff + over_Seff;
SEM = SEM*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff/2; 

j = 0;
S = cell(2*length(kb)/3,1);
sig0 = zeros(length(kb)/3,1); 
sig0_phi = zeros(length(kb)/3,1); 

for i=1:3:length(kb)
    
    S{2*j+1} = ['\Phi = ', num2str(volfrac(i),'%.3f')];
    j = j+1; 
    
end

regressX = linspace(0,max(over_Seff),400);
j = 1; 

for i=1:3:length(kb)
    
    fit = polyfit(over_Seff(i:i+2),sigxz_total(i:i+2),1);
    sigxz_regress = fit(1)*regressX+fit(2); 
    sig0(j) = fit(2); 
    sig0_phi(j) = volfrac(i); 
    sigxz_total(i:i+2);
    S{2*j} = ['\sigma_{xz} = ', num2str(fit(1),3), ' / S_{eff} + ', num2str(fit(2),3)];
    j = j+1; 
    
end

sig0err = [39.96 0.146 0.504];

figure(2)
hold on;
E4 = errorbar(sig0_phi(1:3),sig0(1:3),sig0err/2,'v','MarkerEdgeColor',color{4},'Linewidth',2.5,'color',color{4});
sig0_phi_log = log(sig0_phi(1:3));
sig0_log = log(sig0(1:3));
corr = polyfit(sig0_phi_log,sig0_log,1);
sigxz_regressX = linspace(sig0_phi(1), sig0_phi(3),150);
alpha = exp(corr(2));
beta = corr(1);
sigxz_regress = alpha*sigxz_regressX.^beta; 
h2(4) = plot(sigxz_regressX,sigxz_regress,'linewidth',3,'color',color{4});


%% 6400 
File = fopen('6400.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f ';
sizeA = [10 Inf];
stress = flipud(fscanf(File,formatSpec,sizeA)');
fclose(File);

% nfib nL3 kb volfrac sigxz stdev sem contact 1 2

sigxz_p = stress(:,5);
kb = stress(:,3); 
volfrac = stress(:,4);
nL3 = stress(:,2); 
SEM = stress(:,7);

Seff = pi*kb/nseg^4; 
over_Seff = 1./Seff;

sigxz_total = sigxz_p*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff + over_Seff;
SEM = SEM*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff/2; 

j = 0;
S = cell(2*length(kb)/3,1);
sig0 = zeros(length(kb)/3,1); 
sig0_phi = zeros(length(kb)/3,1); 

for i=1:3:length(kb)
    
    S{2*j+1} = ['\Phi = ', num2str(volfrac(i),'%.3f')];
    j = j+1; 
    
end

regressX = linspace(0,max(over_Seff),400);
j = 1; 

for i=1:3:length(kb)
    
    fit = polyfit(over_Seff(i:i+2),sigxz_total(i:i+2),1);
    sigxz_regress = fit(1)*regressX+fit(2); 
    sig0(j) = fit(2); 
    sig0_phi(j) = volfrac(i); 
    sigxz_total(i:i+2);
    S{2*j} = ['\sigma_{xz} = ', num2str(fit(1),3), ' / S_{eff} + ', num2str(fit(2),3)];
    j = j+1; 
    
end

sig0err = [51.71 2.34 1.34 0.14];    

figure(2)
hold on;
E5 = errorbar(sig0_phi,sig0,sig0err/2,'v','MarkerEdgeColor',color{5},'Linewidth',2.5,'color',color{5});
sig0_phi_log = log(sig0_phi);
sig0_log = log(sig0);
corr = polyfit(sig0_phi_log,sig0_log,1);
sigxz_regressX = linspace(sig0_phi(1), sig0_phi(4),200);
alpha = exp(corr(2));
beta = corr(1);
sigxz_regress = alpha*sigxz_regressX.^beta; 
h2(5) = plot(sigxz_regressX,sigxz_regress,'linewidth',3,'color',color{5});

%% 12800
File = fopen('12800.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f ';
sizeA = [10 Inf];
stress = flipud(fscanf(File,formatSpec,sizeA)');
fclose(File);

% nfib nL3 kb volfrac sigxz stdev sem contact 1 2

sigxz_p = stress(:,5);
kb = stress(:,3); 
volfrac = stress(:,4);
nL3 = stress(:,2); 
SEM = stress(:,7);

Seff = pi*kb/nseg^4; 
over_Seff = 1./Seff;

sigxz_total = sigxz_p*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff + over_Seff;
SEM = SEM*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff/2; 

j = 0;
S = cell(2*length(kb)/3,1);
sig0 = zeros(length(kb)/3,1); 
sig0_phi = zeros(length(kb)/3,1); 

for i=1:3:length(kb)
    
    S{2*j+1} = ['\Phi = ', num2str(volfrac(i),'%.3f')];
    j = j+1; 
    
end

regressX = linspace(0,max(over_Seff),400);
j = 1; 

for i=1:3:length(kb)
    
  
%     s(j) = scatter(over_Seff(i:i+2),sigxz_total(i:i+2),100,'MarkerEdgeColor',color{j},'Linewidth',2.0); 
    over_Seff(i:i+2);
    sigxz_total(i:i+2);
    fit = polyfit(over_Seff(i:i+2),sigxz_total(i:i+2),1);
    sigxz_regress = fit(1)*regressX+fit(2); 
    sig0(j) = fit(2); 
    sig0_phi(j) = volfrac(i); 
    S{2*j} = ['\sigma_{xz} = ', num2str(fit(1),3), ' / S_{eff} + ', num2str(fit(2),'%.2f')];
    j = j+1; 
    
end

sig0err = [15.92 3.15 0.303];

figure(2)
hold on;
E6 = errorbar(sig0_phi(1:3),sig0(1:3),sig0err/2,'v','MarkerEdgeColor',color{6},'Linewidth',2.5,'color',color{6});
sig0_phi_log = log(sig0_phi(1:3));
sig0_log = log(sig0(1:3));
corr = polyfit(sig0_phi_log,sig0_log,1);
sigxz_regressX = linspace(sig0_phi(1), sig0_phi(3),150);
alpha = exp(corr(2));
beta = corr(1);
sigxz_regress = alpha*sigxz_regressX.^beta; 
h2(6) = plot(sigxz_regressX,sigxz_regress,'linewidth',3,'color',color{6});


figure(2)
hold on;
box on; 
set(gca,'xscale','log','yscale','log')
xlabel('\it{\Phi}')
ylabel('\it{\sigma_{0} L^4 / E_Y I}')
set(gca,'fontsize',29,'fontName','Times New Roman','linewidth',2)
set(gcf,'Color','white');

set(E1,'MarkerSize',10,'Marker','o');
set(E2,'MarkerSize',10,'Marker','X');
set(E3,'MarkerSize',10,'Marker','d');
set(E4,'MarkerSize',10,'Marker','+');
set(E5,'MarkerSize',10,'Marker','s');
set(E6,'MarkerSize',10,'Marker','*');
h2(1).LineStyle = '-';
h2(2).LineStyle = '--';
h2(3).LineStyle = ':';
h2(4).LineStyle = '-.';
h2(5).LineStyle = ':';
h2(6).LineStyle = '-.';
legend({'\it{}N_{fib}\rm{} = 320',...
'\it{}\sigma_0\rm{} ~ \it{}\Phi\rm{}^{4.8\pm0.2}',...
'\it{}N_{fib}\rm{} = 640',...
'\it{}\sigma_0\rm{} ~ \it{}\Phi\rm{}^{4.1\pm0.1}',...
'\it{}N_{fib}\rm{} = 1280',...
'\it{}\sigma_0\rm{} ~ \it{}\Phi\rm{}^{5.0\pm0.1}',...
'\it{}N_{fib}\rm{} = 2560',...
'\it{}\sigma_0\rm{} ~ \it{}\Phi\rm{}^{4.5\pm0.1}',...
'\it{}N_{fib}\rm{} = 6400',...
'\it{}\sigma_0\rm{} ~ \it{}\Phi\rm{}^{4.2\pm0.1}',...
'\it{}N_{fib}\rm{} = 12800',...
'\it{}\sigma_0\rm{} ~ \it{}\Phi\rm{}^{4.8\pm0.2}'},...
'Location','northwest')

xlim([0.001 0.02])
