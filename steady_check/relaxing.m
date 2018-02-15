% function [] = relaxing()
%% Read files
clc; clear; close all;

%% Parameters
File = fopen('Parameters.in','r');
Parameters = textscan(File,'%f %*s',Inf,'Delimiter','\n');
fclose(File);
Parameters = cell2mat(Parameters);

%% Stress tensor
File = fopen('Stress_tensor.txt','r');
Stress = fscanf(File,'%f',[7 Inf])';
fclose(File);

%% Number of Contacts
File = fopen('Number_of_Contacts.txt','r');
NumberofContacts = fscanf(File,'%f',[5 Inf])';
fclose(File);

%% Number of Contacts
File = fopen('Eelastic.txt','r');
Eelastic = fscanf(File,'%f',[4 Inf])';
fclose(File);
% Eelastic(1,:) = []; 

%% Define parameters
nfib = Parameters(1);
nseg = Parameters(2);
rps  = Parameters(3);
kb  = Parameters(4);
dt = Parameters(10);
strain = Parameters(11);
side = Parameters(12); % cubic
contact_write = Parameters(15);
stress_write = Parameters(49);
stop_strain = Parameters(68);
elastic_write = Parameters(69);

sigmapstar = Stress(:,4); 
NC = NumberofContacts(:,4);
Energy = Eelastic(:,4); 
tstress = dt*stress_write:dt*stress_write:strain;
tcontact = dt*contact_write:dt*contact_write:strain;
telastic = dt*elastic_write:dt*elastic_write:strain;
telastic = telastic';
%% fiber dimensions
a = 1.60E-05;       % radius (m)
Imom = pi*a^4/4;    % area moment (m^4)
EY = 1e9;           % Young's modulus (Pa m^2)
eta0 = 1;           % fluid viscosity (Pa s)

%% calculate relevant parameters
[Seff,rp,nL3,L,gamma] = pCalc(nfib, nseg, rps, kb, side, a, EY, Imom, eta0);

%% Dimensionalize stress
sigmap = sigmapstar*pi.*nL3./(6*nseg.^3.*log(2*rps))*eta0.*gamma;

%% Particle stress contribution6
figure('units','normalized','outerposition',[0.2 0.1 0.6 0.9])
subplot(3,1,1)
hold on;
xlabel('\gamma')
ylabel('\it{\sigma_{xz,p} L^4/ E_Y I}')
scatter(tstress,sigmap.*L.^4/EY/Imom,10,'MarkerEdgeColor',rgb('MediumBlue'),'MarkerFaceColor',rgb('MediumBlue'));
line([stop_strain stop_strain], [min(sigmap.*L.^4/EY/Imom) max(sigmap.*L.^4/EY/Imom)],...
    'linewidth',2,'color',rgb('OrangeRed')); 
title('Particle stress contribution')

%% Number of contacts
subplot(3,1,2)
hold on;
xlabel('\gamma')
ylabel('<\it{N_C}\rm>')
scatter(tcontact,NC,10,'MarkerEdgeColor',rgb('MediumBlue'),'MarkerFaceColor',rgb('MediumBlue'));
line([stop_strain stop_strain], [min(NC) max(NC)],...
    'linewidth',2,'color',rgb('OrangeRed')); 
title('Number of contacts')
ylim([0 inf])

%% Elastic energy storage

subplot(3,1,3)
hold on
xlabel('\gamma')
ylabel('<\it{}E_{elastic}\rm{}>')
scatter(telastic,Energy,10,'MarkerEdgeColor',rgb('MediumBlue'),'MarkerFaceColor',rgb('MediumBlue'));
line([stop_strain stop_strain], [min(Energy) max(Energy)],...
    'linewidth',2,'color',rgb('OrangeRed')); 
title('Elastic energy storage')
ylim([0 inf])

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 24;
linewidth = 1.5;

%% Figure formatting
for i=1:3
    subplot(3,1,i)
    hold on
    box on
    xlabel('\gamma')
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    xlim([0 inf])
end


% spreadfigures();

% end