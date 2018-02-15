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

%% Intensity
File = fopen('ISV.txt','r');
ISV = fscanf(File,'%f',[8 Inf])';
fclose(File);

%% Number of Contacts
File = fopen('Number_of_Contacts.txt','r');
NC = fscanf(File,'%f',[5 Inf])';
fclose(File);

%% Elastic storage
File = fopen('Eelastic.txt','r');
elastic = fscanf(File,'%f',[4 Inf])';
fclose(File);
elastic(1,:) = []; 

%% Define parameters
nfib = Parameters(1); 
nseg = Parameters(2); 
rps  = Parameters(3); 
kb  = Parameters(4); 
dt = Parameters(10);
strain = Parameters(11); 
side = Parameters(12); 
config_write = Parameters(14); 
contact_write = Parameters(15); 
stress_write = Parameters(49); 
strain_stop = Parameters(68); 
elastic_write = Parameters(69);

%% Calculate parameters
L = 2*rps*nseg; 
nL3 = nfib*L^3/side^3;
Seff = kb*pi/nseg^4;
over_Seff = 1/Seff; 

%% Time 
tStress = dt*stress_write:dt*stress_write:strain;
tCon = dt*contact_write:dt*contact_write:strain;
tISV = 0:dt*config_write:strain;
telastic =  dt*elastic_write:dt*elastic_write:strain;

tStressStop = tStress-strain_stop;
tConStop = tCon-strain_stop;
tISVStop = tISV-strain_stop;
telasticStop = telastic-strain_stop; 

%% Delete unwanted data
% NC(length(tCon)+1:end,:) = [];

%% Stress processing
sigxzStar = Stress(:,4)*pi*nL3/(6*nseg^3*log(2*rps))*over_Seff + over_Seff;

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 24;
linewidth = 1.5;

%% Plotting
figure()
subplot(2,2,1)
scatter(tStressStop,sigxzStar,...
        markersize,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'))
ylabel('\it{\sigma_{xz} L^4/ E_Y I}')

subplot(2,2,2)
scatter(tConStop,NC(:,4), ...
        markersize,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'))
ylabel('N_C / N_{fib}')

subplot(2,2,3)
scatter(tISVStop,ISV(:,2),...
        markersize,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'))
ylabel('\it{}I')

subplot(2,2,4)
scatter(telasticStop,elastic(:,4),...
        markersize,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'))
ylabel('\it{}E_{el}*')


%% Figure formatting
for i=1:4
    subplot(2,2,i)
    hold on
    box on
    xlabel('\gamma')
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
end
