close all;
clc
clear

%% parameters
nfib = 384;
a = 1.60E-05;
Imom = pi*a^4/4;
EY = 1e9; 
eta0 = 1; 
side = 246.5;
kb = 10;
nseg = 8;
rps = 10;

Seff = pi.*kb./nseg.^4;
rp = rps.*nseg;
nL3 = nfib .* (2.*rp./side).^3;

%% dimensional parameters
L = 2*rps.*nseg*a;
gamma = EY*Imom ./ (Seff.*L.^4*eta0);
ts = 1/(Seff*gamma);

%% declare number of files
%casename = {'1a','1b','2a','2b','3a','3b','4a','4b'};
casename = {'4a','4b'};
ncase = length(casename);
filename = cell(4,1);
for i=1:ncase    
   filename{i} = ['Eelastic',casename{i},'.txt'];    
end

%% obtain dimensions and key variables
sizeA = [4 Inf];
File = fopen(filename{1},'r');
data = fscanf(File,'%f',sizeA)';
fclose(File);

strain_stop = data(1,1); 
dt = data(1,2); 
elastic_write = data(1,3);
stress_write = data(1,4);
data(1,:) = []; 

%% time scales
t = data(:,1);
tAfterStop = t-strain_stop; 
ndata = length(t); 
% 
% telastic = t;
% tAfterStopelastic = tAfterStop;

% telastic = t*Seff;
% tAfterStopelastic = tAfterStop*Seff;
% 
telastic = t*Seff*(410/pi);
tAfterStopelastic = tAfterStop*Seff*(410/pi);


nstop = strain_stop / elastic_write /dt; 

%% allocate space for arrays
E = zeros(ndata,ncase); 
Ebend = zeros(ndata,ncase); 
Etwist = zeros(ndata,ncase); 
Enorm = zeros(ndata,ncase); 

%% read all files
for i=1:ncase
    File = fopen(filename{i},'r');
    data = fscanf(File,'%f',sizeA)';
    fclose(File);
    data(1,:) = [];
    E(:,i) = data(:,4);
    Ebend(:,i) = data(:,2);
    Etwist(:,i) = data(:,3);
    Enorm(:,i) = E(:,i)./mean(data(1:nstop,4));
end

%% plot

% 2*E(:,1)*8*pi/Seff/nseg^3

% E = E*8*pi/Seff/nseg^3 *EY*Imom/L / (0.002*pi*eta0*gamma*L^3);
% 
% E = E*8/(0.002*nseg^3);

 E = E*4*pi/nseg^3/Seff;

for j=1:2:ncase
 
    % E_elastic
    figure
    hold on
    scatter(tAfterStopelastic, 2*E(:,j),10,'MarkerEdgeColor',rgb('Crimson'),...
        'MarkerFaceColor',rgb('Crimson'));
    scatter(tAfterStopelastic, 2*E(:,j+1),10,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'));
    ylabel('\it{}E_{elastic}*')
    
    % normalized elastic
    figure
    hold on
    scatter(tAfterStopelastic, Enorm(:,j),10,'MarkerEdgeColor',rgb('Crimson'),...
        'MarkerFaceColor',rgb('Crimson'));
    scatter(tAfterStopelastic, Enorm(:,j+1),10,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'));
    ylabel('\it{}E_{elastic} / E_{elastic,0}')
    
end

%% figure formatting
tickx = 1.5;

for i=1:ncase 
    
    figure(i)
    hold on;
    box on;    
    set(gca,'fontsize',24,'linewidth',2,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
    xlim([tAfterStopelastic(1) inf])
    ylim([0 inf])
    legend('nL3 = 105','nL3 = 52')
    xlabel('dimensionless time')
end