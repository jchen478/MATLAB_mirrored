close all;

%% Parameters
File = fopen('Parameters.in','r');
Parameters = textscan(File,'%f %*s',Inf,'Delimiter','\n');
fclose(File);
Parameters = cell2mat(Parameters);

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

%% Simulation Paramters
nL3 = nfib*L^3/side^3;
Seff = pi*kb/nseg^4; 
over_Seff = 1/Seff;

%% Time 
tINS = 0:dt*config_write:strain;
tStress = dt*stress_write:dt*stress_write:strain;

%% Read file
File = fopen('ISV.txt','r');
data = fscanf(File,'%f',[8 Inf])';
fclose(File);

File = fopen('Stress_tensor.txt','r');
s = fscanf(File,'%f',[7 Inf])';
fclose(File);

%% Define variables
% tINS = data(:,1);
I = data(:,2);
Sx = data(:,3);
Sy = data(:,4);
Sz = data(:,5);
Vx = data(:,6);
Vy = data(:,7);
Vz = data(:,8);
tauxz = s(:,4);

%% dimensionalize
sigxz = tauxz*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff + over_Seff;

%% Plotting
markersize = 10;
fontsize = 20; 
linewidth = 2; 
tickx = 3; 

figure('color','white')
filename = ('nfib12800_mu23');

subplot(2,2,1)
hold on
scatter(tINS,I,...
        markersize,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'))
ylabel('\it{}I')


subplot(2,2,2)
hold on
box on
scatter(tINS,Sx,...
    markersize,'MarkerEdgeColor',rgb('Crimson'),...
    'MarkerFaceColor',rgb('Crimson'))
scatter(tINS,Sy,...
    markersize,'MarkerEdgeColor',rgb('Chocolate'),...
    'MarkerFaceColor',rgb('Chocolate'))
scatter(tINS,Sz,...
    markersize,'MarkerEdgeColor',rgb('MediumBlue'),...
    'MarkerFaceColor',rgb('MediumBlue'))
ylabel('\it{}S')
legend('x','y','z','location','northwest')

subplot(2,2,3)
hold on
scatter(tINS,Vx,...
    markersize,'MarkerEdgeColor',rgb('Crimson'),...
    'MarkerFaceColor',rgb('Crimson'))
scatter(tINS,Vy,...
    markersize,'MarkerEdgeColor',rgb('Chocolate'),...
    'MarkerFaceColor',rgb('Chocolate'))
scatter(tINS,Vz,...
    markersize,'MarkerEdgeColor',rgb('MediumBlue'),...
    'MarkerFaceColor',rgb('MediumBlue'))
ylabel('\it{}V')
legend('x','y','z','location','northwest')

subplot(2,2,4)
hold on
scatter(tStress,sigxz,...
        markersize,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'))
ylabel('\it{\sigma_{xz} L^4/ E_Y I}')

for i=1:4
    
    subplot(2,2,i)
    hold on
    box on
    xlabel('\gamma')
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    xlim([0 inf])
end

set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
print(filename,'-dpng');
savefig(filename);