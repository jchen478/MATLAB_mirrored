%% load files 
close all;

load('result64.mat')
load('nc12800.mat')
load('stress12800.mat')

mu = result64(:,1); 
I = result64(:,2);
Sx = result64(:,3);
Sy = result64(:,4);
Sz = result64(:,5);
Vx = result64(:,6);
Vy = result64(:,7);
Vz = result64(:,8);

%% Intensity and scale
figure('color','white')
hold on
box on;
scatter(Sx,I,40,'MarkerEdgeColor',rgb('Crimson'),...
        'MarkerFaceColor',rgb('Crimson'))
scatter(Sy,I,40,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'))
scatter(Sz,I,40,'MarkerEdgeColor',rgb('Green'),...
        'MarkerFaceColor',rgb('Green'))
xlabel('\it{S}')
ylabel('\it{I}')
legend('x (flow direction)','y (spanwise direction)','z (gradient direction)','location','northwest')
set(gca,'fontsize',30,'fontName','Times New Roman','linewidth',2)
set(gca,'YMinorTick','on','XMinorTick','on')
dx = -0.004;
dy = 0.000;
text(Sx+dx,I+dy,num2str(mu),'fontsize',16,'fontName','Times New Roman','color',rgb('Crimson'))
text(Sy+dx,I+dy,num2str(mu),'fontsize',16,'fontName','Times New Roman','color',rgb('MediumBlue'))
text(Sz+dx,I+dy,num2str(mu),'fontsize',16,'fontName','Times New Roman','color',rgb('Green'))

%% Intensity and scale
figure('color','white')
hold on
box on;
plot(Sx,I,'--o','MarkerSize',10,'Color',rgb('Crimson'),...
    'MarkerEdgeColor',rgb('Crimson'),...
    'MarkerFaceColor',rgb('Crimson'),'LineWidth',2)
plot(Sy,I,'--o','MarkerSize',10,'Color',rgb('MediumBlue'),...
    'MarkerEdgeColor',rgb('MediumBlue'),...
    'MarkerFaceColor',rgb('MediumBlue'),'LineWidth',2)
plot(Sz,I,'--o','MarkerSize',10,'Color',rgb('Green'),...
    'MarkerEdgeColor',rgb('Green'),...
    'MarkerFaceColor',rgb('Green'),'LineWidth',2)
xlabel('\it{S}')
ylabel('\it{I}')
legend('x (flow direction)','y (spanwise direction)','z (gradient direction)','location','northwest')
set(gca,'fontsize',30,'fontName','Times New Roman','linewidth',2)
set(gca,'YMinorTick','on','XMinorTick','on')
dx = -0.004;
dy = 0.000;
text(Sx+dx,I+dy,num2str(mu),'fontsize',16,'fontName','Times New Roman','color',rgb('Crimson'))
text(Sy+dx,I+dy,num2str(mu),'fontsize',16,'fontName','Times New Roman','color',rgb('MediumBlue'))
text(Sz+dx,I+dy,num2str(mu),'fontsize',16,'fontName','Times New Roman','color',rgb('Green'))

%% Scale
figure('color','white')
hold on
box on;
% scatter(mu,Sx,40,'MarkerEdgeColor',rgb('Crimson'),...
%         'MarkerFaceColor',rgb('Crimson'))
% scatter(mu,Sy,40,'MarkerEdgeColor',rgb('MediumBlue'),...
%         'MarkerFaceColor',rgb('MediumBlue'))
% scatter(mu,Sz,40,'MarkerEdgeColor',rgb('Green'),...
%         'MarkerFaceColor',rgb('Green'))
plot(mu,Sx,'--o','MarkerSize',10,'Color',rgb('Crimson'),...
    'MarkerEdgeColor',rgb('Crimson'),...
    'MarkerFaceColor',rgb('Crimson'),'LineWidth',2)
plot(mu,Sy,'--o','MarkerSize',10,'Color',rgb('MediumBlue'),...
    'MarkerEdgeColor',rgb('MediumBlue'),...
    'MarkerFaceColor',rgb('MediumBlue'),'LineWidth',2)
plot(mu,Sz,'--o','MarkerSize',10,'Color',rgb('Green'),...
    'MarkerEdgeColor',rgb('Green'),...
    'MarkerFaceColor',rgb('Green'),'LineWidth',2)
xlabel('\mu')
ylabel('\it{S}')
legend('x (flow direction)','y (spanwise direction)','z (gradient direction)','location','northwest')
set(gca,'fontsize',30,'fontName','Times New Roman','linewidth',2)
set(gca,'YMinorTick','on','XMinorTick','on')

%% Intensity
figure('color','white')
hold on
box on;
% scatter(mu,I,40,'MarkerEdgeColor',rgb('MediumBlue'),...
%         'MarkerFaceColor',rgb('MediumBlue'))
plot(mu,I,'--o','markerSize',5,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'), 'color',rgb('MediumBlue'), 'linewidth',2.0)
    
xlabel('\mu')
ylabel('\it{I}')
set(gca,'fontsize',30,'fontName','Times New Roman','linewidth',2)
set(gca,'YMinorTick','on','XMinorTick','on')

%% V/S3
figure('color','white')

subplot(3,1,1)
hold on
box on;
scatter(mu,Vx./Sx.^3,40,'MarkerEdgeColor',rgb('Crimson'),...
        'MarkerFaceColor',rgb('Crimson'))
xlabel('\mu')
ylabel('V_x/S_x^3')
set(gca,'fontsize',30,'fontName','Times New Roman','linewidth',2)
set(gca,'YMinorTick','on','XMinorTick','on')
% set(gca,'yscale','log')

subplot(3,1,2)
hold on
box on;
scatter(mu,Vy./Sy.^3,40,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'))
xlabel('\mu')
ylabel('V_y/S_y^3')
set(gca,'fontsize',30,'fontName','Times New Roman','linewidth',2)
set(gca,'YMinorTick','on','XMinorTick','on')
% set(gca,'yscale','log')

subplot(3,1,3)
hold on
box on;
scatter(mu,Vz./Sz.^3,40,'MarkerEdgeColor',rgb('Green'),...
        'MarkerFaceColor',rgb('Green'))
xlabel('\mu')
ylabel('V_z/S_z^3')
set(gca,'fontsize',30,'fontName','Times New Roman','linewidth',2)
set(gca,'YMinorTick','on','XMinorTick','on')
% set(gca,'yscale','log')

%% V/S3
figure('color','white')
hold on
box on;
scatter(mu,Vx./(4*Sx.^3),40,'MarkerEdgeColor',rgb('Crimson'),...
        'MarkerFaceColor',rgb('Crimson'))
scatter(mu,Vy./(4*Sy.^3),40,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'))
scatter(mu,Vz./(4*Sz.^3),40,'MarkerEdgeColor',rgb('Green'),...
        'MarkerFaceColor',rgb('Green'))
% plot(mu,Vx./(4*Sx.^3),'--o','MarkerSize',10,'Color',rgb('Crimson'),...
%     'MarkerEdgeColor',rgb('Crimson'),...
%     'MarkerFaceColor',rgb('Crimson'),'LineWidth',2)
% plot(mu,Vy./(4*Sy.^3),'--o','MarkerSize',10,'Color',rgb('MediumBlue'),...
%     'MarkerEdgeColor',rgb('MediumBlue'),...
%     'MarkerFaceColor',rgb('MediumBlue'),'LineWidth',2)
% plot(mu,Vz./(4*Sz.^3),'--o','MarkerSize',10,'Color',rgb('Green'),...
%     'MarkerEdgeColor',rgb('Green'),...
%     'MarkerFaceColor',rgb('Green'),'LineWidth',2)    
xlabel('\mu')
ylabel('\it{}V\rm{}/(4\it{}S^3\rm{})')
legend('x (flow direction)','y (spanwise direction)','z (gradient direction)','location','northeast')
set(gca,'fontsize',30,'fontName','Times New Roman','linewidth',2)
set(gca,'YMinorTick','on','XMinorTick','on')

%% Number of contacts
figure('color','white')
hold on
box on;
% scatter(nc12800(:,1),nc12800(:,2),40,'MarkerEdgeColor',rgb('MediumBlue'),...
%         'MarkerFaceColor',rgb('MediumBlue'))
plot(nc12800(:,1),nc12800(:,2),'--o','markerSize',5,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'), 'color',rgb('MediumBlue'), 'linewidth',2.0)
xlabel('\mu')
ylabel('<N_c>')
set(gca,'fontsize',30,'fontName','Times New Roman','linewidth',2)
set(gca,'YMinorTick','on','XMinorTick','on')

%% Stress
figure('color','white')
hold on
box on;
% scatter(stress12800(:,1),stress12800(:,2),40,'MarkerEdgeColor',rgb('MediumBlue'),...
%         'MarkerFaceColor',rgb('MediumBlue'))
plot(stress12800(:,1),stress12800(:,2),'--o','markerSize',5,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'), 'color',rgb('MediumBlue'), 'linewidth',2.0)
xlabel('\mu')
ylabel('\it{\sigma_{xz} L^4/ E_Y I}')
set(gca,'fontsize',30,'fontName','Times New Roman','linewidth',2)
set(gca,'YMinorTick','on','XMinorTick','on')