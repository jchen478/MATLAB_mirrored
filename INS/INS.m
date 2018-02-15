clc
clear
close all;

%% declare variables
% theta = [0 0.1 0.2 0.3 0.4 0.5 0.6];
% mu = [0 1 3 5 10 15 20];
theta = [0 0.1 0.2 0.3 0.4];
mu = [0 1 3 5 10];
[X,Y] = meshgrid(theta,mu);
format = [length(mu) length(theta)];

%% Read expansion with box results
File = fopen('I_expansionWithBox.txt','r');
I = fscanf(File,'%f',format);
fclose(File);

File = fopen('Sx_expansionWithBox.txt','r');
Sx = fscanf(File,'%f',format)';
fclose(File);

File = fopen('Sy_expansionWithBox.txt','r');
Sy = fscanf(File,'%f',format)';
fclose(File);

File = fopen('Sz_expansionWithBox.txt','r');
Sz = fscanf(File,'%f',format)';
fclose(File);

File = fopen('Vx_expansionWithBox.txt','r');
Vx = fscanf(File,'%f',format)';
fclose(File);

File = fopen('Vy_expansionWithBox.txt','r');
Vy = fscanf(File,'%f',format)';
fclose(File);

File = fopen('Vz_expansionWithBox.txt','r');
Vz = fscanf(File,'%f',format)';
fclose(File);

%% Read low_nL3 results
File = fopen('I_low_nL3.txt','r');
I_low_nL3 = fscanf(File,'%f',format);
fclose(File);

File = fopen('Sx_low_nL3.txt','r');
Sx_low_nL3 = fscanf(File,'%f',format)';
fclose(File);

File = fopen('Sy_low_nL3.txt','r');
Sy_low_nL3 = fscanf(File,'%f',format)';
fclose(File);

File = fopen('Sz_low_nL3.txt','r');
Sz_low_nL3 = fscanf(File,'%f',format)';
fclose(File);

File = fopen('Vx_low_nL3.txt','r');
Vx_low_nL3 = fscanf(File,'%f',format)';
fclose(File);

File = fopen('Vy_low_nL3.txt','r');
Vy_low_nL3 = fscanf(File,'%f',format)';
fclose(File);

File = fopen('Vz_low_nL3.txt','r');
Vz_low_nL3 = fscanf(File,'%f',format)';
fclose(File);

%% Find ratio
I_ratio = abs(I-I_low_nL3)./I_low_nL3*100; 
Sx_ratio = abs(Sx-Sx_low_nL3)./Sx_low_nL3*100; 
Sy_ratio = Sy./Sy_low_nL3; 
Sz_ratio = Sz./Sz_low_nL3; 
Vx_ratio = Vx./Vx_low_nL3; 
Vy_ratio = Vy./Vy_low_nL3; 
Vz_ratio = Vz./Vz_low_nL3; 

%% Plotting

viewPoint = [-12 16];
figure('color','white')

%% intensity
subplot(2,3,1)
hold on;
box on;
rotate3d on;
surf(X,Y,I,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}I')
title('After redispersion')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
colorbar 
view(viewPoint);

subplot(2,3,2)
hold on;
surf(X,Y,I_low_nL3,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}I')
title('Without redispersion')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
colorbar 
view(viewPoint);

subplot(2,3,3)
hold on;
surf(X,Y,I_ratio,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('%\it{}I')
title('Percent different')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
colorbar 
view(viewPoint);

%% Sx
subplot(2,3,4)
hold on;
box on;
rotate3d on;
surf(X,Y,Sx,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}Sx')
title('After redispersion')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
colorbar 
view(viewPoint);

subplot(2,3,5)
hold on;
surf(X,Y,Sx_low_nL3,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}Sx')
title('Without redispersion')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
colorbar 
view(viewPoint);

subplot(2,3,6)
hold on;
surf(X,Y,Sx_ratio,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}Sx')
title('Percent different')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
colorbar 
view(viewPoint);

%% Never concentrated
figure('color','white')
hold on;
box on;
surf(X,Y,I_low_nL3,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}I')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
colorbar 
view(viewPoint);


figure('color','white')
subplot(2,3,1)
hold on;
box on;
surf(X,Y,Sx_low_nL3,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}Sx')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
xlim([0 0.6])
ylim([0 20])
colorbar 
view(viewPoint);

subplot(2,3,2)
hold on;
box on;
surf(X,Y,Sy_low_nL3,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}Sy')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
xlim([0 0.6])
ylim([0 20])
colorbar 
view(viewPoint);

subplot(2,3,3)
hold on;
box on;
surf(X,Y,Sz_low_nL3,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}Sz')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
xlim([0 0.6])
ylim([0 20])
colorbar 
view(viewPoint);

subplot(2,3,4)
hold on;
box on;
surf(X,Y,Vx_low_nL3,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}Vx')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
xlim([0 0.6])
ylim([0 20])
colorbar 
view(viewPoint);

subplot(2,3,5)
hold on;
box on;
surf(X,Y,Vy_low_nL3,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}Vy')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
xlim([0 0.6])
ylim([0 20])
colorbar 
view(viewPoint);

subplot(2,3,6)
hold on;
box on;
surf(X,Y,Vz_low_nL3,'linestyle','none')
xlabel('\theta_{eq}')
ylabel('\mu_{stat}')
zlabel('\it{}Vz')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
xlim([0 0.6])
ylim([0 20])
colorbar 
view(viewPoint);
