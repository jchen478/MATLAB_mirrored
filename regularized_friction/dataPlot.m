close all;
clc; 
clear all;

%% Info
%  data taken from flocculated input with nfib = 160 

%% intact 3x3
File = fopen('data_intact.txt','r');
formatSpec = '%d %d %f %f %f %f';
sizeA = [6 Inf];
data = fscanf(File,formatSpec,sizeA)';
fclose(File);

dUx = data(:,3); 
dUy = data(:,4); 
dUz = data(:,5); 
dUnorm_prev = data(:,6); 

mu_stat = ones(length(dUnorm_prev),1); 
mu_stat = 20*mu_stat; 

%% broken 3x3
File = fopen('data_broken.txt','r');
formatSpec = '%d %d %f %f %f %f';
sizeA = [6 Inf];
data_broken = fscanf(File,formatSpec,sizeA)';
fclose(File);

dUnorm_prev_broken = data_broken(:,6); 

mu_kin = ones(length(dUnorm_prev_broken),1); 
mu_kin = 0.0*mu_kin; 
 
 %% plotting

figure('color','white');
hold on;
scatter(dUnorm_prev, mu_stat, 'filled');
scatter(dUnorm_prev_broken, mu_kin, 'filled');
legend('||dU_{slide}|| intact','||dU_{slide}|| broken');
ylabel('\mu')
xlabel('||dU||')