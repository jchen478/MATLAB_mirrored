clc; 
close all; 
clear;

%% parameters
Lx = 1293;

%% with friction
File = fopen('dx_friction640.txt','r');
formatSpec = '%f';
sizeA = [Inf];
dx = fscanf(File,formatSpec,sizeA); 
fclose(File);

dx = sort(dx);
nbinx = round(Lx./dx);
nbin = nbinx.^3;
ncases = length(dx);
S_friction = cell(ncases,1); 
F_friction = cell(ncases,1);
S_frictionless = cell(ncases,1); 
F_frictionless = cell(ncases,1);

time = 0:0.5:640;
ndata = length(time);

H_friction = zeros(ndata,ncases); 
H_frictionless = zeros(ndata,ncases); 


for i=1:ncases   
    
    % create file name and legend name
    S_friction{i} = ['dx = ',num2str(dx(i),'%.2f\n')];  
    F_friction{i} = ['Shannon_entropy_friction640_',num2str(dx(i),'%.2f\n'),'.txt'];
    
    % read file
    File = fopen(F_friction{i},'r');
    formatSpec = '%f %f';
    sizeA = [2 Inf];
    data = fscanf(File,formatSpec,sizeA)';
    fclose(File);
    H_friction(:,i) = data(2:end,2);
    
end


%% frictionless

for i=1:ncases   
    
    % create file name and legend name
    S_frictionless{i} = ['dx = ',num2str(dx(i),'%.2f\n')];  
    F_frictionless{i} = ['Shannon_entropy_frictionless640_',num2str(dx(i),'%.2f\n'),'.txt'];
    
    % read file
    File = fopen(F_frictionless{i},'r');
    formatSpec = '%f %f';
    sizeA = [2 Inf];
    data = fscanf(File,formatSpec,sizeA)';
    fclose(File);
    H_frictionless(:,i) = data(2:end,2);
    
end

%% normalize
for i=1:ncases
    
    H_friction(:,i) = H_friction(:,i)./log(nbin(i));
    H_frictionless(:,i) = H_frictionless(:,i)./log(nbin(i));
end

%% plotting
figure('color','white')
hold on

for i=1:ncases
    
    subplot(3,3,i)
    hold on
    title(S_friction{i})
    scatter(time,H_friction(:,i),10,'MarkerEdgeColor',rgb('Crimson'),...
        'MarkerFaceColor',rgb('Crimson'));
    scatter(time,H_frictionless(:,i),10,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'));
    
end








