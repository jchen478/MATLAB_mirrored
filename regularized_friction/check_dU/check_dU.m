clc;
close all;
clear

File = fopen('helical.txt','r');
formatSpec = '%f';
sizeA = [1 Inf];
floc = fscanf(File,formatSpec,sizeA)';
fclose(File);

mean(floc)

N = length(floc);
ind = zeros(N,1);

for i=1:N
    ind(i) = i;
end


%% Info
% input is a U-shaped floc 
% nfib = 320; nseg = 5; nL3 = 40; Seff = 0.05;
% dU_slide of broken contacts are plotted as a function of mu_kin
% strain = 2.0, dU_slide printed every 100 steps

mukin = [0 0.2 0.4 0.6 0.8]; 
nfric = [50 100];
U_scale = [50 100 150 200];
ncases = length(mukin)*length(nfric)*length(U_scale); 
S = cell(ncases+1,1); 
F = cell(ncases,1);
avg = zeros(ncases,1);
avgLinex = [1 ncases];
avgLiney = [mean(floc) mean(floc)];

for i=1:length(mukin)  
    for j=1:length(nfric)
           for k=1:length(U_scale)
               n = (i-1)*length(nfric)*length(U_scale)+(j-1)*length(U_scale)+k;
               F{n} = ['helical_new3_',num2str(mukin(i)),'_',num2str(nfric(j)),'_',num2str(U_scale(k)),'.txt'];
               S{n} = ['\mu_{kin} = ',num2str(mukin(i)),' n = ',num2str(nfric(j)), ' U_{scale} = ',num2str(U_scale(k))]; 
           end
    end
end

S{end} = ['Before applied friction'];

figure('color','white')
hold on;

for i=1:ncases
    F{i}
    File = fopen(F{i},'r');
    formatSpec = '%f';
    sizeA = [1 Inf];
    dU = fscanf(File,formatSpec,sizeA)';
    avg(i) = mean(dU);
    fclose(File);
    scatter(ind, dU, 'filled')
end


plot(ind,floc,'linewidth',1.1,'color','black')
% legend(S{:})
title('Helical: Sliding Velocity After Friction Applied \mu_{stat} = 1')
xlabel('contact pair index')
ylabel('||dU_{slide}||')
xlim([0 N])
% ylim([0 60000])
set(gca,'fontsize',16)

figure('color','white')
hold on;
for i=1:ncases
    scatter(ind(i), avg(i), 'filled')
end
plot(avgLinex,avgLiney,'linewidth',1.2)

title('Helical: Average Sliding Velocity After Friction Applied \mu_{stat} = 1')
xlabel('case index')
ylabel('||dU_{slide}||')
xlim([0 ncases])
set(gca,'fontsize',16)
legend(S{:})


File = fopen('output.txt','r');
formatSpec = '%f';
sizeA = [1 Inf];
floc = fscanf(File,formatSpec,sizeA)';
fclose(File);
