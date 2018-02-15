%% cleaning
clc;
close all;

%% Info
% input is a helical floc 
% nfib = 320; nseg = 5; nL3 = 40; Seff = 0.05;
% dU_slide of broken contacts are plotted as a function of mu_kin
% strain = 2.0, dU_slide printed every 100 steps

mukin = [0 5 10 15]; 
ncases = length(mukin); 
S = cell(ncases,1); 
F = cell(ncases,1);

for i=1:ncases   
    S{i} = ['\mu_{kin} = ',num2str(mukin(i))];  
    F{i} = ['mu_kin',num2str(mukin(i)),'.txt'];
end

%% helical floc 
figure('color','white')
hold on;
for i=1:ncases
    File = fopen(F{i},'r');
    formatSpec = '%d %d %f %f %f %f %f';
    sizeA = [7 Inf];
    f = fscanf(File,formatSpec,sizeA)';
    fclose(File);
    dU_slide = f(:,6);
    mu_kinPlot = mukin(i)*ones(length(dU_slide),1);
    scatter(dU_slide, mu_kinPlot, 'filled')
end

% legend(S{:})
title('Helical floc \mu_{stat} = 20')
xlabel('||dU_{slide}||')
ylabel('\mu_{kin}')
set(gca,'fontsize',16)

%% histogram
for i=1:ncases
    File = fopen(F{i},'r');
    formatSpec = '%d %d %f %f %f %f %f';
    sizeA = [7 Inf];
    f = fscanf(File,formatSpec,sizeA)';
    fclose(File);
    dU_slide = f(:,6);   
    figure('color','white')
    histogram(dU_slide,'BinWidth',5, 'BinLimits',[0 100],'Normalization','probability')
    xlabel('||dU_{slide}||')
    title(['Helical floc \mu_{stat} = 20 ',S{i}])
    set(gca,'fontsize',16)
end
