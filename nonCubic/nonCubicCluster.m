close all;
clear;
clc;

%% Cluster results
screen = [250 200 75 50 30 20 10];
nscreen = length(screen);
screenLegend = cell(1,nscreen);

for n=1:nscreen
    
    screenLegend{n} = ['min \it{}N_{fib}\rm{} per cluster = ',num2str(screen(n))];
    
end

legendArr = cell(1);
titleArr = cell(1);
i = 1;

filename={'nfib1280_cluster.txt','nonCubic_run001_cluster.txt','nonCubic_run002_cluster.txt'}; 
color = {rgb('DarkGreen'),rgb('DarkOrange'),rgb('MediumBlue')};

cluster_charFunc(filename{1},'cubic',color{i},i,screen)
titleArr{i} = ('(Lx,Ly,Lz) = (600,600,600)');
legendArr{i} = ('(Lx,Ly,Lz) = (600,600,600)'); i= i+1;

cluster_charFunc(filename{2},'run001',color{i},i,screen)
titleArr{i} = ('(Lx,Ly,Lz) = (450,1066,450)');
legendArr{i} = ('(Lx,Ly,Lz) = (450,1066,450)'); i= i+1;

cluster_charFunc(filename{3},'run002',color{i},i,screen)
titleArr{i} = ('(Lx,Ly,Lz) = (1066,450,450)');
legendArr{i} = ('(Lx,Ly,Lz) = (1066,450,450)'); i= i+1;

%% figure formatting
tickx = 1.5;

nfig = 2;

for j=1:nfig
    
    figure(j)
    hold on;
    box on;
    set(gca,'fontsize',24,'linewidth',2,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
    ylim([0 inf])
    xlim([2 inf])
    legend(legendArr{:},'location','best')
    xlabel('\mu')
    
end

for j=nfig+1:nfig+i-1
    figure(j)
    hold on
    box on;
    set(gca,'fontsize',24,'linewidth',2,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
    ylim([0 inf])
    xlim([2 inf])
    legend(screenLegend{:},'location','best')
    xlabel('\mu')
    title(titleArr{j-nfig})
    
end