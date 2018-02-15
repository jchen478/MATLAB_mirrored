clc;
clear;
close all;

%% Case parameters
floc_cutoff = 30;
maxGap = [2 4];
nMaxGap = length(maxGap);

filename = cell(nMaxGap,1);
legendArr = cell(nMaxGap,1);

%% Figure formatting parameters
nfig = 3;
markersize = 12;
tickx = 1.5;
fontsize = 24;
linewidth = 2.5;
colorArr = {rgb('MediumBlue'),rgb('Orange'),rgb('DarkGreen'),rgb('Crimson')};
runArr = {'run001_','run002_','run003_','run004_','run005_','run006_'};

%% Read file and plot
for i=1:6
    for j=1:nMaxGap

        filename{j} = ['nonCubic_Cluster_size_',runArr{i},num2str(maxGap(j)),'.txt'];
        sizeCharFun(filename{j},markersize,linewidth,colorArr{j},floc_cutoff);
        legendArr{j} = ['maxGap = ',num2str(maxGap(j))];
    end
end

%% Format
for f1=1:2
    figure(f1)
    for f=1:nfig
        set(gcf, 'color','white')
        subplot(1,3,f)
        hold on
        box on
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(gca,'ticklength',tickx*get(gca,'ticklength'))
        set(gcf, 'color','white')
        xlabel('\mu')
        xlim([1 5])
        legend(legendArr,'location','best')
    end
end
for f1=3:3
    figure(f1)
    for f=1:nfig
        set(gcf, 'color','white')
        subplot(1,3,f)
        hold on
        box on
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(gca,'ticklength',tickx*get(gca,'ticklength'))
        set(gcf, 'color','white')
%         legend(legendArr,'location','best')
    end
end