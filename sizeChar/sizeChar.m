clc;
clear;
close all;

%% Case parameters
floc_cutoff = 50;
maxGap = 1;
nfib = [12800 1280 160];

nMaxGap = length(maxGap);
nNfib = length(nfib);
ncase = nMaxGap*nNfib;
filename = cell(ncase,1);
legendNfib = cell(nNfib,1);

%% Figure formatting parameters
nfig = 3;
markersize = 12;
tickx = 1.5;
fontsize = 24;
linewidth = 2.5;
colorArr = {rgb('MediumBlue'),rgb('Orange'),rgb('DarkGreen')};

%% Read file and plot
for k=1:nNfib
    legendNfib{k} = ['$N_{fib} =$ ',num2str(nfib(k))];
    for j=1:nMaxGap
        ind = (k-1)*nMaxGap+j;
        filename{ind} = ['Cluster_size_',num2str(nfib(k)),'_',num2str(maxGap(j)),'.txt'];
        sizeCharFun(filename{ind},markersize,linewidth,colorArr{k},floc_cutoff);
    end
end

%% Format
for i=1:2
    figure(i)
    for f=1:nfig
        
        set(gcf, 'color','white')
        subplot(1,3,f)
        hold on
        box on
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(gca,'ticklength',tickx*get(gca,'ticklength'))
        set(gcf, 'color','white')
        xlabel('$\mu$')
        xlim([-1 24])
        if i == 1 
            
            line([0 25], [1293 1293], 'linewidth',2,'color',rgb('MediumBlue'));
            line([0 25], [600 600], 'linewidth',2,'color',rgb('Orange'));
            line([0 25], [300 300], 'linewidth',2,'color',rgb('DarkGreen'));
%             legend(legendNfib{:},'\it{}L_{12800}','\it{}L_{1280}','\it{}L_{160}','location','bestoutside')
        end
        if i == 2
            legend(legendNfib{:},'location','bestoutside')
        end
    end
end
for i=3:3
    figure(i)
    for f=1:nfig
        set(gcf, 'color','white')
        subplot(1,3,f)
        hold on
        box on
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(gca,'ticklength',tickx*get(gca,'ticklength'))
        set(gcf, 'color','white')
    
    end
end