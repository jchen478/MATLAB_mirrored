
%% Read files
clc; clear; close all;

Seff = 0.0077;
elasticFun('Eelastic4a.txt',Seff)
elasticFun('Eelastic4b.txt',Seff)


%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 24;
linewidth = 1.5;


%% Figure formatting
for i=1:1
    figure(i)
    hold on
    box on
    xlabel('\gamma')
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    xlim([0 500])
end