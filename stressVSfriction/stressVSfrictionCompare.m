close all;
clc;
clear;

dataPath = '../data_stressVSfriction/averagedValues/theta6_';

filename={'nfib160.txt', 'nfib1280.txt' ,'nfib12800.txt'}; 
color = {rgb('DarkGreen'),rgb('DarkOrange'),rgb('MediumBlue')};

legendArr = cell(1);
i = 1;

stressVSfrictionFun([dataPath,filename{1}],color{i})
legendArr{i} = ('$N_{fib} = 160$'); i= i+1;

stressVSfrictionFun([dataPath,filename{2}],color{i})
legendArr{i} = ('$N_{fib} = 1280$'); i= i+1;

stressVSfrictionFun([dataPath,filename{3}],color{i})
legendArr{i} = ('$N_{fib} = 12800$'); i= i+1;

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 24;
linewidth = 3;

%% Figure formatting
for i=1:8
    figure(i)
    hold on
    box on   
    xlim([0 inf])
    xlabel('$\mu$')
    legend(legendArr{:},'location','best','orientation','horizontal')
end

spreadfigures;