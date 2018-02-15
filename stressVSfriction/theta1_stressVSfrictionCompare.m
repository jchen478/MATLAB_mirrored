close all;
clc;
clear;

dataPath = '../data_stressVSfriction/averagedValues/';

filename={'theta1_nfib160.txt', 'theta1_nfib1280.txt' ,'theta1_nfib12800.txt'}; 
color = {rgb('Purple'),rgb('OrangeRed'),rgb('DeepSkyBlue')};

legendArr = cell(1);
i = 1;

stressVSfrictionFun([dataPath,filename{1}],color{i})
legendArr{i} = ('$N_{fib} = 160$'); i= i+1;

stressVSfrictionFun([dataPath,filename{2}],color{i})
legendArr{i} = ('$N_{fib} = 1280$'); i= i+1;

stressVSfrictionFun([dataPath,filename{3}],color{i})
legendArr{i} = ('$N_{fib} = 12800$'); i= i+1;

%% Figure formatting
for i=1:8
    figure(i)
    hold on
    box on   
    xlim([0 inf])
    xlabel('$\mu$')
    legend(legendArr{:},'location','best')
end

spreadfigures;