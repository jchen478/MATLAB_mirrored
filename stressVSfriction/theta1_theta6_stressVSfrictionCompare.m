close all;
clc;
clear;

dataPath = '../data_stressVSfriction/averagedValues/';

filename={'nfib160.txt', 'nfib1280_averaged.txt' ,'nfib12800.txt','theta1_nfib160.txt', 'theta1_nfib1280.txt' ,'theta1_nfib12800.txt'}; 

color = {rgb('DarkGreen'),rgb('DarkOrange'),rgb('MediumBlue'),rgb('Purple'),rgb('OrangeRed'),rgb('DeepSkyBlue')};

legendArr = cell(1);
i = 1;

stressVSfrictionFun([dataPath,filename{1}],color{i})
legendArr{i} = ('$(\theta_{eq}, N_{fib}) = (0.6,160)$'); i= i+1;

stressVSfrictionFun([dataPath,filename{2}],color{i})
legendArr{i} = ('$(\theta_{eq}, N_{fib}) = (0.6,1280)$'); i= i+1;

stressVSfrictionFun([dataPath,filename{3}],color{i})
legendArr{i} = ('$(\theta_{eq}, N_{fib}) = (0.6,12800)$'); i= i+1;

stressVSfrictionFun([dataPath,filename{4}],color{i})
legendArr{i} = ('$(\theta_{eq}, N_{fib}) = (0.1,160)$'); i= i+1;

stressVSfrictionFun([dataPath,filename{5}],color{i})
legendArr{i} = ('$(\theta_{eq}, N_{fib}) = (0.1,1280)$'); i= i+1;

stressVSfrictionFun([dataPath,filename{6}],color{i})
legendArr{i} = ('$(\theta_{eq}, N_{fib}) = (0.1,12800)$'); i= i+1;

%% Figure formatting
for i=1:8
    figure(i)
    hold on
    box on
    xlim([0 inf])
    xlabel('$\mu$')
    legend(legendArr{:},'location','bestoutside')
end

