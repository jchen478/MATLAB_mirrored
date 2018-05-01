close all; 
clc;
clear;


set ( 0 , 'DefaultAxesLineWidth' ,2 ,'DefaultAxesFontName' , 'CMU Serif',...
'DefaultAxesFontSize' ,12 , 'DefaultAxesTickLength' , [ 0.025 0.015],...
'DefaultLineLineWidth' ,2 , 'DefaultTextInterpreter' , 'Latex',...
'DefaultLegendInterpreter' , 'Latex','defaultFigureColor','w')
set ( 0 , 'DefaultAxesTickLabelInterpreter' , 'Latex' )

colorOrder = [  rgb('DarkRed'); rgb('Crimson'); rgb('OrangeRed');
    rgb('Orange'); rgb('Gold'); rgb('Lime');
    rgb('Olive'); rgb('DarkGreen'); rgb('LightSkyBlue');
    rgb('MediumBlue'); rgb('Plum'); rgb('Purple')] 

set(groot,'DefaultAxesColorOrder',colorOrder, ...  
    'DefaultAxesLineStyleOrder','-|--|:|-.');


t=0:0.1:2*pi;
x = sin(t);
y = cos(t);
z = sin(t)+2*cos(t);
figure()
hold on
plot(t,x)
plot(t,z)
xlabel('$t$')
ylabel('$x$')
title('Figure1')
figure()
plot(t,y)
xlabel('$t$')
ylabel('$y$')
title('Figure2')

nfig = 2;
for i=1:nfig
    h = figure(i)
    box on;
    set(gcf, 'Units', 'Inches','Position',[0.5 4 3.5 2.5]);
    
end
