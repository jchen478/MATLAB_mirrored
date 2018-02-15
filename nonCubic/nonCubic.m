close all;
clc;
clear;

filename={'cubic_1280.txt','nonCubic_run001.txt','nonCubic_run002.txt',...
    'nonCubic_run006.txt','nonCubic_run003.txt','nonCubic_run004.txt',...
    'nonCubic_run005.txt'}; 
color = {rgb('MediumBlue'),rgb('OrangeRed'),rgb('Purple'),...
    rgb('DarkGreen'),rgb('Crimson'),rgb('Orange'),rgb('Black')};

legendArr = cell(1);
i = 1;

nonCubicFun('nonCubic_run006.txt',color{i},450,450,1066)
legendArr{i} = ('(Lx,Ly,Lz) = (450,450,1066)'); i= i+1;


% 
% nonCubicFun(filename{1},color{i}, 600, 600, 600)
% legendArr{i} = ('(Lx,Ly,Lz) = (600,600,600)'); i= i+1;
% 
% nonCubicFun(filename{2},color{i},450,1066,450)
% legendArr{i} = ('(Lx,Ly,Lz) = (450,1066,450)'); i= i+1;
% 
% nonCubicFun(filename{3},color{i},1066,450,450)
% legendArr{i} = ('(Lx,Ly,Lz) = (1066,450,450)'); i= i+1;
% 
% nonCubicFun(filename{4},color{i},450,450,1066)
% legendArr{i} = ('(Lx,Ly,Lz) = (450,450,1066)'); i= i+1;

% nonCubicFun(filename{5},color{i},692.6,692.6,450)
% legendArr{i} = ('(Lx,Ly,Lz) = (692.6,692.6,450)'); i= i+1;
% 
% nonCubicFun(filename{6},color{i},692.6,450,692.6)
% legendArr{i} = ('(Lx,Ly,Lz) = (692.6,450,692.6)'); i= i+1;
% 
nonCubicFun(filename{7},color{i},692.6,692.6,450)
legendArr{i} = ('(Lx,Ly,Lz) = (692.6,692.6,450)'); i= i+1;

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
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    xlim([2 inf])
    xlabel('\mu')
    legend(legendArr{:},'location','best')
end

spreadfigures;
