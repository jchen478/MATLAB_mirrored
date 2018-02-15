clc;
clear;
close all;
file = {'320.txt','640.txt','1280.txt','2560.txt','6400.txt','12800.txt'};

% file = {'1280.txt'};

nconc = 4;
yieldStress = zeros(length(file)*nconc,3);

for i=1:length(file)
    
    yieldStress = yieldFunc(file{i},i,yieldStress,nconc);

end

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 20;
linewidth = 2;

%% figure formatting
for i=1:length(file)
    
    figure(i)
    hold on
    box on
    xlabel('1/S_{eff}')
    ylabel('\it{\sigma_{xz} L^4/ E_Y I}')
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    
end

yieldStress