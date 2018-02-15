close all;
clc;
clear

color = {rgb('DarkGreen'),rgb('DarkOrange'),rgb('MediumBlue'),...
    rgb('Purple'),rgb('Crimson'),rgb('DarkGreen'),...
   rgb('Lime'),rgb('Cyan')};

legendArr = cell(1);
i = 1;
%% cases for comparison

% nfib = 3840, nseg = 5, rps = 20, L = 200, R = 200
AR_function('3840run001.txt',color{i})
legendArr{i} = ('nfib = 3840, nseg = 5, rps = 20, L = 200, R = 200'); i= i+1;
% nfib = 1920, nseg = 10, rps = 20, L = 400, R = 200 
AR_function('1920run001.txt',color{i})
legendArr{i} = ('nfib = 1280, nseg = 10, rps = 20, L = 400, R = 200'); i= i+1;
% nfib = 1280, nseg = 15, rps = 20, L = 600, R = 200 
AR_function('1280run001.txt',color{i})
legendArr{i} = ('nfib = 1280, nseg = 15, rps = 20, L = 600, R = 200'); i= i+1;

% % nfib = 3840, nseg = 10, rps = 10, L = 200, R = 225 
% AR_function('3840run002.txt',color{i})
% legendArr{i} = ('nfib = 3840, nseg = 10, rps = 10, L = 200, R = 225'); i= i+1;
% % nfib = 1920, nseg = 20, rps = 10, L = 400, R = 211
% AR_function('1920run002.txt',color{i})
% legendArr{i} = ('nfib = 1920, nseg = 20, rps = 10, L = 400, R = 211'); i= i+1;
% % nfib = 1280, nseg = 30, rps = 10, L = 600, R = 206 
% AR_function('1280run002.txt',color{i})
% legendArr{i} = ('nfib = 1280, nseg = 30, rps = 10, L = 600, R = 206'); i= i+1;

% nfib = 3840, nseg = 10, rps = 10, L = 200, R = 200 
AR_function('3840run003.txt',color{i})
legendArr{i} = ('nfib = 3840, nseg = 10, rps = 10, L = 200, R = 200'); i= i+1;
% nfib = 1280, nseg = 30, rps = 10, L = 600, R = 200 
AR_function('1280run003.txt',color{i})
legendArr{i} = ('nfib = 1280, nseg = 30, rps = 10, L = 600, R = 200'); i= i+1;


%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 24;
linewidth = 3;

%% Figure formatting
for i=1:14
    figure(i)
    hold on
    box on
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    xlim([0 inf])
    legend(legendArr{:},'location','best')
end