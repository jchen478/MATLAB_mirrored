close all;
clc;
clear

color = {rgb('OrangeRed'),rgb('Gold'),rgb('Black'),...
    rgb('HotPink'),rgb('Brown'),rgb('MediumBlue')};

color = {rgb('MediumBlue'),rgb('OrangeRed'),rgb('Purple'),...
    rgb('DarkGreen'),rgb('Crimson'),rgb('Orange')};

legendArr = cell(1);
i = 1;

%% cases for comparison
% % L = 400, N_{fib} = 1920, N_{seg} = 5, rps = 40
% AR_elastic('run007_elastic.txt',color{4})
% legendArr{i} = ('\it{}L\rm{} = 400, \it{}N_{fib}\rm{} = 1920, \it{}N_{seg}\rm{} = 5, \it{}rps\rm{} = 40, \it{}R_U\rm{} = 200'); i= i+1;
%
% % L = 800, N_{fib} = 2048, N_{seg} = 10, rps = 40
% AR_elastic('run008_elastic.txt',color{5})
% legendArr{i} = ('\it{}L\rm{} = 800, \it{}N_{fib}\rm{} = 2048, \it{}N_{seg}\rm{} = 10, \it{}rps\rm{} = 40, \it{}R_U\rm{} = 200'); i= i+1;

%% Comparison of Fiber Curvature
% L = 400, N_{fib} = 1920, N_{seg} = 20, rps = 10
AR_elastic('run004_elastic.txt',rgb('Crimson'))
legendArr{i} = ('\it{}L\rm{} = 400, \it{}N_{fib}\rm{} = 1920, \it{}N_{seg}\rm{} = 20, \it{}rps\rm{} = 10, \it{}R_U\rm{} = 200'); i= i+1;

% % L = 400, N_{fib} = 1920, N_{seg} = 10, rps = 20
% AR_elastic('run006_elastic.txt',rgb('Purple'))
% legendArr{i} = ('\it{}L\rm{} = 400, \it{}N_{fib}\rm{} = 1920, \it{}N_{seg}\rm{} = 10, \it{}rps\rm{} = 20, \it{}R_U\rm{} = 200'); i= i+1;
%
% % L = 400, N_{fib} = 1920, N_{seg} = 20, rps = 10, RU = 400
% AR_elastic('run010_elastic.txt',rgb('Crimson'))
% legendArr{i} = ('\it{}L\rm{} = 400, \it{}N_{fib}\rm{} = 1920, \it{}N_{seg}\rm{} = 20, \it{}rps\rm{} = 10, \it{}R_U\rm = 400'); i= i+1;
%
% % L = 400, N_{fib} = 1920, N_{seg} = 10, rps = 20, RU = 400
% AR_elastic('run011_elastic.txt',rgb('DarkGreen'))
% legendArr{i} = ('\it{}L\rm{} = 400, \it{}N_{fib}\rm{} = 1920, \it{}N_{seg}\rm{} = 10, \it{}rps\rm{} = 20, \it{}R_U\rm = 400'); i= i+1;

%% Comparison of Fiber Aspect Ratio
% % L = 800, N_{fib} = 2048, N_{seg} = 20, rps = 20
% AR_elastic('run009_elastic.txt',rgb('Purple'))
% legendArr{i} = ('\it{}L\rm{} = 800, \it{}N_{fib}\rm{} = 2048, \it{}N_{seg}\rm{} = 20, \it{}rps\rm{} = 20, \it{}R_U\rm{} = 200'); i= i+1;
%
% % L = 400, N_{fib} = 1920, N_{seg} = 20, rps = 10
% AR_elastic('run004_elastic.txt',rgb('MediumBlue'))
% legendArr{i} = ('\it{}L\rm{} = 400, \it{}N_{fib}\rm{} = 1920, \it{}N_{seg}\rm{} = 20, \it{}rps\rm{} = 10, \it{}R_U\rm{} = 200'); i= i+1;
%
% % L = 400, N_{fib} = 1920, N_{seg} = 20, rps = 10, run 2
% AR_elastic('run005_elastic.txt',rgb('OrangeRed'))
% legendArr{i} = ('\it{}L\rm{} = 400, \it{}N_{fib}\rm{} = 1920, \it{}N_{seg}\rm{} = 20, \it{}rps\rm{} = 10, \it{}R_U\rm{} = 200, run 2'); i= i+1;
%
% % L = 200, N_{fib} = 1024, N_{seg} = 10, rps = 10, RU = 200
% AR_elastic('run012_elastic.txt',rgb('Crimson'))
% legendArr{i} = ('\it{}L\rm{} = 200, \it{}N_{fib}\rm{} = 1024, \it{}N_{seg}\rm{} = 10, \it{}rps\rm{} = 10, \it{}R_U\rm{} = 200'); i= i+1;
%
% % L = 100, N_{fib} = 320, N_{seg} = 5, rps = 10, RU = 200
% AR_elastic('run013_elastic.txt',rgb('DarkGreen'))
% legendArr{i} = ('\it{}L\rm{} = 100, \it{}N_{fib}\rm{} = 320, \it{}N_{seg}\rm{} = 5, \it{}rps\rm{} = 10, \it{}R_U\rm{} = 200'); i= i+1;

%% Relaxation
% % L = 400, N_{fib} = 1920, N_{seg} = 20, rps = 10
% AR_elastic('run004_elastic.txt',rgb('MediumBlue'))
% legendArr{i} = ('Under shear: \it{}L\rm{} = 400, \it{}N_{fib}\rm{} = 1920, \it{}N_{seg}\rm{} = 20, \it{}rps\rm{} = 10, \it{}R_U\rm{} = 200'); i= i+1;

% % L = 400, N_{fib} = 1920, N_{seg} = 20, rps = 10
% AR_elastic('run004_relaxed_elastic.txt',rgb('Crimson'))
% legendArr{i} = ('Relaxed:      \it{}L\rm{} = 400, \it{}N_{fib}\rm{} = 1920, \it{}N_{seg}\rm{} = 20, \it{}rps\rm{} = 10, \it{}R_U\rm{} = 200'); i= i+1;

%% Correlation
% correlation
% AR_elastic('correlation_elastic.txt',rgb('MediumBlue'))
% legendArr{i} = ('\it{}L\rm{} = 400, \it{}N_{fib}\rm{} = 1920, \it{}N_{seg}\rm{} = 20, \it{}rps\rm{} = 10, \it{}R_U\rm{} = 200'); i= i+1;

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 26;
linewidth = 2.5;

%% Figure formatting
for i=1:5
    figure(i)
    if (i < 3)
        for j = 1:2
            
            subplot(2,1,j)
            hold on
            box on
            set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
            set(gca,'YMinorTick','on','XMinorTick','on')
            set(gca,'ticklength',tickx*get(gca,'ticklength'))
            set(gcf, 'color','white')
            xlim([0 inf])
            legend(legendArr{:},'location','best')
        end
    else
        hold on
        box on
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(gca,'ticklength',tickx*get(gca,'ticklength'))
        set(gcf, 'color','white')
        xlim([0 inf])
        legend(legendArr{:},'location','best')
    end
    
end

%% close figures
% close(1) % 1/Seff vs eta_sp
% close(2) % 1/Seff vs eta_rel
% close(3) % 1/Seff vs stress
% close(4) % 1/Seff vs I
% close(5) % 1/Seff vs Sx
% close(6)% 1/Seff vs Sy
% close(7) % 1/Seff vs Sz
% close(8) % 1/Seff vs N1
% close(9) % 1/Seff vs N2
% close(10) % scale2 vs eta_rel
% close(11) % scale2 vs stress
% close(12) % scale2 vs I
% close(13) % scale2 vs Sy
% close(14) % scale2 vs eta_sp
