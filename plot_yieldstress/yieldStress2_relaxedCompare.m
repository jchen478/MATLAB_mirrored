clc;
clear;
close all;
% file = {'yieldStress2_1280_phi1.txt','yieldStress2_relaxed_1280_phi1.txt'};
% file = {'yieldStress2_1280_phi2.txt','yieldStress2_relaxed_1280_phi2.txt'};
% file = {'yieldStress2_1280_phi1.txt','yieldStress2_relaxed_1280_phi1.txt','yieldStress2_1280_phi2.txt','yieldStress2_relaxed_1280_phi2.txt'};
file = {'yieldStress2_1280_phi1.txt','yieldStress2_1280_phi2.txt','yieldStress2_1280_phi4.txt'};
% color = {rgb('Maroon'),rgb('DeepPink'),rgb('DarkOrange'),rgb('DarkGreen'),rgb('MediumBlue'),rgb('DarkMagenta')};
color = {rgb('MediumBlue'),rgb('DarkOrange'),rgb('DarkGreen'),rgb('Crimson')};
% legendArr = {'Under shear', 'Relaxed'}; 
legendArr = {'Case 1: under shear', 'Case 1: relaxed','Case 2: under shear', 'Case 2: relaxed'}; 

nfit = 3;
nstiff = 6;
yieldStress = zeros(length(file)*nfit);

for i=1:length(file)
    relax_yield( file{i}, color{i})
end

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 20;
linewidth = 2;

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 26;
linewidth = 2.5;

%% Figure formatting
for i=1:7
    figure(i)
    hold on
    box on
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    xlim([0 inf])
    xlabel('1/\it{}S_{eff}')
    legend(legendArr{:},'location','best')
end

figure(6)
hold on
legendArr = {'Case 1: under shear', 'Case 1: under shear fit', ...
    'Case 2: under shear', 'Case 2: under shear fit'}; 
legend(legendArr{:},'location','best')
%     set(gca,'yscale','log')
    title('Bingham model fit')
    
    figure(7)
hold on
legendArr = {'Case 1: under shear', ...
    'Case 2: under shear'}; 
legend(legendArr{:},'location','best')
%     set(gca,'yscale','log')
    title('Bingham model fit')
% end

% yieldStress
spreadfigures();