clc;
clear;
close all;

% non-flocculated case
% file = {'yieldStress2_1280_phi4.txt','yieldStress2_1280_phi2.txt','yieldStress2_1280_phi1.txt'};
% color = {rgb('Maroon'),rgb('Crimson'),rgb('DarkOrange'),rgb('DarkGreen'),rgb('MediumBlue'),rgb('DarkMagenta')};

% flocculated case
file = {'yieldStress2_1280_06_phi3.txt'};
color = {rgb('Maroon'),rgb('Crimson'),rgb('DarkOrange'),rgb('DarkGreen'),rgb('MediumBlue'),rgb('DarkMagenta')};

nfit = 3;
nstiff = 5; 
yieldStress = zeros(length(file)*nfit);

for i=1:length(file)
    
    yieldModel(file{i},i,yieldStress,nfit,nstiff,color{i});

end

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 20;
linewidth = 2;

%% figure formatting
for i=1:2:6
    figure(i)
    hold on
    box on
    xlabel('1/S_{eff}')
    ylabel('\it{\sigma_{xz} L^4/ E_Y I}')
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
%     legend('\phi = 0.0184','\phi = 0.0184 fit','\phi = 0.0132','\phi = 0.0132 fit','\phi = 0.0073','\phi = 0.0073 fit')
     legend('\phi = 0.003','\phi = 0.003 fit','\phi = 0.0132','\phi = 0.0132 fit','\phi = 0.0073','\phi = 0.0073 fit','location','best')
end
for i=2:2:6
    figure(i)
    hold on
    box on
    xlabel('1/S_{eff}')
    ylabel('RS')
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    legend('\phi = 0.003','\phi = 0.0132','\phi = 0.0073','location','best')
end

