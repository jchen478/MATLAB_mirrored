clc;
close all;

nfib = 1280;
L = 150; 

%% directly read in Lbox and t
File = fopen('Box_size.txt','r');
data = fscanf(File,'%f',[4 Inf])';
fclose(File);

t = data(:,1);
Lbox = data(:,2);

%% calculate Lbox from parameters
% n = 4000; 
% t = linspace(0,2000,n);
% 
% Lx_init = 600;
% fconc = 0.94;
% fexpd = 1.0804;
% nconc = 10;
% nexpd = 8;
% eqconc = 50;
% eqexpd = 100;
% nL3 = zeros(n,1);
% Lbox = zeros(n,1);
% 
% for i=1:n
%     
%     if t(i) < 50
%         Lbox(i) = Lx_init;
%     elseif t(i) >= 50 && t(i) < 550
%         power = floor((t(i)) / eqconc);
%         Lbox(i) = Lx_init*fconc^power;
%     elseif t(i) >= 550 && t(i) < 850
%         Lbox(i) = Lx_init*fconc^nconc;
%     elseif t(i) >= 850 && t(i) < 1650
%         power = ceil((t(i) - 850) / eqexpd);
%         Lbox(i) = Lx_init*fconc^nconc*fexpd^power;
%     elseif t(i) >= 1650
%         Lbox(i) = Lx_init*fconc^nconc*fexpd^nexpd;
%     end
%     
% end


%% nL3 calculation and plot
nL3 = nfib./Lbox.^3*L^3;

figure('color','white')
plot(t,nL3,'color',rgb('MediumBlue'),'linewidth',2);
xlabel('\gamma')
ylabel('nL3')
set(gca,'FontName','Times New Roman','fontsize',36)
box on
set(gca,'YMinorTick','on','XMinorTick','on','linewidth',2)