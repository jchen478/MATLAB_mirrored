clc;
close all;

filename = {'nonCubicShape.txt'};
filename = {'nonCubicShape_run001.txt','nonCubicShape_run002.txt',...
    'nonCubicShape_run003.txt','nonCubicShape_run004.txt',...
    'nonCubicShape_run005.txt','nonCubicShape_run006.txt'};
ncase = length(filename);

colorArr = {rgb('Crimson'),rgb('MediumBlue'),rgb('Orange'),rgb('DarkGreen'),rgb('Purple'),rgb('SkyBlue')};
markersize = 50;
I = zeros(3,3);

fignum = 0;
ncase
%% read data
for i=1:ncase
    
    File = fopen(filename{i},'r');
    data = fscanf(File,'%f',[19 Inf])';
    fclose(File);
    
    %% Define variables
    nfib = data(:,1);
    mu = data(:,2);
    flocID = data(:,3);
    flocNfib = data(:,4);
    rps = data(:,5);
    nseg = data(:,6);
    sidex = data(:,8);
    sidey = data(:,9);
    sidez = data(:,10);
    Rgx = data(:,11);
    Rgy = data(:,12);
    Rgz = data(:,13);
    I1 = data(:,14);
    I2 = data(:,15);
    I3 = data(:,16);
    I5 = data(:,17);
    I6 = data(:,18);
    I9 = data(:,19);
    
    ndata = length(Rgx);
    lambda1 = zeros(ndata,1);
    lambda2 = zeros(ndata,1);
    lambda3 = zeros(ndata,1);
    
    for j=1:ndata
        I(1) = I1(j);
        I(2) = I2(j);
        I(3) = I3(j);
        I(4) = I(2);
        I(5) = I5(j);
        I(6) = I6(j);
        I(7) = I(3);
        I(8) = I(6);
        I(9) = I9(j);
        E = eig(I);
        lambda1(j) = E(1);
        lambda2(j) = E(2);
        lambda3(j) = E(3);
    end
    
    % average and invariants
    lambdaAvg = (lambda1+lambda2+lambda3)/3;
    Ii1 = (lambda1+lambda2+lambda3)./lambda3;
    Ii2 = (lambda1.*lambda2+lambda2.*lambda3+lambda1.*lambda3)./lambda3.^2;
    Ii3 = (lambda1.*lambda2.*lambda3)./lambda3.^3;
    
    % degree of prolateness
    DP = 27*(lambda1-lambdaAvg).*(lambda2-lambdaAvg)...
        .*(lambda3-lambdaAvg)./(lambda1+lambda2+lambda3).^3;
    
    fignum = fignum+1;
    %% Radius of gyration
    figure(fignum)
    subplot(3,2,1)
    hold on
    scatter(mu, Rgx, markersize,'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    line([1 5] , [sidex(1) sidex(1)], 'linewidth',2,'color',colorArr{i});
    ylabel('\it{}R_{g,x}')
    text(0.05,0.2,['\it{}L_x\rm{} = ',num2str(sidex(1))],'Units','normalized','fontsize',24,'fontname','Times New Roman')
    subplot(3,2,3)
    hold on
    scatter(mu, Rgy, markersize,'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    line([1 5] , [sidey(1) sidey(1)], 'linewidth',2,'color',colorArr{i});
    ylabel('\it{}R_{g,y}')
    text(0.05,0.2,['\it{}L_y\rm{} = ',num2str(sidey(1))],'Units','normalized','fontsize',24,'fontname','Times New Roman')
    subplot(3,2,5)
    hold on
    scatter(mu, Rgz, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    line([1 5] , [sidez(1) sidez(1)], 'linewidth',2,'color',colorArr{i});
    ylabel('\it{}R_{g,z}')
    text(0.05,0.2,['\it{}L_z\rm{} = ',num2str(sidez(1))],'Units','normalized','fontsize',24,'fontname','Times New Roman')
    subplot(3,2,2)
    hold on
    scatter(mu, 2*Rgx./sidex, markersize,'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('2\it{}R_{g,x} / L_x')
    subplot(3,2,4)
    hold on
    scatter(mu, 2*Rgy./sidey, markersize,'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('2\it{}R_{g,y} / L_y')
    subplot(3,2,6)
    hold on
    scatter(mu, 2*Rgz./sidez, markersize,'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('2\it{}R_{g,z} / L_z')
    
    fignum = fignum+1;
    %% fiber counting stats
    figure(fignum)
    subplot(2,1,1)
    hold on
    scatter(mu, flocNfib, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('\it{}N_{fib,floc}')
    title('Number of fibers in floc')
    text(0.05,0.2,['\it{}L_x, L_y, L_z \rm{} = ',num2str(sidex(1)), ', ',num2str(sidey(1)),', ',num2str(sidez(1))],'Units','normalized','fontsize',24,'fontname','Times New Roman')
    subplot(2,1,2)
    hold on
    scatter(mu, flocNfib./nfib, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    text(0.05,0.2,['\it{}L_x, L_y, L_z \rm{} = ',num2str(sidex(1)), ', ',num2str(sidey(1)),', ',num2str(sidez(1))],'Units','normalized','fontsize',24,'fontname','Times New Roman')
    ylabel('\it{}N_{fib,floc} / N_{fib,tot}')
    title('Normalized number of fibers in floc')
    
    fignum = fignum+1;
    %% degree of prolateness
    figure(fignum)
    hold on
    scatter(mu, DP, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('DP')
    text(0.05,0.2,['\it{}L_x, L_y, L_z \rm{} = ',num2str(sidex(1)), ', ',num2str(sidey(1)),', ',num2str(sidez(1))],'Units','normalized','fontsize',24,'fontname','Times New Roman')
    title('Degree of prolateness')
    
    fignum = fignum+1;
    %% invariants based on moment of inertia tensor
    figure(fignum)
    subplot(3,1,1)
    hold on
    scatter(mu, Ii1, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    text(0.05,0.2,['\it{}L_x, L_y, L_z \rm{} = ',num2str(sidex(1)), ', ',num2str(sidey(1)),', ',num2str(sidez(1))],'Units','normalized','fontsize',24,'fontname','Times New Roman')
    ylabel('\it{}I_1')
    title('Invariants based on inertia')
    subplot(3,1,2)
    hold on
    scatter(mu, Ii2, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    text(0.05,0.2,['\it{}L_x, L_y, L_z \rm{} = ',num2str(sidex(1)), ', ',num2str(sidey(1)),', ',num2str(sidez(1))],'Units','normalized','fontsize',24,'fontname','Times New Roman')
    ylabel('\it{}I_2')
%     title('Invariants based on inertia')
    subplot(3,1,3)
    hold on
    text(0.05,0.2,['\it{}L_x, L_y, L_z \rm{} = ',num2str(sidex(1)), ', ',num2str(sidey(1)),', ',num2str(sidez(1))],'Units','normalized','fontsize',24,'fontname','Times New Roman')
    scatter(mu, Ii3, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('\it{}I_3')
%     title('Invariants based on inertia')
    
    %% Figure formatting parameters
    nfig = 12;
    tickx = 1.5;
    fontsize = 24;
    linewidth = 3;
    
    %% Figure formatting
    figure(fignum)
    for j=1:3
        subplot(3,1,j)
        hold on
        box on
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(gca,'ticklength',tickx*get(gca,'ticklength'))
        set(gcf, 'color','white')
        xlabel('\mu')
        xlim([1 5])
    end
    figure(fignum-1)
    hold on
    box on
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    xlim([1 5])
    xlabel('\mu')
    line([1 5], [0 0], 'linewidth',2,'color',rgb('OrangeRed'));
    figure(fignum-2)
    for j=1:2
        subplot(2,1,j)
        hold on
        box on
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(gca,'ticklength',tickx*get(gca,'ticklength'))
        set(gcf, 'color','white')
        xlim([1 5])
        xlabel('\mu')
    end
    figure(fignum-3)
    for j=1:6
        subplot(3,2,j)
        hold on
        box on
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(gca,'ticklength',tickx*get(gca,'ticklength'))
        set(gcf, 'color','white')
        xlim([1 5])
        xlabel('\mu')
    end
    
    
end

