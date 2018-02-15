clc;
close all;

filename = {'nfib160_50_Shape.txt','nfib1280_50_Shape.txt'};
filename = {'nfib160_50_Shape.txt','nfib1280_50_Shape.txt','nfib12800_50_Shape.txt'};
legendArrWL = {'\itN_{fib}\rm = 160','Box length (\itN_{fib}\rm = 160)', ...
    '\itN_{fib}\rm = 1280','Box length (\itN_{fib}\rm = 1280)', ...
    '\itN_{fib}\rm = 12800','Box length (\itN_{fib}\rm = 12800)'};
legendArr = {'\itN_{fib}\rm = 160','\itN_{fib}\rm = 1280','\itN_{fib}\rm = 12800'};
ncase = length(filename);

colorArr = {rgb('Crimson'),rgb('MediumBlue'),rgb('Orange')};
markersize = 50;
I = zeros(3,3);

for i=1:ncase
    
    %% read files
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
    
    %% Radius of gyration
    figure(1)
    hold on
    scatter(mu, Rgx, markersize,'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    line([0 25], [sidex(1) sidex(1)], 'linewidth',2,'color',colorArr{i});
    ylabel('\it{}R_{g,x}')
    figure(2)
    hold on
    scatter(mu, Rgy, markersize,'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    line([0 25], [sidex(1) sidex(1)], 'linewidth',2,'color',colorArr{i});
    ylabel('\it{}R_{g,y}')
    figure(3)
    hold on
    scatter(mu, Rgz, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    line([0 25], [sidex(1) sidex(1)], 'linewidth',2,'color',colorArr{i});
    ylabel('\it{}R_{g,z}')
    figure(4)
    hold on
    scatter(mu, 2*Rgx./sidex, markersize,'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('2\it{}R_{g,x} / L_x')
    figure(5)
    hold on
    scatter(mu, 2*Rgy./sidey, markersize,'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('2\it{}R_{g,y} / L_y')
    figure(6)
    hold on
    scatter(mu, 2*Rgz./sidez, markersize,'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('2\it{}R_{g,z} / L_z')
    
    %% fiber counting stats
    figure(7)
    hold on
    scatter(mu, flocNfib, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('\it{}N_{fib,floc}')
    title('Number of fibers in floc')
    figure(8)
    hold on
    scatter(mu, flocNfib./nfib, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('\it{}N_{fib,floc} / N_{fib,tot}')
    title('Normalized number of fibers in floc')
    
    %% degree of prolateness
    figure(9)
    hold on
    scatter(mu, DP, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('DP')
    title('Degree of prolateness')
    
    %% invariants based on moment of inertia tensor
    figure(10)
    hold on
    scatter(mu, Ii1, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('\it{}I_1')
    title('Invariants based on inertia')
    figure(11)
    hold on
    scatter(mu, Ii2, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('\it{}I_2')
    title('Invariants based on inertia')
    figure(12)
    hold on
    scatter(mu, Ii3, markersize, 'filled','markerfacecolor',colorArr{i},'markeredgecolor',colorArr{i})
    ylabel('\it{}I_3')
    title('Invariants based on inertia')
    
end

%% Figure formatting parameters
nfig = 12;
tickx = 1.5;
fontsize = 24;
linewidth = 3;

%% Figure formatting
for i=1:nfig
    figure(i)
    hold on
    box on
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    xlabel('\mu')
    if i < 4
        legend(legendArrWL{:},'location','best')
    else
        legend(legendArr{:},'location','northoutside','orientation','horizontal')
    end
end

figure(9)
hold on
line([0 25], [0 0], 'linewidth',2,'color',rgb('OrangeRed'));