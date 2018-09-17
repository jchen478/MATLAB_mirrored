clc; 
close all;
clear;

%% Interfiber normal forces parameters
fstar = 150; 
fact = 20;
Astar = [0, 9, 20, 30, 35, 50];
decatt = 35;
rep_cut = 0.66; 
overlap = 1.85; 

eta0 = 1; % Pa s
L = 2.50E-3; % m 
b = 1.6E-5; % m
gamma = 10; % 1/s
nseg = 5;
l = L/(2*nseg);

prefactor = -6*pi*l*b*gamma*eta0*10^6 % muN

g = (1.85:0.01:2.66)'; % interfiber distance from center of fiber
h = g-2; % interfiber distance from fiber surface

nH = length(h);
nAstar = length(Astar);

Fatt = zeros(nH,nAstar);
Ftot = zeros(nH,nAstar);
Uval = zeros(nH,nAstar);

maxFtot = zeros(nAstar,2);
zeroFtot = zeros(nAstar,2);

AstarLegend = cell(nAstar,1);
for i=1:nAstar
    AstarLegend{i} = ['$A_N =$ ', num2str(Astar(i))];
end
FtotAstarLegend = cell(nAstar,1);
for i=1:nAstar
    FtotAstarLegend{i} = ['$F_{tot} (A_N =$ ',...
        num2str(Astar(i)),')'];
end
%% Interfiber normal forces calculations
Frep = prefactor*fstar.*exp(-fact*h); 

%% Interfiber maximum forces
syms sij a b c
digits(5)
for i=1:nAstar
    a = prefactor*fstar.*exp(-fact*sij); 
    b = piecewise(sij<0, prefactor*(-1)*Astar(i),sij>=0,prefactor*(-1)*Astar(i)*exp(-decatt*sij.^2));
    Fatt(:,i) = vpa(subs(b,h));
    f = a+b;
    c = a+prefactor*(-1)*Astar(i)*exp(-decatt*sij.^2);
    g = diff(c,sij);
    S = vpasolve(g==0, sij);
    if (i ~= 1)
        maxFtot(i,1) = S;
        maxFtot(i,2) = vpa(subs(f,S));
    end
    S = vpasolve(c==0, sij);
    if (i ~= 1)
        zeroFtot(i,1) = S;
        zeroFtot(i,2) = vpa(subs(f,S));
    end
    Ftot(:,i) = Frep+Fatt(:,i);
    U = -int(f,sij);
    Uval(:,i) = vpa(subs(U,h));
end

%% Plotting Components
AstarColorArr = {rgb('Crimson'), rgb('DarkOrange'), ...
    rgb('DarkGreen'), rgb('MediumBlue'), rgb('Purple'), rgb('Maroon')};

figure('Units','Inches','Position',[1 1 4.8 3.0]);
title('\textbf{Components of normal force}')
box on
hold on
for i=1:nAstar
    plot(h,Fatt(:,i),'linestyle','-','color',AstarColorArr{i})
end
legend(AstarLegend,'Location','best')
ylabel('$F_{att}\ (\mu N)$')
xlabel('$h/b$')
xlim([-0.2 h(end)])
ylim([-0.5 40])

yyaxis right
ylabel('$F_{rep}\ (\mu N)$')
set(gca,'YColor',rgb('MediumVioletRed'))
hold on
plot(h,Frep,'linestyle','-.','color',rgb('MediumVioletRed'))
ylim([-Inf 20])

%% Plotting Total with Maximums 
fig = figure('Units','Inches','Position',[1 1 4.8 3.7]);
title('\textbf{Total normal force}')
box on
hold on
ylabel('$F_{tot}\ (\mu N)$')
xlabel('$h/b$')
xlim([-0.2 h(end)])
ylim([-inf 50])
for i=1:nAstar
    plot(h,Ftot(:,i),'linestyle','-','color',AstarColorArr{i})
end
legend(FtotAstarLegend,'Location','northeast')
for i=2:nAstar
    plot(maxFtot(i,1),maxFtot(i,2),'o','color',AstarColorArr{i})
end


% create a new pair of axes inside current figure
p = get(gca, 'Position');
handle = axes('Parent', gcf,'Units','Inches',...
    'Position', [1.5 1.0 2 1]);
box on 
hold on
for i=1:nAstar
    plot(h,Ftot(:,i),'linestyle','-','color',AstarColorArr{i}) 
end
for i=2:nAstar
    plot(maxFtot(i,1),maxFtot(i,2),'o','color',AstarColorArr{i})
end
set(handle, 'Xlim', [0.1 0.3], 'Ylim', [-1 15]);

%% plot potential energy
figure('Units','Inches','Position',[1 1 5.0 3.0]);
title('\textbf{Interfiber potential}')
box on
hold on
for i=1:nAstar
    plot(h,-Uval(:,i),'linestyle','-','color',AstarColorArr{i})
end
legend(AstarLegend,'Location','best')
ylabel('$-U = \int F d(h/b)$')
xlabel('$h/b$')
xlim([-0.2 h(end)])
set(gca,'yscale','log')

%% Plotting Total with zero total force  
fig = figure('Units','Inches','Position',[1 1 4.8 3.7]);
title('\textbf{Total normal force}')
box on
hold on
ylabel('$F_{tot}\ (\mu N)$')
xlabel('$h/b$')
xlim([-0.2 h(end)])
ylim([-inf 50])
for i=1:nAstar
    plot(h,Ftot(:,i),'linestyle','-','color',AstarColorArr{i})
end
legend(FtotAstarLegend,'Location','northeast')
for i=2:nAstar
    plot(zeroFtot(i,1),zeroFtot(i,2),'o','color',AstarColorArr{i})
end


% create a new pair of axes inside current figure
p = get(gca, 'Position');
handle = axes('Parent', gcf,'Units','Inches',...
    'Position', [1.5 1.0 2 1]);
box on 
hold on
for i=1:nAstar
    plot(h,Ftot(:,i),'linestyle','-','color',AstarColorArr{i}) 
end
for i=2:nAstar
    plot(zeroFtot(i,1),zeroFtot(i,2),'o','color',AstarColorArr{i})
end
set(handle, 'Xlim', [0.04 0.28], 'Ylim', [-5 5]);

%% Plotting Total with two insets total force  

fig = figure('Units','Inches','Position',[1 1 7 4]);
title('\textbf{Total normal force}')
box on
hold on
ylabel('$F_{tot}\ (\mu N)$')
xlabel('$h/b$')
xlim([-0.2 h(end)])
ylim([-inf 50])
for i=1:nAstar
    plot(h,Ftot(:,i),'linestyle','-','color',AstarColorArr{i})
end
legend(FtotAstarLegend,'Location','northeast')
text(1.55,1.6,'\bf{Max} \boldmath$F_{att}$','Units','Inches',...
    'FontSize',14,'FontName','CMU Serif','Interpreter','Latex')
text(3.9,1.6,'\boldmath$F_{tot} = 0$','Units','Inches',...
    'FontSize',14,'FontName','CMU Serif','Interpreter','Latex')

% create a new pair of axes inside current figure
p = get(gca, 'Position');
handle = axes('Parent', gcf,'Units','Inches',...
    'Position', [4.1 1.0 2 1]);
box on 
hold on
for i=1:nAstar
    plot(h,Ftot(:,i),'linestyle','-','color',AstarColorArr{i}) 
end
for i=2:nAstar
    plot(zeroFtot(i,1),zeroFtot(i,2),'o','color',AstarColorArr{i})
end
set(handle, 'Xlim', [0.04 0.28], 'Ylim', [-5 5]);

% create a new pair of axes inside current figure
p = get(gca, 'Position');
handle = axes('Parent', gcf,'Units','Inches',...
    'Position', [1.75 1.0 2 1]);
box on 
hold on
for i=1:nAstar
    plot(h,Ftot(:,i),'linestyle','-','color',AstarColorArr{i}) 
end
for i=2:nAstar
    plot(maxFtot(i,1),maxFtot(i,2),'o','color',AstarColorArr{i})
end
set(handle, 'Xlim', [0.1 0.3], 'Ylim', [-1 15]);

%%
figure('Units','Inches','Position',[1 1 3.5 3.5]);
hold on
box on 
hold on
title('\boldmath$F_{tot} = 0$')
for i=1:nAstar
    plot(h,Ftot(:,i),'linestyle','-','color',AstarColorArr{i}) 
end
legend(FtotAstarLegend,'Location','bestoutside')
for i=2:nAstar
    plot(zeroFtot(i,1),zeroFtot(i,2),'o','color',AstarColorArr{i})
end
xlim([0.04 0.28])
ylim([-5 5])
ylabel('$F_{N}\ (\mu N)$')
xlabel('$h/b$')