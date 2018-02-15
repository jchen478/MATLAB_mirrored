%%%
%%% Plot transient percolation results
%%%

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%%%%%%%%%%%%%%%%%%%% U-shaped fibers %%%%%%%%%%%%%%%%%%%%%
% fileNameArr = {'theta0'};
% fileNameArr = {'theta1'};
% fileNameArr = {'theta3'}; thetaArr = 3;
fileNameArr = {'theta6'}; thetaArr = 6;

nfibArr = [160 240 320 640 1280 3200 6400 10240 12800];
lboxArr = [300 343.4 378 476.2 600 814.3 1026 1200 1293];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];

rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('Olive') rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue') rgb('Plum') rgb('Purple') };
%}
%%{
%%%%%%%%%%%%%%%%%%%%% Helical fibers %%%%%%%%%%%%%%%%%%%%% 
nfibArr = [160 240  320 640 1280 3200 6400];
lboxArr = [300 343.4 378 476.2 600 814.3  1026];
muArr = [0 1 2 3 4 5 10 15 20];

nfibArr = [3200 6400];
lboxArr = [814.3  1026];
muArr = [5 10 20];

thetaArr = [3];
fileNameArr = {'helical'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('DarkGreen') rgb('LightSkyBlue') rgb('Plum')};


colorArr = {rgb('Lime')...
    rgb('DarkGreen') rgb('Plum')};
%}

nTheta = length(thetaArr);
nMu = length(muArr);
nLbox = length(lboxArr);
nNfib = length(nfibArr);
muLegendArr = cell(nMu,1);
thetaNfibLegendArr = cell(nTheta*nNfib,1);
markersize = 50;

for i=1:nMu
    muLegendArr{i} = ['$\mu =$ ',num2str(muArr(i))];
end
if strcmpi(fileNameArr,'helical')
    for i=1:nTheta
        for j=1:nNfib
            thetaNfibLegendArr{(i-1)*nNfib+j} = ['$(\theta_{eq},\phi_{eq},N_{fib}) =$ (0.8, 0.7, ',num2str(nfibArr(j)),')'];
        end
    end
else
    for i=1:nTheta
        for j=1:nNfib
            thetaNfibLegendArr{(i-1)*nNfib+j} = ['$(\theta_{eq},N_{fib}) =$ (0.',num2str(thetaArr(i)),', ',num2str(nfibArr(j)),')'];
        end
    end
end
figStart = 1;

dataFile = cell(2,1);
for ii=1:2
    dataFile{ii} = ['fig14_principal_vs_strain/',fileNameArr{1},'_nfib'];
    for j=1:nNfib
        dataFile{ii} = [ dataFile{ii},num2str(nfibArr(j)),'_'];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot principal directions vs. strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = '../data_stressVSfriction/ClusterP/';
N = zeros(nMu, nNfib); % number of flocs meeting criteria

for j=1:nNfib
%     figure('units','normalized','outerposition',[0.2 0.2 0.25 0.4])
figure('units','normalized','outerposition',[0.1 0.1 0.4 0.5])
    hold on
    title([fileNameArr{1}, ' ', num2str(nfibArr(j))]);
    for i=1:nMu
        
        name = [dataPath,fileNameArr{1},'_clusterp_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File=fopen(name,'r');
        data = fscanf(File,'%f',[35 Inf])';
        fclose(File);
        
        keep = ones(size(data,1),1);
        % condition for elimination
        % 1. Percolation
        %%{
        sidexo = data(:,5);
        mle3 = data(:,35);
        keep = (mle3 >= sqrt(3)*sidexo);
        %}
        % 2. Number of fibers in cluster
        %{
        nC = data(:,10);
        keep = (nC >= 10);
        %}
        % 3. Normalized number of fibers in cluster
        %{
        nfib = data(:,1);
        nC = data(:,10);
        keep = (nC >= 0.01*nfib);
        %}
        % 4. Radius of gyration
        %{
        Rg = data(:,32);
        sidexo = data(:,5);
        keep = (Rg >= 0.6*sidexo);
        %}
        
        % delete empty rows
        data = data.*keep;
        data( ~any(data,2),:) = [];
        
        % define variables
        nfib = data(:,1);
        mu = data(:,2);
        strain = data(:,3);
        fcut = data(:,4);
        sidexo = data(:,5);
        sidexn = data(:,6);
        dx = data(:,7);
        nx = data(:,8);
        cID = data(:,9);
        nC = data(:,10);
        rcmx = data(:,11);
        rcmy = data(:,12);
        rcmz = data(:,13);
        Rxx = data(:,14);
        Rxy = data(:,15);
        Rxz = data(:,16);
        Ryy = data(:,17);
        Ryz = data(:,18);
        Rzz = data(:,19);
        x1 = data(:,20);
        y1 = data(:,21);
        z1 = data(:,22);
        x2 = data(:,23);
        y2 = data(:,24);
        z2 = data(:,25);
        x3 = data(:,26);
        y3 = data(:,27);
        z3 = data(:,28);
        lam1 = data(:,29);
        lam2 = data(:,30);
        lam3 = data(:,31);
        Rg = data(:,32);
        mle1 = data(:,33);
        mle2 = data(:,34);
        mle3 = data(:,35);
        
        % calculate invariants of gyration tensor
        I1 = (lam1+lam2+lam3)./lam3; 
        I2 = (lam1.*lam2+lam2.*lam3+lam1.*lam3)./(lam3.^2);
        I3 = (lam1.*lam2.*lam3)./lam3.^3; 
        
        % shape measurements
        b = lam3.^2 - 0.5*(lam2.^2+lam1.^2); %asphericity 
        c = lam2.^2 - lam1.^2; % acylindricity
        k2 = 3/2*(lam1.^4+lam2.^4+lam3.^4)./(lam1.^2+lam2.^2+lam3.^2).^2-0.5; 
        % relative shape anisotropy
        
        
        N(i,j) = length(lam1);
        
        % different plots
        % 1. Third principal axis
        %{
        scatter3(abs(x3),abs(y3),abs(z3),'filled')
        box on
        xlabel('$\nu_{3,x}$ (flow)')
        ylabel('$\nu_{3,y}$ (vorticity)')
        zlabel('$\nu_{3,z}$ (gradient)')
        %}
        % 2. mle3 * Third principal axis
        %{       
        scatter3(mle3.*abs(x3),mle3.*abs(y3),mle3.*abs(z3),'filled')
        box on
        xlabel('$\nu_{3,x}MLE_3$ (flow)')
        ylabel('$\nu_{3,y}MLE_3$ (vorticity)')
        zlabel('$\nu_{3,z}MLE_3$ (gradient)')
        %}
        % 3. mle3
        %{
        scatter(mu,mle3,'filled')
        xlabel('$\mu$')
        ylabel('$MLE_3$')
        %}
        % 4. mle2
        %{
        scatter(mu,mle2,'filled')
        xlabel('$\mu$')
        ylabel('$MLE_2$')
        %}
        % 5. Rg
        %{
        scatter(mu,Rg,'filled')
        xlabel('$\mu$')
        ylabel('$R_g$')
        %}
        % 6. I1
        %{
        scatter(mu,I1,'filled')
        xlabel('$\mu$')
        ylabel('$I_1$')
        %}
        % 7. I2
        %{
        scatter(mu,I2,'filled')
        xlabel('$\mu$')
        ylabel('$I_2$')
        %}
        % 8. I3
        %{
        scatter(mu,I3,'filled')
        xlabel('$\mu$')
        ylabel('$I_3$')
        %}
        % 9. strain vs. data
        %{
        scatter(strain,c./(3*sidexo).^2,markersize,...
            'MarkerFaceColor',colorArr{i},...
            'MarkerEdgeColor',colorArr{i})
        xlabel('$\gamma$')
        ylabel('$c/(3L_{box})^2$')
        xlim([1000 1500])
        %}
        scatter(strain,nC,markersize,...
            'MarkerFaceColor',colorArr{i},...
            'MarkerEdgeColor',colorArr{i})
        xlabel('$\gamma$')
        ylabel('$N_{fib}$ in floc')
    end
    legend(muLegendArr,'location','bestoutside')
end

se_N = zeros(nMu,nNfib); 
figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('Number of flocs meeting criteria')
ylabel('$N$')
hold on
for j=1:nNfib
    plot(muArr,N(:,j),...
        '-.o','MarkerSize',10, 'linewidth',2.5);
end
xlabel('$\mu$')
xlim([0 inf])
legend(thetaNfibLegendArr,'location','bestoutside')

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
title('Normalized number of flocs meeting criteria')
ylabel('$N/L_{box}$')
hold on
for j=1:nNfib
    plot(muArr,N(:,j)/lboxArr(j),...
        '-.o','MarkerSize',10, 'linewidth',2.5);
end
xlabel('$\mu$')
xlim([0 inf])
legend(thetaNfibLegendArr,'location','bestoutside')

C1 = load('sub1_theta6');