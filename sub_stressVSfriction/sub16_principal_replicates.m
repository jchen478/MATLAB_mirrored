%%%
%%% Plot transient percolation results
%%%

clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
%%%%%%%%%%%%%%%%%%%%% U-shaped fibers %%%%%%%%%%%%%%%%%%%%%
% fileNameArr = {'theta0'}; thetaArr = 0;
% fileNameArr = {'theta1'}; thetaArr = 1;
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
%{
%%%%%%%%%%%%%%%%%%%%% Helical fibers %%%%%%%%%%%%%%%%%%%%%
nfibArr = [160 240  320 640 1280 3200 6400];
lboxArr = [300 343.4 378 476.2 600 814.3  1026];
muArr = [0 1 2 3 4 5 10 15 20];

thetaArr = [3];
fileNameArr = {'helical'};
rpFiber = 75;
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed')...
    rgb('Orange') rgb('Gold') rgb('Lime')...
    rgb('DarkGreen') rgb('LightSkyBlue') rgb('Plum')};
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

% replicate info
nReplicate = 2;
replicate_flag = zeros(nNfib,1);

figStart = 1;

dataFile = cell(2,1);
for ii=1:2
    dataFile{ii} = ['fig16_principal_replicates/',fileNameArr{1},'_nfib'];
    for j=1:nNfib
        dataFile{ii} = [ dataFile{ii},num2str(nfibArr(j)),'_'];
    end
end
dataFile{1} = [ dataFile{1},'principal_direction3'];

% sampling info
sample_strain = 1000:25:1500;
nSample = length(sample_strain); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot principal directions vs. strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = '../data_stressVSfriction/ClusterP/';

Np_avg = zeros(nMu, nNfib);
Np_std = zeros(nMu, nNfib);

for j=1:nNfib
    
    N = zeros(nMu, nReplicate+1); 
    
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_clusterp_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File=fopen(name,'r');
        data = fscanf(File,'%f',[35 Inf])';
        fclose(File);
        
        % eliminate non-percolating structures
        sidexo = data(:,5);
        mle3 = data(:,35);
        keep = (mle3 >= sqrt(3)*sidexo);
        
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
        
        % obtain number of percolating structures
        allN = zeros(nSample,1);
        for ii=1:nSample
            tstrain = sample_strain(ii);
            allN(ii) = sum(strain == tstrain);
        end
        N(i,1) = mean(allN); 
        
        for r=1:nReplicate
            replicate_name = ['../data_stressVSfriction/ClusterP_replicate',...
                num2str(r),'/',fileNameArr{1},'_clusterp_nfib',...
                num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];

            if exist(replicate_name, 'file') == 2
                
                % if replicates exist, open file ...
                replicate_flag(j) = r;
                File=fopen(replicate_name,'r');
                data = fscanf(File,'%f',[35 Inf])';
                fclose(File);             
                
                % eliminate non-percolating structures ...
                sidexo = data(:,5);
                mle3 = data(:,35);
                keep = (mle3 >= sqrt(3)*sidexo);
                
                % delete empty rows ...
                data = data.*keep;
                data( ~any(data,2),:) = [];
                
                % find strain info ...
                strain = data(:,3);
                
                % obtain percolation results
                allN = zeros(nSample,1);
                for ii=1:nSample
                    tstrain = sample_strain(ii);
                    allN(ii) = sum(strain == tstrain);
                end
                N(i,1+r) = mean(allN);
            end
        end
        
        
    end
    % average over replicates
    average_range = replicate_flag(j) + 1; 
    Np_avg(:,j) = mean(N(:,1:average_range),2);
    Np_std(:,j) = std(N(:,1:average_range),0,2);
end

se_N = zeros(nMu,nNfib);
figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])

ylabel('$<N_p>$')
hold on
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,Np_avg(:,j),...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    else
        errorbar(muArr,Np_avg(:,j),...
            Np_std(:,j)/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
    end
end
xlabel('$\mu$')
xlim([0 inf])

legend(thetaNfibLegendArr,'location','bestoutside')

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.6])
ylabel('$<N_p> L/L_{box}$')
hold on
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,Np_avg(:,j)/lboxArr(j)*2*rpFiber,...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    else
        errorbar(muArr,Np_avg(:,j)/lboxArr(j)*2*rpFiber,...
            Np_std(:,j)/2/lboxArr(j)*2*rpFiber,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
    end
end
xlabel('$\mu$')
xlim([0 inf])
legend(thetaNfibLegendArr,'location','bestoutside')
