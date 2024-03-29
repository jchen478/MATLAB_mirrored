%%%
%%% Plot evolution of properties
%%% -- properties include nc, sigP
%%%

clc;
clear;
close all;

% Define common parameters
simulation_cases;

% plotting commands
xLowLim = 0;
xUpLim = 2000;

for j=1:nNfib
    figure('Units','Inches','Position',[3 3 6.5 6.5])
    set(gca,'XMinorTick','on')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NC

dataPath = '../data_stressVSfriction/NC/';
% dataPath = '../data_stressVSfriction/NC_replicate1/';
% dataPath = '../data_stressVSfriction/NC_replicate2/';
%%{
for j=1:nNfib
    figure(j)
    subplot(2,1,1)
    set(gca,'XMinorTick','on')
    hold on;
    ylabel('$N_c$')
    maxNC = 0;
    minNC = 0;
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_NC_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File = fopen(name,'r');
        data = fscanf(File,'%f',[5 Inf])';
        fclose(File);
        if (size(data,1) == 1)
            plot(0,0,'x')
            continue
        end      
        diff = data(2,1) - data(1,1);
        % make sure that strain is increasing
        for ii=2:length(data)
            if (data(ii,1) < data(ii-1,1))
                data(ii,1) = data(ii-1,1) + diff;
            end
        end
        % block average
        block = 25;
        npt = floor(length(data)/block);
        B = zeros(npt,1);
        Bt = zeros(npt,1);
        ind = 1;
        for ii=1:block:length(data)-block+1
            B(ind) = mean(data(ii:ii+block-1,4));
            Bt(ind) = mean(data(ii:ii+block-1,1));
            ind = ind + 1;
        end
        plot(Bt,B,'-o','color',colorArr{i},...
            'MarkerFaceColor',colorArr{i},...
            'MarkerEdgeColor',colorArr{i},...
            'MarkerSize',2)
        % find maximum and minimum
        if min(B) < minNC
            minNC = min(B);
        end
        if max(B) > maxNC
            maxNC = max(B);
        end
    end
    box on
    xlabel('$\gamma$')
    xlim([xLowLim xUpLim])
    legend(muLegendArr,'location','bestoutside')
    title([fileNameArr{1},' $N_{fib} =$ ',num2str(nfibArr(j))])
%     line([1400 1400],[minNC maxNC],'color',rgb('DeepPink'))
end

%}

% % Stress
dataPath = '../data_stressVSfriction/Stress/';
% dataPath = '../data_stressVSfriction/Stress_replicate1/';
% dataPath = '../data_stressVSfriction/Stress_replicate2/';

for j=1:nNfib
    figure(j)
    subplot(2,1,2)
    hold on;
    maxStress = 0;
    minStress = 0;
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_stress_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File = fopen(name,'r');
        data = fscanf(File,'%f',[7 Inf])';
        fclose(File);
        if (size(data,1) == 1)
            plot(0,0,'x')
            continue
        end
        diff = data(2,1) - data(1,1);
        nStep = length(data);
        % make sure that strain is increasing
        for ii=2:nStep
            if (data(ii,1) < data(ii-1,1))
                data(ii,1) = data(ii-1,1) + diff;
            end
        end
        
        nfib = nfibArr(j)*ones(nStep,1);
        nsegArr = nseg*ones(nStep,1);
        side = lboxArr(j)*ones(nStep,1);
        % calculate relevant parameters
        [Seff,rp,nL3,L,gamma] = pCalc(nfib, nsegArr, rps, kb, side, a, EY, Imom, eta0);
        
        % Dimensionalize stresses
        [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
            stressDim(data(:,4), zeros(nStep,1), zeros(nStep,1), ...
            zeros(nStep,1), nseg, rps, eta0, gamma, nL3);
        
        
        % Non-dimensionalize
        sigmap_nondim = sigmap.*L.^4/EY/Imom;
        strain = data(:,1);
        
        % block average
        block = 25;
        npt = floor(length(sigmap_nondim)/block);
        B = zeros(npt,1);
        Bt = zeros(npt,1);
        ind = 1;
        for ii=1:block:length(sigmap_nondim)-block+1
            B(ind) = mean(sigmap_nondim(ii:ii+block-1));
            Bt(ind) = mean(strain(ii:ii+block-1));
            ind = ind + 1;
        end
        plot(Bt,B,'-o','color',colorArr{i},...
            'MarkerFaceColor',colorArr{i},...
            'MarkerEdgeColor',colorArr{i},...
            'MarkerSize',2)
        if min(B) < minStress
            minStress = min(B);
        end
        if max(B) > maxStress
            maxStress = max(B);
        end
    end
    box on
    xlabel('$\gamma$')
    ylabel('$\sigma_{p,xz} L^4/ E_Y I$')
    xlim([xLowLim xUpLim])
    legend(muLegendArr,'location','bestoutside')
    title([fileNameArr{1},' $N_{fib} =$ ',num2str(nfibArr(j))])
%     line([1400 1400],[minStress maxStress],'color',rgb('DeepPink'))
end

