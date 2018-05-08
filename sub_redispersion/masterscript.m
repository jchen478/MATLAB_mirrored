%% analysis for redispersion
clc;
clear;
close all;

%% define simulation cases
simulation_cases;

%% For each parameter set
fig = figure('Units','Inches','Position',[1 1 4.0 6.0]);
%     set(gcf,'name',['mu_',num2str(muArr(i))]);
hold on
for i=1:nMu
    NC_stat_all = zeros(5,2,nAtt,nA);
    sigmap_nondim_stat_all = zeros(5,2,nAtt,nA);
    eta_stat_all = zeros(5,2,nAtt,nA);
    for k=1:nAtt
        for j=1:nA
            display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
            filePrefix = [dataPath,shape,'_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_'];
            %% read input and process data
            [stress_strain, sigxz_out] = read_stress([filePrefix,'Stress_tensor.txt']);
            [box_strain, sidex] = read_box([filePrefix,'Lbox.txt']);
            [NC_strain, NC] = read_NC([filePrefix,'Number_of_Contacts.txt']);
            
            % define lboxArr to calculate stress
            lboxArr = find_lboxArr(stress_strain, box_strain, sidex);
            %
            % find relevant parameters for stress calculations
            nPt = length(stress_strain);
            [Seff,rp,nL3,L,gamma] = pCalc(nfib*ones(nPt,1), nseg*ones(nPt,1), rps, kb, lboxArr, a, EY, Imom, eta0);
            
            % Dimensionalize stresses
            [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
                stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
                zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
            
            etarel = sigma./gamma;
            % Non-dimensionalize
            sigmap_nondim = sigmap.*L.^4/EY/Imom;
            
            %% calculate interval averages
            r = round(intervals(box_strain,sidex),0);
            
            % if simulation has not finished
            if round(stress_strain(end),1) ~= round(box_strain(end),1)
                continue;
            end
            
            NC_stat = interval_average(NC_strain,NC,r');
            sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r');
            eta_stat = interval_average(stress_strain,etarel,r');
            
            eta_stat_all(:,:,k,j) = eta_stat;
            NC_stat_all(:,:,k,j) = NC_stat;
            sigmap_nondim_stat_all(:,:,k,j) = sigmap_nondim_stat;
        end
    end
    stress_plot = squeeze(sigmap_nondim_stat_all(5,1,1,:));
    NC_plot = squeeze(NC_stat_all(5,1,1,:));
    eta_plot = squeeze(eta_stat_all(5,1,1,:));
    
    %% plot data to corresponding plots
    figure(1)
    subplot(2,1,1)
    hold on
    plot(aArr,eta_plot,'-o','Color',muColorArr{i})
    subplot(2,1,2)
    hold on
    plot(aArr,NC_plot,'-o','Color',muColorArr{i})
    
    fig2 = figure('Units','Inches','Position',[1 1 4.0 6.0]);
    for j=1:nA
        stress_plot = squeeze(sigmap_nondim_stat_all(5,1,:,j));
        NC_plot = squeeze(NC_stat_all(5,1,:,j));
        eta_plot = squeeze(eta_stat_all(5,1,:,j));
        subplot(2,1,1)
        hold on
        plot(attArr,eta_plot,'-o','Color',aColorArr{j})
        subplot(2,1,2)
        hold on
        plot(attArr,NC_plot,'-o','Color',aColorArr{j})
    end
end

figure(1)
subplot(2,1,1)
title('\bf{Viscosity}')
hold on
legend(muLegend,'location','best')
xlabel('$a$')
ylabel('$\eta_{app}/\eta_0$')
ylim([1 inf])

subplot(2,1,2)
hold on
title('\bf{Number of contacts}')
legend(muLegend,'location','best')
xlabel('$a$')
ylabel('$N_C$')

for j=1:nMu
    figure(j+1)
    set(gcf,'name',['mu',num2str(muArr(j))])
    subplot(2,1,1)
    title('\bf{Viscosity}')
    hold on
    legend(aLegend,'location','best')
    xlabel('$A_N$')
    ylabel('$\eta_{app}/\eta_0$')
    ylim([1 inf])
    
    subplot(2,1,2)
    hold on
    title('\bf{Number of contacts}')
    legend(aLegend,'location','best')
    xlabel('$A_N$')
    ylabel('$N_C$')
end