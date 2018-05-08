%% analysis for redispersion
clc;
clear;
close all;

%% define simulation cases
simulation_cases;

%% For each parameter set
fig = figure('Units','Inches','Position',[1 1 4.0 6.0]);
hold on
for i=1:nMu
    stress_trend=zeros(nNL3,2);
    NC_trend=zeros(nNL3,2);
    eta_trend=zeros(nNL3,2);
    for j=1:nNL3
        display(['Processing mu',num2str(muArr(i)),'_nL3_',num2str(nL3Arr(j))])
        filePrefix = [dataPath,shape,'_mu',num2str(muArr(i)),'_nL3_',num2str(nL3Arr(j)),'_'];
        %% read input and process data
        [stress_strain, sigxz_out] = read_stress([filePrefix,'Stress_tensor.txt']);
        [NC_strain, NC] = read_NC([filePrefix,'Number_of_Contacts.txt']);
        
        % find relevant parameters for stress calculations
        nPt = length(stress_strain);
        [Seff,rp,nL3,L,gamma] = pCalc(nfib*ones(nPt,1),...
            nseg*ones(nPt,1), rps, kb, sidex(j)*ones(nPt,1), a, EY, Imom, eta0);
        
        % Dimensionalize stresses
        [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
            stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
            zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
        
        % Calculate viscosity
        etarel = sigma./gamma;
%         se_etarel = se_sigma./gamma;
        
        % Non-dimensionalize
        sigmap_nondim = sigmap.*L.^4/EY/Imom;
        
        %% calculate interval averages
        r = [300 NC_strain(end)];
        NC_stat = interval_average(NC_strain,NC,r);
        sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r);
        eta_stat = interval_average(stress_strain,etarel,r);
        stress_trend(j,:) = sigmap_nondim_stat;
        NC_trend(j,:) = NC_stat;
        eta_trend(j,:) = eta_stat;
    end
    
    %% plot data to corresponding plots
    subplot(2,1,1)
    box on
    hold on
    title('\bf{Viscosity}')
    plot(volfracArr,eta_trend(:,1),'-o')
    xlabel('$\phi (\%)$')
    ylabel('$\eta_{app}/\eta_0$')
%     box on
%     hold on
%     title('\bf{Particle stress}')
%     plot(nL3Arr,stress_trend(:,1),'-o')
%     xlabel('$nL^3$')
%     ylabel('$\sigma_{p,xz} L^4/ E_Y I$')
    
    subplot(2,1,2)
    box on
    hold on
    title('\bf{Number of contacts}')
    plot(volfracArr,NC_trend(:,1),'-o')
    xlabel('$\phi (\%)$')
    ylabel('$N_C$')
end

subplot(2,1,1)
hold on
legend(muLegend,'location','northwest')
% set(gca,'yscale','log')

subplot(2,1,2)
hold on
legend(muLegend,'location','southeast')
