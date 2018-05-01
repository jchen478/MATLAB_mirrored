%% analysis for redispersion
clc; 
clear;
close all;

%% define simulation cases
simulation_cases; 

%% For each parameter set
for i=1:nMu
    fig = figure('Units','Inches','Position',[1 1 3.0 3.0]);
    hold on
    for j=1:nA
        
        display(['Processing mu',num2str(muArr(i)),'_a',num2str(aArr(j))])
        filePrefix = [dataPath,shape,'_mu',num2str(muArr(i)),'_a',num2str(aArr(j)),'_'];
        %% read input and process data
        [stress_strain, sigxz_out] = read_stress([filePrefix,'Stress_tensor.txt']);
        [box_strain, sidex] = read_box([filePrefix,'Lbox.txt']);
        [NC_strain, NC] = read_NC([filePrefix,'Number_of_Contacts.txt']);
        
        % define lboxArr to calculate stress
        lboxArr = find_lboxArr(stress_strain, box_strain, sidex);
 
        % find relevant parameters for stress calculations
        nPt = length(stress_strain);
        [Seff,rp,nL3,L,gamma] = pCalc(nfib*ones(nPt,1), nseg*ones(nPt,1), rps, kb, lboxArr, a, EY, Imom, eta0);

        % Dimensionalize stresses
        [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
            stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
            zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
        
        % Non-dimensionalize
        sigmap_nondim = sigmap.*L.^4/EY/Imom;
        
        %% calculate interval averages        
        r = round(intervals(box_strain,sidex),0);
        NC_stat = interval_average(NC_strain,NC,r);
        sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r);
        
        %% plot data to corresponding plots
        subplot(2,1,1)
        hold on
        title('\bf{Particle stress}')
%         scatter(aArr(j),sigmap_nondim_stat(3,1),'filled','markerFaceColor',rgb('MediumBlue'))
        scatter(aArr(j),sigmap_nondim_stat(5,1),'filled','markerFaceColor',rgb('Crimson'))
        xlabel('$a$')
        ylabel('$\sigma_{p,xz} L^4/ E_Y I$')
        
        subplot(2,1,2)
        hold on
        title('\bf{Number of contacts}')
%         scatter(aArr(j),NC_stat(3,1),'filled','markerFaceColor',rgb('MediumBlue'))
        scatter(aArr(j),NC_stat(5,1),'filled','markerFaceColor',rgb('Crimson'))
        xlabel('$a$')
        ylabel('$N_C$')
    end
end