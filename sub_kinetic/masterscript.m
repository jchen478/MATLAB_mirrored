%% analysis for redispersion
clc;
clear;
close all;

%% define simulation cases
simulation_cases;

etaData = zeros(muC.ndata,attC.ndata,mukinC.ndata);
NCData = zeros(muC.ndata,attC.ndata,mukinC.ndata);

%% For each parameter set
for i=1:muC.ndata
    for k=1:attC.ndata
        for j=1:mukinC.ndata
            display(['Processing mu',num2str(muC.value(i)),'_att',...
                num2str(attC.value(k)),'_mukin',num2str(mukinC.value(j),'%.1f')])
            filePrefix = [dataPath,'kinetic_',shape,'_mu',...
                num2str(muC.value(i)),'_att',num2str(attC.value(k)),...
                '_kin',num2str(mukinC.value(j),'%.1f'),'_'];
            %% read input and process data
            [stress_strain, sigxz_out] = read_stress([filePrefix,'Stress_tensor.txt']);
            [NC_strain, NC] = read_NC([filePrefix,'Number_of_Contacts.txt']);
            
            % find relevant parameters for stress calculations
            nPt = length(stress_strain);
            [Seff,rp,nL3,L,gamma] = pCalc(nfib*ones(nPt,1),...
                nseg*ones(nPt,1), rps, kb, sidex*ones(nPt,1),...
                a, EY, Imom, eta0);
            
            % Dimensionalize stresses
            [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
                stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
                zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
            
            etarel = sigma./gamma;
            
            % Non-dimensionalize
            sigmap_nondim = sigmap.*L.^4/EY/Imom;
            
            %% calculate interval averages
            r = [1000 NC_strain(end)];
            NC_stat = interval_average(NC_strain,NC,r);
            sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r);
            eta_stat = interval_average(stress_strain,etarel,r);
            
            etaData(i,k,j) = eta_stat(1,1);
            NCData(i,k,j) = NC_stat(1,1);

        end
    end

end

%% plots
% Differences
% Values
%%{
% plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',attC,muC,mukinC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',mukinC,attC,muC)
% plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',attC,mukinC,muC)
% plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',muC,mukinC,attC)
% plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',muC,attC,mukinC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',mukinC,muC,attC)
%}