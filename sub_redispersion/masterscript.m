%% analysis for redispersion
clc;
clear;
close all;

%% define simulation cases
simulation_cases;

etaData = zeros(nMu,nAtt,nA);
NCData = zeros(nMu,nAtt,nA);
DetaData = zeros(nMu,nAtt,nA);
DNCData = zeros(nMu,nAtt,nA);
NC_total_statData = zeros(nMu,nAtt,nA);
NC_total_no_jointsData = zeros(nMu,nAtt,nA);
overlapData = zeros(nMu,nAtt,nA);
forcData = zeros(nMu,nAtt,nA);
sijData = zeros(nMu,nAtt,nA);
            
% load basis
Basis = load('Basis.mat');
etaDataBasis = Basis.etaData;
NCDataBasis = Basis.NCData;

%% For each parameter set
for i=1:nMu
    for k=1:nAtt
        % obtain basis for under mu and att
        etaBase = etaDataBasis(i,k,1); 
        NCBase = NCDataBasis(i,k,1); 
        for j=1:nA
            display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
            filePrefix = [dataPath,shape,'_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_'];
            %% read input and process data
            if (exist([filePrefix,'Stress_tensor.txt'], 'file') == 0)
                continue;
            end
            [stress_strain, sigxz_out] = read_stress([filePrefix,'Stress_tensor.txt']);
            [box_strain, sidex] = read_box([filePrefix,'Lbox.txt']);
            [NC_strain, NC] = read_NC([filePrefix,'Number_of_Contacts.txt']);
            [Con_strain, contact] = read_contactStat([filePrefix,'ContactStat.txt']);
            
            % define lboxArr to calculate stress
            lboxArr = find_lboxArr(stress_strain, box_strain, sidex);
            %
            % find relevant parameters for stress calculations
            nPt = length(stress_strain);
            [Seff,rp,nL3,L,gamma] = pCalc(nfib*ones(nPt,1),...
                nseg*ones(nPt,1), rps, kb, lboxArr, a, EY, Imom, eta0);
            
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
                stress_strain(end)
                box_strain(end)
                continue;
            end
            
            NC_stat = interval_average(NC_strain,NC,r');
            sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r');
            eta_stat = interval_average(stress_strain,etarel,r');
            
            etaData(i,k,j) = eta_stat(5,1);
            NCData(i,k,j) = NC_stat(5,1);
            
            % contact statistics
            NC_total = contact(:,1);
            NC_total_no_joints=contact(:,2);
            overlap=contact(:,3);
            forc =contact(:,4);
            sij =contact(:,5);
            
            NC_total_stat = interval_average(Con_strain,NC_total,r');
            NC_total_no_joints_stat = interval_average(Con_strain,NC_total_no_joints,r');
            overlap_stat = interval_average(Con_strain,overlap,r');
            forc_stat = interval_average(Con_strain,forc,r');
            sij_stat = interval_average(Con_strain,sij,r');
            
            NC_total_statData(i,k,j) = NC_total_stat(5,1);
            NC_total_no_jointsData(i,k,j) = NC_total_no_joints_stat(5,1);
            overlapData(i,k,j) = overlap_stat(5,1);
            forcData(i,k,j) = forc_stat(5,1);
            sijData(i,k,j) = sij_stat(5,1);
            
            if (etaBase ~= 0)
                DetaData(i,k,j) = etaData(i,k,j)-etaBase;
                DNCData(i,k,j) = NCData(i,k,j)-NCBase;
            end

        end
    end

end

%% Remove unfinished cases
% DetaData(:,4,:) = [];
% DNCData(:,4,:) = [];
% attC = pluck(attC,4);

% % mu0 att30
% DetaData(1,4,:) = NaN;
% DNCData(1,4,:) = NaN;
% 
% % mu15 att20 a5
% DetaData(4,3,4) = NaN;
% DNCData(4,3,4) = NaN;

% DetaData(4,3,4) = 0.5*(DetaData(4,2,4) + DetaData(4,4,4));
% DNCData(4,3,4) = 0.5*(DNCData(4,2,4) + DNCData(4,4,4));

%% plots
% Differences
%%{
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',attC,muC,aC)
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',aC,attC,muC)
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',attC,aC,muC)
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',muC,aC,attC)
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',muC,attC,aC)
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',aC,muC,attC)
%}
% Values
%{
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',attC,muC,aC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',aC,attC,muC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',attC,aC,muC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',muC,aC,attC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',muC,attC,aC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',aC,muC,attC)
%}
% Phase diagram with interfiber forces
%{
muC.name = '$\mu_{stat}$';
plotPhase(DetaData,'$\Delta \eta_{rel}$',attC,muC,aC,0)

%}
% Contact Statistics
%{
plot3dim(NC_total_no_jointsData,forcData,'$<N_C>$ (Unbroken, No Joints)','$<F_{tot}>$',attC,muC,aC)
plot3dim(NC_total_no_jointsData,forcData,'$<N_C>$ (Unbroken, No Joints)','$<F_{tot}>$',aC,attC,muC)
plot3dim(NC_total_no_jointsData,forcData,'$<N_C>$ (Unbroken, No Joints)','$<F_{tot}>$',attC,aC,muC)
plot3dim(NC_total_no_jointsData,forcData,'$<N_C>$ (Unbroken, No Joints)','$<F_{tot}>$',muC,aC,attC)
plot3dim(NC_total_no_jointsData,forcData,'$<N_C>$ (Unbroken, No Joints)','$<F_{tot}>$',muC,attC,aC)
plot3dim(NC_total_no_jointsData,forcData,'$<N_C>$ (Unbroken, No Joints)','$<F_{tot}>$',aC,muC,attC)
%}