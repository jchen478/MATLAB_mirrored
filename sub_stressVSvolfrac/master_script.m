%% analysis for concentration dependence
clc;
clear;
close all;

%% define simulation cases
simulation_cases;

etaData = zeros(nMu,nAtt,nNL3);
NCData = zeros(nMu,nAtt,nNL3);
DyyData = zeros(nMu,nAtt,nNL3);
DzzData = zeros(nMu,nAtt,nNL3);
elasticData = zeros(nMu,nAtt,nNL3);

%% For each parameter set
for i=1:nMu
    for k=1:nAtt
        for j=1:nNL3
            display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_nL3_',num2str(nL3Arr(j))])
            filePrefix = [dataPath,shape,'_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_nL3_',num2str(nL3Arr(j)),'_'];
            %% read input and process data
            if exist([filePrefix,'Stress_tensor.txt'], 'file')
                s = dir([filePrefix,'Stress_tensor.txt']);
                if s.bytes == 0
                    continue;
                end
                [stress_strain, sigxz_out] = read_stress([filePrefix,'Stress_tensor.txt']);
                [NC_strain, NC] = read_NC([filePrefix,'Number_of_Contacts.txt']);
            else
                continue;
            end
            
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
            
            % Non-dimensionalize
            sigmap_nondim = sigmap.*L.^4/EY/Imom;
            
            %% calculate interval averages
            % if simulation has reached strain
%             if NC_strain(end) == 1500
%                 r = [1000 NC_strain(end)]';
%             elseif NC_strain(end) > 350
%                 r = [350 NC_strain(end)]';
%             end
            if nL3Arr(j) == 20
                NC_strain(end)
            end
            if NC_strain(end) >= 350  
                r = [350  NC_strain(end)]';
                NC_stat = interval_average(NC_strain,NC,r);
                sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r);
                eta_stat = interval_average(stress_strain,etarel,r);
                etaData(i,k,j) = eta_stat(2,1);
                NCData(i,k,j) = NC_stat(2,1);
                
                %{
                %% Extra process -- elastic energy
                if exist([filePrefix,'Eelastic.txt'], 'file')
                    s = dir([filePrefix,'Eelastic.txt']);
                    if s.bytes == 0
                        break;
                    end
                    % read input
                    [elastic_strain, Eelastic] = read_elastic([filePrefix,'Eelastic.txt']);
                    elastic_stat = interval_average(elastic_strain,Eelastic,r);
                    
                    elasticData(i,k,j) = elastic_stat(1,1);
                else
                    disp(['Simulation completed but ',...
                        filePrefix,'Eelastic.txt does not exist']);
                end
                
                %% Extra process -- diffusivity
                if exist([filePrefix,'MSD.txt'], 'file')
                    s = dir([filePrefix,'MSD.txt']);
                    if s.bytes == 0
                        break;
                    end
                    % read input
                    [MSD_strain, MSD] = read_MSD([filePrefix,'MSD.txt']);
                    Dyy_stat = interval_slope(MSD_strain,MSD(:,1),r);
                    Dzz_stat = interval_slope(MSD_strain,MSD(:,2),r);
                    if MSD_strain(end) ~= NC_strain(end)
                        MSD_strain(end)
                        NC_strain(end)
                    end
                    DyyData(i,k,j) = Dyy_stat(1,1);
                    DzzData(i,k,j) = Dzz_stat(1,1);
                    if MSD_strain(1) ~= 500
                        disp('MSD start strain is not 500');
                    end
                else
                    disp(['Simulation completed but ',...
                        filePrefix,'MSD.txt does not exist']);
                end
                %}
            end
        end
    end
end

%% save basis
save('../sub_redispersion/Basis.mat','etaData','NCData');

%% plots
% 1. Rheology
%%{
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',attC,muC,volfracC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',volfracC,attC,muC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',attC,volfracC,muC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',muC,volfracC,attC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',muC,attC,volfracC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',volfracC,muC,attC)
%}

%{
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',attC,muC,weightfracArr)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',weightfracArr,attC,muC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',attC,weightfracArr,muC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',muC,weightfracArr,attC)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',muC,attC,weightfracArr)
plot3dim(etaData,NCData,'$\eta_{app}/\eta_0$','$N_C$',weightfracArr,muC,attC)
%}

% 2. Diffusivity
%{
plot3dim(DyyData,DzzData,'$D_{yy}/b^2\dot{\gamma}$','$D_{zz}/b^2\dot{\gamma}$',attC,muC,volfracC)
plot3dim(DyyData,DzzData,'$D_{yy}/b^2\dot{\gamma}$','$D_{zz}/b^2\dot{\gamma}$',volfracC,attC,muC)
plot3dim(DyyData,DzzData,'$D_{yy}/b^2\dot{\gamma}$','$D_{zz}/b^2\dot{\gamma}$',attC,volfracC,muC)
plot3dim(DyyData,DzzData,'$D_{yy}/b^2\dot{\gamma}$','$D_{zz}/b^2\dot{\gamma}$',muC,volfracC,attC)
plot3dim(DyyData,DzzData,'$D_{yy}/b^2\dot{\gamma}$','$D_{zz}/b^2\dot{\gamma}$',muC,attC,volfracC)
plot3dim(DyyData,DzzData,'$D_{yy}/b^2\dot{\gamma}$','$D_{zz}/b^2\dot{\gamma}$',volfracC,muC,attC)
%}

% 3. Elastic energy
%{
plot3dim(elasticData,DzzData,'$E_{elas} / 8\pi\eta_0\dot{\gamma}l^3$','$D_{zz}/b^2\dot{\gamma}$',attC,muC,volfracC)
plot3dim(elasticData,DzzData,'$E_{elas} / 8\pi\eta_0\dot{\gamma}l^3$','$D_{zz}/b^2\dot{\gamma}$',volfracC,attC,muC)
plot3dim(elasticData,DzzData,'$E_{elas} / 8\pi\eta_0\dot{\gamma}l^3$','$D_{zz}/b^2\dot{\gamma}$',attC,volfracC,muC)
plot3dim(elasticData,DzzData,'$E_{elas} / 8\pi\eta_0\dot{\gamma}l^3$','$D_{zz}/b^2\dot{\gamma}$',muC,volfracC,attC)
plot3dim(elasticData,DzzData,'$E_{elas} / 8\pi\eta_0\dot{\gamma}l^3$','$D_{zz}/b^2\dot{\gamma}$',muC,attC,volfracC)
plot3dim(elasticData,DzzData,'$E_{elas} / 8\pi\eta_0\dot{\gamma}l^3$','$D_{zz}/b^2\dot{\gamma}$',volfracC,muC,attC)
%}