%% analysis for redispersion
%  with shear on and off during concentration and rehydration
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;

basisStrain = 1000;

nShear = shear.ndata;
nMu = mu.ndata;
nAtt = att.ndata;

etaData = zeros(nShear,nMu,nAtt);
NCData = zeros(nShear,nMu,nAtt);
etaBasisData = zeros(nMu,nAtt);
NCBasisData = zeros(nMu,nAtt);
DetaData = zeros(nShear,nMu,nAtt);
DNCData = zeros(nShear,nMu,nAtt);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
% all shear conditions share the same basis
for j=1:nMu
    for k=1:nAtt
        display(['Base - Processing mu',...
            num2str(mu.value(j)),...
            ' att',num2str(att.value(k))])
        filePrefix = [dataPath,shape,'_basis_mu',...
            num2str(mu.value(j)),...
            '_att',num2str(att.value(k)),'_'];
        % read input and process data
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
        [Seff,rpdum,nL3,L,gamma] = pCalc(nfib*ones(nPt,1),...
            nseg*ones(nPt,1), rps, kb,...
            sidex*ones(nPt,1), a, EY, Imom, eta0);
        
        % Dimensionalize stresses
        [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
            stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
            zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
        
        % Calculate viscosity
        etarel = sigma./gamma;
        
        % Non-dimensionalize
        sigmap_nondim = sigmap.*L.^4/EY/Imom;
        
        if NC_strain(end) >= basisStrain
            r = [basisStrain  NC_strain(end)]';
            NC_stat = interval_average(NC_strain,NC,r);
            sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r);
            eta_stat = interval_average(stress_strain,etarel,r);
            etaBasisData(j,k) = eta_stat(2,1);
            NCBasisData(j,k) = NC_stat(2,1);
        else
            NC_strain(end)
        end
    end
end

%% plots
% 1. Rheology
%{
plot3dim(etaBasisData,NCBasisData,'$\eta_{app}/\eta_0$','$N_C$',att,mu,rp)
plot3dim(etaBasisData,NCBasisData,'$\eta_{app}/\eta_0$','$N_C$',rp,att,mu)
plot3dim(etaBasisData,NCBasisData,'$\eta_{app}/\eta_0$','$N_C$',att,rp,mu)
plot3dim(etaBasisData,NCBasisData,'$\eta_{app}/\eta_0$','$N_C$',mu,rp,att)
plot3dim(etaBasisData,NCBasisData,'$\eta_{app}/\eta_0$','$N_C$',mu,att,rp)
plot3dim(etaBasisData,NCBasisData,'$\eta_{app}/\eta_0$','$N_C$',rp,mu,att)
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain redispersed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nShear
    for j=1:nMu
        for k=1:nAtt
            display(['Processing ',...
                shear.value{i},...
                ' mu', num2str(mu.value(j)),...
                ' att',num2str(att.value(k))])
            filePrefix = [dataPath,shape,'_',...
                shear.value{i},...
                '_mu',num2str(mu.value(j)),...
                '_att',num2str(att.value(k)),'_a3_'];
            
            %% read input and process data
            if exist([filePrefix,'Stress_tensor.txt'], 'file')
                s = dir([filePrefix,'Stress_tensor.txt']);
                if s.bytes == 0
                    continue;
                end
            else
                continue;
            end
            [stress_strain, sigxz_out] = read_stress([filePrefix,'Stress_tensor.txt']);
            [box_strain, sidex] = read_box([filePrefix,'Lbox.txt']);
            [NC_strain, NC] = read_NC([filePrefix,'Number_of_Contacts.txt']);
            
            % define lboxArr to calculate stress
            lboxArr = find_lboxArr(stress_strain, box_strain, sidex);
            
            % find relevant parameters for stress calculations
            nPt = length(stress_strain);
            [Seff,rpdum,nL3,L,gamma] =...
                pCalc(nfib*ones(nPt,1), nseg*ones(nPt,1), ...
                rps, kb, lboxArr, a, EY, Imom, eta0);
            
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
            
            etaData(i,j,k) = eta_stat(5,1);
            NCData(i,j,k) = NC_stat(5,1);
            
            if (etaBasisData(j,k) ~= 0)
                DetaData(i,j,k) = etaData(i,j,k)-etaBasisData(j,k);
                DNCData(i,j,k) = NCData(i,j,k)-NCBasisData(j,k);
            end
        end
    end
end

%% plots
DetaData (DetaData == 0) = nan;
DNCData (DNCData == 0) = nan;
% Differences
% %{
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',att,mu,shear)
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',att,shear,mu)
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',shearplot,att,mu)
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',shearplot,mu,att)
%}














