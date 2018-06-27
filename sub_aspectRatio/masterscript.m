%% analysis for redispersion with varying aspect ratio
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions 
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;

basisStrain = 1000;

nMu = mu.ndata; 
nRp = rp.ndata;
nAtt = att.ndata;

etaData = zeros(nRp,nMu,nAtt);
NCData = zeros(nRp,nMu,nAtt);
etaBasisData = zeros(nRp,nMu,nAtt);
NCBasisData = zeros(nRp,nMu,nAtt);
DetaData = zeros(nRp,nMu,nAtt);
DNCData = zeros(nRp,nMu,nAtt);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis 
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nRp
    for j=1:nMu
        for k=1:nAtt
             display(['Base - Processing rp',...
                 num2str(rp.value(i)),...
                 ' mu', num2str(mu.value(j)),...
                 ' att',num2str(att.value(k))])
             filePrefix = [dataPath,shape,'_basis_rp',...
                 num2str(rp.value(i)),...
                 '_mu',num2str(mu.value(j)),...
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
                nseg*ones(nPt,1), rps.value(i), kb.value(i),...
                sidex*ones(nPt,1), a, EY, Imom, eta0);
            
            % Dimensionalize stresses
            [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
                stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
                zeros(nPt,1), nseg, rps.value(i), eta0, gamma, nL3);
            
            % Calculate viscosity
            etarel = sigma./gamma;
            
            % Non-dimensionalize
            sigmap_nondim = sigmap.*L.^4/EY/Imom;
            
            if NC_strain(end) >= basisStrain
                r = [basisStrain  NC_strain(end)]';
                NC_stat = interval_average(NC_strain,NC,r);
                sigmap_nondim_stat = interval_average(stress_strain,sigmap_nondim,r);
                eta_stat = interval_average(stress_strain,etarel,r);
                etaBasisData(i,j,k) = eta_stat(2,1);
                NCBasisData(i,j,k) = NC_stat(2,1);
            else
                NC_strain(end)
            end
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
for i=1:nRp
    for j=1:nMu
        for k=1:nAtt
            display(['Processing rp',...
                num2str(rp.value(i)),...
                ' mu', num2str(mu.value(j)),...
                ' att',num2str(att.value(k))])
            filePrefix = [dataPath,shape,'_rp',...
                num2str(rp.value(i)),...
                '_mu',num2str(mu.value(j)),...
                '_att',num2str(att.value(k)),'_a3_'];
            
            %% read input and process data
            if (exist([filePrefix,'Stress_tensor.txt'], 'file') == 0)
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
                rps.value(i), kb.value(i), lboxArr, a, EY, Imom, eta0);
            
            % Dimensionalize stresses
            [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
                stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
                zeros(nPt,1), nseg, rps.value(i), eta0, gamma, nL3);
            
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
            
            if (etaBasisData(i,j,k) ~= 0)
                DetaData(i,j,k) = etaData(i,j,k)-etaBasisData(i,j,k);
                DNCData(i,j,k) = NCData(i,j,k)-NCBasisData(i,j,k);
            end
        end
    end
end

%% plots
DetaData (DetaData == 0) = nan; 
DNCData (DNCData == 0) = nan; 
% Differences
%%{
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',att,mu,rp)
% plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',rp,att,mu)
plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',att,rp,mu)
% plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',mu,rp,att)
% plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',mu,att,rp)
% plot3dim(DetaData,DNCData,'$\Delta (\eta_{app}/\eta_0)$','$\Delta N_C$',rp,mu,att)
%}