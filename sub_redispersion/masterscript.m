%% analysis for redispersion - includes rheology, contacts, elastic energy
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;

basisStrain = 1000;

etaData = zeros(nMu,nAtt,nA);
NCData = zeros(nMu,nAtt,nA);
NC_total_statData = zeros(nMu,nAtt,nA);
NC_total_no_jointsData = zeros(nMu,nAtt,nA);
overlapData = zeros(nMu,nAtt,nA);
forcData = zeros(nMu,nAtt,nA);
sijData = zeros(nMu,nAtt,nA);
EelasData = zeros(nMu,nAtt,nA);
IData = zeros(nMu,nAtt,nA);

etaDataB = zeros(nMu,nAtt);
NCDataB = zeros(nMu,nAtt);
NC_total_statDataB = zeros(nMu,nAtt);
NC_total_no_jointsDataB = zeros(nMu,nAtt);
overlapDataB = zeros(nMu,nAtt);
forcDataB = zeros(nMu,nAtt);
sijDataB = zeros(nMu,nAtt);
EelasDataB = zeros(nMu,nAtt);
IDataB = zeros(nMu,nAtt);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nMu
    for k=1:nAtt
        
        display(['Processing basis mu',num2str(muArr(i)),'_att',num2str(attArr(k))]);
        filePrefix = [dataPath,shape,'_basis_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_'];
        %% read input and process data
        if (exist([filePrefix,'Stress_tensor.txt'], 'file') == 0)
            continue;
        end
        [stress_strain, sigxz_out] = read_stress([filePrefix,'Stress_tensor.txt']);
        [NC_strain, NC] = read_NC([filePrefix,'Number_of_Contacts.txt']);
        
        % find relevant parameters for stress calculations
        nPt = length(stress_strain);
        [Seff,rp,nL3,L,gamma] = pCalc(nfib*ones(nPt,1),...
            nseg*ones(nPt,1), rps, kb, sidex*ones(nPt,1), a, EY, Imom, eta0);
        
        % Dimensionalize stresses
        [sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
            stressDim(sigxz_out, zeros(nPt,1), zeros(nPt,1), ...
            zeros(nPt,1), nseg, rps, eta0, gamma, nL3);
        
        etarel = sigma./gamma;
        
        % Non-dimensionalize
        sigmap_nondim = sigmap.*L.^4/EY/Imom;
        
        %% calculate interval averages
        
        if NC_strain(end) >= basisStrain
            r = [basisStrain  NC_strain(end)]';
            NC_stat = interval_average(NC_strain,NC,r);
            eta_stat = interval_average(stress_strain,etarel,r);
            etaDataB(i,k) = eta_stat(2,1);
            NCDataB(i,k) = NC_stat(2,1);
        else
            NC_strain(end)
        end
        
        % contact statistics
        if (exist([filePrefix,'ContactStat.txt'], 'file') ~= 0)
            [Con_strain, contact] = read_contactStat([filePrefix,'ContactStat.txt']);
            r = [basisStrain  Con_strain(end)]';
            NC_total = contact(:,1);
            NC_total_no_joints = contact(:,2);
            overlap=contact(:,3);
            forc =contact(:,4);
            sij =contact(:,5);
            
            NC_total_stat = interval_average(Con_strain,NC_total,r);
            NC_total_no_joints_stat = interval_average(Con_strain,NC_total_no_joints,r);
            overlap_stat = interval_average(Con_strain,overlap,r);
            forc_stat = interval_average(Con_strain,forc,r);
            sij_stat = interval_average(Con_strain,sij,r);
            
            NC_total_statDataB(i,k) = NC_total_stat(2,1);
            NC_total_no_jointsDataB(i,k) = NC_total_no_joints_stat(2,1);
            overlapDataB(i,k) = overlap_stat(2,1);
            forcDataB(i,k) = forc_stat(2,1);
            sijDataB(i,k) = sij_stat(2,1);
        end
        
        % elastic energy statistics
        if (exist([filePrefix,'Eelastic.txt'], 'file') ~= 0)
            [Eelas_strain, Eelas_nondim] = read_elastic([filePrefix,'Eelastic.txt']);
            r = [basisStrain  Eelas_strain(end)]';
            Eelas_stat = interval_average(Eelas_strain,Eelas_nondim,r);
            EelasDataB(i,k) = Eelas_stat(2,1);
        end
        
        % intensity energy statistics
        if (exist([filePrefix,'Intensity.txt'], 'file') ~= 0)
            [I_strain, I] = read_intensity([filePrefix,'Intensity.txt']);
            r = [basisStrain  I_strain(end)]';
            I_stat = interval_average(I_strain,I,r);
            IDataB(i,k) = I_stat(2,1);
            
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain redispersed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nMu
    for k=1:nAtt
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
            if (exist([filePrefix,'ContactStat.txt'], 'file') ~= 0)
                [Con_strain, contact] = read_contactStat([filePrefix,'ContactStat.txt']);
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
            end
            
            % elastic energy statistics
            if (exist([filePrefix,'Eelastic.txt'], 'file') ~= 0)
                [Eelas_strain, Eelas_nondim] = read_elastic([filePrefix,'Eelastic.txt']);
                Eelas_stat = interval_average(Eelas_strain,Eelas_nondim,r');
                EelasData(i,k,j) = Eelas_stat(5,1);
            end
            
            % intensity energy statistics
            if (exist([filePrefix,'Intensity.txt'], 'file') ~= 0)
                [I_strain, I] = read_intensity([filePrefix,'Intensity.txt']);
                I_stat = interval_average(I_strain,I,r');
                IData(i,k,j) = I_stat(5,1); 
            end
            
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part III - calculate property difference due to redispersion cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DetaData = zeros(nMu,nAtt,nA);
DNCData = zeros(nMu,nAtt,nA);
DNC_total_statData = zeros(nMu,nAtt,nA);
DNC_total_no_jointsData = zeros(nMu,nAtt,nA);
DoverlapData = zeros(nMu,nAtt,nA);
DforcData = zeros(nMu,nAtt,nA);
DsijData = zeros(nMu,nAtt,nA);
DEelasData = zeros(nMu,nAtt,nA);
DIData = zeros(nMu,nAtt,nA);

etaData (etaData == 0) = NaN;
NCData (NCData == 0) = NaN;
NC_total_statData (NC_total_statData == 0) = NaN;
NC_total_no_jointsData (NC_total_no_jointsData == 0) = NaN;
overlapData (overlapData == 0) = NaN;
forcData (forcData == 0) = NaN;
sijData (sijData == 0) = NaN;
EelasData (EelasData == 0) = NaN;
IData (IData == 0) = NaN;

for i=1:nA
    DetaData(:,:,i) = etaData(:,:,i) - etaDataB;
    DNCData(:,:,i) = NCData(:,:,i) - NCDataB;
    DNC_total_statData(:,:,i) = NC_total_statData(:,:,i) - NC_total_statDataB;
    DNC_total_no_jointsData(:,:,i) = NC_total_no_jointsData(:,:,i) - NC_total_no_jointsDataB;
    DoverlapData(:,:,i) = overlapData(:,:,i) - overlapDataB;
    DforcData(:,:,i) = forcData(:,:,i) - forcDataB;
    DsijData(:,:,i) = sijData(:,:,i) - sijDataB;
    DEelasData(:,:,i) = EelasData(:,:,i) - EelasDataB;
    DIData(:,:,i) = IData(:,:,i) - IDataB;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part IV - plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% plots
% Differences
% plot3dim(etaData,DetaData,'$\eta_{app}/\eta_0$','$\Delta (\eta_{app}/\eta_0)$',attC,muC,aC)
% plot3dim(etaData,DetaData,'$\eta_{app}/\eta_0$','$\Delta (\eta_{app}/\eta_0)$',aC,attC,muC)
% plot3dim(etaData,DetaData,'$\eta_{app}/\eta_0$','$\Delta (\eta_{app}/\eta_0)$',attC,aC,muC)
% plot3dim(etaData,DetaData,'$\eta_{app}/\eta_0$','$\Delta (\eta_{app}/\eta_0)$',muC,aC,attC)
% plot3dim(etaData,DetaData,'$\eta_{app}/\eta_0$','$\Delta (\eta_{app}/\eta_0)$',muC,attC,aC)
% plot3dim(etaData,DetaData,'$\eta_{app}/\eta_0$','$\Delta (\eta_{app}/\eta_0)$',aC,muC,attC)
%{
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

%{
plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',attC,muC,aC)
plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',aC,attC,muC)
plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',attC,aC,muC)
plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',muC,aC,attC)
plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',muC,attC,aC)
plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',aC,muC,attC)
%}
% plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',attC,muC,aC)
% plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',aC,attC,muC)
% plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',attC,aC,muC)
% plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',muC,aC,attC)
% plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',muC,attC,aC)
% plot3dim(IData,EelasData,'$I$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',aC,muC,attC)

plot3dim(DIData,IData,'$\Delta I$','$I$',attC,muC,aC)
plot3dim(DIData,IData,'$\Delta I$','$I$',aC,attC,muC)
plot3dim(DIData,IData,'$\Delta I$','$I$',attC,aC,muC)
plot3dim(DIData,IData,'$\Delta I$','$I$',muC,aC,attC)
plot3dim(DIData,IData,'$\Delta I$','$I$',muC,attC,aC)
plot3dim(DIData,IData,'$\Delta I$','$I$',aC,muC,attC)

% plot3dim(DfData,fData,'$\Delta f$','$f$',attC,muC,aC)
% plot3dim(DfData,fData,'$\Delta f$','$f$',aC,attC,muC)
% plot3dim(DfData,fData,'$\Delta f$','$f$',attC,aC,muC)
% plot3dim(DfData,fData,'$\Delta f$','$f$',muC,aC,attC)
% plot3dim(DfData,fData,'$\Delta f$','$f$',muC,attC,aC)
% plot3dim(DfData,fData,'$\Delta f$','$f$',aC,muC,attC)

% plot3dim(DEelasData,EelasData,'$\Delta E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',attC,muC,aC)
% plot3dim(DEelasData,EelasData,'$\Delta E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',aC,attC,muC)
% plot3dim(DEelasData,EelasData,'$\Delta E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',attC,aC,muC)
% plot3dim(DEelasData,EelasData,'$\Delta E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',muC,aC,attC)
% plot3dim(DEelasData,EelasData,'$\Delta E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',muC,attC,aC)
% plot3dim(DEelasData,EelasData,'$\Delta E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$','$E_{elas}/ 8\pi\eta_0\dot{\gamma}l^3$',aC,muC,attC)

% Phase diagram with interfiber forces
%{
muC.name = '$\mu_{stat}$';
plotPhase(DetaData,'$\Delta \eta_{rel}$',attC,muC,aC,0)
%}