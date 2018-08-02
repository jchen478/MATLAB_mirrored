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

etaDataE = zeros(nMu,nAtt,nA);

etaDataB = zeros(nMu,nAtt);
NCDataB = zeros(nMu,nAtt);
NC_total_statDataB = zeros(nMu,nAtt);
NC_total_no_jointsDataB = zeros(nMu,nAtt);
overlapDataB = zeros(nMu,nAtt);
forcDataB = zeros(nMu,nAtt);
sijDataB = zeros(nMu,nAtt);
EelasDataB = zeros(nMu,nAtt);

etaDataBE = zeros(nMu,nAtt);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0;
for i=1:nMu
    for k=1:nAtt
        display(['Processing basis mu',num2str(muArr(i)),'_att',num2str(attArr(k))]);
        filePrefix = [dataPathBasis,shape,'_rp75_basis_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_'];
        %% read input and process data
        [r, etaDataB(i,k), etaDataBE(i,k)] = process_stress_lbox_eta_basis(filePrefix, basisStrain, 2, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0);
        % number of contacts
        NCDataB(i,k) = process_NC(filePrefix, basisStrain, 2, r');
        % contact statistics
        [NC_total_statDataB(i,k), NC_total_no_jointsDataB(i,k), overlapDataB(i,k),forcDataB(i,k),sijDataB(i,k)] = process_contactStat(filePrefix, basisStrain, 2, r');
        % elastic energy statistics
        EelasDataB(i,k) = process_elastic(filePrefix, basisStrain, 2, r'); 
    end
end

%%
etaDataB (etaDataB == 0) = NaN;
etaDataBObj = data3D('Basis $\eta/\eta_0$',etaDataB);
plot2dim1(etaDataBObj,attC,muC); 
plot2dim1EB(etaDataBObj,etaDataBE,attC,muC); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain redispersed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nMu
    for k=1:nAtt
        for j=1:nA
            display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
            filePrefix = [dataPath,shape,'_rp75_rand_shear_shear_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_'];
            %% read input and process dat
            % viscosity
            [r, etaData(i,k,j), etaDataE(i,k,j)] = process_stress_lbox_eta(filePrefix, 5, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0);
            % number of contacts
            NCData(i,k,j) = process_NC(filePrefix, basisStrain, 5, r');
            % contact statistics
            [NC_total_statData(i,k,j), NC_total_no_jointsData(i,k,j), overlapData(i,k,j),forcData(i,k,j),sijData(i,k,j)] = process_contactStat(filePrefix, basisStrain, 5, r');
            % elastic energy statistics
            EelasData(i,k,j) = process_elastic(filePrefix, basisStrain, 5, r');
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

etaData (etaData == 0) = NaN;
NCData (NCData == 0) = NaN;
NC_total_statData (NC_total_statData == 0) = NaN;
NC_total_no_jointsData (NC_total_no_jointsData == 0) = NaN;
overlapData (overlapData == 0) = NaN;
forcData (forcData == 0) = NaN;
sijData (sijData == 0) = NaN;
EelasData (EelasData == 0) = NaN;

EelasDataB (EelasDataB == 0) = NaN;

for i=1:nA
    DetaData(:,:,i) = etaData(:,:,i) - etaDataB;
    DNCData(:,:,i) = NCData(:,:,i) - NCDataB;
    DNC_total_statData(:,:,i) = NC_total_statData(:,:,i) - NC_total_statDataB;
    DNC_total_no_jointsData(:,:,i) = NC_total_no_jointsData(:,:,i) - NC_total_no_jointsDataB;
    DoverlapData(:,:,i) = overlapData(:,:,i) - overlapDataB;
    DforcData(:,:,i) = forcData(:,:,i) - forcDataB;
    DsijData(:,:,i) = sijData(:,:,i) - sijDataB;
    DEelasData(:,:,i) = EelasData(:,:,i) - EelasDataB;
end

%%
etaDataObj = data3D('$\eta/\eta_0$',etaData);
DetaDataObj = data3D('$\Delta \eta/\eta_0$',DetaData);
EelasDataObj = data3D('$E_{el}$',EelasData);
DEelasDataObj = data3D('$\Delta E_{el}$',DEelasData);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part IV - plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot3dim1(etaDataObj,attC,muC,aC); 
plot3dim1(DetaDataObj,attC,muC,aC);
% plot3dim1(EelasDataObj,attC,muC,aC); 
% plot3dim1(DEelasDataObj,attC,muC,aC);
% plot3dim1EB(etaDataObj,etaDataE,attC,muC,aC); 
% Phase diagram with interfiber forces
%%{
muC.name = '$\mu_{stat}$';
plotPhase(DetaData,'$\Delta \eta/\eta_0$',attC,muC,aC,-0.05)
%}