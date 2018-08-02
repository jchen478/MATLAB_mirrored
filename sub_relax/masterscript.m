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
rbStrain = 100;

nShear = shear.ndata;
nMu = mu.ndata;
nAtt = att.ndata;

etaData = zeros(nShear,nMu,nAtt);
NCData = zeros(nShear,nMu,nAtt);
NC_total_statData = zeros(nShear,nMu,nAtt);
NC_total_no_jointsData = zeros(nShear,nMu,nAtt);
overlapData = zeros(nShear,nMu,nAtt);
forcData = zeros(nShear,nMu,nAtt);
sijData = zeros(nShear,nMu,nAtt);
EelasData = zeros(nShear,nMu,nAtt);
IData = zeros(nShear,nMu,nAtt);


etaDataB = zeros(nMu,nAtt);
NCDataB = zeros(nMu,nAtt);
NC_total_statDataB = zeros(nMu,nAtt);
NC_total_no_jointsDataB = zeros(nMu,nAtt);
overlapDataB = zeros(nMu,nAtt);
forcDataB = zeros(nMu,nAtt);
sijDataB = zeros(nMu,nAtt);
EelasDataB = zeros(nMu,nAtt);
IDataB = zeros(nMu,nAtt);
etaDataRB = zeros(nMu,nAtt);
NCDataRB = zeros(nMu,nAtt);
NC_total_statDataRB = zeros(nMu,nAtt);
NC_total_no_jointsDataRB = zeros(nMu,nAtt);
overlapDataRB = zeros(nMu,nAtt);
forcDataRB = zeros(nMu,nAtt);
sijDataRB = zeros(nMu,nAtt);
EelasDataRB = zeros(nMu,nAtt);
IDataRB = zeros(nMu,nAtt);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sheared basis
r = 0;
for j=1:nMu
    for k=1:nAtt
        display(['Base - Processing mu',...
            num2str(mu.value(j)),...
            ' att',num2str(att.value(k))])
        filePrefix = [dataPathBasis,shape,'_rp75_basis_mu',...
            num2str(mu.value(j)),...
            '_att',num2str(att.value(k)),'_'];
        
        %% read input and process data
        [r, etaDataB(j,k)] = process_stress_lbox_eta_basis(filePrefix, basisStrain, 2, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0);
        % number of contacts
        NCDataB(j,k) = process_NC(filePrefix, basisStrain, 2, r');
        % contact statistics
        [NC_total_statDataB(j,k), NC_total_no_jointsDataB(j,k), overlapDataB(j,k),forcDataB(j,k),sijDataB(j,k)] = process_contactStat(filePrefix, basisStrain, 2, r');
        % elastic energy statistics
        EelasDataB(j,k) = process_elastic(filePrefix, basisStrain, 2, r');
    end
end
%% relaxed basis
for j=1:nMu
    for k=1:nAtt
        display(['Relaxed Base - Processing mu',...
            num2str(mu.value(j)),...
            ' att',num2str(att.value(k))])
        filePrefix = [dataPathBasis,shape,'_rp75_relaxed_basis_mu',...
            num2str(mu.value(j)),...
            '_att',num2str(att.value(k)),'_'];
        %% read input and process data
        [r, etaDataRB(j,k)] = process_stress_lbox_eta_basis(filePrefix, rbStrain, 2, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0);
        % number of contacts
        NCDataRB(j,k) = process_NC(filePrefix, rbStrain, 2, r');
        % contact statistics
        [NC_total_statDataRB(j,k), NC_total_no_jointsDataRB(j,k), overlapDataRB(j,k),forcDataRB(j,k),sijDataRB(j,k)] = process_contactStat(filePrefix, rbStrain, 2, r');
        % elastic energy statistics
        EelasDataRB(j,k) = process_elastic(filePrefix, rbStrain, 2, r');
        
    end
end

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
            filePrefix = [dataPath,shape,'_rp75_',...
                shear.value{i},...
                '_mu',num2str(mu.value(j)),...
                '_att',num2str(att.value(k)),'_a3_'];
            
            %% read input and process dat
            % viscosity
            [r, etaData(i,j,k)] = process_stress_lbox_eta(filePrefix, 5, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0);
            % number of contacts
            NCData(i,j,k) = process_NC(filePrefix, basisStrain, 5, r');
            % contact statistics
            [NC_total_statData(i,j,k), NC_total_no_jointsData(i,j,k), overlapData(i,j,k),forcData(i,j,k),sijData(i,j,k)] = process_contactStat(filePrefix, basisStrain, 5, r');
            % elastic energy statistics
            EelasData(i,j,k) = process_elastic(filePrefix, basisStrain, 5, r');
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part III - obtain differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DetaData = zeros(nShear,nMu,nAtt);
DNCData = zeros(nShear,nMu,nAtt);
DNC_total_statData = zeros(nShear,nMu,nAtt);
DNC_total_no_jointsData = zeros(nShear,nMu,nAtt);
DoverlapData = zeros(nShear,nMu,nAtt);
DforcData = zeros(nShear,nMu,nAtt);
DsijData = zeros(nShear,nMu,nAtt);
DEelasData = zeros(nShear,nMu,nAtt);
DIData = zeros(nShear,nMu,nAtt);

etaData (etaData == 0) = NaN;
NCData (NCData == 0) = NaN;
NC_total_statData (NC_total_statData == 0) = NaN;
NC_total_no_jointsData (NC_total_no_jointsData == 0) = NaN;
overlapData (overlapData == 0) = NaN;
forcData (forcData == 0) = NaN;
sijData (sijData == 0) = NaN;
EelasData (EelasData == 0) = NaN;
IData (IData == 0) = NaN;

for i=1:nShear
    % determine basis depending on shear
    if (strcmp(shear.value{i},'shear_noshear') || strcmp(shear.value{i},'noshear_noshear'))
        DetaData(i,:,:) = squeeze(etaData(i,:,:))-etaDataRB;
        DNCData(i,:,:) = squeeze(NCData(i,:,:))-NCDataRB;
        DNC_total_statData(i,:,:) = squeeze(NC_total_statData(i,:,:)) - NC_total_statDataRB;
        DNC_total_no_jointsData(i,:,:) = squeeze(NC_total_no_jointsData(i,:,:)) - NC_total_no_jointsDataRB;
        DoverlapData(i,:,:) = squeeze(overlapData(i,:,:)) - overlapDataRB;
        DforcData(i,:,:) = squeeze(forcData(i,:,:)) - forcDataRB;
        DsijData(i,:,:) = squeeze(sijData(i,:,:)) - sijDataRB;
        DEelasData(i,:,:) = squeeze(EelasData(i,:,:)) - EelasDataRB;
    else
        DetaData(i,:,:) = squeeze(etaData(i,:,:)) - etaDataB;
        DNCData(i,:,:) = squeeze(NCData(i,:,:)) - NCDataB;
        DNC_total_statData(i,:,:) = squeeze(NC_total_statData(i,:,:)) - NC_total_statDataB;
        DNC_total_no_jointsData(i,:,:) = squeeze(NC_total_no_jointsData(i,:,:)) - NC_total_no_jointsDataB;
        DoverlapData(i,:,:) = squeeze(overlapData(i,:,:)) - overlapDataB;
        DforcData(i,:,:) = squeeze(forcData(i,:,:)) - forcDataB;
        DsijData(i,:,:) = squeeze(sijData(i,:,:)) - sijDataB;
        DEelasData(i,:,:) = squeeze(EelasData(i,:,:)) - EelasDataB;
    end
end

%%
etaDataObj = data3D('$\eta/\eta_0$',etaData);
DetaDataObj = data3D('$\Delta \eta/\eta_0$',DetaData);
EelasDataObj = data3D('$E_{elastic}^*$',EelasData);
DEelasDataObj = data3D('$\Delta E_{elastic}^*$',DEelasData);

%% plots
plot3dim1(etaDataObj,att,shear,mu); 
plot3dim1(DetaDataObj,att,shear,mu);
% plot3dim1(EelasDataObj,att,shear,mu); 
% plot3dim1(DEelasDataObj,att,shear,mu);