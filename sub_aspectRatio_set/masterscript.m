%% analysis for redispersion with varying aspect ratio
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions 
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;

basisStrain = 800;

nRp = rp.ndata;
nCond = cond.ndata;

etaData = zeros(nRp,nCond);
NCData = zeros(nRp,nCond);
NC_total_statData =  zeros(nRp,nCond);
NC_total_no_jointsData =  zeros(nRp,nCond);
overlapData = zeros(nRp,nCond);
forcData = zeros(nRp,nCond);
sijData =  zeros(nRp,nCond);
EelasData =  zeros(nRp,nCond);

etaDataE = zeros(nRp,nCond);

etaDataB = zeros(nRp,nCond);
NCDataB = zeros(nRp,nCond);
NC_total_statDataB = zeros(nRp,nCond);
NC_total_no_jointsDataB = zeros(nRp,nCond);
overlapDataB = zeros(nRp,nCond);
forcDataB = zeros(nRp,nCond);
sijDataB = zeros(nRp,nCond);
EelasDataB = zeros(nRp,nCond);

etaDataBE = zeros(nRp,nCond);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis 
%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0;
for i=1:nRp
    for j=1:nCond
        
        display(['Base - Processing rp',...
            num2str(rp.value(i))]);
        filePrefix = [dataPathBasis,shape,'_rp',...
            num2str(rp.value(i)),...
            '_basis_',cond.value{j},'_'];
        
        %% read input and process data
        [r, etaDataB(i,j), etaDataBE(i,j)] = process_stress_lbox_eta_basis(filePrefix, basisStrain, 2, sidex, nfib.value(i), nseg, rps.value(i), kb.value(i), a, EY, Imom, eta0);
        % number of contacts
        NCDataB(i,j) = process_NC(filePrefix, basisStrain, 2, r');
        % contact statistics
        [NC_total_statDataB(i,j), NC_total_no_jointsDataB(i,j), overlapDataB(i,j),forcDataB(i,j),sijDataB(i,j)] = process_contactStat(filePrefix, basisStrain, 2, r');
        % elastic energy statistics
        EelasDataB(i,j) = process_elastic(filePrefix, basisStrain, 2, r');
        
    end
end
etaDataB (etaDataB == 0) = nan;
etaDataBObj = data3D('ND $\eta/\eta_0$',etaDataB);
NCDataBObj = data3D('ND $N_C$',NCDataB);

plot2dim1(etaDataBObj,rp,cond); 
% plot2dim1EB(etaDataBObj,etaDataBE,rp,cond); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain redispersed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nRp
    for j=1:nCond
            display(['Processing rp',...
                num2str(rp.value(i)),...
                ' mu', condMu.value{j},...
                ' att',condAtt.value{j}]);
            filePrefix = [dataPath,shape,'_rp',...
                num2str(rp.value(i)),...
                '_rand_shear_shear_mu',condMu.value{j},...
                '_att',condAtt.value{j},'_a4_'];
            
            %% read input and process dat
            % viscosity
            [r, etaData(i,j), etaDataE(i,j)] = process_stress_lbox_eta(filePrefix, 5, sidex, nfib.value(i), nseg, rps.value(i), kb.value(i), a, EY, Imom, eta0);
            % number of contacts
            NCData(i,j) = process_NC(filePrefix, basisStrain, 5, r');
            % contact statistics
            [NC_total_statData(i,j), NC_total_no_jointsData(i,j), overlapData(i,j),forcData(i,j),sijData(i,j)] = process_contactStat(filePrefix, basisStrain, 5, r');
            % elastic energy statistics
            EelasData(i,j) = process_elastic(filePrefix, basisStrain, 5, r');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part III - calculate property difference due to redispersion cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etaData (etaData == 0) = NaN;
NCData (NCData == 0) = NaN;
NC_total_statData (NC_total_statData == 0) = NaN;
NC_total_no_jointsData (NC_total_no_jointsData == 0) = NaN;
overlapData (overlapData == 0) = NaN;
forcData (forcData == 0) = NaN;
sijData (sijData == 0) = NaN;
EelasData (EelasData == 0) = NaN;

EelasDataB (EelasDataB == 0) = NaN;

DetaData = etaData - etaDataB;
DNCData = NCData - NCDataB;
DNC_total_statData = NC_total_statData - NC_total_statDataB;
DNC_total_no_jointsData = NC_total_no_jointsData - NC_total_no_jointsDataB;
DoverlapData = overlapData - overlapDataB;
DforcData = forcData - forcDataB;
DsijData = sijData - sijDataB;
DEelasData = EelasData - EelasDataB;

%% plots
etaData (etaData == 0) = nan; 
DetaData (DetaData == 0) = nan; 
DNCData (DNCData == 0) = nan; 

etaDataObj = data3D('$\eta/\eta_0$',etaData);
DetaDataObj = data3D('$\Delta \eta/\eta_0$',DetaData);

plot2dim1(etaDataObj,rp,cond); 
plot2dim1(DetaDataObj,rp,cond); 

