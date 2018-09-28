%% analysis for redispersion with varying aspect ratio
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;

basisStrain = 1000;
rbStrain = 200;

nMu = mu.ndata;
nRp = rp.ndata;
nAtt = att.ndata;

etaData = zeros(nRp,nMu,nAtt);
NCData = zeros(nRp,nMu,nAtt);
NC_total_statData = zeros(nRp,nMu,nAtt);
NC_total_no_jointsData = zeros(nRp,nMu,nAtt);
overlapData = zeros(nRp,nMu,nAtt);
forcData = zeros(nRp,nMu,nAtt);
sijData = zeros(nRp,nMu,nAtt);
EelasData = zeros(nRp,nMu,nAtt);

etaDataB = zeros(nRp,nMu,nAtt);
NCDataB = zeros(nRp,nMu,nAtt);
NC_total_statDataB = zeros(nRp,nMu,nAtt);
NC_total_no_jointsDataB = zeros(nRp,nMu,nAtt);
overlapDataB = zeros(nRp,nMu,nAtt);
forcDataB = zeros(nRp,nMu,nAtt);
sijDataB = zeros(nRp,nMu,nAtt);
EelasDataB = zeros(nRp,nMu,nAtt);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0;
for i=1:nRp
    for j=1:nMu
        for k=1:nAtt
            display(['Base - Processing rp',...
                num2str(rp.value(i)),...
                ' mu', num2str(mu.value(j)),...
                ' att',num2str(att.value(k))])
            filePrefix = [dataPathBasis,shape,'_rp',...
                num2str(rp.value(i)),...
                '_basis_mu',num2str(mu.value(j)),...
                '_att',num2str(att.value(k)),'_'];
            
            %% read input and process data
            [r, etaDataB(i,j,k)] = process_stress_lbox_eta_basis(filePrefix, basisStrain, 2, sidex, nfib.value(i), nseg, rps.value(i), kb.value(i), a, EY, Imom, eta0);
            % number of contacts
            NCDataB(i,j,k) = process_NC(filePrefix, basisStrain, 2, r');
            % contact statistics
            [NC_total_statDataB(i,j,k), NC_total_no_jointsDataB(i,j,k), overlapDataB(i,j,k),forcDataB(i,j,k),sijDataB(i,j,k)] = process_contactStat(filePrefix, basisStrain, 2, r');
            % elastic energy statistics
            EelasDataB(i,j,k) = process_elastic([filePrefix,'relaxed_'], rbStrain, 2, [rbStrain 500]');
            EelasDataB(i,j,k) = EelasDataB(i,j,k)*(rp.value(i)/75)^3; 
        end
    end
end

etaDataBObj = data3D('ND $\eta/\eta_0$',etaDataB);
plot3dim1(etaDataBObj,att,rp,mu);

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
                '_rand_shear_shear_mu',num2str(mu.value(j)),...
                '_att',num2str(att.value(k)),'_a3_'];
            
            %% read input and process dat
            % viscosity
            [r, etaData(i,j,k)] = process_stress_lbox_eta(filePrefix, 5, sidex, nfib.value(i), nseg, rps.value(i), kb.value(i), a, EY, Imom, eta0);
            % number of contacts
            NCData(i,j,k) = process_NC(filePrefix, basisStrain, 5, r');
            % contact statistics
            [NC_total_statData(i,j,k), NC_total_no_jointsData(i,j,k), overlapData(i,j,k),forcData(i,j,k),sijData(i,j,k)] = process_contactStat(filePrefix, basisStrain, 5, r');
            % elastic energy statistics
            EelasData(i,j,k) = process_elastic([filePrefix,'relaxed_'], rbStrain, 2, [rbStrain 500]');
            EelasData(i,j,k) = EelasData(i,j,k)*(rp.value(i)/75)^3; 
        end
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
DetaData (DetaData == 0) = nan;
DNCData (DNCData == 0) = nan;
etaDataObj = data3D('$\eta/\eta_0$',etaData);
DetaDataObj = data3D('$\Delta \eta/\eta_0$',DetaData);
EelasDataObj = data3D('DR $E_{elas}$',EelasData);
DEelasDataObj = data3D('$\Delta E_{elas}/E_{scale}$',DEelasData);

plot3dim1(EelasDataObj,att,mu,rp);
plot3dim1(DEelasDataObj,att,mu,rp);
plot3dim1(etaDataObj,att,rp,mu);
plot3dim1(DetaDataObj,att,rp,mu);

%% paper figure
figure('Units','Inches','Position',[1 1 3.0 2.5]);
hold on
box on
xlabel('$A_N$')
ylabel('$\Delta \eta/\eta_0$')
legendArr = cell(nRp*nMu,1);

colorOrder = [ rgb('Orange'); rgb('MediumBlue'); rgb('DarkGreen')];
defaultFormat
for i=1:2
    figure('Units','Inches','Position',[1 1 3.0 2.5]);
    hold on
    box on
    xlabel('$A_N$')
    ylabel('$\Delta \eta/\eta_0$')
    title(['$r_p = $',num2str(rp.value(i))])
    for j=1:nMu
        %         ind = (i-1)*nMu+j;
        plot(att.value,squeeze(DetaData(i,j,:)),'-.o','linewidth',2.0,'color',colorOrder(j,:))
        %         legendArr{ind} = ['$(r_p,\mu)=$(',num2str(rp.value(i)),',',num2str(mu.value(j)),')'];
    end
%     legend(mu.legend,'location','northoutside','orientation','horizontal')
     legend(mu.legend,'location','best')
end
