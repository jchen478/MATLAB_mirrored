%% analysis for redispersion - includes rheology, contacts, elastic energy
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;
defaultFormat
basisStrain = 1000;
rxStrain = 300;
Fzero = [9 20 30 35 50 ; 0.2504 0.1306 0.09691 0.08558 0.06155];
Fmax = [9 20 30 35 50 ; 0.2747 0.1808 0.1526 0.1432 0.1236];

var = {'eta','NC','NC_tot','NC_tot_no_joints','overs','forc','sij'};
nVar = length(var);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
etaB = zeros(nMu,nAtt,nDup);
NCB = zeros(nMu,nAtt,nDup);
NC_totB = zeros(nMu,nAtt,nDup);
NC_tot_no_jointsB = zeros(nMu,nAtt,nDup);
oversB = zeros(nMu,nAtt,nDup);
forcB = zeros(nMu,nAtt,nDup);
sijB = zeros(nMu,nAtt,nDup);

for ii=1:nDup
    for i=1:nMu
        for k=1:nAtt
            display(['Processing basis mu',num2str(muArr(i)),'_att',num2str(attArr(k))]);
            filePrefix = [dataPathBasis,shape,'_rp75_basis_',dup.value{ii},'_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_'];
            %display(filePrefix)
            %% read input and process data
            [r, etaB(i,k,ii)] = process_stress_lbox_eta_basis(filePrefix, basisStrain, 2, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0);
            % number of contacts
            NCB(i,k,ii) = process_NC(filePrefix, basisStrain, 2, r');
            % contact statistics
            [NC_totB(i,k,ii), NC_tot_no_jointsB(i,k,ii), oversB(i,k,ii),forcB(i,k,ii),sijB(i,k,ii)] = process_contactStat(filePrefix, basisStrain, 2, r');
            % normalize
            oversB(i,k) = oversB(i,k)/nfib;
        end
    end
end

%%
etaB (etaB == 0) = nan;

% average over duplicate
etaB_mean = nanmean(etaB,3);
etaB_std = nanstd(etaB,0,3);

% replace data with mean
etaB = etaB_mean;

%%{
etaBObj = data3D('ND $\eta/\eta_0$',etaB);
NCBObj = data3D('ND Unbroken $N_C$',NCB);
NC_totBObj = data3D('ND $N_C$',NC_totB);
NC_tot_no_jointsBObj = data3D('ND $N_C$ no joints',NC_tot_no_jointsB);
oversBObj = data3D('ND Overlaps',oversB);
forcBObj = data3D('ND $F^N$',forcB);
sijBObj = data3D('ND $h_{ij}$',sijB);
NC_brokenBObj = data3D('ND Broken $N_C$',NC_tot_no_jointsB-NCB);



%%
% plot2dim1(etaBObj,attC,muC);
plot2dim1EB(etaBObj, etaB_std,attC,muC);
% plot2dim1(NCBObj,attC,muC);
% plot2dim1(NC_totBObj,attC,muC);
% plot2dim1(NC_tot_no_jointsBObj,attC,muC);
% plot2dim1(NC_brokenBObj,attC,muC);
% plot2dim1(oversBObj,attC,muC);
% plot2dim1(forcBObj,attC,muC);
% plot2dim1(sijBObj,attC,muC);
% scatter(Fzero(1,:),Fzero(2,:),60,'x','MarkerEdgeColor',rgb('DodgerBlue'),'linewidth',2)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain relaxed basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
NCRxB = zeros(nMu,nAtt,nDup);
NC_totRxB = zeros(nMu,nAtt,nDup);
NC_tot_no_jointsRxB = zeros(nMu,nAtt,nDup);
oversRxB = zeros(nMu,nAtt,nDup);
forcRxB = zeros(nMu,nAtt,nDup);
sijRxB = zeros(nMu,nAtt,nDup);
EelasRxB = zeros(nMu,nAtt,nDup);

r = [rxStrain 500];
for ii=1:nDup
    for i=1:nMu
        for k=1:nAtt
            display(['Processing basis mu',num2str(muArr(i)),'_att',num2str(attArr(k))]);
            filePrefix = [dataPathBasis,shape,'_rp75_basis_',dup.value{ii},'_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_relaxed_'];
            % number of contacts
            NCRxB(i,k,ii) = process_NC(filePrefix, rxStrain, 2, r');
            % contact statistics
            [NC_totRxB(i,k,ii), NC_tot_no_jointsRxB(i,k,ii), oversRxB(i,k,ii),forcRxB(i,k,ii),sijRxB(i,k,ii)] = process_contactStat(filePrefix, rxStrain, 2, r');
            % elastic energy statistics
            EelasRxB(i,k,ii) = process_elastic(filePrefix, rxStrain, 2, r');
            % normalize
            oversRxB(i,k,ii) = oversRxB(i,k,ii)/nfib;
        end
    end
end

%%
A = EelasRxB(2:end,:,:);

A (A == 0) = nan;
EelasRxB(2:end,:,:) = A;

% average over duplicate
EelasRxB_mean = nanmean(EelasRxB,3);
EelasRxB_std = nanstd(EelasRxB,0,3);

% replace data with mean
EelasRxB = EelasRxB_mean;

sijRxB(sijRxB == 0) = NaN;
NC_tot_no_jointsDB = NC_tot_no_jointsRxB-NC_tot_no_jointsB;
%%{
EelasRxBObj = data3D('Rx-ND $E_{el}$',EelasRxB);
NCRxBObj = data3D('Rx-ND Unbroken $N_C$',NCRxB);
NC_totRxBObj = data3D('Rx-ND $N_C$',NC_totRxB);
NC_tot_no_jointsRxBObj = data3D('Rx-ND $N_C$ no joints',NC_tot_no_jointsRxB);
oversRxBObj = data3D('Rx-ND Overlaps',oversRxB);
forcRxBObj = data3D('Rx-ND $F^N$',forcRxB);
sijRxBObj = data3D('Rx-ND $h_{ij}$',sijRxB);
NC_brokenRxBObj = data3D('Rx-ND Broken $N_C$',NC_tot_no_jointsRxB-NCRxB);
NC_tot_no_jointsDBObj = data3D('Rx-ND - ND $N_C$ no joints',NC_tot_no_jointsDB);

plot2dim1EB(EelasRxBObj,EelasRxB_std,attC,muC);
% plot2dim1(NCRxBObj,attC,muC);
% plot2dim1(NC_totRxBObj,attC,muC);
% plot2dim1(NC_tot_no_jointsRxBObj,attC,muC);
% plot2dim1(NC_brokenRxBObj,attC,muC);
% plot2dim1(oversRxBObj,attC,muC);
% plot2dim1(forcRxBObj,attC,muC);
% plot2dim1(sijRxBObj,attC,muC);
% scatter(Fzero(1,:),Fzero(2,:),30,'x','MarkerEdgeColor',rgb('Lime'),'linewidth',2)

% plot2dim1(NC_tot_no_jointsDBObj,attC,muC);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part III - obtain redispersed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta = zeros(nMu,nAtt,nA,nDup);
NC = zeros(nMu,nAtt,nA,nDup);
NC_tot = zeros(nMu,nAtt,nA,nDup);
NC_tot_no_joints = zeros(nMu,nAtt,nA,nDup);
overs = zeros(nMu,nAtt,nA,nDup);
forc = zeros(nMu,nAtt,nA,nDup);
sij = zeros(nMu,nAtt,nA,nDup);

for ii=1:nDup
    for i=1:nMu
        for k=1:nAtt
            for j=1:nA
                display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
                filePrefix = [dataPath,shape,'_rp75_rand_shear_shear_',dup.value{ii},'_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_'];
                %% read input and process dat
                % viscosity
                [r, eta(i,k,j,ii)] = process_stress_lbox_eta(filePrefix, 5, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0);
                % number of contacts
                NC(i,k,j,ii) = process_NC(filePrefix, basisStrain, 5, r');
                % contact statistics
                [NC_tot(i,k,j,ii), NC_tot_no_joints(i,k,j,ii), overs(i,k,j,ii),forc(i,k,j,ii),sij(i,k,j,ii)] = process_contactStat(filePrefix, basisStrain, 5, r');
                % normalize
                overs(i,k,j,ii) = overs(i,k,j,ii)/nfib;
            end
        end
    end
end

eta(eta == 0) = NaN;
NC_tot(NC_tot == 0) = NaN;
NC_tot_no_joints(NC_tot_no_joints == 0) = NaN;
sij(sij == 0) = NaN;

eta_mean = nanmean(eta,4);
eta_std = nanstd(eta,0,4);

% replace data with mean
eta = eta_mean;

%%{
etaObj = data3D('DR $\eta/\eta_0$',eta);
NCObj = data3D('DR Unbroken $N_C$',NC);
NC_totObj = data3D('DR $N_C$',NC_tot);
NC_tot_no_jointsObj = data3D('DR $N_C$ no joints',NC_tot_no_joints);
oversObj = data3D('DR Overlaps',overs);
forcObj = data3D('DR $F^N$',forc);
sijObj = data3D('DR $h_{ij}$',sij);
NC_brokenObj = data3D('DR Broken $N_C$',NC_tot_no_joints-NC);


plot3dim1E(etaObj,eta_std,attC,muC,aC);
% plot3dim1(etaObj,attC,muC,aC);
% plot3dim1(NCObj,attC,muC,aC);
% plot3dim1(NC_totObj,attC,muC,aC);
% plot3dim1(NC_tot_no_jointsObj,attC,muC,aC);
% plot3dim1(NC_brokenObj,attC,muC,aC);
% plot3dim1(oversObj,attC,muC,aC);
% plot3dim1(forcObj,attC,muC,aC);
% plot3dim1(sijObj,attC,muC,aC);
% scatter(Fzero(1,:),Fzero(2,:),60,'x','MarkerEdgeColor',rgb('DodgerBlue'),'linewidth',2)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part IV - obtain relaxed redispersion value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NCRx = zeros(nMu,nAtt,nA,nDup);
NC_totRx = zeros(nMu,nAtt,nA,nDup);
NC_tot_no_jointsRx = zeros(nMu,nAtt,nA,nDup);
oversRx = zeros(nMu,nAtt,nA,nDup);
forcRx = zeros(nMu,nAtt,nA,nDup);
sijRx = zeros(nMu,nAtt,nA,nDup);
EelasRx = zeros(nMu,nAtt,nA,nDup);

r = [rxStrain 500];
for ii=1:nDup
    for i=1:nMu
        for k=1:nAtt
            for j=1:nA
                display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
                filePrefix = [dataPath,shape,'_rp75_rand_shear_shear_',dup.value{ii},'_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_relaxed_'];
                % number of contacts
                NCRx(i,k,j,ii) = process_NC(filePrefix, rxStrain, 2, r');
                % contact statistics
                [NC_totRx(i,k,j,ii), NC_tot_no_jointsRx(i,k,j,ii), oversRx(i,k,j,ii),forcRx(i,k,j,ii),sijRx(i,k,j,ii)] = process_contactStat(filePrefix, rxStrain, 2, r');
                % elastic energy statistics
                EelasRx(i,k,j,ii) = process_elastic(filePrefix, rxStrain, 2, [rxStrain 500]');
                % normalize
                oversRx(i,k,j,ii) = oversRx(i,k,j,ii)/nfib;
            end
        end
    end
end

%%
sijRx(sijRx == 0) = NaN;

A = EelasRx(2:end,:,:,:);
A (A == 0) = nan;
EelasRx(2:end,:,:,:) = A;

% average over duplicate
EelasRx_mean = nanmean(EelasRx,4);
EelasRx_std = nanstd(EelasRx,0,4);

% replace data with mean
EelasRx = EelasRx_mean;

NC_tot_no_jointsD = NC_tot_no_jointsRx-NC_tot_no_joints;
%%{
NCRxObj = data3D('Rx-DR Unbroken $N_C$',NCRx);
NC_totRxObj = data3D('Rx-DR $N_C$',NC_totRx);
NC_tot_no_jointsRxObj = data3D('Rx-DR $N_C$ no joints',NC_tot_no_jointsRx);
oversRxObj = data3D('Rx-DR Overlaps',oversRx);
forcRxObj = data3D('Rx-DR $F^N$',forcRx);
sijRxObj = data3D('Rx-DR $h_{ij}$',sijRx);
NC_brokenRxObj = data3D('Rx-DR Broken $N_C$',NC_tot_no_jointsRx-NCRx);
EelasRxObj = data3D('Rx-DR $E_{el}$',EelasRx_mean);
NC_tot_no_jointsDObj = data3D('Rx-DR - DR $N_C$ no joints',NC_tot_no_jointsD);

plot3dim1E(EelasRxObj,EelasRx_std,attC,muC,aC);
% plot3dim1(NCRxObj,attC,muC,aC);
% plot3dim1(NC_totRxObj,attC,muC,aC);
% plot3dim1(NC_tot_no_jointsRxObj,attC,muC,aC);
% plot3dim1(NC_brokenRxObj,attC,muC,aC);
% plot3dim1(oversRxObj,attC,muC,aC);
% plot3dim1(forcRxObj,attC,muC,aC);
% plot3dim1(sijRxObj,attC,muC,aC);
% scatter(Fzero(1,:),Fzero(2,:),30,'x','MarkerEdgeColor',rgb('Lime'),'linewidth',2)
% plot3dim1(NC_tot_no_jointsDObj,attC,muC,aC);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part V - calculate property difference due to redispersion cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Deta = zeros(nMu,nAtt,nA);
DNC = zeros(nMu,nAtt,nA);
DNC_tot = zeros(nMu,nAtt,nA);
DNC_tot_no_joints = zeros(nMu,nAtt,nA);
Dovers = zeros(nMu,nAtt,nA);
Dforc = zeros(nMu,nAtt,nA);
Dsij = zeros(nMu,nAtt,nA);

DNCRx = zeros(nMu,nAtt,nA);
DNC_totRx = zeros(nMu,nAtt,nA);
DNC_tot_no_jointsRx = zeros(nMu,nAtt,nA);
DoversRx = zeros(nMu,nAtt,nA);
DforcRx = zeros(nMu,nAtt,nA);
DsijRx = zeros(nMu,nAtt,nA);
DEelasRx = zeros(nMu,nAtt,nA);
DNC_tot_no_jointsD = zeros(nMu,nAtt,nA);

Deta_std = zeros(nMu,nAtt,nA);
DEelasRx_std = zeros(nMu,nAtt,nA);
for i=1:nA
    Deta(:,:,i) = eta(:,:,i) - etaB;
    %     DNC(:,:,i) = NC(:,:,i) - NCB;
    %     DNC_tot(:,:,i) = NC_tot(:,:,i) - NC_totB;
    %     DNC_tot_no_joints(:,:,i) = NC_tot_no_joints(:,:,i) - NC_tot_no_jointsB;
    %     Dovers(:,:,i) = overs(:,:,i) - oversB;
    %     Dforc(:,:,i) = forc(:,:,i) - forcB;
    %     Dsij(:,:,i) = sij(:,:,i) - sijB;
    %     DNC_tot_no_jointsD(:,:,i) = NC_tot_no_jointsD(:,:,i) - NC_tot_no_jointsDB;
    
    Deta_std(:,:,i) = (eta_std(:,:,i).^2+etaB_std.^2).^0.5;
    DEelasRx_std(:,:,i) = (EelasRx_std(:,:,i).^2+EelasRxB_std.^2).^0.5;
    
    %     DNCRx(:,:,i) = NCRx(:,:,i) - NCRxB;
    %     DNC_totRx(:,:,i) = NC_totRx(:,:,i) - NC_totRxB;
    %     DNC_tot_no_jointsRx(:,:,i) = NC_tot_no_jointsRx(:,:,i) - NC_tot_no_jointsRxB;
    %     DoversRx(:,:,i) = oversRx(:,:,i) - oversRxB;
    %     DforcRx(:,:,i) = forcRx(:,:,i) - forcRxB;
    %     DsijRx(:,:,i) = sijRx(:,:,i) - sijRxB;
    DEelasRx(:,:,i) = EelasRx(:,:,i) - EelasRxB;
end

DetaObj = data3D('$\Delta \eta/\eta_0$',Deta);
DNCObj = data3D('$\Delta$ Unbroken $N_C$',DNC);
DNC_totObj = data3D('$\Delta N_C$',DNC_tot);
DNC_tot_no_jointsObj = data3D('$\Delta N_C$ no joints',DNC_tot_no_joints);
DoversObj = data3D('$\Delta$ Overlaps',Dovers);
DforcObj = data3D('$\Delta F_N$',Dforc);
DsijObj = data3D('$\Delta h_{ij}$',Dsij);
DNC_tot_no_jointsDObj = data3D('$\Delta \delta N_C$ no joints',DNC_tot_no_jointsD);


DEelasRxObj = data3D('$\Delta E_{elas} / E_{scale}$',DEelasRx);
DNCRxObj = data3D('$\Delta_{Rx}$ Unbroken $N_C$',DNCRx);
DNC_totRxObj = data3D('$\Delta_{Rx} N_C$',DNC_totRx);
DNC_tot_no_jointsRxObj = data3D('$\Delta_{Rx} N_C$ no joints',DNC_tot_no_jointsRx);
DoversRxObj = data3D('$\Delta_{Rx}$ Overlaps',DoversRx);
DforcRxObj = data3D('$\Delta_{Rx} F_N$',DforcRx);
DsijRxObj = data3D('$\Delta_{Rx} h_{ij}$',DsijRx);


%%{
% close all
plot3dim1E(DetaObj,Deta_std,attC,muC,aC);
% plot3dim1(DNCObj,attC,muC,aC);
% plot3dim1(DNC_totObj,attC,muC,aC);
% plot3dim1(DNC_tot_no_jointsObj,attC,muC,aC);
% plot3dim1(DoversObj,attC,muC,aC);
% plot3dim1(DforcObj,attC,muC,aC);
% plot3dim1(DsijObj,attC,muC,aC);
% plot3dim1(DNC_tot_no_jointsDObj,attC,muC,aC);
%%
% close all
plot3dim1E(DEelasRxObj,DEelasRx_std,attC,muC,aC);
% plot3dim1(DEelasRxObj,attC,muC,aC);
% plot3dim1(DNCRxObj,attC,muC,aC);
% plot3dim1(DNC_totRxObj,attC,muC,aC);
% plot3dim1(DNC_tot_no_jointsRxObj,attC,muC,aC);
% plot3dim1(DoversRxObj,attC,muC,aC);
% plot3dim1(DforcRxObj,attC,muC,aC);
% plot3dim1(DsijRxObj,attC,muC,aC);

%}

% Phase diagram with interfiber forces
%{
muC.name = '$\mu$';
plotPhase(Deta,'$\Delta \eta/\eta_0$',attC,muC,aC,-0.05)
%}