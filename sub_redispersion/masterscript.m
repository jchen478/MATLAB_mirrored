%% analysis for redispersion - includes rheology, contacts, elastic energy
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;
defaultFormat
basisStrain = 1000;
rxStrain = 300;
Fzero = [9 20 30 35 50 ; 0.2504 0.1306 0.09691 0.08558 0.06155];
Fmax = [9 20 30 35 50 ; 0.2747 0.1808 0.1526 0.1432 0.1236];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
etaB = zeros(nMu,nAtt);
NCB = zeros(nMu,nAtt);
NC_totB = zeros(nMu,nAtt);
NC_tot_no_jointsB = zeros(nMu,nAtt);
oversB = zeros(nMu,nAtt);
forcB = zeros(nMu,nAtt);
sijB = zeros(nMu,nAtt);
etaBSE = zeros(nMu,nAtt);

for i=1:nMu
    for k=1:nAtt
        display(['Processing basis mu',num2str(muArr(i)),'_att',num2str(attArr(k))]);
        filePrefix = [dataPathBasis,shape,'_rp75_basis_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_'];
        %% read input and process data
        [r, etaB(i,k), etaBSE(i,k)] = process_stress_lbox_eta_basis(filePrefix, basisStrain, 2, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0);
        % number of contacts
        NCB(i,k) = process_NC(filePrefix, basisStrain, 2, r');
        % contact statistics
        [NC_totB(i,k), NC_tot_no_jointsB(i,k), oversB(i,k),forcB(i,k),sijB(i,k)] = process_contactStat(filePrefix, basisStrain, 2, r');
        % normalize
        oversB(i,k) = oversB(i,k)/nfib;
    end
end

%%{
etaBObj = data3D('ND $\eta/\eta_0$',etaB);
NCBObj = data3D('ND Unbroken $N_C$',NCB);
NC_totBObj = data3D('ND $N_C$',NC_totB);
NC_tot_no_jointsBObj = data3D('ND $N_C$ no joints',NC_tot_no_jointsB);
oversBObj = data3D('ND Overlaps',oversB);
forcBObj = data3D('ND $F^N$',forcB);
sijBObj = data3D('ND $s_{ij}$',sijB);
NC_brokenBObj = data3D('ND Broken $N_C$',NC_tot_no_jointsB-NCB);

% plot2dim1(etaBObj,attC,muC);
% plot2dim1(NCBObj,attC,muC);
% plot2dim1(NC_totBObj,attC,muC);
% plot2dim1(NC_tot_no_jointsBObj,attC,muC);
% plot2dim1(NC_brokenBObj,attC,muC);
% plot2dim1(oversBObj,attC,muC);
% plot2dim1(forcBObj,attC,muC);
% plot2dim1(sijBObj,attC,muC);
% scatter(Fzero(1,:),Fzero(2,:),30,'x','MarkerEdgeColor',rgb('Lime'),'linewidth',2)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain relaxed basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
NCRxB = zeros(nMu,nAtt);
NC_totRxB = zeros(nMu,nAtt);
NC_tot_no_jointsRxB = zeros(nMu,nAtt);
oversRxB = zeros(nMu,nAtt);
forcRxB = zeros(nMu,nAtt);
sijRxB = zeros(nMu,nAtt);
EelasRxB = zeros(nMu,nAtt);

r = [rxStrain 500];
for i=1:nMu
    for k=1:nAtt
        display(['Processing basis mu',num2str(muArr(i)),'_att',num2str(attArr(k))]);
        filePrefix = [dataPathBasis,shape,'_rp75_basis_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_relaxed_'];
        % number of contacts
        NCRxB(i,k) = process_NC(filePrefix, rxStrain, 2, r');
        % contact statistics
        [NC_totRxB(i,k), NC_tot_no_jointsRxB(i,k), oversRxB(i,k),forcRxB(i,k),sijRxB(i,k)] = process_contactStat(filePrefix, rxStrain, 2, r');
        % elastic energy statistics
        EelasRxB(i,k) = process_elastic(filePrefix, rxStrain, 2, r');
        % normalize
        oversRxB(i,k) = oversRxB(i,k)/nfib;
    end
end
NC_tot_no_jointsDB = NC_tot_no_jointsRxB-NC_tot_no_jointsB;

sijRxB(sijRxB == 0) = NaN;

%%{
EelasRxBObj = data3D('Rx-ND $E_{el}$',EelasRxB);
NCRxBObj = data3D('Rx-ND Unbroken $N_C$',NCRxB);
NC_totRxBObj = data3D('Rx-ND $N_C$',NC_totRxB);
NC_tot_no_jointsRxBObj = data3D('Rx-ND $N_C$ no joints',NC_tot_no_jointsRxB);
oversRxBObj = data3D('Rx-ND Overlaps',oversRxB);
forcRxBObj = data3D('Rx-ND $F^N$',forcRxB);
sijRxBObj = data3D('Rx-ND $s_{ij}$',sijRxB);
NC_brokenRxBObj = data3D('Rx-ND Broken $N_C$',NC_tot_no_jointsRxB-NCRxB);
NC_tot_no_jointsDBObj = data3D('Rx-ND - ND $N_C$ no joints',NC_tot_no_jointsDB);

% plot2dim1(EelasRxBObj,attC,muC);
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
eta = zeros(nMu,nAtt,nA);
NC = zeros(nMu,nAtt,nA);
NC_tot = zeros(nMu,nAtt,nA);
NC_tot_no_joints = zeros(nMu,nAtt,nA);
overs = zeros(nMu,nAtt,nA);
forc = zeros(nMu,nAtt,nA);
sij = zeros(nMu,nAtt,nA);
etaSE = zeros(nMu,nAtt,nA);

for i=1:nMu
    for k=1:nAtt
        for j=1:nA
            display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
            filePrefix = [dataPath,shape,'_rp75_rand_shear_shear_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_'];
            %% read input and process dat
            % viscosity
            [r, eta(i,k,j), etaSE(i,k,j)] = process_stress_lbox_eta(filePrefix, 5, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0);
            % number of contacts
            NC(i,k,j) = process_NC(filePrefix, basisStrain, 5, r');
            % contact statistics
            [NC_tot(i,k,j), NC_tot_no_joints(i,k,j), overs(i,k,j),forc(i,k,j),sij(i,k,j)] = process_contactStat(filePrefix, basisStrain, 5, r');
            % normalize
            overs(i,k,j) = overs(i,k,j)/nfib;
        end
    end
end

eta(eta == 0) = NaN;
NC_tot(NC_tot == 0) = NaN;
NC_tot_no_joints(NC_tot_no_joints == 0) = NaN;
sij(sij == 0) = NaN;

%%{
etaObj = data3D('DR $\eta/\eta_0$',eta);
NCObj = data3D('DR Unbroken $N_C$',NC);
NC_totObj = data3D('DR $N_C$',NC_tot);
NC_tot_no_jointsObj = data3D('DR $N_C$ no joints',NC_tot_no_joints);
oversObj = data3D('DR Overlaps',overs);
forcObj = data3D('DR $F^N$',forc);
sijObj = data3D('DR $s_{ij}$',sij);
NC_brokenObj = data3D('DR Broken $N_C$',NC_tot_no_joints-NC);

% plot3dim1(etaObj,attC,muC,aC);
% plot3dim1(NCObj,attC,muC,aC);
% plot3dim1(NC_totObj,attC,muC,aC);
% plot3dim1(NC_tot_no_jointsObj,attC,muC,aC);
% plot3dim1(NC_brokenObj,attC,muC,aC);
% plot3dim1(oversObj,attC,muC,aC);
% plot3dim1(forcObj,attC,muC,aC);
% plot3dim1(sijObj,attC,muC,aC);
% scatter(Fzero(1,:),Fzero(2,:),30,'x','MarkerEdgeColor',rgb('Lime'),'linewidth',2)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part IV - obtain relaxed redispersion value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NCRx = zeros(nMu,nAtt,nA);
NC_totRx = zeros(nMu,nAtt,nA);
NC_tot_no_jointsRx = zeros(nMu,nAtt,nA);
oversRx = zeros(nMu,nAtt,nA);
forcRx = zeros(nMu,nAtt,nA);
sijRx = zeros(nMu,nAtt,nA);
EelasRx = zeros(nMu,nAtt,nA);

r = [rxStrain 500];
for i=1:nMu
    for k=1:nAtt
        for j=1:nA
            display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
            filePrefix = [dataPath,shape,'_rp75_rand_shear_shear_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_relaxed_'];
            % number of contacts
            NCRx(i,k,j) = process_NC(filePrefix, rxStrain, 2, r');
            % contact statistics
            [NC_totRx(i,k,j), NC_tot_no_jointsRx(i,k,j), oversRx(i,k,j),forcRx(i,k,j),sijRx(i,k,j)] = process_contactStat(filePrefix, rxStrain, 2, r');
            % elastic energy statistics
            EelasRx(i,k,j) = process_elastic(filePrefix, rxStrain, 2, [rxStrain 500]');
            % normalize
            oversRx(i,k,j) = oversRx(i,k,j)/nfib;
        end
    end
end
NC_tot_no_jointsD = NC_tot_no_jointsRx-NC_tot_no_joints;  

sijRx(sijRx == 0) = NaN;

%%{
NCRxObj = data3D('Rx-DR Unbroken $N_C$',NCRx);
NC_totRxObj = data3D('Rx-DR $N_C$',NC_totRx);
NC_tot_no_jointsRxObj = data3D('Rx-DR $N_C$ no joints',NC_tot_no_jointsRx);
oversRxObj = data3D('Rx-DR Overlaps',oversRx);
forcRxObj = data3D('Rx-DR $F^N$',forcRx);
sijRxObj = data3D('Rx-DR $s_{ij}$',sijRx);
NC_brokenRxObj = data3D('Rx-DR Broken $N_C$',NC_tot_no_jointsRx-NCRx);
EelasRxObj = data3D('Rx-DR $E_{el}$',EelasRx);
NC_tot_no_jointsDObj = data3D('Rx-DR - DR $N_C$ no joints',NC_tot_no_jointsD);

% plot3dim1(EelasRxObj,attC,muC,aC);
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

for i=1:nA
    Deta(:,:,i) = eta(:,:,i) - etaB;
    DNC(:,:,i) = NC(:,:,i) - NCB;
    DNC_tot(:,:,i) = NC_tot(:,:,i) - NC_totB;
    DNC_tot_no_joints(:,:,i) = NC_tot_no_joints(:,:,i) - NC_tot_no_jointsB;
    Dovers(:,:,i) = overs(:,:,i) - oversB;
    Dforc(:,:,i) = forc(:,:,i) - forcB;
    Dsij(:,:,i) = sij(:,:,i) - sijB;
    DNC_tot_no_jointsD(:,:,i) = NC_tot_no_jointsD(:,:,i) - NC_tot_no_jointsDB;
        
        
    DNCRx(:,:,i) = NCRx(:,:,i) - NCRxB;
    DNC_totRx(:,:,i) = NC_totRx(:,:,i) - NC_totRxB;
    DNC_tot_no_jointsRx(:,:,i) = NC_tot_no_jointsRx(:,:,i) - NC_tot_no_jointsRxB;
    DoversRx(:,:,i) = oversRx(:,:,i) - oversRxB;
    DforcRx(:,:,i) = forcRx(:,:,i) - forcRxB;
    DsijRx(:,:,i) = sijRx(:,:,i) - sijRxB;
    DEelasRx(:,:,i) = EelasRx(:,:,i) - EelasRxB;
end

DetaObj = data3D('$\Delta \eta/\eta_0$',Deta);
DNCObj = data3D('$\Delta$ Unbroken $N_C$',DNC);
DNC_totObj = data3D('$\Delta N_C$',DNC_tot);
DNC_tot_no_jointsObj = data3D('$\Delta N_C$ no joints',DNC_tot_no_joints);
DoversObj = data3D('$\Delta$ Overlaps',Dovers);
DforcObj = data3D('$\Delta F_N$',Dforc);
DsijObj = data3D('$\Delta s_{ij}$',Dsij);
DNC_tot_no_jointsDObj = data3D('$\Delta \delta N_C$ no joints',DNC_tot_no_jointsD);
   

DEelasRxObj = data3D('$\Delta E_{elas} / E_{scale}$',DEelasRx);
DNCRxObj = data3D('$\Delta_{Rx}$ Unbroken $N_C$',DNCRx);
DNC_totRxObj = data3D('$\Delta_{Rx} N_C$',DNC_totRx);
DNC_tot_no_jointsRxObj = data3D('$\Delta_{Rx} N_C$ no joints',DNC_tot_no_jointsRx);
DoversRxObj = data3D('$\Delta_{Rx}$ Overlaps',DoversRx);
DforcRxObj = data3D('$\Delta_{Rx} F_N$',DforcRx);
DsijRxObj = data3D('$\Delta_{Rx} s_{ij}$',DsijRx);

%%{
% close all
plot3dim1(DetaObj,attC,muC,aC);
% plot3dim1(DNCObj,attC,muC,aC);
% plot3dim1(DNC_totObj,attC,muC,aC);
% plot3dim1(DNC_tot_no_jointsObj,attC,muC,aC);
% plot3dim1(DoversObj,attC,muC,aC);
% plot3dim1(DforcObj,attC,muC,aC);
% plot3dim1(DsijObj,attC,muC,aC);
% plot3dim1(DNC_tot_no_jointsDObj,attC,muC,aC);

% close all
plot3dim1(DEelasRxObj,attC,muC,aC);
% plot3dim1(DNCRxObj,attC,muC,aC);
% plot3dim1(DNC_totRxObj,attC,muC,aC);
% plot3dim1(DNC_tot_no_jointsRxObj,attC,muC,aC);
% plot3dim1(DoversRxObj,attC,muC,aC);
% plot3dim1(DforcRxObj,attC,muC,aC);
% plot3dim1(DsijRxObj,attC,muC,aC);

%}

% Phase diagram with interfiber forces
%%{
muC.name = '$\mu$';
plotPhase(Deta,'$\Delta \eta/\eta_0$',attC,muC,aC,-0.05)
%}