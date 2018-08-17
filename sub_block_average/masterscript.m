%% analysis for redispersion - includes rheology, contacts, elastic energy
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;
blockSize = 500;
basisStrain = 1000;

blockData(nMu,nAtt,nA) = strainValue();

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
%         figure('Units','Inches','Position',[1 1 3.0 2.5]);
%         hold on
%         title(['Basis: $\mu$ ',num2str(muC.value(i))])
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
%         legend(attC.legend)
%         ylabel('$\eta/\eta_0$')
%         set(gca,'yscale','log')
%         xlabel('$\gamma$')
end

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
            [r, block_strain, block_etarel] = block_stress_lbox_eta(filePrefix, blockSize, 5, sidex, nfib, nseg, rps, kb, a, EY, Imom, eta0);
            blockData(i,k,j).name = ['$(\mu,A_N)=($',num2str(muArr(i)),',',num2str(attArr(k)),')'];
            blockData(i,k,j).block_etarel = block_etarel;
            blockData(i,k,j).block_strain = block_strain;
            blockData(i,k,j).block_ndata = length(block_strain);
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
for i=1:nMu
    for k=1:nAtt
        for j=1:nA
            blockData(i,k,j).block_detarel = blockData(i,k,j).block_etarel - etaDataB(i,k);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part IV - plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legendArr = cell(nMu*nAtt*nA,1);
ind = 1; 


for i=1:nMu
    figure('Units','Inches','Position',[1 1 4.5 3.0]);
    hold on
    box on
    title(muC.legend{i})
    for k=1:nAtt
%             if blockData(i,k,j).block_ndata == 0
        for j=1:nA
%                 continue; 
%             end
            plot(blockData(i,k,j).block_strain, blockData(i,k,j).block_detarel,'o')
%             legendArr{ind} = blockData(i,k,j).name;
%             ind = ind + 1;
        end
    end
    legend(attC.legend,'location','best')
    xlabel('$\gamma$')
    ylabel('$\Delta \eta/\eta_0$')
end


% legend(legendArr{1:ind-1})