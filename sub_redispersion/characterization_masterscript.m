%% analysis for redispersion - include intensity and cluster analysis
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;
basisStrain = 1000;
IData = zeros(nMu,nAtt,nA);
nClusterData = zeros(nMu,nAtt,nA);
SData = zeros(nMu,nAtt,nA);
flocEData = zeros(nMu,nAtt,nA);

IDataB = zeros(nMu,nAtt);

IData_high = zeros(nMu,nAtt,nA);
nClusterData_high = zeros(nMu,nAtt,nA);
SData_high = zeros(nMu,nAtt,nA);
flocEData_high = zeros(nMu,nAtt,nA);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0; 
for i=1:nMu
    for k=1:nAtt
        
        display(['Processing basis mu',num2str(muArr(i)),'_att',num2str(attArr(k))]);
        filePrefix = [dataPathBasis,shape,'_rp75_basis_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_'];

        %% read input and process data
        IDataB(i,k) = process_intensity(filePrefix, basisStrain, 2, r);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain redispersed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nMu
    for k=1:nAtt
        for j=1:nA
            display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
            filePrefix = [dataPath,shape,'_rp75_rand_shear_shear_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_'];
            
            %% read input and process data
            [box_strain, sidex] = read_box([filePrefix,'Lbox.txt']);
            % define intervals
            r = round(intervals(box_strain,sidex),0)';
            % intensity energy statistics
            IData(i,k,j) = process_intensity(filePrefix, basisStrain, 5, r);
            % number of clusters
            [nClusterData(i,k,j), SData(i,k,j), flocEData(i,k,j)] = process_cluster([filePrefix, '']);
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain redispersed value at high conc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nMu
    for k=1:nAtt
        for j=1:nA
            display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
            filePrefix = [dataPath,shape,'_rp75_rand_shear_shear_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_'];
            
            %% read input and process data
            [box_strain, sidex] = read_box([filePrefix,'Lbox.txt']);
            % define intervals
            r = round(intervals(box_strain,sidex),0)';
            % intensity energy statistics
            IData_high(i,k,j) = process_intensity([filePrefix,'high_'], basisStrain, 5, r);
            % number of clusters
            [nClusterData_high(i,k,j), SData_high(i,k,j), flocEData_high(i,k,j)] = process_cluster([[filePrefix,'high_'], '']);
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part III - calculate property difference due to redispersion cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DIData = zeros(nMu,nAtt,nA);
IData (IData == 0) = NaN;
IDataB (IDataB == 0) = NaN;
flocEData(flocEData == 0) = NaN;

SData_high(1,4:5,1) = NaN;
SData(1,4:5,1) = NaN;
nClusterData_high(1,4:5,1) = NaN;
nClusterData(1,4:5,1) = NaN;

for i=1:nA
    DIData(:,:,i) = IData(:,:,i) - IDataB;
end

IDataObj = data3D('$I$',IData); 
IDataBObj = data3D('$I_{ND}$',IDataB); 
SFracDataObj = data3D('Low $\bar{S}^*$',SData/nfib); 
flocEDataObj = data3D('Low $E_{floc}$',flocEData);
DIDataObj = data3D('$\Delta I$',DIData);
nClusterDataObj = data3D('Low $N_{cluster}$',nClusterData);

SFracDataHighObj = data3D('High $\bar{S}^*$',SData_high/nfib); 
flocEDataHighObj = data3D('High $E_{floc}$',flocEData_high);
nClusterDataHighObj = data3D('High $N_{cluster}$',nClusterData_high);

nClusterDiffLowHighObj = data3D('$\Delta N_{cluster}$',nClusterData-nClusterData_high);
nClusterDiffLowHighNormObj = data3D('$N_{cluster} \%$ Difference',100*(nClusterData-nClusterData_high)./nClusterData);
nClusterDiffLowHighRatioObj = data3D('High $N_{cluster} /$ Low $N_{cluster}$',nClusterData_high./nClusterData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part IV - plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot3dim1(SFracDataObj,attC,muC,aC);
plot3dim1(flocEDataObj,attC,muC,aC); 
% plot3dim1(nClusterDataObj,attC,muC,aC); 
% set(gca,'yscale','log')
% plot3dim1(DIDataObj,attC,muC,aC); 
% plot3dim1(IDataObj,attC,muC,aC); 
% plot2dim1(IDataBObj,attC,muC); 

% plot3dim1(SFracDataHighObj,attC,muC,aC); 
plot3dim1(flocEDataHighObj,attC,muC,aC); 
% plot3dim1(nClusterDataHighObj,attC,muC,aC); 
% set(gca,'yscale','log')

% plot3dim1(nClusterDiffLowHighObj,attC,muC,aC); 
% plot3dim1(nClusterDiffLowHighNormObj,attC,muC,aC); 
% plot3dim1(nClusterDiffLowHighRatioObj,attC,muC,aC);
% set(gca,'Ytick',0:1:6)