%% analysis for redispersion - include intensity and cluster analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Notations
%%%%%%%%%%%%%%%%%%%%%%%%%%
% A can be any below
% p - DR property
% pB - ND property
% pRx - DR relaxed property
% pRxB - ND relaxed property
% Dp - DR-ND
% Dp_std - std of Dp

% AObj - object that stores data name and value
% A_mean - mean of A
% A_std - std of A

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;
basisStrain = 1000;
flocCutoff = 0.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
IB = zeros(nMu,nAtt,nDup);
nClusterB = zeros(nMu,nAtt,nDup);
maxClusterB = zeros(nMu,nAtt,nDup);
RgB = zeros(nMu,nAtt,nDup);
volB = zeros(nMu,nAtt,nDup);

r = 0;
for ii=1:nDup
    for i=1:nMu
        for k=1:nAtt
            
            display(['Processing basis mu',num2str(muArr(i)),'_att',num2str(attArr(k))]);
            filePrefix = [dataPathBasis,shape,'_rp75_basis_',dup.value{ii},'_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_'];
            
            %% read input and process data
            IB(i,k,ii) = process_intensity(filePrefix, basisStrain, 2, r);
            
            % number of clusters
            [nClusterB(i,k,ii), maxClusterB(i,k,ii), RgB(i,k,ii), volB(i,k,ii)] = process_cluster(filePrefix, nfib*flocCutoff);
            
            % normalize
            nClusterB(i,k,ii) = nClusterB(i,k,ii)/nfib;
            maxClusterB(i,k,ii) = maxClusterB(i,k,ii)/nfib;
            
            if isfinite(RgB(i,k,ii))
                RgB(i,k,ii)= RgB(i,k,ii)/sidex;
                volB(i,k,ii) =  volB(i,k,ii)/sidex^3;
            end
            
        end
    end
end

IB (IB == 0) = NaN;
nClusterB (nClusterB == 0) = NaN;
maxClusterB (maxClusterB == 0) = NaN;

% average over duplicate
IB_mean = nanmean(IB,3);
IB_std = nanstd(IB,0,3);
nClusterB_mean = nanmean(nClusterB,3);
nClusterB_std = nanstd(nClusterB,0,3);
maxClusterB_mean = nanmean(maxClusterB,3);
maxClusterB_std = nanstd(maxClusterB,0,3);

% replace data with mean
IB = IB_mean;
nClusterB = nClusterB_mean;
maxClusterB = maxClusterB_mean;

IBObj = data3D('ND $I$',IB);
nClusterBObj = data3D('ND $N_{cluster}/N_{fib}$',nClusterB);
maxClusterBObj = data3D('ND max $S/N_{fib}$',maxClusterB);
% RgBObj = data3D('ND $R_g/L_{box}$ for max $S$',RgB);
% volBObj = data3D('ND $8\lambda_1\lambda_2\lambda_3/L_{box}^3$',volB);

plot2dim1EB(nClusterBObj,nClusterB_std,attC,muC);
% set(gca, 'yscale','log')
plot2dim1EB(maxClusterBObj,maxClusterB_std,attC,muC);
plot2dim1EB(IBObj,IB_std,attC,muC);
% plot2dim1(RgBObj,attC,muC);
% plot2dim1(volBObj,attC,muC);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain redispersed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = zeros(nMu,nAtt,nA,nDup);
nCluster = zeros(nMu,nAtt,nA,nDup);
maxCluster = zeros(nMu,nAtt,nA,nDup);
Rg = zeros(nMu,nAtt,nA,nDup);
vol = zeros(nMu,nAtt,nA,nDup);

for ii=1:nDup
    for i=1:nMu
        for k=1:nAtt
            for j=1:nA
                display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
                filePrefix = [dataPath,shape,'_rp75_rand_shear_shear_',dup.value{ii},'_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_'];
                
                if (exist([filePrefix,'Stress_tensor.txt'], 'file') == 0)
                    continue;
                end
                
                %% read input and process data
                [box_strain, sidex] = read_box([filePrefix,'Lbox.txt']);
                % define intervals
                r = round(intervals(box_strain,sidex),0)';
                % intensity energy statistics
                I(i,k,j,ii) = process_intensity(filePrefix, basisStrain, 5, r);
                % number of clusters
                [nCluster(i,k,j,ii), maxCluster(i,k,j,ii), Rg(i,k,j,ii), vol(i,k,j,ii)] = process_cluster(filePrefix, nfib*flocCutoff);
                
                % normalize
                nCluster(i,k,j,ii) = nCluster(i,k,j,ii)/nfib;
                maxCluster(i,k,j,ii) = maxCluster(i,k,j,ii)/nfib;
                sidex = sidex(end);
                
                if isfinite(Rg(i,k,j,ii))
                    Rg(i,k,j,ii)= Rg(i,k,j,ii)/sidex;
                    vol(i,k,j,ii) =  vol(i,k,j,ii)/sidex^3;
                end
            end
        end
    end
end

I (I == 0) = NaN;
nCluster (nCluster == 0) = NaN;
maxCluster (maxCluster == 0) = NaN;
% Rg (Rg == 0) = NaN;
% vol (vol == 0) = NaN;

I_mean = nanmean(I,4);
I_std = nanstd(I,0,4);
nCluster_mean = nanmean(nCluster,4);
nCluster_std = nanstd(nCluster,0,4);
maxCluster_mean = nanmean(maxCluster,4);
maxCluster_std = nanstd(maxCluster,0,4);

% replace data with mean
I = I_mean;
nCluster = nCluster_mean;
maxCluster = maxCluster_mean;

IObj = data3D('DR $I$',I);
nClusterObj = data3D('DR $N_{cluster}/N_{fib}$',nCluster);
maxClusterObj = data3D('DR max $S/N_{fib}$',maxCluster);
% RgObj = data3D('DR $R_g/L_{box}$ for max $S$',Rg);
% volObj = data3D('DR $8\lambda_1\lambda_2\lambda_3/L_{box}^3$',vol);

plot3dim1E(nClusterObj,nCluster_std,attC,muC,aC);
% set(gca, 'yscale','log')
plot3dim1E(maxClusterObj,maxCluster_std,attC,muC,aC);
plot3dim1E(IObj,I_std,attC,muC,aC);
% plot3dim1(RgObj,attC,muC,aC);
% plot3dim1(volObj,attC,muC,aC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part III - calculate property difference due to redispersion cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DI = zeros(nMu,nAtt,nA);
DnCluster = zeros(nMu,nAtt,nA);
DmaxCluster = zeros(nMu,nAtt,nA);
% DRg = zeros(nMu,nAtt,nA);
% Dvol = zeros(nMu,nAtt,nA);

DI_std = zeros(nMu,nAtt,nA);
DnCluster_std = zeros(nMu,nAtt,nA);
DmaxCluster_std = zeros(nMu,nAtt,nA);

for i=1:nA
    DI(:,:,i) = I(:,:,i) - IB;
    DnCluster(:,:,i) = nCluster(:,:,i) - nClusterB;
    DmaxCluster(:,:,i) = maxCluster(:,:,i) - maxClusterB;
    
    DI_std(:,:,i) = (I_std(:,:,i).^2+IB_std.^2).^0.5;
    DnCluster_std(:,:,i) = (nCluster_std(:,:,i).^2+nClusterB_std.^2).^0.5;
    DmaxCluster_std(:,:,i) = (maxCluster_std(:,:,i).^2+maxClusterB_std.^2).^0.5;
    
    %     DRg(:,:,i) = Rg(:,:,i) - RgB;
    %     Dvol(:,:,i) = vol(:,:,i) - volB;
end
DIObj = data3D('$\Delta I$',DI);
DnClusterObj = data3D('$\Delta N_{cluster}/N_{fib}$',DnCluster);
DmaxClusterObj = data3D('$\Delta S/N_{fib}$',DmaxCluster);
% DRgObj = data3D('$\Delta R_g/L_{box}$',DRg);
% DvolObj = data3D('$\Delta 8\lambda_1\lambda_2\lambda_3/L_{box}^3$',Dvol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part IV - plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot3dim1E(DnClusterObj,DnCluster_std,attC,muC,aC);
plot3dim1E(DmaxClusterObj,DmaxCluster_std,attC,muC,aC);
plot3dim1E(DIObj,DI_std,attC,muC,aC);
% plot3dim1(DRgObj,attC,muC,aC);
% plot3dim1(DvolObj,attC,muC,aC);