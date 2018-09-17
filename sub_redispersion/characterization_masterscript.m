%% analysis for redispersion - include intensity and cluster analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;
basisStrain = 1000;
flocCutoff = 0.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
IB = zeros(nMu,nAtt);
nClusterB = zeros(nMu,nAtt);
maxClusterB = zeros(nMu,nAtt);
RgB = zeros(nMu,nAtt);
volB = zeros(nMu,nAtt);

r = 0;
for i=1:nMu
    for k=1:nAtt
        
        display(['Processing basis mu',num2str(muArr(i)),'_att',num2str(attArr(k))]);
        filePrefix = [dataPathBasis,shape,'_rp75_basis_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_'];
        
        %% read input and process data
        IB(i,k) = process_intensity(filePrefix, basisStrain, 2, r);
        % number of clusters
        [nClusterB(i,k), maxClusterB(i,k), RgB(i,k), volB(i,k)] = process_cluster(filePrefix, nfib*flocCutoff);
        
        % normalize
        nClusterB(i,k) = nClusterB(i,k)/nfib;
        maxClusterB(i,k) = maxClusterB(i,k)/nfib;
        
        if isfinite(RgB(i,k))
            RgB(i,k) = RgB(i,k)/sidex;
            volB(i,k) =  volB(i,k)/sidex^3;
        end
        
    end
end

IB (IB == 0) = NaN;
nClusterB (nClusterB == 0) = NaN;
maxClusterB (maxClusterB == 0) = NaN;

IBObj = data3D('ND $I$',IB);
nClusterBObj = data3D('ND $N_{cluster}/N_{fib}$',nClusterB);
maxClusterBObj = data3D('ND max $S/N_{fib}$',maxClusterB);
RgBObj = data3D('ND $R_g/L_{box}$ for max $S$',RgB);
volBObj = data3D('ND $8\lambda_1\lambda_2\lambda_3/L_{box}^3$',volB);

plot2dim1(nClusterBObj,attC,muC);
set(gca, 'yscale','log')
plot2dim1(maxClusterBObj,attC,muC);
plot2dim1(IBObj,attC,muC);
plot2dim1(RgBObj,attC,muC);
plot2dim1(volBObj,attC,muC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain redispersed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = zeros(nMu,nAtt,nA);
nCluster = zeros(nMu,nAtt,nA);
maxCluster = zeros(nMu,nAtt,nA);
Rg = zeros(nMu,nAtt,nA);
vol = zeros(nMu,nAtt,nA);

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
            I(i,k,j) = process_intensity(filePrefix, basisStrain, 5, r);
            % number of clusters
            [nCluster(i,k,j), maxCluster(i,k,j), Rg(i,k,j), vol(i,k,j)] = process_cluster(filePrefix, nfib*flocCutoff);
            
            % normalize
            nCluster(i,k,j) = nCluster(i,k,j)/nfib;
            maxCluster(i,k,j) = maxCluster(i,k,j)/nfib;
            sidex = sidex(end);
            
            if isfinite(Rg(i,k,j))
                Rg(i,k,j) = Rg(i,k,j)/sidex;
                vol(i,k,j) =  vol(i,k,j)/sidex^3;
            end
        end
    end
end

I (I == 0) = NaN;
nCluster (nCluster == 0) = NaN;
maxCluster (maxCluster == 0) = NaN;
Rg (Rg == 0) = NaN;
vol (vol == 0) = NaN;

IObj = data3D('DR $I$',I);
nClusterObj = data3D('DR $N_{cluster}/N_{fib}$',nCluster);
maxClusterObj = data3D('DR max $S/N_{fib}$',maxCluster);
RgObj = data3D('DR $R_g/L_{box}$ for max $S$',Rg);
volObj = data3D('DR $8\lambda_1\lambda_2\lambda_3/L_{box}^3$',vol);

plot3dim1(nClusterObj,attC,muC,aC);
set(gca, 'yscale','log')
plot3dim1(maxClusterObj,attC,muC,aC);
plot3dim1(IObj,attC,muC,aC);
plot3dim1(RgObj,attC,muC,aC);
plot3dim1(volObj,attC,muC,aC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part III - calculate property difference due to redispersion cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DI = zeros(nMu,nAtt,nA);
DnCluster = zeros(nMu,nAtt,nA);
DmaxCluster = zeros(nMu,nAtt,nA);
DRg = zeros(nMu,nAtt,nA);
Dvol = zeros(nMu,nAtt,nA);

for i=1:nA
    DI(:,:,i) = I(:,:,i) - IB;
    DnCluster(:,:,i) = nCluster(:,:,i) - nClusterB;
    DmaxCluster(:,:,i) = maxCluster(:,:,i) - maxClusterB;
    DRg(:,:,i) = Rg(:,:,i) - RgB;
    Dvol(:,:,i) = vol(:,:,i) - volB;
end
DIObj = data3D('$\Delta I$',DI);
DnClusterObj = data3D('$\Delta N_{cluster}/N_{fib}$',DnCluster);
DmaxClusterObj = data3D('$\Delta S/N_{fib}$',DmaxCluster);
DRgObj = data3D('$\Delta R_g/L_{box}$',DRg);
DvolObj = data3D('$\Delta 8\lambda_1\lambda_2\lambda_3/L_{box}^3$',Dvol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part IV - plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot3dim1(DnClusterObj,attC,muC,aC);
plot3dim1(DmaxClusterObj,attC,muC,aC);
plot3dim1(DIObj,attC,muC,aC);
plot3dim1(DRgObj,attC,muC,aC);
plot3dim1(DvolObj,attC,muC,aC);