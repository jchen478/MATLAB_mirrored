function [ nCluster, S, FlocE] = process_cluster( filePrefix)
%PROCESS_INTENSITY Find stats for intensity
nCluster = 0;
S = 0;
FlocE = 0; 
if (exist([filePrefix,'nCluster.txt'], 'file') ~= 0)
    [nCluster] = read_ncluster([filePrefix,'nCluster.txt']);
end
if (exist([filePrefix,'Cluster_results.txt'], 'file') ~= 0)
    [Sarr] = read_cluster_results([filePrefix,'Cluster_results.txt'],nCluster);
    Ssum = sum(Sarr.^2);
    S = sqrt(Ssum)/nCluster;
end
if (exist([filePrefix,'flocEnergy.txt'], 'file') ~= 0)
    [FlocEarr] = read_flocEnergy([filePrefix,'flocEnergy.txt'],nCluster);
    FlocESum = sum(FlocEarr.^2);
    FlocE = sqrt(FlocESum)/nCluster;
end
end

