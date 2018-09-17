function [ nCluster, maxCluster, Rg, vol] = process_cluster(filePrefix, cutoff)
%PROCESS_INTENSITY Find stats for intensity
nCluster = 0;
maxCluster = 0;
Rg = 0;
vol = 0;
if (exist([filePrefix,'nCluster.txt'], 'file') ~= 0)
    [strain, nCluster, maxCluster, Rg, lam] = read_ncluster([filePrefix,'nCluster.txt']);
    
    match = maxCluster >= cutoff;
    % if cluster at any time is smaller than threshold
    if sum(match) ~= length(match)
        Rg = NaN;
        vol = NaN;
    else
        Rg = mean(Rg);
        vol = mean(lam(:,1).*lam(:,2).*lam(:,3));
    end
    maxCluster = mean(maxCluster);
    nCluster = mean(nCluster);
end

end