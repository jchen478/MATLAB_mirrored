clear
clc
load('nfib160_17_DC.mat')

%% check if matrix is symmetric
symmetric = max(max(A - A'));
if symmetric > 0
    disp('matrix not symmetric')
end

%% define variables
num = sum(A,2);
nfib = length(A);
cluster_access = (sum(A,2) ~= 1); 

%% loop through and combine
combine = false;
for m=1:nfib
    if(cluster_access(m) == 1)
        for n=1:nfib
            if ((cluster_access(n) == 1) && (m ~= n))
                for ind=1:nfib
                    if((A(n,ind) == 1) && (A(m,ind) == 1))
                        combine = true;
                        cluster_access(n) = 0;
                       % ind_start = ind;
                        break;
                    end
                end
                if(combine == true)
                   % disp(['m n pair ind_start ',num2str(m-1),' ', num2str(n-1),'  ',num2str(ind_start-1)])
                    A(m,n) = 1;
                    for ind=1:nfib
                        if (A(n,ind) == 1)
                            A(m,ind) = 1;
                        end
                    end
                    n = 1;
                    combine = false;
                end
            end
        end
    end
end



%% view cluster count
nCluster = 0;
clusterInfo = [];
cluster_count = cluster_access .* sum(A,2);
for m=1:nfib
    if (cluster_count(m) > 1)
        nCluster = nCluster + 1;
        clusterInfo(nCluster,1) = m;
        clusterInfo(nCluster,2) = cluster_count(m);
    end
end

for i=1:nCluster
    
    m = clusterInfo(i,1);
    pos = 3;
    for ind=1:nfib
        %A(1,11)
        % disp(['m ind all ',num2str(m-1),' ', num2str(ind-1),'  ',num2str(A(m,ind))])
        if A(m,ind) == 1%
            %disp('in cluster')
            clusterInfo(i,pos) = ind-1;
            pos = pos + 1;
        end
    end
    
end
C = A(1,:)';
B = clusterInfo(1,3:end);
B = B';
