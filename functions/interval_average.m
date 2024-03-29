function [ interval_stat ] = interval_average(strain,A,r)
% Final average of A in each region of redispersion

N = size(r,1); % number of intervals
interval_stat = zeros(N,2); % average and standard deviation
index = zeros(N,2); % bracketing index

% find bracketing index
index(1,1) = 1;
n = 1;
for i=1:length(strain)
    if strain(i) == r(n)
        index(n,2) = i;
        if n ~= N 
            index(n+1,1)=i+1;
        end
        n = n+1;
        if n > N
            break;
        end
    end
end
% index
% find averages and standard error of the mean
for n=1:N
    interval_stat(n,1) = mean(A(index(n,1):index(n,2)));
    % standard error of the mean
    interval_stat(n,2) = std(A(index(n,1):index(n,2)))/sqrt(length(A(index(n,1):index(n,2))));
    % standard deviation
%     interval_stat(n,2) = std(A(index(n,1):index(n,2)));
end

end