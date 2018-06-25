function [ interval_stat ] = interval_slope(strain,A,r)
% Final average of A in each region

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
    end
end

% find averages and standard deviations
for n=1:N
    [interval_stat(n,1), conf] = regress(strain(index(n,1):index(n,2)),...
        A(index(n,1):index(n,2)));
    interval_stat(n,2) = conf(2)-conf(1);
end

end