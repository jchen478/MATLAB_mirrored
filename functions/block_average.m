function [ result ] = block_average( data, blockSize )
%BLOCK_AVERAGE Summary of this function goes here
%   Detailed explanation goes here
nblocks = floor(length(data) / blockSize); 
result = zeros(nblocks, 1);
for i=1:nblocks
    result(i) = mean(data((i-1)*blockSize+1:i*blockSize));
end
end

