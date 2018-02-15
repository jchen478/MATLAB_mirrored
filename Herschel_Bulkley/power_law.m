function [ k, n ] = power_law( x, y, kinit, ninit )
%POWER_LAW Summary of this function goes here
%   Detailed explanation goes here

%% initialize parameters
iter = 100000;
theta = [0; kinit; ninit];
dJ = zeros(3,1);
J = zeros(iter,1);

%% parameters
alpha = 0.00000001; % learning rate
m = length(x); % number of samples

h = theta(1)+ theta(2)*x.^(theta(3));
%% gradient descent
for i = 1:iter
    
    h = theta(1)+ theta(2)*x.^(theta(3));
    dJ(2) = 1/m*sum(x.^theta(3).*(h-y));
    dJ(3) = 1/m*sum(theta(2)*theta(3)*x.^(theta(3)-1).*(h-y));
    
    theta = theta - alpha*dJ;
    
    h = theta(1)+ theta(2)*x.^(theta(3));
    J(i) = 1/(2*m)*sum((h-y).^2);
    
end

k = theta(2);
n = theta(3);

figure()
plot(J)

end

