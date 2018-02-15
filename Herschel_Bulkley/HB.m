clc;
clear;
close all;

%% Training set

x = [ 39.7887
   19.8944
    9.9472
    6.6315
    4.9736
    3.9789];

y = [ 208.6316
  159.0290
   84.5832
   82.1127
   62.2669
   50.0672];

m = length(x); 

%% initialize parameters 
iter = 50000;
theta = [0; 23; 0.8];
dJ = zeros(3,1); 
J = zeros(iter,1); 

%% initial guess plot
h = theta(1)+ theta(2)*x.^(theta(3)); 
figure()
hold on
scatter(x,y)
plot(x,h)

%% learning rate
alpha = 0.0000001; 

%% gradient descent
for i = 1:iter 

%    dJ(1) = 1/m*sum(h-y); 
   dJ(2) = 1/m*sum(x.^theta(3).*(h-y));
   dJ(3) = 1/m*sum(theta(2)*theta(3)*x.^(theta(3)-1).*(h-y)); 
%    dJ(2) = 1/m*sum(x.^theta(3).*(h-y)) + 2*theta(2);
%    dJ(3) = 1/m*sum(theta(2)*theta(3)*x.^(theta(3)-1).*(h-y)) + 2*theta(3); 
   
   theta = theta - alpha*dJ; 
   
   h = theta(1)+ theta(2)*x.^(theta(3)); 
   J(i) = 1/(2*m)*sum((h-y).^2);
    
end

figure()
plot(J)

x2 = linspace(0,max(x),500);
y2 = theta(1)+ theta(2)*x2.^(theta(3));

figure()
hold on
scatter(x,y)
plot(x2,y2)


figure()
hold on
plot(x,(h-y).^2)