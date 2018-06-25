clc; 
close all;
clear;

%% Interfiber normal forces parameters
fstar = 150; 
fact = 20;
Astar = 100;
decatt = 35;
rep_cut = 0.66; 
overlap = 1.85; 

eta0 = 1; % Pa s
L = 2.50E-3; % m 
b = 1.6E-5; % m
gamma = 10; % 1/s
nseg = 5;
l = L/(2*nseg);

prefactor = -6*pi*l*b*gamma*eta0*10^6; % muN

g = (1.85:0.005:2.66)'; % interfiber distance from center of fiber
h = g-2; % interfiber distance from fiber surface

nH = length(h);
nAstar = length(Astar);

%% Interfiber maximum forces

syms sij Frep(sij) Fatt(sij) U(sij) Ftot(sij) f g 

digits(5)
Frep = prefactor*fstar.*exp(-fact*sij);
figure(1) 
hold on
title('Frep')
plot(h,vpa(subs(Frep,h)))

figure(2)
hold on
title('Fatt')

figure(3)
hold on
title('U')

for i=1:nAstar
    f = prefactor*(-1)*Astar(i)*exp(-decatt*sij.^2);
    g = prefactor*(-1)*Astar(i);
    Fatt(sij) = piecewise(sij<0,g,sij>=0,f);
    U = int(Frep+Fatt,sij);
    
    figure(2) 
    hold on
    plot(h,vpa(subs(Fatt,h)))
    
    figure(3)
    hold on
    plot(h,vpa(subs(U,h)))
    
end
