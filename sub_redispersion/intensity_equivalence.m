I = 0.02; 
volfrac = 0.01;

syms R

expr = I*(1-volfrac) - R^3*(1/R^3-1)^2 + (1-R^3)*volfrac;
Rsol = vpasolve(expr==0, R)