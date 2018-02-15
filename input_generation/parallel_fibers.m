%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generates input for parallel fibers of three segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all; 

rp = 15;
%% Fiber 1

r2 = [0;0;0];
r1 = [-sqrt(2)*rp;-sqrt(2)*rp;0];
r3 = [sqrt(2)*rp;sqrt(2)*rp;0];
rcm1 = [1/3*(r3(1)+r2(1)+r1(1));1/3*(r3(2)+r2(2)+r1(2));1/3*(r3(3)+r2(3)+r1(3))];

p1 = [1/sqrt(2);1/sqrt(2);0];
p2 = p1;
p3 = p1; 

px = p1(1);
py = p1(2);
pz = p1(3);

% Euler parameters
cc1 = sqrt(1-py*py);
ai11 = pz/cc1;
ai12 = 0.0;
ai13 = -px/cc1;
ai21 = -px*py/cc1;
ai22 = cc1;
ai23 = -py*pz/cc1;

q0 = 0.5*sqrt(1+ai11+ai22+pz);
q1 = 0.25*(ai23-py)/q0;
q2 = 0.25*(px-ai13)/q0;
q3 = 0.25*(ai12-ai21)/q0;

format long E
Eu = [q0; q1; q2; q3]

