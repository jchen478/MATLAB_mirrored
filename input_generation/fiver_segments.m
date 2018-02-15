%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generates input for a five segment fiber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all; 

%% Parameters
rp = 15;

%% Orientation vector
p1 = [0.6 0.8 0]
p2 = p1; 
p3 = p1;
p4 = p1;
p5 = p1; 

%% Position
r3 = [0;0;1]';
r4 = r3+rp*p3+rp*p4;
r5 = r4+rp*p4+rp*p5;
r2 = r3-rp*p3-rp*p2;
r1 = r2-rp*p2-rp*p1;

r1
r2
r3
r4
r5

%% Euler Parameters
px = p1(1);
py = p1(2);
pz = p1(3);

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

Eu = [q0 q1 q2 q3]


%% Center of Mas
rcm = 1/5*(r1+r2+r3+r4+r5)

%% Output
format long E

Eu