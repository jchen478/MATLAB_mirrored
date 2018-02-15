%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generates input for stacked fibers of three segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all; 

%% Parameters
rp = 15;
theta = 0.6;

%% Orientation vector
Ry = [ cos(theta) 0 sin(theta);
     0 1 0 ;
      -sin(theta) 0 cos(theta)];
  
p2 = [1;0;0];
p3 = Ry*p2;

p1 = p3;
p1(3) = -p1(3);

p1 
p2 
p3

theta = 0.3;
Rz = [ cos(theta) sin(theta) 0;
      -sin(theta) cos(theta) 0;
      0 0 1];
  
p2 = Rz*p2
p3 = Rz*p3
p1 = p3;
p1(3) = -p1(3)


%% Position
r2 = [0;0;0];
r3 = r2+rp*p2+rp*p3;
r1 = r2-rp*p2-rp*p1;

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

Eu1 = [q0; q1; q2; q3];

px = p2(1);
py = p2(2);
pz = p2(3);

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

Eu2 = [q0; q1; q2; q3];

px = p3(1);
py = p3(2);
pz = p3(3);

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

Eu3 = [q0; q1; q2; q3];

%% Center of Mas
rcm = [0.0;0.0;(r1(3)+r2(3)+r3(3))/3];

rcm = [1/3*(r3(1)+r2(1)+r1(1));1/3*(r3(2)+r2(2)+r1(2));1/3*(r3(3)+r2(3)+r1(3))]

%% Output
format long E

r = zeros(3,3);
r(1,:) = r1;
r(2,:) = r2;
r(3,:) = r3;
scatter3(r(:,1),r(:,2),r(:,3))

r1
r2
r3

xlabel('x')
ylabel('y')
zlabel('z')

Eu1
Eu2
Eu3
