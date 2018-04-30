%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generates input fibers and the corresponding curl index of specified shape 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all; 

%% Function Inputt
rps = 15;
nseg = 5; 
theta = 0.4;
phi = 0.3; 
p1 = [1 1 1]; 
p1 = p1 / norm(p1);


%% Function return values
p = zeros(nseg,3); 
r = zeros(nseg,3);  

%% Orientation vectors
% for each segment 
% determine equilibrium orientation from previous segment
% find orthonormal basis (x,y)
basis = zeros(3,3,nseg);
R_phi = [cos(phi) -sin(phi) 0 ; sin(phi) cos(phi) 0; 0 0 1]; 
R_theta = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];

x1 = [0 -p1(3) p1(2)];
y1 = cross(p1,x1);
x1 = x1/norm(x1);
y1 = y1/norm(y1);
basis(:,:,1) = [x1 ; y1 ; p1]; % first segment
p(1,:) = p1; 

for i=2:nseg
    basis(:,:,i) = R_theta*R_phi*basis(:,:,i-1);
    p(i,:) = basis(3,:,i);
end

%% Position vectors
r(1,:) = rps*p(1,:);
for i=2:nseg
    r(i,:) = r(i-1,:)+rps*p(i-1,:)+rps*p(i,:);
end

%% curl index
rend = r(end,:) + rps*p(end,:);
e = norm(rend) 
curl = (2*rps*nseg)/e -1

%% plotting
figure
hold on
quiver3(r(:,1),r(:,2),r(:,3),...
    rps*p(:,1),rps*p(:,2),rps*p(:,3),...
    0,'showArrowHead','off',...
    'linewidth',3, 'color',rgb('MediumBlue'))
quiver3(r(:,1),r(:,2),r(:,3),...
    -rps*p(:,1),-rps*p(:,2),-rps*p(:,3),...
    0,'showArrowHead','off',...
    'linewidth',3)