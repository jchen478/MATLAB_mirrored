clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Generate random positions
%%% (later read from file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%{
nP = 2095*27;
File = fopen('r_1280_20_1500_1.txt','r');
data = fscanf(File,'%f',[5 nP])';
fclose(File);
rx = data(:,3);
ry = data(:,4);
rz = data(:,5); 
figure()
hold on
% box on
scatter3(rx(1:2095),ry(1:2095),rz(1:2095),'filled','markerFaceColor',rgb('Crimson'))
scatter3(rx(2096:end),ry(2096:end),rz(2096:end),'filled','markerFaceColor',rgb('MediumBlue'))
xlabel('x')
ylabel('y')
zlabel('z')
%}

%%{
nP = 540*27;
File = fopen('rAscii_extend.in','r');
data = fscanf(File,'%f',[5 nP])';
fclose(File);
rx = data(:,3);
ry = data(:,4);
rz = data(:,5);   
figure()
hold on 
% box on
scatter3(rx(1:540),ry(1:540),rz(1:540),'filled','markerFaceColor',rgb('Crimson'))
scatter3(rx(541:end),ry(541:end),rz(541:end),'filled','markerFaceColor',rgb('MediumBlue'))
xlabel('x')
ylabel('y')
zlabel('z')
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Extend in all 26 boxes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rx_extend = zeros(27*nP,1); 
ry_extend = zeros(27*nP,1); 
rz_extend = zeros(27*nP,1); 

%%% 1. Central box
range = 1:nP; 
rx_extend(range) = rx; 
ry_extend(range) = ry; 
rz_extend(range) = rz; 

%%% 2. Right box
range = nP+1:2*nP; 
rx_extend(range) = rx+side; 
ry_extend(range) = ry; 
rz_extend(range) = rz; 

%%% 3. Left box
range = 2*nP+1:3*nP; 
rx_extend(range) = rx-side; 
ry_extend(range) = ry; 
rz_extend(range) = rz; 

%%% 4. Top box
range = 3*nP+1:4*nP; 
rx_extend(range) = rx; 
ry_extend(range) = ry+side; 
rz_extend(range) = rz; 

%%% 5. Bot box
range = 4*nP+1:5*nP; 
rx_extend(range) = rx; 
ry_extend(range) = ry-side; 
rz_extend(range) = rz; 

%%% 6. Front box
range = 5*nP+1:6*nP; 
rx_extend(range) = rx; 
ry_extend(range) = ry; 
rz_extend(range) = rz+side; 

%%% 7. Back box
range = 6*nP+1:7*nP; 
rx_extend(range) = rx; 
ry_extend(range) = ry; 
rz_extend(range) = rz-side; 

%%% 8. Top right box
range = 7*nP+1:8*nP; 
rx_extend(range) = rx+side; 
ry_extend(range) = ry+side; 
rz_extend(range) = rz; 

%%% 9. Top left box
range = 8*nP+1:9*nP; 
rx_extend(range) = rx-side; 
ry_extend(range) = ry+side; 
rz_extend(range) = rz; 

%%% 10. Bottom right box
range = 9*nP+1:10*nP; 
rx_extend(range) = rx+side; 
ry_extend(range) = ry-side; 
rz_extend(range) = rz; 

%%% 11. Bottom left box
range = 10*nP+1:11*nP; 
rx_extend(range) = rx-side; 
ry_extend(range) = ry-side; 
rz_extend(range) = rz; 

%%% 12. Top front box
range = 11*nP+1:12*nP; 
rx_extend(range) = rx; 
ry_extend(range) = ry+side; 
rz_extend(range) = rz+side; 

%%% 13. Top back box
range = 12*nP+1:13*nP; 
rx_extend(range) = rx; 
ry_extend(range) = ry+side; 
rz_extend(range) = rz-side; 

%%% 14. Bottom front box
range = 13*nP+1:14*nP; 
rx_extend(range) = rx; 
ry_extend(range) = ry-side; 
rz_extend(range) = rz+side; 

%%% 15. Bottom back box
range = 14*nP+1:15*nP; 
rx_extend(range) = rx; 
ry_extend(range) = ry-side; 
rz_extend(range) = rz-side; 

%%% 16. Right front box
range = 15*nP+1:16*nP; 
rx_extend(range) = rx+side; 
ry_extend(range) = ry; 
rz_extend(range) = rz+side; 

%%% 17. Right back box
range = 16*nP+1:17*nP; 
rx_extend(range) = rx+side; 
ry_extend(range) = ry; 
rz_extend(range) = rz-side; 

%%% 18. Left front box
range = 17*nP+1:18*nP; 
rx_extend(range) = rx-side; 
ry_extend(range) = ry; 
rz_extend(range) = rz+side; 

%%% 19. Left back box
range = 18*nP+1:19*nP; 
rx_extend(range) = rx-side; 
ry_extend(range) = ry; 
rz_extend(range) = rz-side; 

%%% 20. Corner 1
range = 19*nP+1:20*nP; 
rx_extend(range) = rx+side; 
ry_extend(range) = ry+side; 
rz_extend(range) = rz+side; 

%%% 21. Corner 2
range = 20*nP+1:21*nP; 
rx_extend(range) = rx+side; 
ry_extend(range) = ry+side; 
rz_extend(range) = rz-side; 

%%% 22. Corner 3
range = 21*nP+1:22*nP; 
rx_extend(range) = rx+side; 
ry_extend(range) = ry-side; 
rz_extend(range) = rz+side; 

%%% 23. Corner 4
range = 22*nP+1:23*nP; 
rx_extend(range) = rx+side; 
ry_extend(range) = ry-side; 
rz_extend(range) = rz-side; 

%%% 24. Corner 5
range = 23*nP+1:24*nP; 
rx_extend(range) = rx-side; 
ry_extend(range) = ry+side; 
rz_extend(range) = rz+side; 

%%% 25. Corner 6
range = 24*nP+1:25*nP; 
rx_extend(range) = rx-side; 
ry_extend(range) = ry+side; 
rz_extend(range) = rz-side; 

%%% 26. Corner 7
range = 25*nP+1:26*nP; 
rx_extend(range) = rx-side; 
ry_extend(range) = ry-side; 
rz_extend(range) = rz+side; 

%%% 27. Corner 8
range = 26*nP+1:27*nP; 
rx_extend(range) = rx-side; 
ry_extend(range) = ry-side; 
rz_extend(range) = rz-side; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Plot and compare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure()
hold on
% box on
scatter3(rx,ry,rz,'filled','markerFaceColor',rgb('Crimson'))
xlabel('x')
ylabel('y')
zlabel('z')
xlim([-side/2 side/2])
ylim([-side/2 side/2])
zlim([-side/2 side/2])

figure()
hold on
% box on
scatter3(rx_extend(1:nP),ry_extend(1:nP),rz_extend(1:nP),'filled','markerFaceColor',rgb('Crimson'))
scatter3(rx_extend(nP+1:end),ry_extend(nP+1:end),rz_extend(nP+1:end),'filled','markerFaceColor',rgb('MediumBlue'))
xlabel('x')
ylabel('y')
zlabel('z')
xlim([-3*side/2 3*side/2])
ylim([-3*side/2 3*side/2])
zlim([-3*side/2 3*side/2])

%}