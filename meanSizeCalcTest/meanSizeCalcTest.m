close all
clc
clear

%% read file
File = fopen('rAscii.in','r');
data = fscanf(File,'%f',[5 Inf])';
fclose(File);

%% parameters
side = 1293;
N = length(data);
dr = 0.05;
nBin = side/dr/2+1;

%% allocate pdf arrays
gx = zeros(nBin,1);
gy = zeros(nBin,1);
gz = zeros(nBin,1);
gr = zeros(2*nBin,1);

%% assign variables
rx = data(:,3) + side/2;
ry = data(:,4) + side/2;
rz = data(:,5) + side/2;

%% shift coordinates so that all is in box
shift = rx > side;
rx = rx - side*shift; 
shift = ry > side;
ry = ry - side*shift; 
shift = rz > side;
rz = rz - side*shift; 

%% visualize particle position
figure()
hold on
box on
scatter3(rx,ry,rz,10,'filled','markerfacecolor',rgb('Crimson'))
xlim([0 side])
ylim([0 side])
zlim([0 side])
xlabel('x')
ylabel('y')
zlabel('z')

% %% Calculate P(r) with PBC
% for i=1:N-1
%     for j=i+1:N
%         % distance
%         dx = abs(rx(i)-rx(j));
%         dy = abs(ry(i)-ry(j));
%         dz = abs(rz(i)-rz(j));
%         % correct for PBC
%         if dx > side/2
%             dx = abs(dx - side);
%         end
%         if dy > side/2
%             dy = abs(dy - side);
%         end
%         if dz > side/2
%             dz = abs(dz - side);
%         end
%         r = sqrt(dx*dx+dy*dy+dz*dz); 
%         if r > side /2 
%             r = abs(r - side);
%         end
%         
%         binLocx = ceil (dx / dr) + 1 ;
%         binLocy = ceil (dy / dr) + 1 ;
%         binLocz = ceil (dz / dr) + 1 ;
%         binLocr = ceil (r / dr) + 1 ;
%         %disp([num2str(rx(i)),' ',num2str(rx(j)),' ',num2str(dx),' ',num2str(binLoc)])
%         gx(binLocx) = gx(binLocx) + 2;
%         gy(binLocy) = gy(binLocy) + 2;
%         gz(binLocz) = gz(binLocz) + 2;
%         gr(binLocr) = gr(binLocr) + 2;
%     end
% end

%% check resutls
figure()
plot(1:nBin,gx)

figure()
plot(1:nBin,gy)

figure()
plot(1:nBin,gz)

figure()
plot((1:2*nBin)*dr,gr)

Ax = gx(1:nBin);
Ay = gy(1:nBin);
Az = gz(1:nBin);
Ar = gr; 
rho = N / side^3;
nid = 2*dr*side^2*rho;
Ax = Ax / (N*nid);
Ay = Ay / (N*nid);
Az = Az / (N*nid);

for i=1:length(gr)
   
    r = i*dr; 
    nid = 4*pi*r^2*dr*rho; 
    Ar(i) = Ar(i) / (N*nid); 
    
end

figure()
subplot(2,1,1)
hold on
scatter(1:nBin,gx,10,'filled','MarkerEdgeColor',rgb('Crimson'),'MarkerFaceColor',rgb('Crimson'))
scatter(1:nBin,gy,10,'filled','MarkerEdgeColor',rgb('Orange'),'MarkerFaceColor',rgb('Orange'))
scatter(1:nBin,gz,10,'filled','MarkerEdgeColor',rgb('Green'),'MarkerFaceColor',rgb('Green'))
legend('x','y','z')

subplot(2,1,2)
hold on
box on
scatter(0:dr:side/2,Ax,10,'filled','MarkerEdgeColor',rgb('Crimson'),'MarkerFaceColor',rgb('Crimson'));
scatter(0:dr:side/2,Ay,10,'filled','MarkerEdgeColor',rgb('Orange'),'MarkerFaceColor',rgb('Orange'));
scatter(0:dr:side/2,Az,10,'filled','MarkerEdgeColor',rgb('Green'),'MarkerFaceColor',rgb('Green'));
legend('$i = x$','$i = y$','$i = z$','location','best')
ylabel('$g_i$')
xlabel('Distance in $i$ direction')
xlim([0 inf])

figure()
plot((1:nBin)*dr,Ar(1:nBin),'color',rgb('MediumBlue'))
