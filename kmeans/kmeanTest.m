close all
clc
clear


File = fopen('rAscii.in','r');
data = fscanf(File,'%f',[5 Inf])';
fclose(File);

rx = data(:,3); 
ry = data(:,4); 
rz = data(:,5); 

figure()
hold on
scatter3(rx,ry,rz,10,'filled')


X = [rx ry rz];
opts = statset('Display','final');
[cidx, ctrs] = kmeans(X, 50, 'Distance','city', ...
    'Replicates',5, 'Options',opts);


figure()
box on
hold on
for id=1:20
    scatter3(X(cidx==id,1),X(cidx==id,2),X(cidx==id,3),10,'filled');
end
xlabel('x')
ylabel('y')
zlabel('z')
side = 1293;

xlim([-side/2 side/2])
ylim([-side/2 side/2])
zlim([-side/2 side/2])

