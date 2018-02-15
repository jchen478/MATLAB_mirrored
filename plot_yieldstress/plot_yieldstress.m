close all;
clc;

File = fopen('1280.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f ';
sizeA = [10 Inf];
stress = flipud(fscanf(File,formatSpec,sizeA)');
fclose(File);

yieldStress(stress,'N_{fib} = 1280',1)

