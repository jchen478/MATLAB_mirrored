function [strain, ncluster, maxCluster, Rg, lam] = read_ncluster( filename )
%READ_NC read Number_of_Contacts.txt

File = fopen(filename,'r');
data = fscanf(File,'%f',[7 Inf])';
fclose(File);
strain = data(:,1);
ncluster = data(:,2);
maxCluster = data(:,3);
Rg = data(:,4);
lam = data(:,5:7); 
end

