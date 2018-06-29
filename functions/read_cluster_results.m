function [ Sarr] = read_cluster_results( filename, nCluster )
%READ_NC read Number_of_Contacts.txt
File = fopen(filename,'r');
data = fscanf(File,'%f',[2 nCluster])';
fclose(File);
Sarr = data(:,2);
end

