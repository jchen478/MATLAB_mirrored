function [ FlocEarr] = read_flocEnergy( filename, nCluster )
%READ_NC read Number_of_Contacts.txt
File = fopen(filename,'r');
data = fscanf(File,'%f',[3 nCluster])';
fclose(File);
FlocEarr = data(:,3);
end

