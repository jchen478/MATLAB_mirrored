function [ data] = read_ncluster( filename )
%READ_NC read Number_of_Contacts.txt

File = fopen(filename,'r');
data = fscanf(File,'%f',[1 Inf])';
fclose(File);

end

