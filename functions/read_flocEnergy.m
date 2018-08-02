function [ NfibInFloc, FlocEarr] = read_flocEnergy( filename, nCluster )
%READ_NC read Number_of_Contacts.txt
File = fopen(filename,'r');
data = fscanf(File,'%f',[3 nCluster])';
fclose(File);
NfibInFloc = data(:,2);
% FlocEarr = data(:,3);
FlocEarr = 0;
for i=1:length(NfibInFloc)
    if NfibInFloc(i) >= 1280*0.9
        FlocEarr = data(i,3);
    end
end

end

