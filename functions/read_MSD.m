function [ strain, MSD ] = read_MSD( filename )
%READ_NC read Number_of_Contacts.txt

File = fopen(filename,'r');
data = fscanf(File,'%f',[3 Inf])';
fclose(File);

strain = data(:,1);
MSD(:,1) = data(:,2);
MSD(:,2) = data(:,3);

diff = data(2,1) - data(1,1);
nStep = length(strain);
% make sure that strain is increasing
for ii=2:nStep
    if (strain(ii) < strain(ii-1))
        strain(ii) = strain(ii-1) + diff;
    end
end

% round strain to nearest decimal point for later comparison
strain = round(strain,1); 

end

