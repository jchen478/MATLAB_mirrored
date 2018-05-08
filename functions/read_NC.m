function [ strain, NC ] = read_NC( filename )
%READ_NC read Number_of_Contacts.txt

File = fopen(filename,'r');
data = fscanf(File,'%f',[5 Inf])';
fclose(File);

strain = data(:,1);
num_groups = data(:,2);
total_contacts = data(:,3);
NC = data(:,4);
overlap = data(:,5);

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

