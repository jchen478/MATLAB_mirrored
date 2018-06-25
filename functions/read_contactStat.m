function [ strain, contact ] = read_contactStat( filename )
%READ_NC read Number_of_Contacts.txt

File = fopen(filename,'r');
data = fscanf(File,'%f',[9 Inf])';
fclose(File);

strain = data(:,1);
num_groups = data(:,2);
total_contacts = data(:,3);
NC_total = data(:,4);
total_contacts_no_joints = data(:,5);
NC_total_no_joints = data(:,6);
overlap = data(:,7);
forc = data(:,8);
sij = data(:,9);

contact = [NC_total NC_total_no_joints overlap forc sij];

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

