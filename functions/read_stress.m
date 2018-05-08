function [ strain, sigxz_out ] = read_stress( filename )
%READ_STRESS read Stress_tensor.txt

File = fopen(filename,'r');
data = fscanf(File,'%f',[7 Inf])';
fclose(File);

strain = data(:,1);
sigxx_out = data(:,2);
sigyx_out = data(:,3);
sigxz_out = data(:,4);
sigyy_out = data(:,5);
sigyz_out = data(:,6);
sigzz_out = data(:,7);
N1_out = sigxx_out-sigzz_out;
N2_out = sigzz_out-sigyy_out;

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

