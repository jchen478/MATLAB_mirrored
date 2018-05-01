function [ strain, sidex, sidey, sidez ] = read_box( filename )
%READ_BOX read Lbox.txt

File = fopen(filename,'r');
data = fscanf(File,'%f',[7 Inf])';
fclose(File);
strain = data(:,1);
sidex = data(:,2);
sidey = data(:,3);
sidez = data(:,4);

strain = round(strain,1); 

end

