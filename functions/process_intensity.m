function [ I ] = process_intensity( filePrefix, basisStrain, ithStat, rspec )
%PROCESS_INTENSITY Find stats for intensity
I = 0; 
if (exist([filePrefix,'Intensity.txt'], 'file') ~= 0)
    [I_strain, I] = read_intensity([filePrefix,'Intensity.txt']);
    if I_strain(end) < basisStrain
        I = 0; 
        return
    end
    if (length(rspec) == 1)
        r = [basisStrain  I_strain(end)]';
    else
        r = rspec;
    end
    I_stat = interval_average(I_strain,I,r);
    I = I_stat(ithStat,1);
end
end

