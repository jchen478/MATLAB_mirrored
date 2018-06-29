function [ Eelas ] = process_elastic( filePrefix, basisStrain, ithStat, rspec )
%PROCESS_INTENSITY Find stats for intensity
Eelas = 0;
if (exist([filePrefix,'Eelastic.txt'], 'file') ~= 0)
    [Eelas_strain, Eelas_nondim] = read_elastic([filePrefix,'Eelastic.txt']);
    if (length(rspec) == 1)
        r = [basisStrain  Eelas_strain(end)]';
    else
        r = rspec;
    end
    if (Eelas_strain(end) <= basisStrain)
        return
    end
    Eelas_stat = interval_average(Eelas_strain,Eelas_nondim,r);
    Eelas = Eelas_stat(ithStat,1);
end
end

