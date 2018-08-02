function [ NC ] = process_NC( filePrefix, basisStrain, ithStat, rspec )
%PROCESS_INTENSITY Find stats for intensity
NC = 0;
if rspec == 0
    return
end
if (exist([filePrefix,'Number_of_Contacts.txt'], 'file') ~= 0)
    [NC_strain, NCread] = read_NC([filePrefix,'Number_of_Contacts.txt']);
    if (max(size(rspec)) == 1)
        r = [basisStrain  NC_strain(end)]';
    else
        r = rspec;
    end
    if (NC_strain(end) <= basisStrain)
        return
    end
    if NC_strain(end) < r(end)
        return
    end
    NC_stat = interval_average(NC_strain,NCread,r);
    NC = NC_stat(ithStat,1);
end
end

