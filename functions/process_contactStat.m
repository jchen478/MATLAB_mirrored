function [ NC_total_statData, NC_total_no_jointsData,overlapData, forcData, sijData ] = process_contactStat( filePrefix, basisStrain, ithStat, rspec )
%PROCESS_INTENSITY Find stats for intensity
NC_total_statData = 0;
NC_total_no_jointsData =  0;
overlapData =  0;
forcData =  0;
sijData =  0;
if (exist([filePrefix,'ContactStat.txt'], 'file') ~= 0)
    [Con_strain, contact] = read_contactStat([filePrefix,'ContactStat.txt']);
    NC_total = contact(:,1);
    NC_total_no_joints = contact(:,2);
    overlap=contact(:,3);
    forc =contact(:,4);
    sij =contact(:,5);
    
    if (length(rspec) == 1)
        r = [basisStrain  Con_strain(end)]';
    else
        r = rspec;
    end
    if (Con_strain(end) <= basisStrain)
        return
    end
    NC_total_stat = interval_average(Con_strain,NC_total,r);
    NC_total_no_joints_stat = interval_average(Con_strain,NC_total_no_joints,r);
    overlap_stat = interval_average(Con_strain,overlap,r);
    forc_stat = interval_average(Con_strain,forc,r);
    sij_stat = interval_average(Con_strain,sij,r);
    
    NC_total_statData = NC_total_stat(ithStat,1);
    NC_total_no_jointsData = NC_total_no_joints_stat(ithStat,1);
    overlapData = overlap_stat(ithStat,1);
    forcData = forc_stat(ithStat,1);
    sijData = sij_stat(ithStat,1);
end

end

