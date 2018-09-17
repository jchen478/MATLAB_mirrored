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
        plot(0,0)
        return
    end
    Eelas_stat = interval_average(Eelas_strain,Eelas_nondim,r);
    Eelas = Eelas_stat(ithStat,1);
%     plot(Eelas_strain, Eelas_nondim);
%     Eelas_nondim
    %% discard regions not needed
    %     if (length(r) > 0)
    %         Eelas_strain (Eelas_strain <= r(4)) = 0;
    %         ind = 0;
    %         for i=length(Eelas_strain):-1:0
    %             if Eelas_strain(i) == 0
    %                 ind = i;
    %                 break;
    %             end
    %         end
    %         Eelas_strain(1:ind) = [];
    %         Eelas(1:ind) = [];
    %         Eelas_strain = Eelas_strain - Eelas_strain(1);
%     subplot(2,1,1)
%     hold on
%     plot(Eelas_strain, Eelas_nondim)
%     subplot(2,1,2)
%     hold on
%     plot(Eelas_strain, Eelas_nondim/Eelas_nondim(1))
    
    Eelas_stat = interval_average(Eelas_strain,Eelas_nondim,r);
    Eelas = Eelas_stat(ithStat,1);
    
    %     end
end
end

