function [ lboxArr ] = find_lboxArr( stress_strain, box_strain, sidex )
%FIND_LBOXARR Define lboxArr for stress calculations

nStep = length(stress_strain);
ind = 1;
lboxArr = zeros(nStep,1);

for ii=1:length(box_strain)
    if box_strain(ii) == stress_strain(ind,1)
        lboxArr(ind) = sidex(ii);
        if ind == length(stress_strain(:,1))
            break;
        end
        ind = ind+1;
    end
end

end

