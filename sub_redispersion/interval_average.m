function [ avg, std ] = interval_average(boxstrain,Lbox)
% Final average of A in each region of redispersion

% initialization
ind = 1;
comp = Lbox(1);

% region1 - shearing at same strain
for i=2:length(Lbox)
    if Lbox(i) ~= comp
        region1=boxstrain(i-1);
        ind = i;
        break
    end
end

% region2 - shearing while concentrating
for i=ind:length(Lbox)
    if Lbox(i) == Lbox(i-1)
        region2=boxstrain(i-1);
        ind = i;
        break
    end
end

% region3 - shearing at high concentration
for i=ind:length(Lbox)
    if Lbox(i) = Lbox(i-1)
        region2=boxstrain(i-1);
        ind = i;
        break
    end
end


region1
region2