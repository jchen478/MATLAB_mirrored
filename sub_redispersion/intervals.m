function [ r ] = intervals( boxstrain,Lbox )
% Determine averaging intervals
% initialization
ind = 1;
comp = Lbox(1);

% region1 - shearing at same strain
for i=2:length(Lbox)
    if Lbox(i) ~= comp
        r(1)=boxstrain(i-1);
        ind = i;
        break
    end
end

% region2 - shearing while concentrating
for i=ind:length(Lbox)
    if Lbox(i) == Lbox(i-1)
        r(2)=boxstrain(i-1);
        ind = i;
        break
    end
end

% region3 - shearing at high concentration
comp = Lbox(ind);
for i=ind:length(Lbox)
    if Lbox(i) ~= comp
        r(3)=boxstrain(i-1);
        ind = i;
        break
    end
end

% region4 - shearing while expanding
for i=ind:length(Lbox)
    if Lbox(i) == Lbox(i-1)
        r(4)=boxstrain(i-1);
        ind = i;
        break
    end
end

% region5 - shearing at low concentration
r(5)=boxstrain(end);

end

