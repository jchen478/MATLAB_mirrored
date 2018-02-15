function [] = saveFigure( f )
nFig = length(f);
for i=1:nFig
    figure(f(i))
    print(['f',num2str(i)],'-dpng')
end
end

