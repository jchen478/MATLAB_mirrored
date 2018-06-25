function [] = saveAsPDF( f )
nFig = length(f);
for i=1:nFig
    figure(f(i))
    h = figure(f(i));
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,['f',num2str(i)],'-dpdf','-r0')
end
end

