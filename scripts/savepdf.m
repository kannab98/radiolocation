function savepdf = savepdf(name)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
print(gcf,name,'-dpdf','-r0');
savepfg=1;
end
