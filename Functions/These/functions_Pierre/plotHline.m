function h=plotHline(Y,colors)

XX=xlim;
h=plot([XX(1),XX(2)],[Y,Y],'Color',colors,'linewidth',1.5);


end