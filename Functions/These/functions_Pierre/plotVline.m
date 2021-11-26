function h=plotVline(X,colors)

YY=ylim;
h=plot([X,X],[YY(1),YY(2)],'Color',colors,'linewidth',1.5);


end