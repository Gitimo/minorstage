function octaveplot (datafile)

m=load(datafile);
row=size(m)(1);
col=size(m)(2);
xmin=m(1,1);
xmax=m(row,1);
hold off
ymax=1.25*max(m(:,col));
ymin=-ymax;
plot (m(:,1),m(:,2))
axis ([xmin,xmax,ymin,ymax])
hold on
