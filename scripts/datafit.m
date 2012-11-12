clf;
 
# Function that will be fit
function [y]=double_exp(x,par)
  y=par(1)*exp(-par(2)*x)-par(1)*exp(-par(3)*x);
end
 
# Generate a double exponent
j=5e-02;
t1= 2.6032e+02;
t2= 4.2686e+02;
m=load("../data/2012-11-12/noARC/testing/IV.data");
y=m(:,2);
x=m(:,1);

xmin=0;
xmax=max(x);
ymin=1.2*min(y);
ymax=1.2*max(y);
# Add some noise to the line

weights=ones(size(y));%.*abs((y/0.0008));
 
# Perform the fit
pin=[j,t1,t2];
[f,p,cvg,iter,corp,covp]=leasqr(x,y,pin,"double_exp",.0001,20,weights);

diff=double_exp(x,p)-y;

# Print out the results
p
covp
sum(abs((y-double_exp(x,p))))
 
# Plots
subplot(2,1,1)
hold on;
%errorbar(x,y,1./weights);
plot(x,double_exp(x,p));
plot(m(:,1),m(:,2),"4");
title "Data and fit";
axis([xmin,xmax,ymin,ymax]);
 
subplot(2,1,2);
plot(m(:,1),diff(:));
title "Fitted values minus measured values";
axis([xmin,xmax,1.2*min(diff),1.2*max(diff)]);
hold off;

print -deps dummy.eps
