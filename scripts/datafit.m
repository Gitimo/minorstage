clear all; clf;
 
# Function that will be fit
function [y]=double_exp(x,par)
  y=par(1)*exp(-par(2)*x)-par(3)*exp(-par(4)*x);
end
 
# Generate a double exponent
J0b=0.02;
t1=1/0.0001;
J0p=0.02;
t2=1/0.005;
m=load("../data/2012-10-16/flash_middle_light_off_cut");
y=m(:,2);
x=m(:,1);

xmin=0;
xmax=max(x);
ymin=1.2*min(y);
ymax=1.2*max(y);
# Add some noise to the line

weights=ones(size(x));

 
# Perform the fit
pin=[J0b,t1,J0p,t2];
[f,p,cvg,iter,corp,covp]=leasqr(x,y,pin,"double_exp",.0001,20,weights);

diff=double_exp(x,p)-y;

# Print out the results
p
covp
sum((y-double_exp(x,p)).^2.*weights.^2)
 
# Plots
subplot(2,1,1)
hold on;
%errorbar(x,y,1./weights);
plot(x,double_exp(x,p));
plot(m(:,1),m(:,2),"4");
title "Data and fit";
axis([xmin,xmax,ymin,ymax]);
 
subplot(2,1,2);
plot(m(:,1),diff);
title "Fitted values minus measured values";
axis([xmin,xmax,1.2*min(diff),1.2*max(diff)]);
hold off;

print -deps dummy.eps
