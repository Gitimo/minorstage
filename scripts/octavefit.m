clear all; clf;
 
# Function that will be fit
function [y]=double_exp(x,par)
  y=par(1)*exp(-par(2)*x)-par(3)*exp(-par(4)*x);
end
 
# Generate a double exponent
J0b=1.25;
t1=1;
J0p=1.25
t2=5;
x=[0:0.1:10]';
y=J0b*exp(-t1*x)-J0p*exp(-t2*x);
 
# Add some noise to the line
sigma=0.1;
weights=ones(size(x))/sigma;
y=y+randn(size(x))*sigma;
 
# Perform the fit
pin=[J0b,t1,J0p,t2];
[f,p,cvg,iter,corp,covp]=leasqr(x,y,pin,"double_exp",.0001,20,weights);
 
# Print out the results
p
covp
sum((y-double_exp(x,p)).^2.*weights.^2)
 
# Plots
subplot(2,1,1)
hold on;
errorbar(x,y,1./weights);
plot(x,double_exp(x,p));
 
subplot(2,1,2);
errorbar(x,y-double_exp(x,p),1./weights);
