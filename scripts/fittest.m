clear all; clf;

# Define some physical constants
k=1.3806488e-23;
q=1.60217646e-19;
T=293;

# Function that will be fit
function [y]=double_exp(x,par)
  y=par(1)*exp(-par(2)*x)-par(1)*exp(-par(3)*x);
end

# Function to save
function savefig()
s=input("Save the plot?[y/n]","s");
if (s=="y")
	clear plotname;
	plotname=input("Filename for plot?","s");
	print(gcf,plotname,'-deps');
endif
clf;
end

# Generate a double exponent
j=5e-02;
t1= 2.6032e+02;
t2= 4.2686e+02;

#Load data and assign vectors
fhalf="../data/" ; shalf=input("Filename?","s") ; fname=strcat(fhalf,shalf);
m=load(fname);
t=m(:,1) ; ir1=m(:,2) ; ir2=m(:,3) ; it1=m(:,4) ; vt=m(:,5) ;

# Generic weights matrix
weights=ones(size(ir1));
ir1(7764)=(ir1(7763)+ir1(7765))/2;
# Perform the fits
pin=[j,t1,t2];
[fir1norm,pir1norm,cvg,iter,corp,covp]=leasqr(t,ir1,pin,"double_exp",.0001,20,weights);
(sum((fir1norm-ir1).^2)).^0.5
s=input("wait","s");
#hold on;
#hg = errorbar(tcut,ir1cut,0.0005);
#set(hg,"marker","+");
#set(hg,"linestyle","none");
#axis([0,max(t),min(fir1),max(fir1)]);
#plot(t,fir1);
#hold off;
#savefig();

[fir1small,pir1small,cvg,iter,corp,covp]=leasqr(t,ir1,pin,"double_exp",.0001,20,weights./ir1*0.0008);
(sum((fir1small-ir1).^2)).^0.5
s=input("wait","s");
#hold on;
#hg = errorbar(tcut,ir1cut,0.0005);
#set(hg,"marker","+");
#set(hg,"linestyle","none");
#axis([0,max(t),min(fir1),max(fir1)]);
#plot(t,fir1);
#hold off;
#savefig();

[fir1big,pir1big,cvg,iter,corp,covp]=leasqr(t,ir1,pin,"double_exp",.0001,20,weights.*ir1/0.0008);
(sum((fir1big-ir1).^2)).^0.5
s=input("wait","s");
#hold on;
#hg = errorbar(tcut,ir1cut,0.0005);
#set(hg,"marker","+");
#set(hg,"linestyle","none");
#axis([0,max(t),min(fir1),max(fir1)]);
#plot(t,fir1);
#hold off;
#savefig();

