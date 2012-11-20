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
 
# Perform the fits
pin=[j,t1,t2];
[fir1,pir1,cvg,iter,corp,covp]=leasqr(t,ir1,pin,"double_exp",.0001,20,weights);
[fir2,pir2,cvg,iter,corp,covp]=leasqr(t,ir2,pin,"double_exp",.0001,20,weights);
[fit1,pit1,cvg,iter,corp,covp]=leasqr(t,it1,pin,"double_exp",.0001,20,weights);

# Calculate reference cell's I_sc - is that correct?
it2=fit1.*fir2./fir1;

# Calculate illumination in "suns" and convert voltages to currents
r=0.04248 ; equiv=13.7e-3 ; illu=fir2/(r*equiv);
ir1=ir1/r; ir2=ir2/r;
it1=it1/0.05; it2=it2/0.05; #Have to find out Resistance precisely!

# Function to fit J-v

function [J]=currentvoltage(V,par)
k=1.3806488e-23;
q=1.60217646e-19;
T=293;
	J=par(1)*(exp(q*V/(k*T))-1)+par(2)*(exp(q*V/(2*k*T))-1);
end
weight=ones(size(vt(100:8000)));

[J,param,cvg,iter,corp,covp]=leasqr(vt(100:8000),it2(100:8000),[1e-21,1e-11],"currentvoltage",.0001,20,weight);

#Plotting ... saving does not work properly
tcut(1:floor(8000/100)+1)=0;
Jcut(1:floor(8000/100)+1)=0;

for m=1:8000
	if (m==30)
		vtcut(m)=vt(m); it2cut(m)=it2(m);
	elseif (mod(m,100)==0)
		vtcut(m)=vt(m); it2cut(m)=it2(m);
	end
end

hold on;
hg = errorbar(vtcut,it2cut,0.0005);
set(hg,"marker","+");
set(hg,"linestyle","none");
axis([min(vt(100:8000)),max(vt(100:8000)),min(J),max(J)]);
plot(vt(100:8000),J);
xlabel "Voc in V";
ylabel "J sc";
secondj=num2str(param(2));
firstj=num2str(param(1));
fhalf=strcat(firstj,"*(exp(qV/kT)-1) + ");
shalf=strcat(fhalf,secondj); 
tot=strcat(shalf,"*(exp(qV/2kT)-1)");
tit=strcat("Fit: ",tot);
title(tit);
savefig();
hold off;

